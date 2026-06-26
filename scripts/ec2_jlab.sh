#!/usr/bin/env bash
set -euo pipefail

APP_NAME="retrotransposon-miner"
INSTANCE_NAME="${APP_NAME}-wgs"
INSTANCE_TYPE="${INSTANCE_TYPE:-r6i.4xlarge}"
HOST_ALIAS="retro-ec2"
JLAB_ALIAS="jlab"
JLAB_LOCAL_PORT="${JLAB_LOCAL_PORT:-8890}"
REGION="$(aws configure get region || true)"
REGION="${REGION:-us-east-1}"

mkdir -p "${HOME}/.ssh"
KEY_BASENAME="${APP_NAME}-${REGION}"
KEY_PATH="${HOME}/.ssh/${KEY_BASENAME}.pem"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "Missing command: $1"; exit 1; }
}

awsq() {
  aws --region "${REGION}" "$@"
}

imds_get() {
  local path="$1"
  local token
  token="$(curl -fsS -X PUT "http://169.254.169.254/latest/api/token" \
    -H "X-aws-ec2-metadata-token-ttl-seconds: 60" 2>/dev/null || true)"
  [[ -n "${token}" ]] || return 1
  curl -fsS "http://169.254.169.254/latest/${path}" \
    -H "X-aws-ec2-metadata-token: ${token}" 2>/dev/null
}

get_latest_al2023_ami() {
  awsq ssm get-parameter \
    --name /aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64 \
    --query 'Parameter.Value' --output text
}

get_default_vpc() {
  awsq ec2 describe-vpcs \
    --filters Name=isDefault,Values=true \
    --query 'Vpcs[0].VpcId' --output text
}

get_default_subnet() {
  local vpc_id
  vpc_id="$(get_default_vpc)"
  awsq ec2 describe-subnets \
    --filters Name=vpc-id,Values="${vpc_id}" Name=default-for-az,Values=true \
    --query 'Subnets[0].SubnetId' --output text
}

ensure_key_pair() {
  if [[ -f "${KEY_PATH}" ]]; then
    chmod 400 "${KEY_PATH}"
    log "Using existing key file: ${KEY_PATH}"
    return
  fi

  if awsq ec2 describe-key-pairs --key-names "${KEY_BASENAME}" >/dev/null 2>&1; then
    log "KeyPair '${KEY_BASENAME}' exists in AWS, but local PEM not found."
    log "Create/import a new key pair manually or delete old key pair and rerun."
    exit 1
  fi

  log "Creating key pair: ${KEY_BASENAME}"
  awsq ec2 create-key-pair \
    --key-name "${KEY_BASENAME}" \
    --query 'KeyMaterial' --output text > "${KEY_PATH}"
  chmod 400 "${KEY_PATH}"
}

ensure_security_group() {
  local vpc_id sg_name sg_id my_ip
  vpc_id="$(get_default_vpc)"
  sg_name="${APP_NAME}-ssh-sg"

  sg_id="$(awsq ec2 describe-security-groups \
    --filters Name=group-name,Values="${sg_name}" Name=vpc-id,Values="${vpc_id}" \
    --query 'SecurityGroups[0].GroupId' --output text 2>/dev/null || true)"

  if [[ -z "${sg_id}" || "${sg_id}" == "None" ]]; then
    log "Creating security group: ${sg_name}"
    sg_id="$(awsq ec2 create-security-group \
      --group-name "${sg_name}" \
      --description "SSH access for ${APP_NAME}" \
      --vpc-id "${vpc_id}" \
      --query 'GroupId' --output text)"
  fi

  my_ip="$(curl -s https://checkip.amazonaws.com)/32"
  log "Ensuring SSH ingress from ${my_ip} on ${sg_id}"
  awsq ec2 authorize-security-group-ingress \
    --group-id "${sg_id}" \
    --ip-permissions "IpProtocol=tcp,FromPort=22,ToPort=22,IpRanges=[{CidrIp=${my_ip},Description=local-ssh}]" \
    >/dev/null 2>&1 || true

  echo "${sg_id}"
}

get_instance_id() {
  # Primary: explicit Name-tag match used by bootstrap lifecycle.
  local iid
  iid="$(
    awsq ec2 describe-instances \
      --filters Name=tag:Name,Values="${INSTANCE_NAME}" Name=instance-state-name,Values=pending,running,stopping,stopped \
      --query 'Reservations[].Instances[].InstanceId' \
      --output text 2>/dev/null | awk 'NF{print $1; exit}'
  )"
  if [[ -n "${iid}" && "${iid}" != "None" ]]; then
    echo "${iid}"
    return 0
  fi

  # Fallback 1: resolve the configured SSH alias host IP and map it to an instance.
  local host_ip=""
  if command -v ssh >/dev/null 2>&1; then
    host_ip="$(ssh -G "${HOST_ALIAS}" 2>/dev/null | awk '$1=="hostname"{print $2; exit}' || true)"
  fi
  if [[ -n "${host_ip}" ]]; then
    iid="$(
      awsq ec2 describe-instances \
        --filters Name=ip-address,Values="${host_ip}" Name=instance-state-name,Values=pending,running,stopping,stopped \
        --query 'Reservations[].Instances[].InstanceId' \
        --output text 2>/dev/null | awk 'NF{print $1; exit}'
    )"
    if [[ -n "${iid}" && "${iid}" != "None" ]]; then
      echo "${iid}"
      return 0
    fi
  fi

  # Fallback 2: when running on EC2, use current instance metadata.
  iid="$(imds_get "meta-data/instance-id" || true)"
  if [[ -n "${iid}" ]]; then
    echo "${iid}"
    return 0
  fi

  return 1
}

ensure_instance() {
  local instance_id ami subnet sg_id
  instance_id="$(get_instance_id || true)"

  if [[ -n "${instance_id}" ]]; then
    log "Found existing instance: ${instance_id}"
    echo "${instance_id}"
    return
  fi

  ami="$(get_latest_al2023_ami)"
  subnet="$(get_default_subnet)"
  sg_id="$(ensure_security_group)"
  ensure_key_pair

  log "Creating instance ${INSTANCE_NAME} (${INSTANCE_TYPE})"
  instance_id="$(awsq ec2 run-instances \
    --image-id "${ami}" \
    --instance-type "${INSTANCE_TYPE}" \
    --key-name "${KEY_BASENAME}" \
    --subnet-id "${subnet}" \
    --security-group-ids "${sg_id}" \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=${INSTANCE_NAME}},{Key=App,Value=${APP_NAME}}]" \
    --count 1 \
    --query 'Instances[0].InstanceId' --output text)"

  echo "${instance_id}"
}

start_instance() {
  local iid
  iid="$(ensure_instance)"
  log "Starting instance ${iid}"
  awsq ec2 start-instances --instance-ids "${iid}" >/dev/null || true
  awsq ec2 wait instance-running --instance-ids "${iid}"
  awsq ec2 wait instance-status-ok --instance-ids "${iid}"
  log "Instance is running and healthy."
}

stop_instance() {
  local iid
  iid="$(get_instance_id || true)"
  [[ -n "${iid}" ]] || { log "No instance found."; exit 1; }
  log "Stopping ${iid}"
  if ! awsq ec2 stop-instances --instance-ids "${iid}" >/dev/null 2>&1; then
    log "Unable to stop via EC2 API with current credentials."
    log "If you're currently on this EC2 host, use: sudo shutdown -h now"
    exit 1
  fi
  awsq ec2 wait instance-stopped --instance-ids "${iid}"
  log "Instance stopped."
}

reboot_instance() {
  local iid
  iid="$(get_instance_id)"
  [[ -n "${iid}" ]] || { log "No instance found."; exit 1; }
  log "Rebooting ${iid}"
  awsq ec2 reboot-instances --instance-ids "${iid}"
  awsq ec2 wait instance-status-ok --instance-ids "${iid}"
  log "Instance healthy after reboot."
}

ensure_eip() {
  local iid alloc_id assoc_id
  iid="$(get_instance_id)"
  [[ -n "${iid}" ]] || { log "No instance found."; exit 1; }

  alloc_id="$(awsq ec2 describe-addresses \
    --filters Name=tag:Name,Values="${INSTANCE_NAME}-eip" \
    --query 'Addresses[0].AllocationId' --output text 2>/dev/null || true)"

  if [[ -z "${alloc_id}" || "${alloc_id}" == "None" ]]; then
    log "Allocating Elastic IP"
    alloc_id="$(awsq ec2 allocate-address --domain vpc --query 'AllocationId' --output text)"
    awsq ec2 create-tags --resources "${alloc_id}" --tags "Key=Name,Value=${INSTANCE_NAME}-eip" >/dev/null
  fi

  log "Associating EIP ${alloc_id} to ${iid}"
  assoc_id="$(awsq ec2 associate-address --instance-id "${iid}" --allocation-id "${alloc_id}" --allow-reassociation --query 'AssociationId' --output text)"
  log "EIP association: ${assoc_id}"
}

get_public_ip() {
  local iid
  iid="$(get_instance_id)"
  awsq ec2 describe-instances \
    --instance-ids "${iid}" \
    --query 'Reservations[0].Instances[0].PublicIpAddress' \
    --output text
}

write_ssh_config() {
  local ip cfg
  ip="$(get_public_ip)"
  cfg="${HOME}/.ssh/config"
  touch "${cfg}"

  awk '
    BEGIN{skip=0}
    /^Host retro-ec2$/ {skip=1; next}
    /^Host jlab$/ {skip=1; next}
    /^Host / {skip=0}
    skip==0 {print}
  ' "${cfg}" > "${cfg}.tmp"
  mv "${cfg}.tmp" "${cfg}"

  cat >> "${cfg}" <<EOF

Host ${HOST_ALIAS}
  HostName ${ip}
  User ec2-user
  IdentityFile ${KEY_PATH}
  IdentitiesOnly yes
  ServerAliveInterval 30
  ServerAliveCountMax 6
  TCPKeepAlive yes
  StrictHostKeyChecking accept-new

Host ${JLAB_ALIAS}
  HostName ${ip}
  User ec2-user
  IdentityFile ${KEY_PATH}
  IdentitiesOnly yes
  LocalForward ${JLAB_LOCAL_PORT} 127.0.0.1:8888
  ServerAliveInterval 30
  ServerAliveCountMax 6
  TCPKeepAlive yes
  StrictHostKeyChecking accept-new
EOF

  log "Updated ${cfg} with ${HOST_ALIAS} and ${JLAB_ALIAS}."
}

start_jlab() {
  ssh "${HOST_ALIAS}" "python3 -m pip -q install --user jupyterlab >/dev/null 2>&1 || true; mkdir -p ~/.jlab; [ -f ~/.jlab/token ] || openssl rand -hex 24 > ~/.jlab/token; TOKEN=\$(cat ~/.jlab/token); nohup jupyter lab --no-browser --ip 127.0.0.1 --port 8888 --ServerApp.token=\"\$TOKEN\" > ~/jupyter.log 2>&1 & sleep 3; echo \"TOKEN=\$TOKEN\"; jupyter server list"
}

stop_jlab() {
  ssh "${HOST_ALIAS}" "pkill -f 'jupyter-lab|jupyter lab' || true"
  log "JupyterLab stopped."
}

start_tunnel() {
  pkill -f "ssh -N ${JLAB_ALIAS}" >/dev/null 2>&1 || true
  nohup ssh -N "${JLAB_ALIAS}" >/tmp/jlab-tunnel.log 2>&1 &
  sleep 1
  log "Tunnel running on http://127.0.0.1:${JLAB_LOCAL_PORT}"
}

status() {
  local iid
  iid="$(get_instance_id || true)"
  if [[ -z "${iid}" ]]; then
    log "No instance found for Name=${INSTANCE_NAME}"
    exit 0
  fi
  awsq ec2 describe-instances \
    --instance-ids "${iid}" \
    --query 'Reservations[0].Instances[0].{InstanceId:InstanceId,State:State.Name,PublicIp:PublicIpAddress,Type:InstanceType,Name:Tags[?Key==`Name`]|[0].Value}' \
    --output table
}

bootstrap() {
  require_cmd aws
  require_cmd ssh
  require_cmd curl

  ensure_key_pair
  ensure_security_group >/dev/null
  start_instance
  ensure_eip
  write_ssh_config
  status
  start_jlab
  start_tunnel

  local token
  token="$(ssh "${HOST_ALIAS}" "cat ~/.jlab/token")"
  echo
  echo "Open JupyterLab:"
  echo "http://127.0.0.1:${JLAB_LOCAL_PORT}/lab?token=${token}"
}

case "${1:-}" in
  bootstrap) bootstrap ;;
  start-instance) start_instance ;;
  stop-instance) stop_instance ;;
  reboot-instance) reboot_instance ;;
  ensure-eip) ensure_eip ;;
  ssh-config) write_ssh_config ;;
  start-jlab) start_jlab ;;
  stop-jlab) stop_jlab ;;
  start-tunnel) start_tunnel ;;
  status) status ;;
  *)
    echo "Usage: $0 {bootstrap|start-instance|stop-instance|reboot-instance|ensure-eip|ssh-config|start-jlab|stop-jlab|start-tunnel|status}"
    exit 1
    ;;
esac
