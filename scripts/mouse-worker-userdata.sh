#!/bin/bash
# Cloud-init for one extra mouse aligner. The accession is the instance tag.
set -euxo pipefail
exec > /var/log/mouse-worker-setup.log 2>&1
dnf install -y xorg-x11-server-Xvfb xorg-x11-xauth || true
aws s3 cp s3://l1tx-data/mouse/align/bootstrap/setup_mouse_align_worker.sh /tmp/setup_mouse_align_worker.sh
chmod 755 /tmp/setup_mouse_align_worker.sh
TOKEN="$(curl -fsS -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 300")"
ACC="$(curl -fsS -H "X-aws-ec2-metadata-token: ${TOKEN}" "http://169.254.169.254/latest/meta-data/tags/instance/Accession")"
install -d -o ec2-user -g ec2-user /home/ec2-user
printf '%s\n' "${ACC}" > /home/ec2-user/accession.txt
chown ec2-user:ec2-user /home/ec2-user/accession.txt
sudo -u ec2-user bash /tmp/setup_mouse_align_worker.sh
