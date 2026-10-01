#!/usr/bin/env bash
# Install snpEff 5.4c and the GRCh38.115 database.
#
# Order:
#   1. Already unpacked under the prefix: write the launcher and stop.
#   2. A local tarball (for example after the public-data S3 sync): unpack it.
#   3. The project S3 object, if it exists: download that, do not hit upstream.
#   4. Otherwise build the bundle from bioconda + snpeff-public and upload it.
#
#   s3://l1tx-data/public/tools/snpeff/snpeff-5.4.0c-GRCh38.115.tar.gz
#
# Requires Java 21. The rtm-miner conda env provides it.
#   bash scripts/install_snpeff.sh
#   bash scripts/install_snpeff.sh /path/to/prefix
set -euo pipefail

PREFIX="${1:-${SNPEFF_PREFIX:-${HOME}/snpEff}}"
S3_URI="${SNPEFF_S3_URI:-s3://l1tx-data/public/tools/snpeff/snpeff-5.4.0c-GRCh38.115.tar.gz}"
ARCHIVE_NAME="snpeff-5.4.0c-GRCh38.115.tar.gz"
LOCAL_ARCHIVE="${PREFIX}/${ARCHIVE_NAME}"
BIOCONDA_URL="${SNPEFF_BIOCONDA_URL:-https://conda.anaconda.org/bioconda/noarch/snpeff-5.4.0c-hdfd78af_0.conda}"
DATABASE_URL="${SNPEFF_DATABASE_URL:-https://snpeff-public.s3.amazonaws.com/databases/v5_4/snpEff_v5_4_GRCh38.115.zip}"

require_java21() {
  if ! command -v java >/dev/null 2>&1; then
    echo "WARN: Java 21 is not on PATH. snpEff files were installed; annotate-genes needs Java 21." >&2
    return 0
  fi
  local spec major
  spec="$(java -XshowSettings:properties -version 2>&1 | awk -F= '/java.specification.version/{gsub(/ /,"",$2); print $2}')"
  major="${spec%%.*}"
  if [[ -z "${major}" || "${major}" -lt 21 ]]; then
    echo "WARN: Java ${spec:-unknown} is too old for snpEff 5.4c. Need Java 21 or newer." >&2
  fi
}

s3_object_exists() {
  command -v aws >/dev/null 2>&1 && aws s3 ls "${S3_URI}" >/dev/null 2>&1
}

write_launcher() {
  mkdir -p "${PREFIX}/bin"
  cat > "${PREFIX}/bin/snpEff" <<EOF
#!/bin/bash
xmx=4g
args=()
for a in "\$@"; do
  case "\$a" in
    -Xmx*) xmx="\${a#-Xmx}" ;;
    *) args+=("\$a") ;;
  esac
done
cd "${PREFIX}"
exec java -Xmx"\$xmx" -jar "${PREFIX}/snpEff.jar" "\${args[@]}"
EOF
  chmod +x "${PREFIX}/bin/snpEff"
}

unpack_archive() {
  echo "[snpeff] unpacking ${1}"
  tar -C "${PREFIX}" -xzf "${1}"
}

build_upstream_bundle() {
  echo "[snpeff] S3 cache miss ${S3_URI}; building from bioconda and snpeff-public"
  local work conda_pkg db_zip
  work="$(mktemp -d)"
  conda_pkg="${work}/snpeff.conda"
  db_zip="${work}/GRCh38.115.zip"
  curl -fsSL -o "${conda_pkg}" "${BIOCONDA_URL}"
  curl -fsSL -o "${db_zip}" "${DATABASE_URL}"
  if ! command -v zstd >/dev/null 2>&1; then
    echo "ERROR: zstd is required to unpack the bioconda snpEff package on a cache miss." >&2
    exit 1
  fi
  python3 - "${conda_pkg}" "${work}/pkg.tar.zst" <<'PY'
import sys, zipfile
from pathlib import Path
conda_pkg, dest = Path(sys.argv[1]), Path(sys.argv[2])
with zipfile.ZipFile(conda_pkg) as zf:
    name = next(n for n in zf.namelist() if n.startswith("pkg-") and n.endswith(".tar.zst"))
    dest.write_bytes(zf.read(name))
PY
  zstd -q -d "${work}/pkg.tar.zst" -o "${work}/pkg.tar"
  mkdir -p "${work}/pkg"
  tar -C "${work}/pkg" -xf "${work}/pkg.tar"
  mkdir -p "${PREFIX}"
  cp "${work}/pkg/share/snpeff-5.4.0c-0/snpEff.jar" "${work}/pkg/share/snpeff-5.4.0c-0/snpEff.config" "${PREFIX}/"
  mkdir -p "${PREFIX}/data"
  unzip -q -o "${db_zip}" -d "${PREFIX}"
  tar -C "${PREFIX}" -czf "${LOCAL_ARCHIVE}" snpEff.jar snpEff.config data
  rm -rf "${work}"
  if command -v aws >/dev/null 2>&1; then
    echo "[snpeff] uploading ${S3_URI}"
    aws s3 cp --only-show-errors "${LOCAL_ARCHIVE}" "${S3_URI}"
  else
    echo "[snpeff] aws CLI missing; left the bundle at ${LOCAL_ARCHIVE} and did not upload" >&2
  fi
}

require_java21
mkdir -p "${PREFIX}"

if [[ -f "${PREFIX}/snpEff.jar" && -d "${PREFIX}/data/GRCh38.115" ]]; then
  echo "[snpeff] already installed at ${PREFIX}"
elif [[ -s "${LOCAL_ARCHIVE}" ]]; then
  unpack_archive "${LOCAL_ARCHIVE}"
elif s3_object_exists; then
  echo "[snpeff] found ${S3_URI}"
  aws s3 cp --only-show-errors "${S3_URI}" "${LOCAL_ARCHIVE}"
  unpack_archive "${LOCAL_ARCHIVE}"
else
  build_upstream_bundle
fi

if [[ ! -f "${PREFIX}/snpEff.jar" || ! -d "${PREFIX}/data/GRCh38.115" ]]; then
  echo "ERROR: snpEff install incomplete under ${PREFIX}" >&2
  exit 1
fi
write_launcher
echo "[snpeff] installed ${PREFIX}/bin/snpEff"
