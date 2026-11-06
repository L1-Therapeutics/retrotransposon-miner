#!/usr/bin/env bash
# Stop this instance after align_mouse_sra.py has uploaded one BAM.
# A missing "finished" line leaves the instance up.
set -eu
ACC="${1:?accession}"
LOG="${HOME}/logs/mouse-align.log"
PIDFILE="${HOME}/logs/mouse-align.pid"
NOTE="${HOME}/logs/mouse-align-shutdown.log"

while true; do
  pid=""
  if [[ -f "${PIDFILE}" ]]; then
    pid="$(tr -d '[:space:]' < "${PIDFILE}" || true)"
  fi
  if [[ -z "${pid}" ]] || ! kill -0 "${pid}" 2>/dev/null; then
    break
  fi
  sleep 30
done

if [[ -f "${LOG}" ]] && grep -q "finished ${ACC} ->" "${LOG}"; then
  echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) ${ACC} uploaded; stopping instance" >> "${NOTE}"
  sudo shutdown -h now
else
  echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) ${ACC} did not finish an upload; leaving instance running" >> "${NOTE}"
fi
