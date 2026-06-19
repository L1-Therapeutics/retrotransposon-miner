#!/usr/bin/env bash
set -euo pipefail

# Install system packages needed for IGV batch snapshots on headless Linux.
# Conda does not ship a supported end-user Xvfb package; use the OS package manager.

if [[ "$(uname -s)" != "Linux" ]]; then
  exit 0
fi

if command -v xvfb-run >/dev/null 2>&1 || command -v Xvfb >/dev/null 2>&1; then
  echo "[headless-igv] Xvfb already available"
  exit 0
fi

install_with_dnf() {
  sudo dnf install -y xorg-x11-server-Xvfb xorg-x11-xauth
}

install_with_apt() {
  sudo apt-get update
  sudo apt-get install -y xvfb
}

echo "[headless-igv] Xvfb not found; installing system packages for headless IGV..."

if command -v dnf >/dev/null 2>&1; then
  install_with_dnf
elif command -v yum >/dev/null 2>&1; then
  sudo yum install -y xorg-x11-server-Xvfb xorg-x11-xauth
elif command -v apt-get >/dev/null 2>&1; then
  install_with_apt
else
  echo "WARN: no supported package manager found (dnf/yum/apt-get)." >&2
  echo "Install Xvfb manually, e.g.: sudo dnf install -y xorg-x11-server-Xvfb xorg-x11-xauth" >&2
  exit 0
fi

if command -v xvfb-run >/dev/null 2>&1 || command -v Xvfb >/dev/null 2>&1; then
  echo "[headless-igv] Xvfb installed successfully"
else
  echo "WARN: package install finished but Xvfb/xvfb-run still not on PATH." >&2
fi
