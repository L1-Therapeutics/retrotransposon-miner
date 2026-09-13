#!/usr/bin/env bash
set -e

echo "[+] Running Ruff linter..."
python3 -m pip install -q ruff black mypy
ruff check src/ scripts/ benchmarks/ tests/

echo "[+] Running Black format check..."
black --check src/ scripts/ benchmarks/ tests/

echo "[+] Running Mypy type checker on core package..."
mypy src/retro_miner
