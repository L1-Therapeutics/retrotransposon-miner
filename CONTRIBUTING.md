# Contributing to retrotransposon-miner

Thanks for your interest in contributing.

## Ways to Contribute

- Report bugs with clear reproduction steps and expected behavior.
- Propose features with concrete use cases.
- Improve documentation, examples, and reproducibility notes.
- Submit code fixes and enhancements via pull requests.
- For questions/support, open an issue (or discussion if enabled) or email `william [at] l1tx [dot] com`.

## Before You Start

1. Open an issue for substantial changes so scope can be discussed early.
2. Keep changes focused and reviewable.
3. Prefer tests or validation notes alongside behavioral changes.

## Development Setup

From the repository root:

```bash
bash scripts/bootstrap_env.sh
bash scripts/install_ucsc_tools.sh
conda activate rtm-miner || micromamba activate rtm-miner
pip install -e ".[dev]"
bash scripts/validate_environment.sh
```

## Branching and Pull Requests

- Branch from `main`.
- Use descriptive branch names (for example: `fix/igv-timeout-handling`).
- Keep pull requests small and self-contained where possible.
- In the PR description include:
  - problem statement,
  - summary of approach,
  - test/validation steps,
  - any limitations or follow-up work.

## Coding Guidelines

- Preserve existing CLI contracts and output column names unless explicitly changing API behavior.
- Favor explicit, reproducible pipeline behavior over implicit defaults.
- Add concise comments only where logic is non-obvious.
- Update README/docs when user-facing behavior changes.

## Testing and Validation

At minimum, run:

```bash
bash scripts/validate_environment.sh
pip install -e ".[dev]"
ruff check src tests
python -m pytest -q
```

If your change affects calling behavior, include a small reproducible run and output summary in the PR.

### Focused Test Runs

For specific test files (for example, PR #32 scope):

```bash
python -m pytest -q tests/test_extract_one_pass.py tests/test_mate_fetch_sweep.py
```

For the EC2 helper script (`scripts/ec2_jlab.sh`) or environment validation changes,
the shell-facing tests are self-contained (they run the real scripts against a fake
`aws`/tool PATH, never real credentials or network):

```bash
python -m pytest -q tests/test_ec2_jlab_ami.py \
  tests/test_ec2_jlab_bind_resolve_ip.py \
  tests/test_ec2_jlab_ensure_eip.py \
  tests/test_ec2_jlab_start.py \
  tests/test_ec2_jlab_status.py \
  tests/test_ec2_jlab_stop_reboot.py \
  tests/test_validate_environment.py
```

CI (`.github/workflows/ci.yml`) runs the full `pytest` suite plus `ruff check` and
`bash -n` syntax checks on every push/PR, so local failures map 1:1 to CI.

## License

By submitting a contribution, you agree that your contributions are licensed under the Apache License 2.0 in this repository.
