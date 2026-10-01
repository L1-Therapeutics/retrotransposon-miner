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

### Type checking

`mypy` runs over the whole package using the `[tool.mypy]` config in `pyproject.toml`:

```bash
mypy
```

The policy is:

- **Every module under `src/retro_miner` is type-checked by default.** New modules must be clean; do not add entries to the grandfather list.
- Four legacy pipeline modules (`mei_support`, `local_assembly`, `candidate_loci`, `read_architecture`) are temporarily exempt via `ignore_errors` while their dynamic-pandas usage is cleaned up incrementally.
- The literature-anchored marker modules (`mei_markers`, `transduction_tags`) are held to a strict flag set (`disallow_untyped_defs`, `disallow_incomplete_defs`, `warn_return_any`). The flags are spelled out individually rather than via `strict = true`, which leaks across modules on some mypy builds.
- Type stubs (`pandas-stubs`, `boto3-stubs`) are listed in `environment.yml`; install them in your env or mypy will report missing-import noise.
- When a legacy module is fully clean, remove it from the `ignore_errors` override in the same PR. When a new module earns strict status, add it to the strict override.

## Testing and Validation

At minimum, run:

```bash
bash scripts/validate_environment.sh
```

If your change affects calling behavior, include a small reproducible run and output summary in the PR.

## License

By submitting a contribution, you agree that your contributions are licensed under the Apache License 2.0 in this repository.
