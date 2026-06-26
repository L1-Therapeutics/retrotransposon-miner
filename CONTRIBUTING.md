# Contributing to retrotransposon-miner

Thanks for your interest in contributing.

## Ways to Contribute

- Report bugs with clear reproduction steps and expected behavior.
- Propose features with concrete use cases.
- Improve documentation, examples, and reproducibility notes.
- Submit code fixes and enhancements via pull requests.
- For questions/support, open an issue (or discussion if enabled) or email `william@l1tx.com`.

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

## Testing and Validation

At minimum, run:

```bash
bash scripts/validate_environment.sh
```

If your change affects calling behavior, include a small reproducible run and output summary in the PR.

## License

By submitting a contribution, you agree that your contributions are licensed under the Apache License 2.0 in this repository.
