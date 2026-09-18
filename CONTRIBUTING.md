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

## Testing and Validation

At minimum, run:

```bash
bash scripts/validate_environment.sh
```

If your change affects calling behavior, include a small reproducible run and output summary in the PR.

## Working on `scripts/ec2_jlab.sh`

`scripts/ec2_jlab.sh` provisions and manages the EC2 instance used for larger
runs. Because its failure modes cost money (a wrong AMI, an orphaned Elastic
IP, an instance left running), its functions are unit-tested rather than
exercised against live AWS.

### Test harness

Tests drive individual shell functions through an internal hook instead of
running the CLI end to end:

```bash
EC2_JLAB_TEST_FN=resolve_ami bash scripts/ec2_jlab.sh
```

The hook dispatches to one named function, then exits with that function's
status. If you add a function that needs direct coverage, extend the `case`
block guarded by `EC2_JLAB_TEST_FN` near the bottom of the script.

Each test places a **fake `aws` executable first on `PATH`** and asserts that
`PATH` resolution picks the fake rather than a real binary. AWS calls are
matched by argument string and answered with canned stdout/exit codes,
including real-world error text (`AccessDeniedException`, `ParameterNotFound`,
endpoint-connection failures) so error paths are covered rather than assumed.

Run them with:

```bash
python -m pytest tests/test_ec2_jlab_ami.py \
  tests/test_ec2_jlab_bind_resolve_ip.py tests/test_ec2_jlab_ensure_eip.py \
  tests/test_ec2_jlab_start.py tests/test_ec2_jlab_status.py \
  tests/test_ec2_jlab_stop_reboot.py
```

### Rules for these tests

- **No real AWS.** No credentials, no instance metadata, no network. Tests
  scrub AWS environment variables before launching any subprocess, so a
  contributor with live credentials cannot accidentally bill the account.
- **Never reach `run-instances`** (or any other state-changing call) from a
  test. Assert on the arguments that *would* have been sent instead.
- **No `shell=True`.** Subprocesses are invoked with an argument list.
- Prefer asserting the **sequence of AWS calls** made, not just the return
  value: the ordering is the contract (for example, an explicit `AMI_ID` must
  short-circuit before any SSM or `describe-images` lookup).

### AMI resolution contract

`resolve_ami` tries three sources in order, and **validates the result of every
one** before returning it:

1. `AMI_ID` environment override — checked with `describe-images`, never
   trusted blindly, and never falls through to a lookup.
2. The public AL2023 SSM parameter — the returned value is syntax-checked
   rather than assumed well-formed.
3. `ec2 describe-images` fallback — filtered to Amazon-owned, AL2023
   kernel 6.1, `x86_64`, EBS-backed, `available`, newest by `CreationDate`.

Validation rejects an AMI that does not exist, is not `available`, is not
`x86_64`, or is not EBS-backed. Only the resolved AMI id goes to stdout;
diagnostics go to stderr, so callers can safely capture the id by command
substitution.

## License

By submitting a contribution, you agree that your contributions are licensed under the Apache License 2.0 in this repository.
