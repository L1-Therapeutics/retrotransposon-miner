"""Focused tests for scripts/ec2_jlab.sh start-instance reliability.

The real shell script is executed against a fake `aws` executable placed
first on PATH. No real AWS credentials, metadata, or network are used, and
no subprocess is ever launched with shell=True.
"""

from __future__ import annotations

import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "ec2_jlab.sh"
IID = "i-0b6effaaa91380898"
MULTILINE_START_ERR = (
    "An error occurred (IncorrectInstanceState) when calling the StartInstances operation.\n"
    'Instance i-0b6effaaa91380898 is in state "stopping"; retry after the stop completes.'
)

FAKE_AWS = """\
#!/usr/bin/env python3
import os
import sys


def _log(args):
    path = os.environ.get("FAKE_AWS_LOG", "")
    if path:
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(" ".join(args) + "\\n")


def _main():
    args = sys.argv[1:]
    _log(args)

    # aws configure get region  (script top-level)
    if args[:3] == ["configure", "get", "region"]:
        print(os.environ.get("FAKE_AWS_REGION", "us-east-1"))
        return 0

    # strip --region <region>
    if args and args[0] == "--region":
        args = args[2:]

    if args[:2] == ["ec2", "describe-instances"]:
        describe_err = os.environ.get("FAKE_AWS_DESCRIBE_ERR", "")
        if describe_err:
            print(describe_err, file=sys.stderr)
            return 255
        flipped = (
            os.environ.get("FAKE_AWS_MISSING_AFTER_WAIT", "0") == "1"
            and os.path.exists(os.environ.get("FAKE_AWS_STATE", ""))
        )
        missing = os.environ.get("FAKE_AWS_MISSING", "0") == "1" or flipped
        if missing:
            print(
                "An error occurred (InvalidInstanceID.NotFound) when calling "
                "the DescribeInstances operation: The instance ID does not exist",
                file=sys.stderr,
            )
            return 255
        return 0

    if args[:2] == ["ec2", "start-instances"]:
        rc = int(os.environ.get("FAKE_AWS_START_RC", "0"))
        start_err = os.environ.get("FAKE_AWS_START_ERR", "")
        if rc and start_err:
            print(start_err, file=sys.stderr)
        elif rc:
            print(
                "An error occurred (UnauthorizedOperation) when calling the "
                "StartInstances operation: You are not authorized",
                file=sys.stderr,
            )
        return rc

    if args[:2] == ["ec2", "wait"]:
        sub = args[2]
        key = (
            "FAKE_AWS_WAIT_RC_RUNNING"
            if sub == "instance-running"
            else "FAKE_AWS_WAIT_RC_STATUS"
        )
        rc = int(os.environ.get(key, "0"))
        if rc:
            print(f"Waiter {sub} failed: Max attempts exceeded", file=sys.stderr)
            with open(os.environ.get("FAKE_AWS_STATE", ""), "w", encoding="utf-8") as fh:
                fh.write("flipped\\n")
        return rc

    return 0


if __name__ == "__main__":
    sys.exit(_main())
"""


class ScriptRunner:
    def __init__(self, root: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        bindir = root / "bin"
        bindir.mkdir()
        aws = bindir / "aws"
        aws.write_text(FAKE_AWS, encoding="utf-8")
        aws.chmod(aws.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        self.bindir = bindir
        self.home = root / "home"
        self.home.mkdir()
        self.tmpdir = root / "tmpdir"
        self.tmpdir.mkdir()
        self.state_file = root / "instance-state.env"
        self.log_file = root / "aws-calls.log"
        self.state_marker = root / "fake-aws.state"

        monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
        monkeypatch.setenv("HOME", str(self.home))
        monkeypatch.setenv("INSTANCE_STATE_FILE", str(self.state_file))
        monkeypatch.setenv("INSTANCE_ID", IID)
        monkeypatch.delenv("INSTANCE_NAME", raising=False)
        monkeypatch.setenv("FAKE_AWS_LOG", str(self.log_file))
        monkeypatch.setenv("FAKE_AWS_STATE", str(self.state_marker))
        monkeypatch.setenv("FAKE_AWS_REGION", "us-east-1")

        self.returncode: int | None = None
        self.stdout = ""
        self.stderr = ""
        self._calls: list[str] = []

    def run(self, *args: str, **env_overrides: str) -> None:
        self._calls = []
        self.log_file.write_text("", encoding="utf-8")
        self.state_marker.unlink(missing_ok=True)

        env = dict(os.environ)
        env.update(env_overrides)
        # Defense in depth: never hand real AWS credentials to any subprocess.
        for var in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN", "AWS_PROFILE"):
            env.pop(var, None)
        env["AWS_EC2_METADATA_DISABLED"] = "true"
        env["TMPDIR"] = str(self.tmpdir)

        if "PATH" not in env_overrides:
            assert shutil.which("aws", path=env["PATH"]) == str(self.bindir / "aws"), (
                "the fake aws must be what PATH resolves, never the real binary"
            )

        proc = subprocess.run(
            [str(SCRIPT), *args],
            capture_output=True,
            text=True,
            env=env,
            cwd=str(REPO_ROOT),
        )
        self.returncode = proc.returncode
        self.stdout = proc.stdout
        self.stderr = proc.stderr
        self._calls = [
            line for line in self.log_file.read_text(encoding="utf-8").splitlines() if line
        ]

    @property
    def calls(self) -> list[str]:
        return list(self._calls)

    def called_with(self, *tokens: str) -> bool:
        for call in self._calls:
            pos = 0
            for token in tokens:
                idx = call.find(token, pos)
                if idx < 0:
                    break
                pos = idx + len(token)
            else:
                return True
        return False

    def assert_no_temp_leftovers(self) -> None:
        assert list(self.tmpdir.iterdir()) == [], "temp files leaked into TMPDIR"


@pytest.fixture()
def runner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> ScriptRunner:
    return ScriptRunner(tmp_path, monkeypatch)


def test_start_succeeds_when_instance_is_healthy(runner: ScriptRunner) -> None:
    runner.run("start-instance")

    assert runner.returncode == 0
    assert "Starting instance " + IID in runner.stderr
    assert "Instance is running and healthy." in runner.stderr
    # A/B: healthy path runs the exact expected AWS sequence, in order.
    assert runner.calls == [
        "configure get region",
        f"--region us-east-1 ec2 describe-instances --instance-ids {IID}",
        f"--region us-east-1 ec2 start-instances --instance-ids {IID}",
        f"--region us-east-1 ec2 wait instance-running --instance-ids {IID}",
        f"--region us-east-1 ec2 wait instance-status-ok --instance-ids {IID}",
    ]
    runner.assert_no_temp_leftovers()


def test_start_failure_exits_nonzero_with_diagnostic(runner: ScriptRunner) -> None:
    # B/E: start-instances fails (default UnauthorizedOperation); error is preserved.
    runner.run("start-instance", FAKE_AWS_START_RC="255")

    assert runner.returncode != 0
    assert "Failed to start instance " + IID in runner.stderr
    assert "UnauthorizedOperation" in runner.stderr
    assert "Instance is running and healthy." not in runner.stderr
    runner.assert_no_temp_leftovers()


def test_multiline_start_error_is_preserved(runner: ScriptRunner) -> None:
    # N: AWS stderr with newlines and quotes survives word-for-word, indented.
    runner.run("start-instance", FAKE_AWS_START_RC="255", FAKE_AWS_START_ERR=MULTILINE_START_ERR)

    assert runner.returncode != 0
    assert "Failed to start instance " + IID in runner.stderr
    assert "  An error occurred (IncorrectInstanceState)" in runner.stderr
    assert '  Instance i-0b6effaaa91380898 is in state "stopping"' in runner.stderr
    runner.assert_no_temp_leftovers()


def test_missing_instance_before_start_reports_and_suggests_workflow(runner: ScriptRunner) -> None:
    # C: instance gone before start; NotFound must yield the missing diagnostic.
    runner.run("start-instance", FAKE_AWS_MISSING="1")

    assert runner.returncode != 0
    assert "Selected instance no longer exists: " + IID in runner.stderr
    assert "list-instances" in runner.stderr
    assert "use <instance-id-or-name>" in runner.stderr
    assert "bootstrap" in runner.stderr
    assert "Instance is running and healthy." not in runner.stderr
    assert runner.calls == ["configure get region", f"--region us-east-1 ec2 describe-instances --instance-ids {IID}"]
    runner.assert_no_temp_leftovers()


def test_describe_auth_failure_is_not_labeled_missing(runner: ScriptRunner) -> None:
    # Objective 6: an auth/network describe failure must not read as "no longer exists".
    runner.run(
        "start-instance",
        FAKE_AWS_DESCRIBE_ERR=(
            "An error occurred (UnauthorizedOperation) when calling the "
            "DescribeInstances operation: You are not authorized"
        ),
    )

    assert runner.returncode != 0
    assert "Unable to query instance " + IID in runner.stderr
    assert "UnauthorizedOperation" in runner.stderr
    assert "Selected instance no longer exists" not in runner.stderr
    assert runner.called_with("ec2", "describe-instances", "--instance-ids", IID), \
        "preflight describe must have run exactly once"
    assert not runner.called_with("ec2", "start-instances"), \
        "no start-instances should run when the preflight query fails"


def test_instance_disappearing_during_running_wait_is_reported(runner: ScriptRunner) -> None:
    # H: start succeeds, running waiter fails, describe then says the instance is gone.
    runner.run("start-instance", FAKE_AWS_WAIT_RC_RUNNING="255", FAKE_AWS_MISSING_AFTER_WAIT="1")

    assert runner.returncode != 0
    assert "Selected instance no longer exists: " + IID in runner.stderr
    assert "list-instances" in runner.stderr
    assert runner.called_with("ec2", "wait", "instance-running", "--instance-ids", IID)
    assert not runner.called_with("ec2", "wait", "instance-status-ok", "--instance-ids", IID)
    runner.assert_no_temp_leftovers()


def test_running_waiter_timeout_with_instance_present_fails_clearly(runner: ScriptRunner) -> None:
    # G: waiter times out but instance exists: distinct, actionable message.
    runner.run("start-instance", FAKE_AWS_WAIT_RC_RUNNING="255")

    assert runner.returncode != 0
    assert "did not reach the running state" in runner.stderr
    assert "Selected instance no longer exists" not in runner.stderr
    runner.assert_no_temp_leftovers()


def test_status_ok_waiter_timeout_with_instance_present_fails_clearly(runner: ScriptRunner) -> None:
    # I: status-ok waiter times out but the instance still exists.
    runner.run("start-instance", FAKE_AWS_WAIT_RC_STATUS="255")

    assert runner.returncode != 0
    assert "is running but not status-ok yet" in runner.stderr
    assert "did not reach the running state" not in runner.stderr
    assert runner.called_with("ec2", "wait", "instance-running", "--instance-ids", IID)
    assert runner.called_with("ec2", "wait", "instance-status-ok", "--instance-ids", IID)
    runner.assert_no_temp_leftovers()


def test_instance_disappearing_during_status_ok_wait_is_reported(runner: ScriptRunner) -> None:
    # J: running waiter passes, status-ok waiter fails, then describe says it is gone.
    runner.run("start-instance", FAKE_AWS_WAIT_RC_STATUS="255", FAKE_AWS_MISSING_AFTER_WAIT="1")

    assert runner.returncode != 0
    assert "Selected instance no longer exists: " + IID in runner.stderr
    assert "list-instances" in runner.stderr


def test_aws_cli_missing_aborts_with_clear_message(runner: ScriptRunner) -> None:
    # K: an absent aws CLI must fail clearly and never touch a real binary.
    minibin = runner.bindir.parent / "minibin"
    minibin.mkdir()
    for cmd in ("bash", "env", "mkdir", "dirname", "basename", "date", "sed"):
        src = shutil.which(cmd)
        if src is None:
            pytest.skip(f"{cmd} not available for the minimal PATH")
        os.symlink(src, minibin / cmd)

    runner.run("start-instance", PATH=str(minibin))

    assert runner.returncode != 0
    assert "Missing command: aws" in runner.stdout
    assert "Selected instance no longer exists" not in runner.stderr


def test_no_real_aws_command_is_invoked(runner: ScriptRunner) -> None:
    # D: every aws-shaped call in a failing start went through the fake, in order.
    runner.run("start-instance", FAKE_AWS_START_RC="255")

    expected = [
        "configure get region",
        f"--region us-east-1 ec2 describe-instances --instance-ids {IID}",
        f"--region us-east-1 ec2 start-instances --instance-ids {IID}",
    ]
    assert runner.calls == expected
    assert not runner.called_with("ec2", "wait")


def test_help_output_is_unaffected(runner: ScriptRunner) -> None:
    # E: help/usage output is unchanged and requires no AWS calls beyond the fake.
    for flag in ("help", "--help"):
        runner.run(flag)
        assert runner.returncode == 0
        assert "Usage:" in runner.stdout
        assert "list-instances" in runner.stdout
        assert "start-instance" in runner.stdout
        assert "Instance is running and healthy." not in runner.stderr