"""Focused tests for scripts/ec2_jlab.sh stop-instance / reboot-instance reliability.

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
GREP_UNAUTHORIZED_ESCAPED = "UnauthorizedOperation"
MULTILINE_STOP_ERR = (
    "An error occurred (IncorrectInstanceState) when calling the StopInstances operation.\n"
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

    if args[:2] == ["ec2", "stop-instances"]:
        rc = int(os.environ.get("FAKE_AWS_STOP_RC", "0"))
        stop_err = os.environ.get("FAKE_AWS_STOP_ERR", "")
        if stop_err:
            print(stop_err, file=sys.stderr)
        elif rc:
            print(
                "An error occurred (UnauthorizedOperation) when calling the "
                "StopInstances operation: You are not authorized",
                file=sys.stderr,
            )
        return rc

    if args[:2] == ["ec2", "reboot-instances"]:
        rc = int(os.environ.get("FAKE_AWS_REBOOT_RC", "0"))
        reboot_err = os.environ.get("FAKE_AWS_REBOOT_ERR", "")
        if reboot_err:
            print(reboot_err, file=sys.stderr)
        elif rc:
            print(
                "An error occurred (UnauthorizedOperation) when calling the "
                "RebootInstances operation: You are not authorized",
                file=sys.stderr,
            )
        return rc

    if args[:2] == ["ec2", "wait"]:
        sub = args[2]
        key = {
            "instance-stopped": "FAKE_AWS_WAIT_RC_STOPPED",
            "instance-status-ok": "FAKE_AWS_WAIT_RC_STATUS",
            "instance-running": "FAKE_AWS_WAIT_RC_RUNNING",
        }.get(sub)
        rc = int(os.environ.get(key or "FAKE_AWS_WAIT_RC_STATUS", "0"))
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


def test_stop_succeeds_when_instance_running(runner: ScriptRunner) -> None:
    runner.run("stop-instance")

    assert runner.returncode == 0
    assert "Stopping " + IID in runner.stderr
    assert "Instance stopped." in runner.stderr
    # A: healthy path runs the exact expected AWS sequence, in order.
    assert runner.calls == [
        "configure get region",
        f"--region us-east-1 ec2 describe-instances --instance-ids {IID}",
        f"--region us-east-1 ec2 stop-instances --instance-ids {IID}",
        f"--region us-east-1 ec2 wait instance-stopped --instance-ids {IID}",
    ]
    runner.assert_no_temp_leftovers()


def test_stop_missing_instance_fails_fast(runner: ScriptRunner) -> None:
    runner.run("stop-instance", FAKE_AWS_MISSING="1")

    assert runner.returncode != 0
    assert "Selected instance no longer exists: " + IID in runner.stderr
    assert "list-instances" in runner.stderr
    assert "use <instance-id-or-name>" in runner.stderr
    assert "Instance stopped." not in runner.stderr
    assert not runner.called_with("ec2", "stop-instances"), \
        "no stop-instances should run when the preflight query fails"
    assert not runner.called_with("ec2", "wait")
    runner.assert_no_temp_leftovers()


def test_stop_api_failure_preserves_error_and_sudo_hint(runner: ScriptRunner) -> None:
    runner.run("stop-instance", FAKE_AWS_STOP_RC="255")

    assert runner.returncode != 0
    assert "Failed to stop instance " + IID in runner.stderr
    assert GREP_UNAUTHORIZED_ESCAPED in runner.stderr
    assert "sudo shutdown -h now" in runner.stderr
    assert not runner.called_with("ec2", "wait")
    runner.assert_no_temp_leftovers()


def test_stop_multiline_error_is_preserved(runner: ScriptRunner) -> None:
    runner.run("stop-instance", FAKE_AWS_STOP_RC="255", FAKE_AWS_STOP_ERR=MULTILINE_STOP_ERR)

    assert runner.returncode != 0
    assert "  An error occurred (IncorrectInstanceState)" in runner.stderr
    assert "  Instance i-0b6effaaa91380898 is in state \"stopping\"" in runner.stderr
    runner.assert_no_temp_leftovers()


def test_stop_waiter_timeout_with_instance_present(runner: ScriptRunner) -> None:
    runner.run("stop-instance", FAKE_AWS_WAIT_RC_STOPPED="255")

    assert runner.returncode != 0
    assert "did not reach the stopped state" in runner.stderr
    assert "Selected instance no longer exists" not in runner.stderr
    assert "sudo shutdown -h now" not in runner.stderr
    runner.assert_no_temp_leftovers()


def test_stop_waiter_timeout_with_instance_vanished(runner: ScriptRunner) -> None:
    runner.run("stop-instance", FAKE_AWS_WAIT_RC_STOPPED="255", FAKE_AWS_MISSING_AFTER_WAIT="1")

    assert runner.returncode != 0
    assert "Selected instance no longer exists: " + IID in runner.stderr
    assert "did not reach the stopped state" not in runner.stderr
    runner.assert_no_temp_leftovers()


def test_stop_describe_auth_failure_is_not_labeled_missing(runner: ScriptRunner) -> None:
    runner.run(
        "stop-instance",
        FAKE_AWS_DESCRIBE_ERR=(
            "An error occurred (UnauthorizedOperation) when calling the "
            "DescribeInstances operation: You are not authorized"
        ),
    )

    assert runner.returncode != 0
    assert "Unable to query instance " + IID in runner.stderr
    assert "UnauthorizedOperation" in runner.stderr
    assert "Selected instance no longer exists" not in runner.stderr
    assert not runner.called_with("ec2", "stop-instances")
    runner.assert_no_temp_leftovers()


def test_reboot_succeeds_when_instance_healthy(runner: ScriptRunner) -> None:
    runner.run("reboot-instance")

    assert runner.returncode == 0
    assert "Rebooting " + IID in runner.stderr
    assert "Instance healthy after reboot." in runner.stderr
    # A: healthy path runs the exact expected AWS sequence, in order.
    assert runner.calls == [
        "configure get region",
        f"--region us-east-1 ec2 describe-instances --instance-ids {IID}",
        f"--region us-east-1 ec2 reboot-instances --instance-ids {IID}",
        f"--region us-east-1 ec2 wait instance-status-ok --instance-ids {IID}",
    ]
    runner.assert_no_temp_leftovers()


def test_reboot_missing_instance_fails_fast(runner: ScriptRunner) -> None:
    runner.run("reboot-instance", FAKE_AWS_MISSING="1")

    assert runner.returncode != 0
    assert "Selected instance no longer exists: " + IID in runner.stderr
    assert "Instance healthy after reboot." not in runner.stderr
    assert not runner.called_with("ec2", "reboot-instances")
    assert not runner.called_with("ec2", "wait")
    runner.assert_no_temp_leftovers()


def test_reboot_api_failure_preserves_error(runner: ScriptRunner) -> None:
    runner.run("reboot-instance", FAKE_AWS_REBOOT_RC="255")

    assert runner.returncode != 0
    assert "Failed to reboot instance " + IID in runner.stderr
    assert GREP_UNAUTHORIZED_ESCAPED in runner.stderr
    assert not runner.called_with("ec2", "wait")
    runner.assert_no_temp_leftovers()


def test_reboot_waiter_timeout_with_instance_present(runner: ScriptRunner) -> None:
    runner.run("reboot-instance", FAKE_AWS_WAIT_RC_STATUS="255")

    assert runner.returncode != 0
    assert "did not become status-ok after reboot" in runner.stderr
    assert "Selected instance no longer exists" not in runner.stderr
    assert "Instance healthy after reboot." not in runner.stderr
    runner.assert_no_temp_leftovers()


def test_reboot_waiter_timeout_with_instance_vanished(runner: ScriptRunner) -> None:
    runner.run("reboot-instance", FAKE_AWS_WAIT_RC_STATUS="255", FAKE_AWS_MISSING_AFTER_WAIT="1")

    assert runner.returncode != 0
    assert "Selected instance no longer exists: " + IID in runner.stderr
    assert "Instance healthy after reboot." not in runner.stderr
    runner.assert_no_temp_leftovers()


def test_aws_cli_missing_aborts_with_clear_message(runner: ScriptRunner) -> None:
    minibin = runner.bindir.parent / "minibin"
    minibin.mkdir()
    for cmd in ("bash", "env", "mkdir", "dirname", "basename", "date", "sed"):
        src = shutil.which(cmd)
        if src is None:
            pytest.skip(f"{cmd} not available for the minimal PATH")
        os.symlink(src, minibin / cmd)

    for subcommand in ("stop-instance", "reboot-instance"):
        runner.run(subcommand, PATH=str(minibin))

        assert runner.returncode != 0
        assert "Missing command: aws" in runner.stdout
        assert "Selected instance no longer exists" not in runner.stderr


def test_help_output_is_unaffected(runner: ScriptRunner) -> None:
    runner.run("help")
    assert runner.returncode == 0
    assert "Usage:" in runner.stdout
    assert "stop-instance" in runner.stdout
    assert "reboot-instance" in runner.stdout
    assert "Instance stopped." not in runner.stderr
