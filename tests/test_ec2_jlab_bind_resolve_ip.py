"""Focused tests for scripts/ec2_jlab.sh bind-instance and resolve-ip reliability.

The real shell script is executed against a fake `aws` executable placed
first on PATH. A fake `curl` (which exits 1) is also placed first on PATH so
the IMDS fallback inside get_instance_id can never touch a real network. No
real AWS credentials, metadata, or network are used, and no subprocess is
ever launched with shell=True.
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
PROBE_CALL = f"--region us-east-1 ec2 describe-instances --instance-ids {IID}"
NAME_QUERY = "Reservations[0].Instances[0].Tags[?Key==`Name`]|[0].Value"
IMAGE_QUERY = "Reservations[0].Instances[0].ImageId"
KEY_QUERY = "Reservations[0].Instances[0].KeyName"
IP_QUERY = "Reservations[0].Instances[0].PublicIpAddress"
FAKE_IP = "203.0.113.10"

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
        if os.environ.get("FAKE_AWS_MISSING", "0") == "1":
            print(
                "An error occurred (InvalidInstanceID.NotFound) when calling "
                "the DescribeInstances operation: The instance ID does not exist",
                file=sys.stderr,
            )
            return 255
        query = ""
        if "--query" in args:
            query = args[args.index("--query") + 1]
        if "PublicIpAddress" in query:
            print(os.environ.get("FAKE_IP", "203.0.113.10"))
        elif "KeyName" in query:
            print("fake-key")
        elif "ImageId" in query:
            print("ami-0f111111111111111")
        elif "Tags[?Key==`Name`]" in query:
            print("fake-name")
        return 0

    if args[:2] == ["ec2", "describe-images"]:
        print("al2023-ami-2023.1.1.0-kernel-6.1-x86_64")
        return 0

    return 0


if __name__ == "__main__":
    sys.exit(_main())
"""

FAKE_CURL = """\
#!/usr/bin/env python3
import sys
sys.exit(1)
"""


class BindRunner:
    def __init__(self, root: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        bindir = root / "bin"
        bindir.mkdir()
        for name, body in (("aws", FAKE_AWS), ("curl", FAKE_CURL)):
            exe = bindir / name
            exe.write_text(body, encoding="utf-8")
            exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        self.bindir = bindir
        self.home = root / "home"
        self.home.mkdir()
        self.tmpdir = root / "tmpdir"
        self.tmpdir.mkdir()
        self.state_file = root / "instance-state.env"
        self.log_file = root / "aws-calls.log"

        monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
        monkeypatch.setenv("HOME", str(self.home))
        monkeypatch.setenv("INSTANCE_STATE_FILE", str(self.state_file))
        monkeypatch.delenv("INSTANCE_ID", raising=False)
        monkeypatch.delenv("INSTANCE_NAME", raising=False)
        monkeypatch.delenv("AMI_ID", raising=False)
        monkeypatch.setenv("FAKE_AWS_LOG", str(self.log_file))
        monkeypatch.setenv("FAKE_AWS_REGION", "us-east-1")

        self.returncode: int | None = None
        self.stdout = ""
        self.stderr = ""
        self._calls: list[str] = []

    def run(self, *args: str, **env_overrides: str) -> None:
        self._calls = []
        self.log_file.write_text("", encoding="utf-8")

        env = dict(os.environ)
        env.update(env_overrides)
        # Defense in depth: never hand real AWS credentials to any subprocess.
        for var in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN", "AWS_PROFILE"):
            env.pop(var, None)
        env["AWS_EC2_METADATA_DISABLED"] = "true"
        env["TMPDIR"] = str(self.tmpdir)

        assert shutil.which("aws", path=env["PATH"]) == str(self.bindir / "aws"), (
            "the fake aws must be what PATH resolves, never the real binary"
        )
        assert shutil.which("curl", path=env["PATH"]) == str(self.bindir / "curl"), (
            "the fake curl must be what PATH resolves, never the real network tool"
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
def runner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> BindRunner:
    return BindRunner(tmp_path, monkeypatch)


def test_bind_healthy_binds_instance_and_writes_ssh_config(runner: BindRunner) -> None:
    runner.run("use", IID)

    assert runner.returncode == 0
    assert f"Bound {IID}. Start with:" in runner.stderr
    assert str(SCRIPT) in runner.stderr
    # Probe, then write_ssh_config's name/image/key lookups, in order.
    assert runner.calls[0] == "configure get region"
    assert runner.calls[1] == PROBE_CALL
    assert f"--query {NAME_QUERY} --output text" in runner.calls[2]
    assert f"--query {IMAGE_QUERY} --output text" in runner.calls[3]
    assert runner.called_with("ec2", "describe-images")
    assert f"--query {KEY_QUERY} --output text" in runner.calls[-1]
    assert (runner.home / ".ssh" / "config").exists(), "ssh config should be written"
    runner.assert_no_temp_leftovers()


def test_bind_terminated_instance_reported_not_found(runner: BindRunner) -> None:
    runner.run("use", IID, FAKE_AWS_MISSING="1")

    assert runner.returncode != 0
    assert f"Instance not found in us-east-1: {IID}" in runner.stderr
    assert "Selected instance no longer exists" not in runner.stderr
    # Probe fails first; no ssh-config writes or later lookups follow.
    assert runner.calls == ["configure get region", PROBE_CALL]
    assert not runner.called_with("ec2", "describe-images")
    assert not (runner.home / ".ssh" / "config").exists()
    runner.assert_no_temp_leftovers()


def test_bind_auth_failure_is_not_labeled_not_found(runner: BindRunner) -> None:
    runner.run(
        "use", IID,
        FAKE_AWS_DESCRIBE_ERR=(
            "An error occurred (UnauthorizedOperation) when calling the "
            "DescribeInstances operation: You are not authorized"
        ),
    )

    assert runner.returncode != 0
    assert "Unable to query instance " + IID in runner.stderr
    assert "UnauthorizedOperation" in runner.stderr
    assert "Instance not found" not in runner.stderr
    assert runner.calls == ["configure get region", PROBE_CALL]
    runner.assert_no_temp_leftovers()


def test_resolve_ip_returns_public_ip(runner: BindRunner) -> None:
    runner.run("resolve-ip", INSTANCE_ID=IID)

    assert runner.returncode == 0
    assert runner.stdout == FAKE_IP + "\n"
    assert runner.calls == ["configure get region", PROBE_CALL, f"--region us-east-1 ec2 describe-instances --instance-ids {IID} --query {IP_QUERY} --output text"]
    runner.assert_no_temp_leftovers()


def test_resolve_ip_missing_instance_gives_diagnostic(runner: BindRunner) -> None:
    runner.run("resolve-ip", INSTANCE_ID=IID, FAKE_AWS_MISSING="1")

    assert runner.returncode != 0
    assert "Selected instance no longer exists: " + IID in runner.stderr
    assert runner.stdout == ""
    assert runner.calls == ["configure get region", PROBE_CALL]
    assert not runner.called_with(IP_QUERY)
    runner.assert_no_temp_leftovers()


def test_resolve_ip_auth_failure_is_not_labeled_missing(runner: BindRunner) -> None:
    runner.run(
        "resolve-ip",
        INSTANCE_ID=IID,
        FAKE_AWS_DESCRIBE_ERR=(
            "An error occurred (UnauthorizedOperation) when calling the "
            "DescribeInstances operation: You are not authorized"
        ),
    )

    assert runner.returncode != 0
    assert "Unable to query instance " + IID in runner.stderr
    assert "Selected instance no longer exists" not in runner.stderr
    assert runner.stdout == ""
    runner.assert_no_temp_leftovers()
