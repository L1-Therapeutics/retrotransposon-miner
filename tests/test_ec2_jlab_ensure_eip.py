"""Focused tests for scripts/ec2_jlab.sh ensure-eip reliability.

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
APP = "retrotransposon-miner"
EIP_TAG = f"{APP}-{IID}-eip"
ALLOC = "eipalloc-1234567890"
ASSOC = "eipassoc-0987654321"
PROBE_CALL = f"--region us-east-1 ec2 describe-instances --instance-ids {IID}"
TAG_DESCRIBE = (
    f"--region us-east-1 ec2 describe-addresses --filters Name=tag:Name,Values={EIP_TAG} "
    "--query Addresses[0].AllocationId --output text"
)
INSTANCE_DESCRIBE = (
    f"--region us-east-1 ec2 describe-addresses --filters Name=instance-id,Values={IID} "
    "--query Addresses[0].AllocationId --output text"
)
ALLOCATE = (
    "--region us-east-1 ec2 allocate-address --domain vpc --query AllocationId --output text"
)
CREATE_TAGS = (
    f"--region us-east-1 ec2 create-tags --resources {ALLOC} "
    f"--tags Key=Name,Value={EIP_TAG} Key=App,Value={APP}"
)
ASSOCIATE = (
    f"--region us-east-1 ec2 associate-address --instance-id {IID} "
    f"--allocation-id {ALLOC} --allow-reassociation --query AssociationId --output text"
)
RELEASE = f"--region us-east-1 ec2 release-address --allocation-id {ALLOC}"

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

    if args[:3] == ["configure", "get", "region"]:
        print(os.environ.get("FAKE_AWS_REGION", "us-east-1"))
        return 0

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
        return 0

    if args[:2] == ["ec2", "describe-addresses"]:
        mode = os.environ.get("FAKE_EIP_MODE", "none")
        this = "tag" if any("Name=tag:Name" in a for a in args) else "instance"
        if this == mode:
            print("eipalloc-1234567890")
        else:
            print("None")
        return 0

    if args[:2] == ["ec2", "allocate-address"]:
        rc = int(os.environ.get("FAKE_ALLOC_RC", "0"))
        if rc:
            print(
                "An error occurred (UnauthorizedOperation) when calling the "
                "AllocateAddress operation: You are not authorized",
                file=sys.stderr,
            )
            return rc
        print("eipalloc-1234567890")
        return 0

    if args[:2] == ["ec2", "associate-address"]:
        rc = int(os.environ.get("FAKE_ASSOC_RC", "0"))
        if rc:
            print(
                "An error occurred (InvalidAllocationID.NotFound) when calling "
                "the AssociateAddress operation: The allocation ID does not exist",
                file=sys.stderr,
            )
            return rc
        print("eipassoc-0987654321")
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


class EipRunner:
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
        monkeypatch.setenv("INSTANCE_ID", IID)
        monkeypatch.delenv("INSTANCE_NAME", raising=False)
        monkeypatch.setenv("FAKE_AWS_LOG", str(self.log_file))
        monkeypatch.setenv("FAKE_AWS_REGION", "us-east-1")

        self.returncode: int | None = None
        self.stdout = ""
        self.stderr = ""
        self._calls: list[str] = []

    def run(self, **env_overrides: str) -> None:
        self._calls = []
        self.log_file.write_text("", encoding="utf-8")

        env = dict(os.environ)
        env.update(env_overrides)
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
            [str(SCRIPT), "ensure-eip"],
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
def runner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> EipRunner:
    return EipRunner(tmp_path, monkeypatch)


def test_eip_healthy_allocates_and_associates(runner: EipRunner) -> None:
    runner.run()

    assert runner.returncode == 0
    assert f"EIP association: {ASSOC}" in runner.stderr
    assert runner.calls == [
        "configure get region",
        PROBE_CALL,
        TAG_DESCRIBE,
        INSTANCE_DESCRIBE,
        ALLOCATE,
        CREATE_TAGS,
        ASSOCIATE,
    ]
    assert not runner.called_with("release-address")
    runner.assert_no_temp_leftovers()


def test_eip_reuses_existing_tagged_address(runner: EipRunner) -> None:
    runner.run(FAKE_EIP_MODE="tag")

    assert runner.returncode == 0
    assert f"EIP association: {ASSOC}" in runner.stderr
    assert runner.calls == [
        "configure get region",
        PROBE_CALL,
        TAG_DESCRIBE,
        ASSOCIATE,
    ]
    assert not runner.called_with("allocate-address")
    assert not runner.called_with("create-tags")
    runner.assert_no_temp_leftovers()


def test_eip_reuses_existing_instance_address(runner: EipRunner) -> None:
    runner.run(FAKE_EIP_MODE="instance")

    assert runner.returncode == 0
    assert f"EIP association: {ASSOC}" in runner.stderr
    assert runner.calls == [
        "configure get region",
        PROBE_CALL,
        TAG_DESCRIBE,
        INSTANCE_DESCRIBE,
        ASSOCIATE,
    ]
    assert not runner.called_with("allocate-address")
    runner.assert_no_temp_leftovers()


def test_eip_missing_instance_fails_before_allocating(runner: EipRunner) -> None:
    runner.run(FAKE_AWS_MISSING="1")

    assert runner.returncode != 0
    assert "Selected instance no longer exists: " + IID in runner.stderr
    assert runner.calls == ["configure get region", PROBE_CALL]
    assert not runner.called_with("describe-addresses")
    assert not runner.called_with("allocate-address")
    assert not runner.called_with("associate-address")
    runner.assert_no_temp_leftovers()


def test_eip_auth_failure_is_not_labeled_missing(runner: EipRunner) -> None:
    runner.run(
        FAKE_AWS_DESCRIBE_ERR=(
            "An error occurred (UnauthorizedOperation) when calling the "
            "DescribeInstances operation: You are not authorized"
        ),
    )

    assert runner.returncode != 0
    assert "Unable to query instance " + IID in runner.stderr
    assert "Selected instance no longer exists" not in runner.stderr
    assert not runner.called_with("allocate-address")
    runner.assert_no_temp_leftovers()


def test_eip_fresh_allocation_failure(runner: EipRunner) -> None:
    runner.run(FAKE_ALLOC_RC="255")

    assert runner.returncode != 0
    assert "Failed to allocate an Elastic IP." in runner.stderr
    assert not runner.called_with("associate-address")
    assert not runner.called_with("release-address")
    runner.assert_no_temp_leftovers()


def test_eip_association_failure_releases_fresh_allocation(runner: EipRunner) -> None:
    runner.run(FAKE_ASSOC_RC="255")

    assert runner.returncode != 0
    assert f"Failed to associate EIP {ALLOC} to {IID}" in runner.stderr
    assert "Released freshly allocated EIP " + ALLOC in runner.stderr
    assert runner.called_with(RELEASE)
    assert runner.calls.index(RELEASE) > runner.calls.index(ALLOCATE)
    runner.assert_no_temp_leftovers()


def test_eip_association_failure_keeps_preexisting_address(runner: EipRunner) -> None:
    runner.run(FAKE_EIP_MODE="tag", FAKE_ASSOC_RC="255")

    assert runner.returncode != 0
    assert f"Failed to associate EIP {ALLOC} to {IID}" in runner.stderr
    assert "Released freshly allocated EIP" not in runner.stderr
    assert not runner.called_with("release-address")
    assert not runner.called_with("allocate-address")
    runner.assert_no_temp_leftovers()


def test_fake_aws_is_the_resolved_binary(runner: EipRunner) -> None:
    assert shutil.which("aws", path=os.environ["PATH"]) == str(runner.bindir / "aws")
    assert shutil.which("curl", path=os.environ["PATH"]) == str(runner.bindir / "curl")
