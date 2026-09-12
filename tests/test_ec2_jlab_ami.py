"""Focused tests for scripts/ec2_jlab.sh AMI selection reliability.

AMI resolution is exercised through the internal `EC2_JLAB_TEST_FN=resolve_ami`
test hook against a fake `aws` executable placed first on PATH. No real AWS
credentials, metadata, or network are used, no subprocess is launched with
shell=True, and no test may ever reach `run-instances`.
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

PARAM = "/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64"
CONFIGURE_CALL = "configure get region"
SSM_CALL = (
    f"--region us-east-1 ssm get-parameter --name {PARAM} "
    "--query Parameter.Value --output text"
)
VALIDATE_QUERY = "Images[0].[ImageId,State,Architecture,RootDeviceType]"
SEARCH_QUERY = "sort_by(Images, &CreationDate)[-1].ImageId"
DEFAULT_AMI = "ami-0c49f867152a8ba12"

SSM_ACCESS_DENIED = (
    "An error occurred (AccessDeniedException) when calling the GetParameter "
    "operation: User: arn:aws:iam::123456789012:user/jane is not authorized "
    "to perform: ssm:GetParameter on resource: TestParameter"
)
SSM_NOT_FOUND = (
    "An error occurred (ParameterNotFound) when calling the GetParameter "
    "operation: Parameter /aws/service/example not found."
)
SSM_NETWORK_ERR = (
    "Could not connect to the endpoint URL: "
    "'https://ssm.us-east-1.amazonaws.com/'"
)
EC2_UNAUTHORIZED = (
    "An error occurred (UnauthorizedOperation) when calling the "
    "DescribeImages operation: You are not authorized to perform this "
    "operation. Encoded authorization failure message: ABCXYZ"
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


def _describe_validation(args):
    err = os.environ.get("FAKE_EC2_ERR", "")
    if err:
        print(err, file=sys.stderr)
        return 255
    ami = args[args.index("--image-ids") + 1]
    mode = os.environ.get("FAKE_VALIDATE_MODE", "ok")
    if mode == "notfound":
        print(
            "An error occurred (InvalidAMIID.NotFound) when calling the "
            "DescribeImages operation: The image id does not exist",
            file=sys.stderr,
        )
        return 255
    if mode == "not-amazon":
        print("None")
        return 0
    state = "unavailable" if mode == "unavailable" else "available"
    arch = "arm64" if mode == "wrong-arch" else "x86_64"
    root = "instance-store" if mode == "not-ebs" else "ebs"
    print("%s\\t%s\\t%s\\t%s" % (ami, state, arch, root))
    return 0


def _describe_search(args):
    rc = int(os.environ.get("FAKE_FALLBACK_RC", "0"))
    if rc:
        print(os.environ.get("FAKE_FALLBACK_ERR", "ec2 describe-images failed"),
              file=sys.stderr)
        return rc
    mode = os.environ.get("FAKE_FALLBACK_MODE", "ok")
    if mode == "none":
        print("None")
        return 0
    if mode == "malformed":
        print("ami-0000-not-hex")
        return 0
    images = [t for t in os.environ.get("FAKE_FALLBACK_IMAGES", "").split("|") if t]
    if not images:
        print("%s" % os.environ.get("FAKE_FALLBACK_AMI", "ami-0c49f867152a8ba12"))
        return 0
    newest = max(images, key=lambda t: t.split(":", 1)[1])
    print(newest.split(":", 1)[0])
    return 0


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

    if args[:2] == ["ec2", "run-instances"]:
        print("run-instances must never be reached by AMI tests", file=sys.stderr)
        return 255

    if args[:2] == ["ssm", "get-parameter"]:
        rc = int(os.environ.get("FAKE_SSM_RC", "0"))
        err = os.environ.get("FAKE_SSM_ERR", "")
        if rc and err:
            print(err, file=sys.stderr)
        elif rc:
            print(
                "An error occurred (ServiceUnavailableException) when calling "
                "the GetParameter operation: The service is unavailable",
                file=sys.stderr,
            )
        out = os.environ.get("FAKE_SSM_OUT", "ami-0c49f867152a8ba12")
        if rc == 0 and out:
            print(out)
        return rc

    if args[:2] == ["ec2", "describe-images"] and "--image-ids" in args:
        return _describe_validation(args)

    if args[:2] == ["ec2", "describe-images"]:
        return _describe_search(args)

    return 0


if __name__ == "__main__":
    sys.exit(_main())
"""


class AmiRunner:
    def __init__(self, root: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        bindir = root / "bin"
        bindir.mkdir()
        aws = bindir / "aws"
        aws.write_text(FAKE_AWS, encoding="utf-8")
        aws.chmod(aws.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        self.bindir = bindir
        self.home = root / "home"
        self.home.mkdir()
        self.state_file = root / "instance-state.env"
        self.log_file = root / "aws-calls.log"

        monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
        monkeypatch.setenv("HOME", str(self.home))
        monkeypatch.setenv("INSTANCE_STATE_FILE", str(self.state_file))
        monkeypatch.delenv("AMI_ID", raising=False)
        monkeypatch.delenv("INSTANCE_ID", raising=False)
        monkeypatch.setenv("FAKE_AWS_LOG", str(self.log_file))
        monkeypatch.setenv("FAKE_AWS_REGION", "us-east-1")
        monkeypatch.setenv("EC2_JLAB_TEST_FN", "resolve_ami")

        self.returncode: int | None = None
        self.stdout = ""
        self.stderr = ""
        self._calls: list[str] = []

    def run(self, **env_overrides: str) -> None:
        self._calls = []
        self.log_file.write_text("", encoding="utf-8")

        env = dict(os.environ)
        env.update(env_overrides)
        # Defense in depth: never hand real AWS credentials to any subprocess.
        for var in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN", "AWS_PROFILE"):
            env.pop(var, None)
        env["AWS_EC2_METADATA_DISABLED"] = "true"

        # Q: the fake must be what PATH resolves, never the real aws binary.
        assert shutil.which("aws", path=env["PATH"]) == str(self.bindir / "aws"), (
            "the fake aws must be what PATH resolves, never the real binary"
        )

        proc = subprocess.run(
            [str(SCRIPT)],
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


@pytest.fixture()
def runner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> AmiRunner:
    return AmiRunner(tmp_path, monkeypatch)


def test_explicit_valid_ami_bypasses_ssm_and_search(runner: AmiRunner) -> None:
    # A: explicit valid AMI_ID succeeds with exactly one validation request.
    runner.run(AMI_ID=DEFAULT_AMI)

    assert runner.returncode == 0
    assert runner.stdout == DEFAULT_AMI + "\n"
    assert runner.calls == [
        CONFIGURE_CALL,
        f"--region us-east-1 ec2 describe-images --owners amazon --image-ids {DEFAULT_AMI} --query {VALIDATE_QUERY} --output text",
    ]
    assert not runner.called_with("ssm", "get-parameter")
    assert not runner.called_with("--filters")
    assert not runner.called_with("ec2", "run-instances")


def test_explicit_malformed_ami_exits_before_any_lookup(runner: AmiRunner) -> None:
    # B: malformed AMI_ID fails before any SSM lookup and before run-instances.
    runner.run(AMI_ID="not-an-ami")

    assert runner.returncode != 0
    assert runner.stdout == ""
    assert "Invalid AMI id not-an-ami" in runner.stderr
    assert runner.calls == [CONFIGURE_CALL]
    assert not runner.called_with("ssm", "get-parameter")
    assert not runner.called_with("ec2", "run-instances")


def test_explicit_ami_not_found(runner: AmiRunner) -> None:
    # C: valid syntax but the image does not exist.
    runner.run(AMI_ID=DEFAULT_AMI, FAKE_VALIDATE_MODE="notfound")

    assert runner.returncode != 0
    assert "does not exist" in runner.stderr
    assert not runner.called_with("ssm", "get-parameter")
    assert not runner.called_with("ec2", "run-instances")


def test_explicit_ami_unavailable(runner: AmiRunner) -> None:
    # D: image exists but is not in the available state.
    runner.run(AMI_ID=DEFAULT_AMI, FAKE_VALIDATE_MODE="unavailable")

    assert runner.returncode != 0
    assert "is not available (state=unavailable)" in runner.stderr
    assert not runner.called_with("ssm", "get-parameter")
    assert not runner.called_with("ec2", "run-instances")


def test_explicit_ami_wrong_architecture(runner: AmiRunner) -> None:
    # E: image is not x86_64.
    runner.run(AMI_ID=DEFAULT_AMI, FAKE_VALIDATE_MODE="wrong-arch")

    assert runner.returncode != 0
    assert "is not x86_64 (architecture=arm64)" in runner.stderr
    assert not runner.called_with("ec2", "run-instances")


def test_explicit_ami_not_ebs_backed(runner: AmiRunner) -> None:
    # F: image is not EBS-backed.
    runner.run(AMI_ID=DEFAULT_AMI, FAKE_VALIDATE_MODE="not-ebs")

    assert runner.returncode != 0
    assert "is not EBS-backed (root-device-type=instance-store)" in runner.stderr
    assert not runner.called_with("ec2", "run-instances")


def test_explicit_ami_not_amazon_owned(runner: AmiRunner) -> None:
    # G: image filtered out by the --owners amazon request.
    runner.run(AMI_ID=DEFAULT_AMI, FAKE_VALIDATE_MODE="not-amazon")

    assert runner.returncode != 0
    assert "does not exist or is not Amazon-owned" in runner.stderr
    assert not runner.called_with("ec2", "run-instances")


def test_ssm_succeeds_with_valid_ami(runner: AmiRunner) -> None:
    # H: SSM returns a valid AMI; no fallback is attempted.
    runner.run()

    assert runner.returncode == 0
    assert runner.stdout == DEFAULT_AMI + "\n"
    assert runner.calls == [CONFIGURE_CALL, SSM_CALL]
    assert not runner.called_with("--filters")


@pytest.mark.parametrize("ssm_out", ["", "None"])
def test_ssm_empty_or_none_falls_back(runner: AmiRunner, ssm_out: str) -> None:
    # I: SSM succeeds (rc 0) but yields no usable value; fallback is used.
    runner.run(FAKE_SSM_RC="0", FAKE_SSM_OUT=ssm_out)

    assert runner.returncode == 0
    assert runner.stdout == DEFAULT_AMI + "\n"
    assert runner.called_with("ssm", "get-parameter")
    assert runner.called_with("ec2", "describe-images", "--owners", "amazon", "--filters")


@pytest.mark.parametrize(
    "ssm_err",
    [SSM_ACCESS_DENIED, SSM_NOT_FOUND, SSM_NETWORK_ERR],
)
def test_ssm_failure_falls_back(runner: AmiRunner, ssm_err: str) -> None:
    # J/K/L: AccessDenied, ParameterNotFound, and network errors all fall back.
    runner.run(FAKE_SSM_RC="255", FAKE_SSM_ERR=ssm_err)

    assert runner.returncode == 0
    assert runner.stdout == DEFAULT_AMI + "\n"
    assert runner.called_with("ssm", "get-parameter")
    assert runner.called_with("ec2", "describe-images", "--owners", "amazon", "--filters")
    assert runner.called_with(
        "Name=name,Values=al2023-ami-2023.*-kernel-6.1-x86_64"
    )
    assert runner.called_with("ec2", "describe-images", "--image-ids", DEFAULT_AMI)
    assert not runner.called_with("ec2", "run-instances")


def test_fallback_no_matching_images(runner: AmiRunner) -> None:
    # M: describe-images returns no images; clear failure, never run-instances.
    runner.run(
        FAKE_SSM_RC="255",
        FAKE_SSM_ERR=SSM_ACCESS_DENIED,
        FAKE_FALLBACK_MODE="none",
    )

    assert runner.returncode != 0
    assert runner.stdout == ""
    assert "No Amazon Linux 2023 image matches the fallback filters" in runner.stderr
    assert "AMI_ID=ami-<id>" in runner.stderr
    assert not runner.called_with("ec2", "run-instances")


def test_fallback_malformed_output(runner: AmiRunner) -> None:
    # N: describe-images returns a non-AMI token.
    runner.run(FAKE_SSM_RC="255", FAKE_FALLBACK_MODE="malformed")

    assert runner.returncode != 0
    assert runner.stdout == ""
    assert "malformed AMI id" in runner.stderr
    assert "AMI_ID=ami-<id>" in runner.stderr
    assert not runner.called_with("ec2", "run-instances")


def test_fallback_api_failure_preserves_both_errors(runner: AmiRunner) -> None:
    # O: SSM and the fallback both fail; both error contexts are preserved.
    runner.run(
        FAKE_SSM_RC="255",
        FAKE_SSM_ERR=SSM_ACCESS_DENIED,
        FAKE_FALLBACK_RC="255",
        FAKE_FALLBACK_ERR=EC2_UNAUTHORIZED,
    )

    assert runner.returncode != 0
    assert runner.stdout == ""
    assert "ssm:GetParameter" in runner.stderr
    assert "UnauthorizedOperation" in runner.stderr
    assert "AMI_ID=ami-<id>" in runner.stderr
    assert "bootstrap" in runner.stderr
    assert not runner.called_with("ec2", "run-instances")


def test_fallback_selects_newest_creation_date(runner: AmiRunner) -> None:
    # P: newest CreationDate wins deterministically via the sort query.
    images = (
        "ami-0aaaa111111111111:2023-01-15T00:00:00.000Z|"
        "ami-0bbb222222222222:2025-06-01T00:00:00.000Z|"
        "ami-0ccc333333333333:2024-03-10T00:00:00.000Z"
    )
    runner.run(FAKE_SSM_RC="255", FAKE_FALLBACK_IMAGES=images)

    assert runner.returncode == 0
    assert runner.stdout == "ami-0bbb222222222222\n"
    assert runner.called_with(SEARCH_QUERY)
    assert not runner.called_with("ec2", "run-instances")


def test_fake_aws_is_the_resolved_binary(runner: AmiRunner) -> None:
    # Q: PATH resolution must always pick the fake (run() also asserts this).
    assert shutil.which("aws", path=os.environ["PATH"]) == str(runner.bindir / "aws")
