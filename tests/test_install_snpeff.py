"""snpEff setup uses the S3 bundle before any upstream download."""

from __future__ import annotations

import os
from pathlib import Path
import stat
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "install_snpeff.sh"


def _exe(path: Path, text: str) -> None:
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


def _java(bin_dir: Path) -> None:
    _exe(
        bin_dir / "java",
        "#!/bin/sh\nprintf '%s\\n' '    java.specification.version = 21' >&2\n",
    )


def _run(prefix: Path, bin_dir: Path, extra_env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env.update(extra_env)
    return subprocess.run(
        ["bash", str(SCRIPT), str(prefix)],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )


def test_installed_tree_does_not_touch_s3(tmp_path: Path) -> None:
    prefix = tmp_path / "snpeff"
    (prefix / "data" / "GRCh38.115").mkdir(parents=True)
    (prefix / "snpEff.jar").write_bytes(b"jar")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _java(bin_dir)
    log = tmp_path / "aws.log"
    _exe(bin_dir / "aws", f"#!/bin/sh\necho \"$@\" >> {log}\nexit 1\n")
    result = _run(prefix, bin_dir, {})
    assert result.returncode == 0, result.stderr
    assert not log.exists()
    assert (prefix / "bin" / "snpEff").is_file()


def test_s3_hit_downloads_the_bundle_and_skips_upstream(tmp_path: Path) -> None:
    prefix = tmp_path / "snpeff"
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _java(bin_dir)
    bundle = tmp_path / "bundle.tar.gz"
    staging = tmp_path / "stage"
    (staging / "data" / "GRCh38.115").mkdir(parents=True)
    (staging / "snpEff.jar").write_bytes(b"jar")
    (staging / "snpEff.config").write_text("db\n")
    with tarfile.open(bundle, "w:gz") as tar:
        for name in ("snpEff.jar", "snpEff.config", "data"):
            tar.add(staging / name, arcname=name)
    log = tmp_path / "aws.log"
    _exe(
        bin_dir / "aws",
        f"""#!/bin/sh
echo "$@" >> {log}
if [ "$1" = s3 ] && [ "$2" = ls ]; then
  exit 0
fi
if [ "$1" = s3 ] && [ "$2" = cp ]; then
  cp {bundle} "$5"
  exit 0
fi
exit 1
""",
    )
    _exe(bin_dir / "curl", f"#!/bin/sh\necho curl >> {log}\nexit 1\n")
    result = _run(prefix, bin_dir, {"SNPEFF_S3_URI": "s3://l1tx-data/public/tools/snpeff/bundle.tar.gz"})
    assert result.returncode == 0, result.stderr + result.stdout
    calls = log.read_text()
    assert "s3 ls" in calls
    assert "s3 cp" in calls
    assert "curl" not in calls
    assert (prefix / "data" / "GRCh38.115").is_dir()
