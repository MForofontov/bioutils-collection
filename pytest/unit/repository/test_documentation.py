"""Public examples must parse and execute outside the checkout."""

import json
import runpy
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit


def test_package_documentation(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[3]
    result = subprocess.run(
        [sys.executable, str(root / "scripts/check_documentation.py")],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=150,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads(result.stdout)
    assert report["attempted"] > 0
    assert not report["parse_errors"]


def test_documentation_scanner_rejects_broken_newline_escaping(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[3]
    check = runpy.run_path(str(root / "scripts/check_documentation.py"))["check"]
    package = tmp_path / "bioutils_collection"
    package.mkdir()
    (package / "broken.py").write_text(
        'def example():\n    """\n    >>> print("first\\nsecond")\n    first\n    second\n    """\n',
        encoding="utf-8",
    )
    report = check(tmp_path, parse_only=True)
    assert report["status"] == "failed"
    assert report["parse_errors"][0]["file"] == "bioutils_collection/broken.py"
