"""Reject stale or incomplete artifacts before any installation is attempted."""

import io
import runpy
import tarfile
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
CHECKER = runpy.run_path(str(ROOT / "scripts/check_distribution.py"))


def artifacts(directory: Path, corruption: str) -> tuple[Path, Path, str]:
    version = CHECKER["package_version"](ROOT)
    sources = {
        path.relative_to(ROOT).as_posix(): path.read_bytes()
        for path in (ROOT / "bioutils_collection").rglob("*.py")
    }
    metadata = f"Name: bioutils-collection\nVersion: {version}\n".encode()
    wheel = directory / "package.whl"
    source = directory / "package.tar.gz"
    with zipfile.ZipFile(wheel, "w") as archive:
        for name, content in sources.items():
            if corruption == "wheel-source" and name.endswith("/_version.py"):
                content += b"\n# stale content\n"
            archive.writestr(name, content)
        archive.writestr(
            "package.dist-info/METADATA",
            metadata if corruption != "wheel-version" else b"Version: 0.0.0\n",
        )
        archive.writestr(
            "package.dist-info/licenses/LICENSE", (ROOT / "LICENSE").read_bytes()
        )
        if corruption != "wheel-marker":
            archive.writestr("bioutils_collection/py.typed", b"")
    with tarfile.open(source, "w:gz") as archive:
        files = {
            **sources,
            "PKG-INFO": metadata,
            "LICENSE": (ROOT / "LICENSE").read_bytes(),
        }
        if corruption != "source-marker":
            files["bioutils_collection/py.typed"] = b""
        if corruption == "source-version":
            files["PKG-INFO"] = b"Version: 0.0.0\n"
        if corruption == "source-content":
            files["bioutils_collection/_version.py"] += b"\n# stale content\n"
        for name, content in files.items():
            member = tarfile.TarInfo("package/" + name)
            member.size = len(content)
            archive.addfile(member, io.BytesIO(content))
    return wheel, source, version


def test_matching_artifact_contents(tmp_path: Path) -> None:
    wheel, source, version = artifacts(tmp_path, "")
    CHECKER["verify_contents"](ROOT, wheel, source, version)


@pytest.mark.parametrize(
    "corruption",
    [
        "wheel-version",
        "wheel-source",
        "wheel-marker",
        "source-version",
        "source-content",
        "source-marker",
    ],
)
def test_invalid_artifact_contents(tmp_path: Path, corruption: str) -> None:
    wheel, source, version = artifacts(tmp_path, corruption)
    with pytest.raises(ValueError):
        CHECKER["verify_contents"](ROOT, wheel, source, version)
