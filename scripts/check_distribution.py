"""Verify artifacts and install wheel/sdist separately outside the checkout.

Default installs use fresh virtual environments and resolve dependencies. Local
offline checks can explicitly reuse installed dependencies/build tools; that
mode is reported and does not establish clean dependency resolution.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import re
import subprocess
import sys
import tarfile
import tempfile
import venv
import zipfile
from email.parser import Parser
from pathlib import Path
from typing import Any

SMOKE = r"""
import importlib, importlib.metadata, json, pkgutil, re, sys
from pathlib import Path
environment, root = map(Path, sys.argv[1:3])
expected = sys.argv[3]
import bioutils_collection as package
assert Path(package.__file__).resolve().is_relative_to(environment.resolve()), package.__file__
assert not Path(package.__file__).resolve().is_relative_to(root.resolve())
assert package.__version__ == expected
assert importlib.metadata.version('bioutils-collection') == expected
modules = [package] + [importlib.import_module(info.name) for info in pkgutil.walk_packages(package.__path__, package.__name__+'.')]
for module in modules:
    for name in getattr(module, '__all__', []):
        assert hasattr(module, name), (module.__name__, name)
text = (root/'README.md').read_text(encoding='utf-8').replace(chr(13), '')
examples = re.findall(chr(96)*3+'python'+chr(10)+'(.*?)'+chr(96)*3, text, re.S)
for example in examples:
    exec(compile(example, '<README example>', 'exec'), {})
print(json.dumps({'version':package.__version__, 'modules':len(modules), 'readme_examples':len(examples), 'package_path':package.__file__, 'dependencies':{name:importlib.metadata.version(name) for name in ('numpy','scipy','scikit-learn','numba')}}))
"""


def package_version(root: Path) -> str:
    tree = ast.parse(
        (root / "bioutils_collection/_version.py").read_text(encoding="utf-8")
    )
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "__version__"
            for target in node.targets
        ):
            value = ast.literal_eval(node.value)
            if isinstance(value, str):
                return value
    raise ValueError("Cannot read the literal package version")


def verify_contents(root: Path, wheel: Path, source: Path, version: str) -> None:
    expected = {
        path.relative_to(root).as_posix(): path.read_text(encoding="utf-8")
        for path in (root / "bioutils_collection").rglob("*.py")
    }
    with zipfile.ZipFile(wheel) as archive:
        metadata_paths = [
            name for name in archive.namelist() if name.endswith(".dist-info/METADATA")
        ]
        if len(metadata_paths) != 1:
            raise ValueError(
                "Wheel must contain exactly one distribution metadata file"
            )
        metadata = Parser().parsestr(archive.read(metadata_paths[0]).decode())
        if metadata["Version"] != version:
            raise ValueError(
                f"Wheel version {metadata['Version']} differs from {version}"
            )
        if (
            re.sub(r"[-_.]+", "-", str(metadata["Name"])).lower()
            != "bioutils-collection"
        ):
            raise ValueError("Wheel distribution name differs from bioutils-collection")
        files = set(archive.namelist())
        if {
            name
            for name in files
            if name.startswith("bioutils_collection/") and name.endswith(".py")
        } != set(expected):
            raise ValueError("Wheel package source inventory differs from the checkout")
        for name, text in expected.items():
            if archive.read(name).decode("utf-8").replace("\r\n", "\n") != text:
                raise ValueError(f"Stale wheel source: {name}")
        if "bioutils_collection/py.typed" not in files or not any(
            name.endswith("/LICENSE") for name in files
        ):
            raise ValueError("Wheel lacks its typing marker or license")
    with tarfile.open(source) as archive:
        members = {
            member.name: member for member in archive.getmembers() if member.isfile()
        }
        package_files = {
            name.split("/", 1)[1]: member
            for name, member in members.items()
            if "/" in name and name.split("/", 1)[1].startswith("bioutils_collection/")
        }
        if {name for name in package_files if name.endswith(".py")} != set(expected):
            raise ValueError(
                "Source distribution package inventory differs from checkout"
            )
        for name, text in expected.items():
            stream = archive.extractfile(package_files[name])
            if (
                stream is None
                or stream.read().decode("utf-8").replace("\r\n", "\n") != text
            ):
                raise ValueError(f"Stale source distribution: {name}")
        metadata_member = next(
            (
                member
                for name, member in members.items()
                if name.count("/") == 1 and name.endswith("/PKG-INFO")
            ),
            None,
        )
        if metadata_member is None:
            raise ValueError("Source distribution lacks PKG-INFO")
        stream = archive.extractfile(metadata_member)
        metadata = Parser().parsestr(stream.read().decode()) if stream else {}
        if metadata.get("Version") != version:
            raise ValueError("Source distribution version differs from checkout")
        if (
            re.sub(r"[-_.]+", "-", str(metadata.get("Name"))).lower()
            != "bioutils-collection"
        ):
            raise ValueError(
                "Source distribution name differs from bioutils-collection"
            )
        if "bioutils_collection/py.typed" not in package_files or not any(
            name.endswith("/LICENSE") for name in members
        ):
            raise ValueError("Source distribution lacks its typing marker or license")


def install_and_check(
    artifact: Path, root: Path, version: str, reuse: bool, no_isolation: bool
) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="bioutils-distribution-") as temporary:
        directory = Path(temporary)
        environment = directory / "environment"
        venv.EnvBuilder(with_pip=True).create(environment)
        python = environment / (
            "Scripts/python.exe" if os.name == "nt" else "bin/python"
        )
        if reuse:
            site = environment / (
                "Lib/site-packages"
                if os.name == "nt"
                else f"lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages"
            )
            paths = [
                str(Path(path).resolve())
                for path in sys.path
                if path and Path(path).name in {"site-packages", "dist-packages"}
            ]
            (site / "audit-dependencies.pth").write_text(
                "\n".join(paths) + "\n", encoding="utf-8"
            )
        args = [str(python), "-m", "pip", "install", "--disable-pip-version-check"]
        if reuse:
            # Force our artifact into this environment despite any reused
            # editable distribution metadata for the same package.
            args += ["--no-deps", "--ignore-installed"]
        if no_isolation and artifact.suffix != ".whl":
            args += ["--no-build-isolation"]
        child_environment = {
            key: value
            for key, value in os.environ.items()
            if key not in {"PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV"}
        }
        subprocess.run(
            [*args, str(artifact)],
            cwd=directory,
            env=child_environment,
            check=True,
            timeout=300,
        )
        subprocess.run(
            [str(python), "-m", "pip", "check"],
            cwd=directory,
            env=child_environment,
            check=True,
            timeout=60,
        )
        smoke = directory / "installed_smoke.py"
        smoke.write_text(SMOKE, encoding="utf-8")
        result = subprocess.run(
            [str(python), "-I", str(smoke), str(environment), str(root), version],
            cwd=directory,
            env=child_environment,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
        if result.returncode:
            sys.stderr.write(result.stdout + result.stderr)
            result.check_returncode()
        return {
            "artifact": artifact.name,
            "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
            **json.loads(result.stdout.splitlines()[-1]),
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--reuse-dependencies", action="store_true")
    parser.add_argument("--no-build-isolation", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    version = package_version(root)
    wheels, sources = (
        list(args.artifact_dir.glob("*.whl")),
        list(args.artifact_dir.glob("*.tar.gz")),
    )
    if len(wheels) != 1 or len(sources) != 1:
        parser.error("artifact directory must contain exactly one wheel and one sdist")
    wheel, source = wheels[0].resolve(), sources[0].resolve()
    report: dict[str, Any] = {
        "version": version,
        "dependency_mode": "reused" if args.reuse_dependencies else "fresh resolution",
        "build_isolation": not args.no_build_isolation,
        "installations": [],
    }
    try:
        verify_contents(root, wheel, source, version)
        for artifact in (wheel, source):
            report["installations"].append(
                install_and_check(
                    artifact,
                    root,
                    version,
                    args.reuse_dependencies,
                    args.no_build_isolation,
                )
            )
        report["status"] = "passed"
    except Exception as error:
        report.update(status="failed", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(
                json.dumps(report, indent=2) + "\n", encoding="utf-8"
            )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
