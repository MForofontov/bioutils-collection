"""Parse and execute package doctests in a bounded, temporary working directory."""

from __future__ import annotations

import argparse
import ast
import doctest
import importlib
import io
import json
import math
import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


def check(root: Path, parse_only: bool = False) -> dict[str, Any]:
    """Reject malformed examples before executing any package code."""
    errors = []
    paths = sorted((root / "bioutils_collection").rglob("*.py"))
    for path in paths:
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(
                node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
            ):
                doc = ast.get_docstring(node)
                if doc:
                    try:
                        doctest.DocTestParser().get_examples(doc)
                    except ValueError as error:
                        errors.append(
                            {
                                "file": path.relative_to(root).as_posix(),
                                "object": getattr(node, "name", "module"),
                                "error": str(error),
                            }
                        )
    report: dict[str, Any] = {
        "parse_errors": errors,
        "modules": [],
        "attempted": 0,
        "failed": 0,
        "python": sys.version,
        "platform": platform.platform(),
    }
    if not errors and not parse_only:
        sys.path.insert(0, str(root))
        for path in paths:
            name = ".".join(path.relative_to(root).with_suffix("").parts)
            if name.endswith(".__init__"):
                name = name.removesuffix(".__init__")
            module = importlib.import_module(name)
            output = io.StringIO()
            runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
            for test in doctest.DocTestFinder().find(module):
                runner.run(test, out=output.write)
            report["modules"].append(
                {
                    "name": name,
                    "attempted": runner.tries,
                    "failed": runner.failures,
                    "output": output.getvalue(),
                }
            )
            report["attempted"] += runner.tries
            report["failed"] += runner.failures
    report["status"] = "failed" if errors or report["failed"] else "passed"
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--parse-only", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be finite and positive")
    root = Path(__file__).resolve().parents[1]
    if args.worker:
        report = check(root, args.parse_only)
    else:
        environment = os.environ.copy()
        for name in (
            "OPENBLAS_NUM_THREADS",
            "OMP_NUM_THREADS",
            "MKL_NUM_THREADS",
            "NUMBA_NUM_THREADS",
        ):
            environment[name] = "1"
        command = [sys.executable, "-I", str(Path(__file__).resolve()), "--worker"]
        if args.parse_only:
            command.append("--parse-only")
        with tempfile.TemporaryDirectory(prefix="bioutils-doctests-") as directory:
            try:
                result = subprocess.run(
                    command,
                    cwd=directory,
                    env=environment,
                    capture_output=True,
                    text=True,
                    timeout=args.timeout,
                )
                try:
                    report = json.loads(result.stdout)
                except json.JSONDecodeError:
                    report = {
                        "status": "failed",
                        "error": result.stderr or result.stdout,
                    }
            except subprocess.TimeoutExpired:
                report = {
                    "status": "failed",
                    "error": f"Documentation exceeded {args.timeout} seconds",
                }
    rendered = json.dumps(report, indent=2) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    if report["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
