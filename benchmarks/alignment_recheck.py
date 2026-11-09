"""Calibrated alignment checks with unchanged controls and frozen sources.

Run alone, without concurrent heavy checks. A control drift above 10% makes the
timings inconclusive; repeat once and retain both reports. Memory is measured
separately, and every complete output must match before any speed claim.
"""

from __future__ import annotations

import argparse
import gc
import importlib.metadata
import json
import os
import platform
import random
import runpy
import statistics
import sys
import tempfile
import time
import tracemalloc
from functools import partial
from pathlib import Path
from typing import Any

PAIRED = runpy.run_path(str(Path(__file__).with_name("paired_audit.py")))


def worker(request: dict[str, Any]) -> dict[str, Any]:
    for name in (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMBA_NUM_THREADS",
    ):
        os.environ[name] = "1"
    sys.path.insert(0, request["source"])
    started = time.perf_counter()
    import bioutils_collection as package

    imported = time.perf_counter() - started
    assert Path(package.__file__).resolve().is_relative_to(Path(request["source"]))
    results = []
    for case in request["cases"]:
        rng = random.Random(20261102 + case["rows"] * 1001 + case["columns"])
        if case["pattern"] == "random":
            left = "".join(rng.choices("ACGT", k=case["rows"]))
            right = "".join(rng.choices("ACGT", k=case["columns"]))
        elif case["pattern"] == "repetitive":
            left = ("ACGT" * case["rows"])[: case["rows"]]
            right = ("TGCA" * case["columns"])[: case["columns"]]
        else:
            left, right = "A" * case["rows"], "C" * case["columns"]
        if case["name"] == "global_alignment":
            operation = partial(package.needleman_wunsch, left, right, **case["scores"])
        else:
            operation = partial(package.levenshtein_distance, left, right)

        gc.collect()
        if request["stage"] == "time":
            operation()  # Specialization/warmup is outside the timed batch.
            loops = 1
            # Calibration is separate from the reported batch.
            while True:
                started = time.perf_counter()
                for _ in range(loops):
                    operation()
                if time.perf_counter() - started >= 0.1:
                    break
                loops *= 2
            while True:
                started = time.perf_counter()
                for _ in range(loops):
                    result = operation()
                elapsed = time.perf_counter() - started
                if elapsed >= 0.1:
                    break
                loops *= 2
            metrics = {
                "seconds": elapsed / loops,
                "batch_seconds": elapsed,
                "loops": loops,
            }
        else:
            before = PAIRED["peak_rss"]()
            result = operation()
            after = PAIRED["peak_rss"]()
            del result
            gc.collect()
            tracemalloc.start()
            result = operation()
            _, allocated = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            metrics = {
                "rss_before_bytes": before,
                "peak_rss_bytes": after,
                "peak_python_bytes": allocated,
            }
        results.append(
            {"id": case["id"], "digest": PAIRED["digest"](result, None), **metrics}
        )
    return {"import_seconds": imported, "results": results}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", default="186a85f")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    cases = []
    for rows, columns in ((100, 50), (300, 300), (1000, 100), (1000, 1000)):
        for pattern in ("random", "repetitive", "tie-heavy"):
            for name in ("global_alignment", "edit_distance"):
                cases.append(
                    {
                        "id": f"{name}/{rows}x{columns}/{pattern}",
                        "name": name,
                        "rows": rows,
                        "columns": columns,
                        "pattern": pattern,
                        "scores": {"match": 0, "mismatch": 0, "gap": 0}
                        if pattern == "tie-heavy"
                        else {"match": 1, "mismatch": -1, "gap": -1},
                    }
                )
    with tempfile.TemporaryDirectory(prefix="bioutils-calibrated-") as directory:
        sources = {role: Path(directory) / role for role in ("baseline", "current")}
        manifests = {
            role: PAIRED["snapshot"](
                root,
                path,
                args.baseline if role == "baseline" else None,
                args.output.resolve(),
            )
            for role, path in sources.items()
        }
        script = sources["current"] / "benchmarks/alignment_recheck.py"
        runs: dict[str, dict[str, list[dict[str, Any]]]] = {
            case["id"]: {role: [] for role in sources} for case in cases
        }
        digests: dict[str, dict[str, str]] = {case["id"]: {} for case in cases}
        startup = []
        for sample in range(5):
            order = (
                ("baseline", "current") if sample % 2 == 0 else ("current", "baseline")
            )
            for role in order:
                result = PAIRED["run_worker"](
                    script,
                    sources[role],
                    script,
                    role,
                    "time",
                    cases if sample % 2 == 0 else list(reversed(cases)),
                )
                startup.append(
                    {
                        "sample": sample,
                        "role": role,
                        "import_seconds": result["import_seconds"],
                    }
                )
                for row in result["results"]:
                    old = digests[row["id"]].setdefault(role, row["digest"])
                    assert old == row["digest"], "Nondeterministic output"
                    runs[row["id"]][role].append(row)
            print(f"Calibrated round {sample + 1}/5 complete", flush=True)
        measurements = []
        for index, case in enumerate(cases):
            memory = {}
            for role in sources:
                row = PAIRED["run_worker"](
                    script, sources[role], script, role, "memory", [case]
                )["results"][0]
                assert row["digest"] == digests[case["id"]][role]
                memory[role] = {
                    key: value
                    for key, value in row.items()
                    if key not in {"id", "digest"}
                }
            medians = {
                role: statistics.median(
                    row["seconds"] for row in runs[case["id"]][role]
                )
                for role in sources
            }
            ratio = medians["current"] / medians["baseline"]
            allocated = memory["current"]["peak_python_bytes"] / max(
                1, memory["baseline"]["peak_python_bytes"]
            )
            equivalent = len(set(digests[case["id"]].values())) == 1
            measurements.append(
                {
                    **case,
                    "equivalent": equivalent,
                    "digests": digests[case["id"]],
                    "runs": runs[case["id"]],
                    "median_seconds": medians,
                    "runtime_ratio": ratio,
                    "python_memory_ratio": allocated,
                    "memory": memory,
                    "meets_gate": equivalent
                    and (ratio <= 0.9 or (allocated <= 0.75 and ratio <= 1.1)),
                }
            )
            print(f"Memory {index + 1}/{len(cases)} complete", flush=True)
        controls = [row for row in measurements if row["name"] == "edit_distance"]
        control_path = "bioutils_collection/alignment_functions/levenshtein_distance.py"
        assert (
            manifests["baseline"]["files"][control_path]
            == manifests["current"]["files"][control_path]
        )
        report = {
            "baseline": args.baseline,
            "snapshot_manifests": manifests,
            "python": sys.version,
            "platform": platform.platform(),
            "dependencies": {
                name: importlib.metadata.version(name)
                for name in ("numpy", "scipy", "scikit-learn", "numba")
            },
            "samples": 5,
            "minimum_batch_seconds": 0.1,
            "seed": 20261102,
            "startup": startup,
            "timing_conclusive": all(
                0.9 <= row["runtime_ratio"] <= 1.1 for row in controls
            ),
            "all_outputs_equivalent": all(row["equivalent"] for row in measurements),
            "measurements": measurements,
            "notes": "Frozen packages and measuring code; actual sequence dimensions recorded. Calibration, imports, warmup, output hashing and memory probes are outside timed batches. Python allocations and absolute native process peak RSS are separate. Controls are unchanged byte-for-byte. No speed claim if control medians drift above 10%; repeat once, retaining both reports. No JIT or process-pool work is timed by these alignment-only cases.",
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"Saved {args.output}; timing conclusive: {report['timing_conclusive']}")


if __name__ == "__main__":
    if "--worker" in sys.argv:
        print(json.dumps(worker(json.load(sys.stdin))))
    else:
        main()
