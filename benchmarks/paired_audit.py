"""Compare a Git baseline with a complete working-tree snapshot.

Five counterbalanced timing rounds run separately from allocation/process-memory
measurements. Closed/streaming selector references are explicitly distinguished
from old implementations whose output was incorrect. No repository checkout or
tracked file is changed. Run with --output benchmarks/paired_audit_results.json.
"""

from __future__ import annotations

import argparse
import ctypes
import gc
import hashlib
import importlib.metadata
import io
import json
import os
import platform
import random
import runpy
import shutil
import statistics
import subprocess
import sys
import tarfile
import tempfile
import time
import tracemalloc
from collections.abc import Callable
from pathlib import Path, PurePosixPath
from typing import Any

NAMES = (
    "open_syncmers",
    "closed_syncmers",
    "syncmer_density",
    "minimizer_sketch",
    "canonical_sketch",
    "streaming_sketch",
    "kmer_frequency",
    "gc_windows",
    "translation_python",
    "translation_numba",
    "fasta_parse",
    "global_alignment",
    "local_alignment",
    "edit_distance",
)
OPTIMIZATIONS = {
    "open_syncmers",
    "closed_syncmers",
    "syncmer_density",
    "minimizer_sketch",
    "canonical_sketch",
    "streaming_sketch",
    "kmer_frequency",
}


def peak_rss() -> int:
    if os.name != "nt":
        import resource

        value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return int(value if sys.platform == "darwin" else value * 1024)
    from ctypes import wintypes

    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("faults", wintypes.DWORD)] + [
            (name, ctypes.c_size_t)
            for name in (
                "peak",
                "working",
                "pool_peak",
                "pool",
                "nonpaged_peak",
                "nonpaged",
                "pagefile",
                "pagefile_peak",
            )
        ]

    counters = Counters()
    counters.cb = ctypes.sizeof(counters)
    process = ctypes.windll.kernel32.GetCurrentProcess
    process.restype = wintypes.HANDLE
    read = ctypes.windll.psapi.GetProcessMemoryInfo
    read.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
    if not read(process(), ctypes.byref(counters), counters.cb):
        raise ctypes.WinError()
    return int(counters.peak)


def make_case(
    ns: dict[str, Any], case: dict[str, Any], role: str, directory: Path
) -> tuple[Callable[[], Any], Path | None]:
    size, pattern, name = case["size"], case["pattern"], case["name"]
    rng = random.Random(20261030 + size)
    seq = (
        "".join(rng.choices("ACGT", k=size))
        if pattern == "random"
        else ("ACGT" * ((size + 3) // 4))[:size]
    )
    k, s, w = 31, 9, 32
    reference = role == "baseline" and case["comparison"] == "corrected-reference"
    if name == "syncmer_density":
        return lambda: ns["syncmer_density"](seq, k, s, "open"), None
    if name in {"open_syncmers", "closed_syncmers"}:
        method = "open" if name == "open_syncmers" else "closed"
        function = ns["naive_syncmers"] if reference else ns["syncmer_positions"]
        return lambda: function(seq, k, s, method), None
    if name == "streaming_sketch":
        if reference:
            return lambda: ns["naive_sketch"](seq, k, w), None
        return lambda: ns["minimizer_sketch_from_generator"](
            ((seq[i : i + k], i) for i in range(size - k + 1)), w
        ), None
    functions = {
        "minimizer_sketch": "minimizer_sketch",
        "canonical_sketch": "canonical_minimizer_sketch",
    }
    if name in functions:
        return lambda: ns[functions[name]](seq, k, w), None
    if name == "kmer_frequency":
        return lambda: ns["kmer_frequency"](seq, k), None
    if name == "gc_windows":
        return lambda: ns["gc_content_windows"](seq, 100, 10), None
    if name.startswith("translation_"):
        seq = seq[: size // 3 * 3]
        function = ns[
            "translate_dna_fast"
            if name.endswith("numba")
            else "translate_dna_to_protein"
        ]
        return lambda: function(seq), None
    if name == "fasta_parse":
        fasta = "".join(f">r{i}\n{seq[i : i + 300]}\n" for i in range(0, size, 300))
        return lambda: list(ns["parse_fasta"](fasta)), None
    if name in {"global_alignment", "local_alignment", "edit_distance"}:
        width = min(300, max(20, size // 100))
        left, right = seq[:width], seq[-width // 2 :]
        function = ns[
            {
                "global_alignment": "needleman_wunsch",
                "local_alignment": "smith_waterman",
                "edit_distance": "levenshtein_distance",
            }[name]
        ]
        return lambda: function(left, right), None
    if name == "batch_process":
        sequences = ["ATG" + "GCC" * i + "TAA" for i in range(12)]
        return lambda: ns["translate_batch"](sequences, n_processes=2), None
    if name == "fasta_stream":
        source, output = directory / "input.fa", directory / "output.fa"
        seq = "ATG" * (size // 3)
        with source.open("w", encoding="utf-8") as stream:
            width = 300 if pattern == "many_records" else len(seq)
            for i in range(0, len(seq), width):
                stream.write(f">r{i}\n{seq[i : i + width]}\n")
        return lambda: ns["translate_fasta_streaming"](source, output), output
    raise ValueError(name)


def digest(result: Any, output: Path | None) -> str:
    value = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    if output is not None:
        value += output.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(value).hexdigest()


def worker(request: dict[str, Any]) -> dict[str, Any]:
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[variable] = "1"
    sys.path.insert(0, request["source"])
    started = time.perf_counter()
    ns = runpy.run_path(request["reference_script"])
    imported = time.perf_counter() - started
    # Ensure editable imports cannot silently bypass the selected snapshot.
    import bioutils_collection

    if (
        not Path(bioutils_collection.__file__)
        .resolve()
        .is_relative_to(Path(request["source"]).resolve())
    ):
        raise RuntimeError("Benchmark imported outside its selected snapshot")
    jit = None
    if any(case["name"] == "translation_numba" for case in request["cases"]):
        started = time.perf_counter()
        ns["translate_dna_fast"]("ATGGCC" * 100)
        jit = time.perf_counter() - started
    results = []
    with tempfile.TemporaryDirectory(prefix="bioutils-paired-worker-") as temporary:
        for case in request["cases"]:
            function, output = make_case(ns, case, request["role"], Path(temporary))
            gc.collect()
            if request["stage"] == "time":
                loops = (
                    20 if case["size"] <= 1000 else (3 if case["size"] <= 10000 else 1)
                )
                if case["name"] in {"batch_process", "fasta_stream"}:
                    loops = 1
                started = time.perf_counter()
                for _ in range(loops):
                    result = function()
                metrics = {
                    "seconds": (time.perf_counter() - started) / loops,
                    "loops": loops,
                }
            else:
                # Native peak is sampled before tracing or hashing results.
                before = peak_rss()
                result = function()
                after = peak_rss()
                del result
                gc.collect()
                tracemalloc.start()
                result = function()
                _, allocated = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                metrics = {
                    "rss_before_bytes": before,
                    "peak_rss_bytes": after,
                    "peak_python_bytes": allocated,
                }
            results.append(
                {"id": case["id"], "digest": digest(result, output), **metrics}
            )
            del result
    return {
        "import_seconds": imported,
        "numba_first_call_seconds": jit,
        "results": results,
    }


def snapshot(
    root: Path, directory: Path, ref: str | None, excluded: Path
) -> dict[str, Any]:
    directory.mkdir()
    if ref is not None:
        archive = subprocess.run(
            ["git", "archive", ref], cwd=root, check=True, capture_output=True
        ).stdout
        with tarfile.open(fileobj=io.BytesIO(archive)) as source:
            for member in source.getmembers():
                parts = PurePosixPath(member.name)
                if parts.is_absolute() or ".." in parts.parts or not member.isfile():
                    continue
                target = directory.joinpath(*parts.parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                stream = source.extractfile(member)
                if stream is not None:
                    target.write_bytes(stream.read())
    else:
        files = subprocess.run(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
            cwd=root,
            check=True,
            capture_output=True,
        ).stdout
        for relative in sorted(set(files.decode().split("\0")) - {""}):
            path = root / relative
            if path.resolve() == excluded or not path.is_file():
                continue
            target = directory / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
    hashes = {
        path.relative_to(directory).as_posix(): hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }
    package_hashes = {
        name: value
        for name, value in hashes.items()
        if name.startswith("bioutils_collection/")
    }
    return {
        "tree_sha256": hashlib.sha256(
            json.dumps(hashes, sort_keys=True).encode()
        ).hexdigest(),
        "package_sha256": hashlib.sha256(
            json.dumps(package_hashes, sort_keys=True).encode()
        ).hexdigest(),
        "files": hashes,
    }


def run_worker(
    script: Path,
    source: Path,
    reference: Path,
    role: str,
    stage: str,
    cases: list[dict[str, Any]],
) -> dict[str, Any]:
    request = {
        "source": str(source),
        "reference_script": str(reference),
        "role": role,
        "stage": stage,
        "cases": cases,
    }
    result = subprocess.run(
        [sys.executable, "-I", str(script), "--worker"],
        input=json.dumps(request),
        capture_output=True,
        text=True,
        check=True,
        timeout=240,
    )
    return json.loads(result.stdout)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", default="9e1cc06")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sizes", nargs="+", type=int, default=[1000, 10000, 100000])
    parser.add_argument("--samples", type=int, default=5)
    parser.add_argument(
        "--cases",
        nargs="+",
        choices=[*NAMES, "fasta_stream", "batch_process"],
        default=[*NAMES, "fasta_stream", "batch_process"],
    )
    args = parser.parse_args()
    if args.samples < 5 or min(args.sizes) < 100:
        parser.error("use at least five samples and sizes >= 100")
    root = Path(__file__).resolve().parents[1]
    script = Path(__file__).resolve()
    cases = []
    for name in args.cases:
        patterns = (
            ("many_records", "one_record")
            if name == "fasta_stream"
            else ("random", "repetitive")
        )
        for size in args.sizes:
            for pattern in patterns:
                if name == "batch_process" and (
                    size != args.sizes[0] or pattern != "random"
                ):
                    continue
                for comparison in (
                    ("implementation", "corrected-reference")
                    if name in {"closed_syncmers", "streaming_sketch"}
                    else ("implementation",)
                ):
                    cases.append(
                        {
                            "id": f"{name}/{size}/{pattern}/{comparison}",
                            "name": name,
                            "size": size,
                            "pattern": pattern,
                            "comparison": comparison,
                        }
                    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="bioutils-paired-snapshots-") as temporary:
        directories = {role: Path(temporary) / role for role in ("baseline", "current")}
        manifests = {
            role: snapshot(
                root,
                directory,
                args.baseline if role == "baseline" else None,
                args.output.resolve(),
            )
            for role, directory in directories.items()
        }
        reference = directories["current"] / "benchmarks/major_audit.py"
        # Freeze the measuring program as well as the packages and references.
        script = directories["current"] / "benchmarks/paired_audit.py"
        timings: dict[str, dict[str, list[float]]] = {
            case["id"]: {"baseline": [], "current": []} for case in cases
        }
        digests: dict[str, dict[str, str]] = {case["id"]: {} for case in cases}
        startup = []
        for sample in range(args.samples):
            order = (
                ("baseline", "current") if sample % 2 == 0 else ("current", "baseline")
            )
            ordered_cases = cases if sample % 2 == 0 else list(reversed(cases))
            for role in order:
                result = run_worker(
                    script, directories[role], reference, role, "time", ordered_cases
                )
                startup.append(
                    {
                        "sample": sample,
                        "role": role,
                        "import_seconds": result["import_seconds"],
                        "numba_first_call_seconds": result["numba_first_call_seconds"],
                    }
                )
                for row in result["results"]:
                    timings[row["id"]][role].append(row["seconds"])
                    previous = digests[row["id"]].setdefault(role, row["digest"])
                    if previous != row["digest"]:
                        raise RuntimeError(f"Nondeterministic output: {row['id']}")
            print(f"Timing round {sample + 1}/{args.samples} complete", flush=True)
        measurements = []
        for index, case in enumerate(cases):
            memory = {}
            for role in ("baseline", "current"):
                row = run_worker(
                    script, directories[role], reference, role, "memory", [case]
                )["results"][0]
                if row["digest"] != digests[case["id"]][role]:
                    raise RuntimeError(f"Memory run changed output: {case['id']}")
                memory[role] = {
                    key: value
                    for key, value in row.items()
                    if key not in {"id", "digest"}
                }
            medians = {
                role: statistics.median(values)
                for role, values in timings[case["id"]].items()
            }
            equivalent = (
                digests[case["id"]]["baseline"] == digests[case["id"]]["current"]
            )
            ratio = medians["current"] / medians["baseline"]
            memory_ratio = memory["current"]["peak_python_bytes"] / max(
                memory["baseline"]["peak_python_bytes"], 1
            )
            verdict = (
                "semantic correction; timing not comparable"
                if not equivalent
                else (
                    "meets balanced gate"
                    if ratio <= 0.9 or (memory_ratio <= 0.75 and ratio <= 1.1)
                    else "runtime regression"
                    if ratio > 1.1
                    else "neutral"
                )
            )
            measurements.append(
                {
                    **case,
                    "equivalent": equivalent,
                    "digests": digests[case["id"]],
                    "seconds": timings[case["id"]],
                    "median_seconds": medians,
                    "runtime_ratio": ratio,
                    "python_memory_ratio": memory_ratio,
                    "memory": memory,
                    "verdict": verdict,
                    "optimization_candidate": case["name"] in OPTIMIZATIONS,
                }
            )
            if (index + 1) % 10 == 0:
                print(f"Memory checks {index + 1}/{len(cases)} complete", flush=True)
        baseline = subprocess.run(
            ["git", "rev-parse", args.baseline],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        report = {
            "baseline_commit": baseline,
            "snapshot_manifests": manifests,
            "python": sys.version,
            "platform": platform.platform(),
            "dependencies": {
                name: importlib.metadata.version(name)
                for name in ("numpy", "scipy", "scikit-learn", "numba")
            },
            "seed": 20261030,
            "samples": args.samples,
            "startup": startup,
            "notes": "Complete Git-tracked/nonignored-untracked snapshots exclude only this generated report. Timing rounds counterbalance version and case order. Operation timings exclude imports, preparation, JIT warmup, output hashing and memory instrumentation. RSS is absolute process peak including imports, preparation and JIT; pre-operation peak is also recorded. Python allocations are separate un-timed calls. Non-equivalent old outputs have no speedup claim. Shared-host load is uncontrolled.",
            "measurements": measurements,
        }
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"Saved {len(measurements)} comparisons to {args.output}")


if __name__ == "__main__":
    if "--worker" in sys.argv:
        print(json.dumps(worker(json.load(sys.stdin))))
    else:
        main()
