"""Reproducible correctness-first benchmarks; no timing gates in functional CI.

Run: python benchmarks/major_audit.py --output benchmarks/major_audit_results.json
Peak bytes are Python allocations measured separately with tracemalloc, not RSS.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import platform
import random
import statistics
import subprocess
import sys
import tempfile
import time
import tracemalloc
from collections import Counter
from collections.abc import Callable
from pathlib import Path
from typing import Any

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

from bioutils_collection import (  # noqa: E402
    canonical_minimizer_sketch,
    gc_content_windows,
    generate_minimizers,
    kmer_frequency,
    levenshtein_distance,
    minimizer_sketch,
    minimizer_sketch_from_generator,
    needleman_wunsch,
    parse_fasta,
    reverse_complement,
    smith_waterman,
    syncmer_density,
    syncmer_positions,
)
from bioutils_collection.translation_functions import (  # noqa: E402
    translate_batch,
    translate_dna_fast,
    translate_dna_to_protein,
    translate_fasta_streaming,
)


def naive_syncmers(seq: str, k: int, s: int, method: str) -> list[int]:
    result = []
    for start in range(len(seq) - k + 1):
        smers = [seq[j : j + s] for j in range(start, start + k - s + 1)]
        offset = smers.index(min(smers))
        if offset == 0 or (method == "closed" and offset == k - s):
            result.append(start)
    return result


def naive_sketch(
    seq: str, k: int, w: int, canonical: bool = False
) -> dict[str, list[int]]:
    pairs = [(seq[i : i + k], i) for i in range(len(seq) - k + 1)]
    if canonical:
        pairs = [(min(value, reverse_complement(value)), i) for value, i in pairs]
    result: dict[str, list[int]] = {}
    for i in range(max(1, len(pairs) - w + 1)):
        value, pos = min(pairs[i : i + w], key=lambda pair: pair[0])
        result.setdefault(value, []).append(pos)
    return result


def measure(function: Callable[[], Any], repeats: int) -> dict[str, Any]:
    timings = []
    for _ in range(repeats):
        start = time.perf_counter()
        result = function()
        timings.append(time.perf_counter() - start)
        del result
    tracemalloc.start()
    result = function()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    del result
    return {
        "seconds": timings,
        "median_seconds": statistics.median(timings),
        "peak_python_bytes": peak,
    }


def cold_numba_process() -> dict[str, float]:
    """Use a new empty disk cache and report import and first-call separately."""
    with tempfile.TemporaryDirectory() as cache:
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                "import json,time; t=time.perf_counter(); from bioutils_collection.translation_functions import translate_dna_fast; imported=time.perf_counter()-t; t=time.perf_counter(); translate_dna_fast('ATGGCC'*100); print(json.dumps({'import_seconds':imported,'first_call_seconds':time.perf_counter()-t}))",
            ],
            env={**os.environ, "NUMBA_CACHE_DIR": cache},
            check=True,
            capture_output=True,
            text=True,
        )
        return json.loads(result.stdout)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sizes", type=int, nargs="+", default=[1000, 10000, 100000])
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1 or min(args.sizes) < 100:
        parser.error("use positive repeats and sizes >= 100")
    rows = []

    def record(
        name: str, bases: int, pattern: str, function: Callable[[], Any]
    ) -> None:
        rows.append(
            {
                "name": name,
                "bases": bases,
                "pattern": pattern,
                **measure(function, args.repeats),
            }
        )

    randomizer = random.Random(20261030)
    for size in args.sizes:
        for pattern in ("random", "repetitive"):
            seq = (
                "".join(randomizer.choices("ACGT", k=size))
                if pattern == "random"
                else ("ACGT" * ((size + 3) // 4))[:size]
            )
            k, s, w = 31, 9, 32
            for method in ("open", "closed"):
                expected = naive_syncmers(seq, k, s, method)
                assert syncmer_positions(seq, k, s, method) == expected
                assert syncmer_density(seq, k, s, method) == len(expected) / (
                    size - k + 1
                )
                record(
                    f"syncmer_{method}_reference",
                    size,
                    pattern,
                    lambda k=k, method=method, s=s, seq=seq: naive_syncmers(
                        seq, k, s, method
                    ),
                )
                record(
                    f"syncmer_{method}",
                    size,
                    pattern,
                    lambda k=k, method=method, s=s, seq=seq: syncmer_positions(
                        seq, k, s, method
                    ),
                )
            record(
                "syncmer_density_materialized",
                size,
                pattern,
                lambda k=k, s=s, seq=seq, size=size: (
                    len(syncmer_positions(seq, k, s, "closed")) / (size - k + 1)
                ),
            )
            record(
                "syncmer_density",
                size,
                pattern,
                lambda k=k, s=s, seq=seq: syncmer_density(seq, k, s, "closed"),
            )
            for canonical in (False, True):
                expected_sketch = naive_sketch(seq, k, w, canonical)
                function = canonical_minimizer_sketch if canonical else minimizer_sketch
                assert function(seq, k, w) == expected_sketch
                label = "canonical_sketch" if canonical else "minimizer_sketch"
                record(
                    f"{label}_reference",
                    size,
                    pattern,
                    lambda canonical=canonical, k=k, seq=seq, w=w: naive_sketch(
                        seq, k, w, canonical
                    ),
                )
                record(
                    label,
                    size,
                    pattern,
                    lambda function=function, k=k, seq=seq, w=w: function(seq, k, w),
                )
            pairs = [(seq[i : i + k], i) for i in range(size - k + 1)]
            assert minimizer_sketch_from_generator(iter(pairs), w) == naive_sketch(
                seq, k, w
            )
            record(
                "streaming_sketch",
                size,
                pattern,
                lambda pairs=pairs, w=w: minimizer_sketch_from_generator(
                    iter(pairs), w
                ),
            )
            kmers = [value for value, _ in pairs]
            assert generate_minimizers(seq, k, w) == [
                min(kmers[i : i + w]) for i in range(len(kmers) - w + 1)
            ]
            record(
                "minimizers",
                size,
                pattern,
                lambda k=k, seq=seq, w=w: generate_minimizers(seq, k, w),
            )
            assert kmer_frequency(seq, k) == dict(Counter(kmers))
            record(
                "kmer_frequency",
                size,
                pattern,
                lambda k=k, seq=seq: kmer_frequency(seq, k),
            )
            record(
                "kmer_list_reference",
                size,
                pattern,
                lambda k=k, seq=seq, size=size: dict(
                    Counter([seq[i : i + k] for i in range(size - k + 1)])
                ),
            )
            record(
                "kmer_generator_candidate",
                size,
                pattern,
                lambda k=k, seq=seq, size=size: dict(
                    Counter(seq[i : i + k] for i in range(size - k + 1))
                ),
            )
            window, step = 100, 10
            expected_gc = [
                (
                    i,
                    i + window,
                    round(
                        (
                            seq[i : i + window].count("G")
                            + seq[i : i + window].count("C")
                        )
                        / window
                        * 100,
                        2,
                    ),
                )
                for i in range(0, size - window + 1, step)
            ]
            assert gc_content_windows(seq, window, step) == expected_gc
            record(
                "gc_windows",
                size,
                pattern,
                lambda seq=seq, window=window, step=step: gc_content_windows(
                    seq, window, step
                ),
            )
            fasta = "".join(f">r{i}\n{seq[i : i + 300]}\n" for i in range(0, size, 300))
            assert "".join(value for _, value in parse_fasta(fasta)) == seq
            record(
                "fasta_parse",
                size,
                pattern,
                lambda fasta=fasta: list(parse_fasta(fasta)),
            )
            coding = seq[: size // 3 * 3]
            assert translate_dna_fast(coding) == translate_dna_to_protein(coding)
            record(
                "translation_python",
                len(coding),
                pattern,
                lambda coding=coding: translate_dna_to_protein(coding),
            )
            record(
                "translation_numba_warm",
                len(coding),
                pattern,
                lambda coding=coding: translate_dna_fast(coding),
            )
        # Alignment has quadratic cost; retain deliberately uneven lengths.
        n = {0: 30, 1: 100}.get(args.sizes.index(size), 300)
        a, b = "ACGT" * (n // 4), "AGGT" * (n // 8)
        for name, aligner in [
            ("global_alignment", needleman_wunsch),
            ("local_alignment", smith_waterman),
            ("edit_distance", levenshtein_distance),
        ]:
            record(
                name,
                len(a) + len(b),
                "uneven",
                lambda a=a, aligner=aligner, b=b: aligner(a, b),
            )

    startup = measure(
        lambda: subprocess.run(
            [sys.executable, "-c", "import bioutils_collection"],
            check=True,
            capture_output=True,
        ),
        args.repeats,
    )
    # New process includes import, Numba cache loading or compilation, and call.
    jit_first = measure(
        lambda: subprocess.run(
            [
                sys.executable,
                "-c",
                "from bioutils_collection.translation_functions import translate_dna_fast; translate_dna_fast('ATGGCC'*100)",
            ],
            check=True,
            capture_output=True,
        ),
        args.repeats,
    )
    seqs = ["ATG" + "GCC" * i + "TAA" for i in range(12)]
    assert translate_batch(seqs, n_processes=1) == translate_batch(seqs, n_processes=2)
    record(
        "batch_serial",
        sum(map(len, seqs)),
        "uneven",
        lambda: translate_batch(seqs, n_processes=1),
    )
    record(
        "batch_process_startup",
        sum(map(len, seqs)),
        "uneven",
        lambda: translate_batch(seqs, n_processes=2),
    )
    with tempfile.TemporaryDirectory() as directory:
        source, output = Path(directory) / "in.fa", Path(directory) / "out.fa"
        for pattern in ("many_records", "one_record"):
            for size in (100000, 1000000):
                seq = "ATG" * (size // 3)
                source.write_text(
                    "".join(
                        f">r{i}\n{seq[i : i + 300]}\n" for i in range(0, len(seq), 300)
                    )
                    if pattern == "many_records"
                    else f">one\n{seq}\n",
                    encoding="utf-8",
                )
                record(
                    "fasta_stream_translate",
                    len(seq),
                    pattern,
                    lambda output=output, source=source: translate_fasta_streaming(
                        source, output
                    ),
                )
    payload = {
        "seed": 20261030,
        "baseline_commit": "9e1cc06420f14940f9e640af554539be29886fb5",
        "python": sys.version,
        "platform": platform.platform(),
        "dependencies": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "scipy", "scikit-learn", "numba"]
        },
        "notes": "Times exclude tracemalloc. Peak Python allocations exclude native arrays, worker RSS and preallocated inputs. Numba first-process timings may use existing disk cache; they are not guaranteed cold compilation. References use corrected selection semantics. Alignment fixtures are capped due to quadratic cost.",
        "import_fresh_process": startup,
        "numba_first_process": jit_first,
        "numba_empty_cache_processes": [
            cold_numba_process() for _ in range(args.repeats)
        ],
        "measurements": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Saved {len(rows)} measurements to {args.output}")


if __name__ == "__main__":
    main()
