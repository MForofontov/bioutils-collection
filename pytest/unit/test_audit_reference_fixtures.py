"""Cross-module biological and mathematical fixtures for the major audit."""

from collections import Counter
from functools import cache
from itertools import product

import pytest
from bioutils_collection import (
    canonical_minimizer_sketch,
    generate_minimizers,
    motif_clustering,
    motif_search,
    needleman_wunsch,
    palindromic_sequence_finder,
    restriction_site_finder,
    reverse_complement,
    smith_waterman,
    tandem_repeat_finder,
    validate_dna_sequence,
)
from bioutils_collection.sequence_statistics.codon_usage_frequency import (
    codon_usage_frequency,
)
from bioutils_collection.sequence_statistics.effective_number_of_codons import (
    effective_number_of_codons,
)
from bioutils_collection.sequence_statistics.kmer_frequency import kmer_frequency
from bioutils_collection.sequence_statistics.relative_synonymous_codon_usage import (
    relative_synonymous_codon_usage,
)
from bioutils_collection.sequence_statistics.sequence_conservation import (
    sequence_conservation,
)
from bioutils_collection.sequence_statistics.sequence_entropy import sequence_entropy
from bioutils_collection.translation_functions import (
    translate_batch,
    translate_dna_to_protein,
    translate_large_sequence,
)


@cache
def path_scores(a: str, b: str) -> frozenset[int]:
    """Enumerate scores of every possible global alignment of tiny strings."""
    if not a:
        return frozenset([-len(b)])
    if not b:
        return frozenset([-len(a)])
    return frozenset(
        [value + (2 if a[0] == b[0] else -2) for value in path_scores(a[1:], b[1:])]
        + [value - 1 for value in path_scores(a[1:], b)]
        + [value - 1 for value in path_scores(a, b[1:])]
    )


def test_alignment_optima_and_tracebacks() -> None:
    strings = ["".join(chars) for n in range(1, 4) for chars in product("AC", repeat=n)]
    for a, b in product(strings, repeat=2):
        score, left, right = needleman_wunsch(a, b, match=2, mismatch=-2)
        assert score == max(path_scores(a, b))
        assert left.replace("-", "") == a
        assert right.replace("-", "") == b
        assert score == sum(
            -1 if "-" in (x, y) else (2 if x == y else -2)
            for x, y in zip(left, right, strict=True)
        )
        expected_local = max(
            max(path_scores(a[i:j], b[k:end]))
            for i in range(len(a) + 1)
            for j in range(i, len(a) + 1)
            for k in range(len(b) + 1)
            for end in range(k, len(b) + 1)
        )
        score, left, right = smith_waterman(a, b, match=2, mismatch=-2)
        assert score == expected_local
        assert left.replace("-", "") in a and right.replace("-", "") in b
        assert score == sum(
            -1 if "-" in (x, y) else (2 if x == y else -2)
            for x, y in zip(left, right, strict=True)
        )


def test_minimizers_and_canonical_sketch_against_direct_windows() -> None:
    seq = "ACGTAACGTACCC"
    for k, w in product([1, 3, 5], [1, 3, 30]):
        kmers = [seq[i : i + k] for i in range(len(seq) - k + 1)]
        expected = [min(kmers[i : i + w]) for i in range(max(1, len(kmers) - w + 1))]
        assert generate_minimizers(seq, k, w) == expected
        canonical = [
            (min(value, reverse_complement(value)), i) for i, value in enumerate(kmers)
        ]
        sketch = {}
        for start in range(max(1, len(canonical) - w + 1)):
            value, position = min(
                canonical[start : start + w], key=lambda pair: pair[0]
            )
            sketch.setdefault(value, []).append(position)
        assert canonical_minimizer_sketch(seq, k, w) == sketch


def test_codon_statistics_and_entropy() -> None:
    assert relative_synonymous_codon_usage("TTTTTTTTCTTT")["TTT"] == 1.5
    assert codon_usage_frequency("ATGATGTTT") == {"ATG": 2 / 3, "TTT": 1 / 3}
    all_codons = ["".join(bases) for bases in product("ACGT", repeat=3)]
    assert effective_number_of_codons("".join(all_codons) * 100) == 61.0
    assert effective_number_of_codons("TTT" * 100) == 20.0
    assert sequence_entropy("ACGT" * 20) == 2.0
    assert sequence_conservation(["ATGC", "ATGT", "ATGA"]) == [1.0, 1.0, 1.0, 1 / 3]
    for k in [1, 3, 6]:
        seq = "AaCGTACGT"
        assert kmer_frequency(seq, k) == dict(
            Counter(seq[i : i + k] for i in range(len(seq) - k + 1))
        )


def test_motifs_repeats_restriction_and_validation() -> None:
    assert motif_search("AAGAGG", "ARG") == [0, 3]
    assert restriction_site_finder("aaaa", ["AA"]) == {"AA": [0, 1, 2]}
    assert palindromic_sequence_finder("GAATTC", min_length=6) == [(0, 6, "GAATTC")]
    assert tandem_repeat_finder("ATATAT") == [(0, 6, "AT"), (1, 5, "TA"), (2, 6, "AT")]
    assert validate_dna_sequence("acgtn", allow_ambiguous=True)["gc_content"] == 50.0
    assert not validate_dna_sequence("ACGTN")["is_valid"]


def test_unambiguous_clustering_partition() -> None:
    labels = motif_clustering(["AAAA", "AAAA", "TTTT", "TTTT"], 2)
    assert labels[0] == labels[1] and labels[2] == labels[3] and labels[0] != labels[2]


@pytest.mark.parametrize("fast", [False, True])
def test_serial_parallel_translation_ordering(fast: bool) -> None:
    seqs = ["ATG" + "GCC" * i + "TAA" for i in range(12)]
    table = {"ATG": "start", "GCC": "A", "TAA": "*"}
    expected = ["start" + "A" * i + "*" for i in range(12)]
    assert translate_batch(seqs, table, n_processes=1, use_fast=fast) == expected
    assert translate_batch(seqs, table, n_processes=2, use_fast=fast) == expected
    sequence = "ATGGCC" * 25
    assert translate_large_sequence(
        sequence, table, chunk_size=13, n_processes=2, use_fast=fast
    ) == translate_dna_to_protein(sequence, table)
