"""Independent boundary and format fixtures for the 0.3.3 maintenance audit."""

import pytest
from bioutils_collection import (
    fasta_stats,
    genbank_to_fasta,
    remove_low_complexity_regions,
    restriction_site_finder,
)

pytestmark = pytest.mark.unit


@pytest.mark.parametrize("newline", ["\n", "\r\n", "\r"])
@pytest.mark.parametrize("separator", [" ", "\t"])
def test_genbank_sequence_whitespace(newline: str, separator: str) -> None:
    record = newline.join(
        ["LOCUS       seq1", "ORIGIN", f"        1 acgt{separator}tgca", "//"]
    )
    assert genbank_to_fasta(record) == ">seq1\nACGTTGCA\n"


def test_fasta_gc_is_case_insensitive_and_record_averaged() -> None:
    # GC fractions are 2/4, 4/4 and 0/0 (defined as zero), not pooled counts.
    assert fasta_stats(">a\naTgC\n>b\nggcc\n>empty\n") == (3, 8 / 3, 0.5)
    # Ambiguous bases remain in the per-record denominator.
    assert fasta_stats(">n\ngcNN\n") == (1, 4.0, 0.5)


def test_restriction_site_cannot_be_empty() -> None:
    with pytest.raises(ValueError, match="site.*empty"):
        restriction_site_finder("ACGT", ["AC", ""])


@pytest.mark.parametrize("site", [None, 3, b"AC"])
def test_restriction_sites_must_be_strings(site: object) -> None:
    with pytest.raises(TypeError, match="site.*string"):
        restriction_site_finder("ACGT", [site])  # type: ignore[list-item]


@pytest.mark.parametrize("sequence", ["", "A", "ACGT"])
def test_masking_rejects_nan_threshold_even_for_short_input(sequence: str) -> None:
    with pytest.raises(ValueError, match="complexity_threshold"):
        remove_low_complexity_regions(sequence, complexity_threshold=float("nan"))


def test_masking_union_and_strict_entropy_boundary() -> None:
    # Only the two AAA windows have entropy < 0.5: their union is [0, 4).
    assert remove_low_complexity_regions("AAAACG", 3, 0.5) == "NNNNCG"
    # A zero-entropy window is equal to a zero threshold, so remains unmasked.
    assert remove_low_complexity_regions("AAAA", 3, 0) == "AAAA"
