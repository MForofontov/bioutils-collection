"""Independent small-input expectations for the second maintenance audit."""

import math
from itertools import product

import pytest
from bioutils_collection import (
    canonical_minimizer_sketch,
    canonical_minimizers,
    generate_minimizers,
    minimizer_positions,
    minimizer_sketch,
    motif_search,
    needleman_wunsch,
    remove_low_complexity_regions,
    sequence_pattern_match,
    smith_waterman,
    unique_minimizer_positions,
)

pytestmark = pytest.mark.unit


def test_masking_against_direct_entropy_and_interval_union() -> None:
    for length in range(7):
        for bases in product("ACN", repeat=length):
            seq = "".join(bases)
            for width, threshold in product((2, 3, 5), (0.0, 1.0, 2.0)):
                covered = set()
                for start in range(length - width + 1):
                    window = seq[start : start + width]
                    probabilities = [window.count(base) / width for base in set(window)]
                    entropy = -math.fsum(p * math.log2(p) for p in probabilities)
                    if entropy < threshold:
                        covered.update(range(start, start + width))
                expected = "".join(
                    "X" if i in covered else base for i, base in enumerate(seq)
                )
                assert (
                    remove_low_complexity_regions(seq, width, threshold, "X")
                    == expected
                )


def test_masking_preserves_case_and_counts_n_as_a_symbol() -> None:
    assert remove_low_complexity_regions("aaaCGT", 3, 0.5) == "NNNCGT"
    assert remove_low_complexity_regions("acgt", 4, 2) == "acgt"
    # Five equally frequent symbols have entropy log2(5) > 2 bits.
    assert remove_low_complexity_regions("aCgTN", 5, 2) == "aCgTN"


@pytest.mark.parametrize("threshold", [float("nan"), float("inf"), -float("inf")])
@pytest.mark.parametrize("seq", ["", "A", "ACGT"])
def test_masking_nonfinite_thresholds(seq: str, threshold: float) -> None:
    with pytest.raises(ValueError, match="complexity_threshold"):
        remove_low_complexity_regions(seq, complexity_threshold=threshold)


def test_iupac_searches_against_allowed_base_products() -> None:
    allowed = {"A": "A", "R": "AG", "Y": "CT", "N": "ACGT"}
    for length in range(2, 5):
        for bases in product("ACGT", repeat=length):
            seq = "".join(bases)
            for motif in ("R", "Y", "N", "AN", "RY"):
                words = {
                    "".join(word) for word in product(*(allowed[c] for c in motif))
                }
                expected = [
                    i
                    for i in range(length - len(motif) + 1)
                    if seq[i : i + len(motif)] in words
                ]
                assert motif_search(seq.lower(), motif.lower()) == expected
                assert (
                    sequence_pattern_match(seq.lower(), motif.lower(), use_iupac=True)
                    == expected
                )


def test_selector_outputs_against_direct_complete_windows() -> None:
    complement = str.maketrans("ACGT", "TGCA")
    for length in range(3, 7):
        for bases in product("AT", repeat=length):
            seq = "".join(bases)
            for k, width in product((1, 2, 3), (1, 3, 8)):
                pairs = [(seq[i : i + k], i) for i in range(length - k + 1)]
                canonical = [
                    (min(word, word.translate(complement)[::-1]), i)
                    for word, i in pairs
                ]
                selected = []
                canonical_selected = []
                for start in range(max(1, len(pairs) - width + 1)):
                    selected.append(
                        min(pairs[start : start + width], key=lambda pair: pair[0])
                    )
                    canonical_selected.append(
                        min(canonical[start : start + width], key=lambda pair: pair[0])
                    )
                sketch: dict[str, list[int]] = {}
                canonical_sketch: dict[str, list[int]] = {}
                for word, pos in selected:
                    sketch.setdefault(word, []).append(pos)
                for word, pos in canonical_selected:
                    canonical_sketch.setdefault(word, []).append(pos)
                assert generate_minimizers(seq.lower(), k, width) == [
                    word for word, _ in selected
                ]
                assert minimizer_positions(seq, k, width) == [
                    pos for _, pos in selected
                ]
                assert minimizer_sketch(seq, k, width) == sketch
                assert unique_minimizer_positions(seq, k, width) == {
                    word: sorted(set(positions)) for word, positions in sketch.items()
                }
                assert canonical_minimizers(seq, k, width) == [
                    word for word, _ in canonical_selected
                ]
                assert canonical_minimizer_sketch(seq, k, width) == canonical_sketch


def test_alignment_ties_have_documented_traceback_order() -> None:
    assert needleman_wunsch("AA", "A") == (0, "AA", "-A")
    assert needleman_wunsch("AC", "CA", 0, 0, 0) == (0, "AC", "CA")
    assert smith_waterman("A", "AA") == (2, "A", "A")
    assert smith_waterman("AC", "CA", 0, 0, 0) == (0, "", "")
