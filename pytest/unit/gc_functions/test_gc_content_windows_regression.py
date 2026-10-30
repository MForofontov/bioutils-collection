"""Compare window boundaries and GC percentages with a direct implementation."""

import random

import pytest
from bioutils_collection.gc_functions.gc_content_windows import gc_content_windows

pytestmark = [pytest.mark.unit, pytest.mark.gc_functions]


@pytest.mark.parametrize("window_size", [1, 2, 7, 25, 101])
@pytest.mark.parametrize("step_size", [1, 2, 7, 25, 101, 200])
def test_gc_windows_match_reference(window_size: int, step_size: int) -> None:
    rng = random.Random(42)
    sequence = "".join(rng.choices("ATGCNatgcn", k=101))
    expected = []
    for start in range(0, len(sequence) - window_size + 1, step_size):
        window = sequence[start : start + window_size].upper()
        count = window.count("G") + window.count("C")
        expected.append(
            (start, start + window_size, round(count / window_size * 100, 2))
        )

    assert gc_content_windows(sequence, window_size, step_size) == expected


def test_gc_windows_skip_incomplete_tail() -> None:
    assert gc_content_windows("GGCCATN", 4, 3) == [
        (0, 4, 100.0),
        (3, 7, 25.0),
    ]


def test_gc_windows_reject_invalid_base_outside_windows() -> None:
    with pytest.raises(ValueError, match="Invalid DNA bases"):
        gc_content_windows("ATGCX", 4, 10)
