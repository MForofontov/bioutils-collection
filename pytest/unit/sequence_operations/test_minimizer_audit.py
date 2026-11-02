import itertools

import pytest
from bioutils_collection import minimizer_sketch, minimizer_sketch_from_generator


@pytest.mark.parametrize("w", [1, 2, 3, 10])
def test_streaming_sketch_uses_only_full_windows_or_one_short_window(w: int) -> None:
    for values in itertools.product("AC", repeat=4):
        items = [(value, 5 * i) for i, value in enumerate(values)]
        selected = [
            min(items[i : i + w], key=lambda pair: pair[0])
            for i in range(max(1, len(items) - w + 1))
        ]
        expected = {}
        for value, position in selected:
            expected.setdefault(value, []).append(position)
        assert minimizer_sketch_from_generator(iter(items), w) == expected


@pytest.mark.parametrize("w", [1, 3, 100])
def test_zero_density_is_empty(w: int) -> None:
    assert minimizer_sketch("ACGTACGT", 3, w, density=0.0) == {}
