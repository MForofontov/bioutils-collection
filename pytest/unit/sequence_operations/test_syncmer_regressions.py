"""Independent selection oracle, including ties and internal minima."""

import itertools

import pytest
from bioutils_collection import generate_syncmers, syncmer_density, syncmer_positions


@pytest.mark.parametrize("method", ["open", "closed"])
def test_syncmers_exhaustive(method: str) -> None:
    for bases in itertools.product("AC", repeat=7):
        seq = "".join(bases)
        for k, s in [(3, 1), (4, 2), (5, 1), (5, 3)]:
            expected = []
            for i in range(len(seq) - k + 1):
                smers = [seq[j : j + s] for j in range(i, i + k - s + 1)]
                offset = smers.index(min(smers))
                if offset == 0 or (method == "closed" and offset == k - s):
                    expected.append(i)
            assert syncmer_positions(seq.lower(), k, s, method) == expected
            assert generate_syncmers(seq, k, s, method) == [
                seq[i : i + k] for i in expected
            ]
            assert syncmer_density(seq, k, s, method) == len(expected) / (
                len(seq) - k + 1
            )


def test_closed_syncmer_internal_minimum_and_tie() -> None:
    assert generate_syncmers("TATT", 4, 1, "closed") == []
    assert syncmer_positions("TAA", 3, 1, "closed") == []
    assert syncmer_positions("TTA", 3, 1, "closed") == [0]
