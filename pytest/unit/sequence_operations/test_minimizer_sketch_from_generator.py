"""Unit tests for minimizer_sketch_from_generator function."""

import pytest

from bioutils_collection import minimizer_sketch_from_generator


def test_minimizer_sketch_from_generator_basic_case() -> None:
    """Test case 1: Basic sketch generation from generator."""
    kmers = [("ATG", 0), ("TGC", 1), ("GCG", 2), ("CGA", 3)]
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    assert isinstance(result, dict)
    assert all(isinstance(k, str) for k in result.keys())
    assert all(isinstance(v, list) for v in result.values())


def test_minimizer_sketch_from_generator_positions() -> None:
    """Test case 2: Sketch contains correct positions."""
    kmers = [("ATG", 0), ("TGC", 1), ("GCG", 2)]
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    for positions in result.values():
        assert all(isinstance(p, int) for p in positions)
        assert all(p >= 0 for p in positions)


def test_minimizer_sketch_from_generator_simple_input() -> None:
    """Test case 3: Simple k-mer list."""
    kmers = [("AA", 0), ("AC", 1), ("CC", 2)]
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    assert isinstance(result, dict)


def test_minimizer_sketch_from_generator_single_window() -> None:
    """Test case 4: Window size equals k-mer count."""
    kmers = [("ATG", 0), ("TGC", 1)]
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    assert len(result) >= 1


def test_minimizer_sketch_from_generator_window_larger_than_input() -> None:
    """Test case 5: Window larger than number of k-mers."""
    kmers = [("ATG", 0), ("TGC", 1)]
    result = minimizer_sketch_from_generator(iter(kmers), 10)
    assert isinstance(result, dict)
    assert len(result) >= 1


def test_minimizer_sketch_from_generator_empty_input() -> None:
    """Test case 6: Empty generator."""
    kmers = []
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    assert isinstance(result, dict)
    assert len(result) == 0


def test_minimizer_sketch_from_generator_single_kmer() -> None:
    """Test case 7: Single k-mer input."""
    kmers = [("ATG", 0)]
    result = minimizer_sketch_from_generator(iter(kmers), 1)
    assert isinstance(result, dict)
    assert "ATG" in result


def test_minimizer_sketch_from_generator_returns_dict() -> None:
    """Test case 8: Verify function returns a dict."""
    kmers = [("ATG", 0), ("TGC", 1), ("GCG", 2)]
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    assert isinstance(result, dict)


def test_minimizer_sketch_from_generator_all_keys_strings() -> None:
    """Test case 9: All dict keys are strings."""
    kmers = [("ATG", 0), ("TGC", 1), ("GCG", 2)]
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    assert all(isinstance(k, str) for k in result.keys())


def test_minimizer_sketch_from_generator_all_values_lists() -> None:
    """Test case 10: All dict values are lists."""
    kmers = [("ATG", 0), ("TGC", 1), ("GCG", 2)]
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    assert all(isinstance(v, list) for v in result.values())


def test_minimizer_sketch_from_generator_all_positions_integers() -> None:
    """Test case 11: All positions are integers."""
    kmers = [("ATG", 0), ("TGC", 1), ("GCG", 2)]
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    for positions in result.values():
        assert all(isinstance(p, int) for p in positions)


def test_minimizer_sketch_from_generator_selects_minimums() -> None:
    """Test case 12: Selects minimum k-mers."""
    kmers = [("CCC", 0), ("AAA", 1), ("GGG", 2)]
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    assert "AAA" in result  # Should select AAA as minimum


def test_minimizer_sketch_from_generator_repeated_kmers() -> None:
    """Test case 13: Handle repeated k-mers."""
    kmers = [("ATG", 0), ("ATG", 1), ("ATG", 2)]
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    assert "ATG" in result


def test_minimizer_sketch_from_generator_different_positions() -> None:
    """Test case 14: Same k-mer at different positions."""
    kmers = [("ATG", 0), ("GGG", 5), ("ATG", 10)]
    result = minimizer_sketch_from_generator(iter(kmers), 2)
    if "ATG" in result:
        assert isinstance(result["ATG"], list)


def test_minimizer_sketch_from_generator_streaming_use() -> None:
    """Test case 15: Simulate streaming usage."""
    def kmer_generator():
        for i, kmer in enumerate(["ATG", "TGC", "GCG", "CGA"]):
            yield (kmer, i)
    
    result = minimizer_sketch_from_generator(kmer_generator(), 2)
    assert isinstance(result, dict)


def test_minimizer_sketch_from_generator_large_window() -> None:
    """Test case 16: Large window size."""
    kmers = [(f"K{i}", i) for i in range(100)]
    result = minimizer_sketch_from_generator(iter(kmers), 10)
    assert isinstance(result, dict)


def test_minimizer_sketch_from_generator_consistency() -> None:
    """Test case 17: Multiple calls with same data produce same result."""
    kmers_data = [("ATG", 0), ("TGC", 1), ("GCG", 2)]
    result1 = minimizer_sketch_from_generator(iter(kmers_data), 2)
    result2 = minimizer_sketch_from_generator(iter(kmers_data), 2)
    assert result1 == result2


def test_minimizer_sketch_from_generator_type_error_w_not_int() -> None:
    """Test case 18: TypeError when w is not an integer."""
    kmers = [("ATG", 0), ("TGC", 1)]
    with pytest.raises(TypeError, match="w must be int, got str"):
        minimizer_sketch_from_generator(iter(kmers), "2")


def test_minimizer_sketch_from_generator_value_error_w_negative() -> None:
    """Test case 19: ValueError when w is negative."""
    kmers = [("ATG", 0), ("TGC", 1)]
    with pytest.raises(ValueError, match="w must be positive"):
        minimizer_sketch_from_generator(iter(kmers), -1)


def test_minimizer_sketch_from_generator_value_error_w_zero() -> None:
    """Test case 20: ValueError when w is zero."""
    kmers = [("ATG", 0), ("TGC", 1)]
    with pytest.raises(ValueError, match="w must be positive"):
        minimizer_sketch_from_generator(iter(kmers), 0)
