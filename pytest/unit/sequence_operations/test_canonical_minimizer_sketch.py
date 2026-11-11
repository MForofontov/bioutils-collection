"""Unit tests for canonical_minimizer_sketch function."""

import pytest

from bioutils_collection import canonical_minimizer_sketch


def test_canonical_minimizer_sketch_basic_case() -> None:
    """Test case 1: Basic sketch generation."""
    result = canonical_minimizer_sketch("ATGCGATCG", 3, 4)
    assert isinstance(result, dict)
    assert all(isinstance(k, str) for k in result.keys())
    assert all(isinstance(v, list) for v in result.values())


def test_canonical_minimizer_sketch_has_positions() -> None:
    """Test case 2: Sketch contains positions."""
    result = canonical_minimizer_sketch("ATGCGATCG", 3, 4)
    for positions in result.values():
        assert all(isinstance(p, int) for p in positions)
        assert all(p >= 0 for p in positions)


def test_canonical_minimizer_sketch_simple_sequence() -> None:
    """Test case 3: Simple repeating sequence."""
    result = canonical_minimizer_sketch("AAACCCGGG", 2, 3)
    assert isinstance(result, dict)
    assert "AA" in result or "TT" in result  # Canonical form


def test_canonical_minimizer_sketch_all_positions_valid() -> None:
    """Test case 4: All positions are within sequence bounds."""
    seq = "ATGCGATCG"
    result = canonical_minimizer_sketch(seq, 3, 4)
    for positions in result.values():
        assert all(0 <= p < len(seq) for p in positions)


def test_canonical_minimizer_sketch_lowercase() -> None:
    """Test case 5: Handle lowercase input."""
    result_upper = canonical_minimizer_sketch("ATGC", 2, 2)
    result_lower = canonical_minimizer_sketch("atgc", 2, 2)
    assert result_upper == result_lower


def test_canonical_minimizer_sketch_returns_dict() -> None:
    """Test case 6: Verify function returns a dict."""
    result = canonical_minimizer_sketch("ATGCGATCG", 3, 4)
    assert isinstance(result, dict)


def test_canonical_minimizer_sketch_all_keys_strings() -> None:
    """Test case 7: All dict keys are strings."""
    result = canonical_minimizer_sketch("ATGCGATCG", 3, 4)
    assert all(isinstance(k, str) for k in result.keys())


def test_canonical_minimizer_sketch_all_values_lists() -> None:
    """Test case 8: All dict values are lists."""
    result = canonical_minimizer_sketch("ATGCGATCG", 3, 4)
    assert all(isinstance(v, list) for v in result.values())


def test_canonical_minimizer_sketch_all_positions_integers() -> None:
    """Test case 9: All positions are integers."""
    result = canonical_minimizer_sketch("ATGCGATCG", 3, 4)
    for positions in result.values():
        assert all(isinstance(p, int) for p in positions)


def test_canonical_minimizer_sketch_fewer_kmers_than_window() -> None:
    """Test case 10: Handle when sequence has fewer k-mers than window size."""
    result = canonical_minimizer_sketch("ATGC", 3, 10)
    assert len(result) == 1


def test_canonical_minimizer_sketch_non_empty_lists() -> None:
    """Test case 11: All position lists are non-empty."""
    result = canonical_minimizer_sketch("ATGCGATCG", 3, 4)
    assert all(len(v) > 0 for v in result.values())


def test_canonical_minimizer_sketch_large_sequence() -> None:
    """Test case 12: Handle large sequence."""
    large_seq = "ATGC" * 1000
    result = canonical_minimizer_sketch(large_seq, 5, 10)
    assert isinstance(result, dict)


def test_canonical_minimizer_sketch_all_minimizers_uppercase() -> None:
    """Test case 13: All minimizers are uppercase."""
    result = canonical_minimizer_sketch("atgc", 2, 2)
    assert all(k.isupper() for k in result.keys())


def test_canonical_minimizer_sketch_type_error_seq_not_string() -> None:
    """Test case 14: TypeError when seq is not a string."""
    with pytest.raises(TypeError, match="seq must be str, got int"):
        canonical_minimizer_sketch(12345, 3, 4)


def test_canonical_minimizer_sketch_type_error_k_not_int() -> None:
    """Test case 15: TypeError when k is not an integer."""
    with pytest.raises(TypeError, match="k must be int, got str"):
        canonical_minimizer_sketch("ATGC", "3", 4)


def test_canonical_minimizer_sketch_type_error_w_not_int() -> None:
    """Test case 16: TypeError when w is not an integer."""
    with pytest.raises(TypeError, match="w must be int, got list"):
        canonical_minimizer_sketch("ATGC", 3, [4])


def test_canonical_minimizer_sketch_value_error_empty() -> None:
    """Test case 17: ValueError when sequence is empty."""
    with pytest.raises(ValueError, match="seq cannot be empty"):
        canonical_minimizer_sketch("", 3, 4)


def test_canonical_minimizer_sketch_value_error_k_negative() -> None:
    """Test case 18: ValueError when k is negative."""
    with pytest.raises(ValueError, match="k must be positive"):
        canonical_minimizer_sketch("ATGC", -1, 4)


def test_canonical_minimizer_sketch_value_error_w_negative() -> None:
    """Test case 19: ValueError when w is negative."""
    with pytest.raises(ValueError, match="w must be positive"):
        canonical_minimizer_sketch("ATGC", 3, -1)


def test_canonical_minimizer_sketch_value_error_k_too_long() -> None:
    """Test case 20: ValueError when k is longer than sequence."""
    with pytest.raises(ValueError, match="k cannot be longer than sequence"):
        canonical_minimizer_sketch("ATG", 5, 4)
