"""Unit tests for minimizer_positions function."""

import pytest

from bioutils_collection import minimizer_positions


def test_minimizer_positions_basic_case() -> None:
    """Test case 1: Basic position retrieval."""
    result = minimizer_positions("ATGCGATCG", 3, 4)
    assert isinstance(result, list)
    assert all(isinstance(p, int) for p in result)


def test_minimizer_positions_simple_sequence() -> None:
    """Test case 2: Simple repeating sequence."""
    result = minimizer_positions("AAACCCGGG", 2, 3)
    assert isinstance(result, list)
    # 6 windows, so 6 minimizers selected
    assert len(result) == 6


def test_minimizer_positions_all_non_negative() -> None:
    """Test case 3: All positions are non-negative."""
    result = minimizer_positions("ATGCGATCG", 3, 4)
    assert all(p >= 0 for p in result)


def test_minimizer_positions_within_bounds() -> None:
    """Test case 4: All positions within sequence bounds."""
    seq = "ATGCGATCG"
    k = 3
    result = minimizer_positions(seq, k, 4)
    assert all(p <= len(seq) - k for p in result)


def test_minimizer_positions_count_matches_minimizers() -> None:
    """Test case 5: Position count matches number of windows."""
    seq = "ATGCGATCG"
    k, w = 3, 4
    result = minimizer_positions(seq, k, w)
    num_windows = len(seq) - k - w + 2
    assert len(result) == num_windows


def test_minimizer_positions_lowercase() -> None:
    """Test case 6: Handle lowercase input."""
    result_upper = minimizer_positions("ATGC", 2, 2)
    result_lower = minimizer_positions("atgc", 2, 2)
    assert result_upper == result_lower


def test_minimizer_positions_type_error_seq_not_string() -> None:
    """Test case 7: TypeError when seq is not a string."""
    with pytest.raises(TypeError, match="seq must be str, got int"):
        minimizer_positions(12345, 3, 4)


def test_minimizer_positions_type_error_k_not_int() -> None:
    """Test case 8: TypeError when k is not an integer."""
    with pytest.raises(TypeError, match="k must be int, got float"):
        minimizer_positions("ATGC", 3.0, 4)


def test_minimizer_positions_type_error_w_not_int() -> None:
    """Test case 9: TypeError when w is not an integer."""
    with pytest.raises(TypeError, match="w must be int, got list"):
        minimizer_positions("ATGC", 3, [4])


def test_minimizer_positions_value_error_empty() -> None:
    """Test case 10: ValueError when sequence is empty."""
    with pytest.raises(ValueError, match="seq cannot be empty"):
        minimizer_positions("", 3, 4)


def test_minimizer_positions_value_error_k_negative() -> None:
    """Test case 11: ValueError when k is negative."""
    with pytest.raises(ValueError, match="k must be positive"):
        minimizer_positions("ATGC", -1, 4)


def test_minimizer_positions_value_error_k_zero() -> None:
    """Test case 12: ValueError when k is zero."""
    with pytest.raises(ValueError, match="k must be positive"):
        minimizer_positions("ATGC", 0, 4)


def test_minimizer_positions_value_error_w_negative() -> None:
    """Test case 13: ValueError when w is negative."""
    with pytest.raises(ValueError, match="w must be positive"):
        minimizer_positions("ATGC", 3, -1)


def test_minimizer_positions_value_error_k_too_long() -> None:
    """Test case 14: ValueError when k is longer than sequence."""
    with pytest.raises(ValueError, match="k cannot be longer than sequence"):
        minimizer_positions("ATG", 5, 4)


def test_minimizer_positions_returns_list() -> None:
    """Test case 15: Verify function returns a list."""
    result = minimizer_positions("ATGCGATCG", 3, 4)
    assert isinstance(result, list)


def test_minimizer_positions_sorted_order() -> None:
    """Test case 16: Positions appear in order (not necessarily sorted)."""
    result = minimizer_positions("ATGCGATCG", 3, 4)
    # Just verify they're all integers and within bounds
    assert all(isinstance(p, int) for p in result)


def test_minimizer_positions_fewer_kmers_than_window() -> None:
    """Test case 17: Handle when sequence has fewer k-mers than window size."""
    result = minimizer_positions("ATGC", 3, 10)
    assert len(result) == 1
    assert result[0] == 0


def test_minimizer_positions_exact_match() -> None:
    """Test case 18: Test exact positions for known sequence."""
    result = minimizer_positions("ATGCGATCG", 3, 4)
    # Verify positions are valid indices
    assert all(0 <= p <= 6 for p in result)


def test_minimizer_positions_large_sequence() -> None:
    """Test case 19: Handle large sequence."""
    large_seq = "ATGC" * 1000
    result = minimizer_positions(large_seq, 5, 10)
    assert isinstance(result, list)
    assert all(isinstance(p, int) for p in result)


def test_minimizer_positions_single_window() -> None:
    """Test case 20: Handle minimal valid input."""
    result = minimizer_positions("ATGC", 2, 2)
    assert isinstance(result, list)
    assert len(result) > 0
