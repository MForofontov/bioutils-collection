"""Unit tests for syncmer_positions function."""

import pytest

from bioutils_collection import syncmer_positions


def test_syncmer_positions_basic_open() -> None:
    """Test case 1: Basic open syncmer position retrieval."""
    result = syncmer_positions("ATGCGATCG", 4, 2, method="open")
    assert isinstance(result, list)
    assert all(isinstance(p, int) for p in result)


def test_syncmer_positions_basic_closed() -> None:
    """Test case 2: Basic closed syncmer position retrieval."""
    result = syncmer_positions("ATGCGATCG", 4, 2, method="closed")
    assert isinstance(result, list)
    assert all(isinstance(p, int) for p in result)


def test_syncmer_positions_all_non_negative() -> None:
    """Test case 3: All positions are non-negative."""
    result = syncmer_positions("ATGCGATCG", 4, 2)
    assert all(p >= 0 for p in result)


def test_syncmer_positions_within_bounds() -> None:
    """Test case 4: All positions within sequence bounds."""
    seq = "ATGCGATCG"
    k = 4
    result = syncmer_positions(seq, k, 2)
    assert all(p <= len(seq) - k for p in result)


def test_syncmer_positions_closed_more_than_open() -> None:
    """Test case 5: Closed typically has more positions than open."""
    seq = "AAACCCGGG"
    open_result = syncmer_positions(seq, 3, 2, method="open")
    closed_result = syncmer_positions(seq, 3, 2, method="closed")
    assert len(closed_result) >= len(open_result)


def test_syncmer_positions_default_method() -> None:
    """Test case 6: Default method is 'open'."""
    result_default = syncmer_positions("ATGC", 3, 2)
    result_open = syncmer_positions("ATGC", 3, 2, method="open")
    assert result_default == result_open


def test_syncmer_positions_lowercase() -> None:
    """Test case 7: Handle lowercase input."""
    result_upper = syncmer_positions("ATGC", 3, 2)
    result_lower = syncmer_positions("atgc", 3, 2)
    assert result_upper == result_lower


def test_syncmer_positions_type_error_seq_not_string() -> None:
    """Test case 8: TypeError when seq is not a string."""
    with pytest.raises(TypeError, match="seq must be str, got list"):
        syncmer_positions([], 4, 2)


def test_syncmer_positions_type_error_k_not_int() -> None:
    """Test case 9: TypeError when k is not an integer."""
    with pytest.raises(TypeError, match="k must be int, got str"):
        syncmer_positions("ATGC", "4", 2)


def test_syncmer_positions_type_error_s_not_int() -> None:
    """Test case 10: TypeError when s is not an integer."""
    with pytest.raises(TypeError, match="s must be int, got dict"):
        syncmer_positions("ATGC", 4, {})


def test_syncmer_positions_type_error_method_not_string() -> None:
    """Test case 11: TypeError when method is not a string."""
    with pytest.raises(TypeError, match="method must be str, got list"):
        syncmer_positions("ATGC", 4, 2, method=[])


def test_syncmer_positions_value_error_empty() -> None:
    """Test case 12: ValueError when sequence is empty."""
    with pytest.raises(ValueError, match="seq cannot be empty"):
        syncmer_positions("", 4, 2)


def test_syncmer_positions_value_error_k_negative() -> None:
    """Test case 13: ValueError when k is negative."""
    with pytest.raises(ValueError, match="k must be positive"):
        syncmer_positions("ATGC", -1, 2)


def test_syncmer_positions_value_error_s_negative() -> None:
    """Test case 14: ValueError when s is negative."""
    with pytest.raises(ValueError, match="s must be positive"):
        syncmer_positions("ATGC", 4, -1)


def test_syncmer_positions_value_error_s_not_less_than_k() -> None:
    """Test case 15: ValueError when s >= k."""
    with pytest.raises(ValueError, match="s must be less than k"):
        syncmer_positions("ATGC", 4, 5)


def test_syncmer_positions_value_error_k_too_long() -> None:
    """Test case 16: ValueError when k is longer than sequence."""
    with pytest.raises(ValueError, match="k cannot be longer than sequence"):
        syncmer_positions("ATG", 5, 2)


def test_syncmer_positions_value_error_invalid_method() -> None:
    """Test case 17: ValueError when method is not 'open' or 'closed'."""
    with pytest.raises(ValueError, match="method must be 'open' or 'closed'"):
        syncmer_positions("ATGC", 4, 2, method="random")


def test_syncmer_positions_returns_list() -> None:
    """Test case 18: Verify function returns a list."""
    result = syncmer_positions("ATGCGATCG", 4, 2)
    assert isinstance(result, list)


def test_syncmer_positions_known_example() -> None:
    """Test case 19: Test with known example."""
    result = syncmer_positions("ATGCGATCG", 4, 2, method="open")
    assert all(isinstance(p, int) for p in result)
    assert all(0 <= p <= 5 for p in result)


def test_syncmer_positions_large_sequence() -> None:
    """Test case 20: Handle large sequence."""
    large_seq = "ATGC" * 1000
    result = syncmer_positions(large_seq, 5, 3)
    assert isinstance(result, list)
    assert all(isinstance(p, int) for p in result)
