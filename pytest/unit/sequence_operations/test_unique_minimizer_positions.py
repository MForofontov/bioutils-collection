"""Unit tests for unique_minimizer_positions function."""

import pytest

from bioutils_collection import unique_minimizer_positions


def test_unique_minimizer_positions_basic_case() -> None:
    """Test case 1: Basic unique positions retrieval."""
    result = unique_minimizer_positions("ATGCGATCG", 3, 4)
    assert isinstance(result, dict)
    assert all(isinstance(k, str) for k in result.keys())
    assert all(isinstance(v, list) for v in result.values())


def test_unique_minimizer_positions_sorted_positions() -> None:
    """Test case 2: Positions are sorted for each minimizer."""
    result = unique_minimizer_positions("AAACCCGGG", 2, 3)
    for positions in result.values():
        assert positions == sorted(positions)


def test_unique_minimizer_positions_no_duplicates() -> None:
    """Test case 3: No duplicate positions for each minimizer."""
    result = unique_minimizer_positions("AAACCCGGG", 2, 3)
    for positions in result.values():
        assert len(positions) == len(set(positions))


def test_unique_minimizer_positions_all_within_bounds() -> None:
    """Test case 4: All positions within sequence bounds."""
    seq = "ATGCGATCG"
    k = 3
    result = unique_minimizer_positions(seq, k, 4)
    for positions in result.values():
        assert all(0 <= p <= len(seq) - k for p in positions)


def test_unique_minimizer_positions_lowercase() -> None:
    """Test case 5: Handle lowercase input."""
    result_upper = unique_minimizer_positions("ATGC", 2, 2)
    result_lower = unique_minimizer_positions("atgc", 2, 2)
    assert result_upper == result_lower


def test_unique_minimizer_positions_returns_dict() -> None:
    """Test case 6: Verify function returns a dict."""
    result = unique_minimizer_positions("ATGCGATCG", 3, 4)
    assert isinstance(result, dict)


def test_unique_minimizer_positions_all_keys_strings() -> None:
    """Test case 7: All dict keys are strings."""
    result = unique_minimizer_positions("ATGCGATCG", 3, 4)
    assert all(isinstance(k, str) for k in result.keys())


def test_unique_minimizer_positions_all_values_lists() -> None:
    """Test case 8: All dict values are lists."""
    result = unique_minimizer_positions("ATGCGATCG", 3, 4)
    assert all(isinstance(v, list) for v in result.values())


def test_unique_minimizer_positions_all_positions_integers() -> None:
    """Test case 9: All positions are integers."""
    result = unique_minimizer_positions("ATGCGATCG", 3, 4)
    for positions in result.values():
        assert all(isinstance(p, int) for p in positions)


def test_unique_minimizer_positions_fewer_kmers_than_window() -> None:
    """Test case 10: Handle when sequence has fewer k-mers than window size."""
    result = unique_minimizer_positions("ATGC", 3, 10)
    assert len(result) == 1
    assert "ATG" in result


def test_unique_minimizer_positions_known_example() -> None:
    """Test case 11: Test with known example."""
    result = unique_minimizer_positions("AAACCCGGG", 2, 3)
    assert "AA" in result
    assert "CC" in result


def test_unique_minimizer_positions_large_sequence() -> None:
    """Test case 12: Handle large sequence."""
    large_seq = "ATGC" * 1000
    result = unique_minimizer_positions(large_seq, 5, 10)
    assert isinstance(result, dict)


def test_unique_minimizer_positions_non_empty_lists() -> None:
    """Test case 13: All position lists are non-empty."""
    result = unique_minimizer_positions("ATGCGATCG", 3, 4)
    assert all(len(v) > 0 for v in result.values())


def test_unique_minimizer_positions_type_error_seq_not_string() -> None:
    """Test case 14: TypeError when seq is not a string."""
    with pytest.raises(TypeError, match="seq must be str, got int"):
        unique_minimizer_positions(12345, 3, 4)


def test_unique_minimizer_positions_type_error_k_not_int() -> None:
    """Test case 15: TypeError when k is not an integer."""
    with pytest.raises(TypeError, match="k must be int, got str"):
        unique_minimizer_positions("ATGC", "3", 4)


def test_unique_minimizer_positions_type_error_w_not_int() -> None:
    """Test case 16: TypeError when w is not an integer."""
    with pytest.raises(TypeError, match="w must be int, got dict"):
        unique_minimizer_positions("ATGC", 3, {})


def test_unique_minimizer_positions_value_error_empty() -> None:
    """Test case 17: ValueError when sequence is empty."""
    with pytest.raises(ValueError, match="seq cannot be empty"):
        unique_minimizer_positions("", 3, 4)


def test_unique_minimizer_positions_value_error_k_negative() -> None:
    """Test case 18: ValueError when k is negative."""
    with pytest.raises(ValueError, match="k must be positive"):
        unique_minimizer_positions("ATGC", -1, 4)


def test_unique_minimizer_positions_value_error_w_negative() -> None:
    """Test case 19: ValueError when w is negative."""
    with pytest.raises(ValueError, match="w must be positive"):
        unique_minimizer_positions("ATGC", 3, -1)


def test_unique_minimizer_positions_value_error_k_too_long() -> None:
    """Test case 20: ValueError when k is longer than sequence."""
    with pytest.raises(ValueError, match="k cannot be longer than sequence"):
        unique_minimizer_positions("ATG", 5, 4)
