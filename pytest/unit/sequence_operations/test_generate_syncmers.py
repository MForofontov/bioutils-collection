"""Unit tests for generate_syncmers function."""

import pytest

from bioutils_collection import generate_syncmers


def test_generate_syncmers_basic_open() -> None:
    """Test case 1: Basic open syncmer generation."""
    result = generate_syncmers("ATGCGATCG", 4, 2, method="open")
    assert isinstance(result, list)
    assert all(len(s) == 4 for s in result)


def test_generate_syncmers_basic_closed() -> None:
    """Test case 2: Basic closed syncmer generation."""
    result = generate_syncmers("ATGCGATCG", 4, 2, method="closed")
    assert isinstance(result, list)
    assert all(len(s) == 4 for s in result)


def test_generate_syncmers_closed_more_than_open() -> None:
    """Test case 3: Closed syncmers typically produce more results than open."""
    seq = "AAACCCGGG"
    open_result = generate_syncmers(seq, 3, 2, method="open")
    closed_result = generate_syncmers(seq, 3, 2, method="closed")
    assert len(closed_result) >= len(open_result)


def test_generate_syncmers_simple_sequence() -> None:
    """Test case 4: Simple repeating sequence."""
    result = generate_syncmers("AAACCCGGG", 3, 2, method="closed")
    assert isinstance(result, list)
    assert len(result) > 0


def test_generate_syncmers_default_method() -> None:
    """Test case 5: Default method is 'open'."""
    result_default = generate_syncmers("ATGC", 3, 2)
    result_open = generate_syncmers("ATGC", 3, 2, method="open")
    assert result_default == result_open


def test_generate_syncmers_lowercase() -> None:
    """Test case 6: Handle lowercase input."""
    result_upper = generate_syncmers("ATGC", 3, 2)
    result_lower = generate_syncmers("atgc", 3, 2)
    assert result_upper == result_lower


def test_generate_syncmers_type_error_seq_not_string() -> None:
    """Test case 7: TypeError when seq is not a string."""
    with pytest.raises(TypeError, match="seq must be str, got int"):
        generate_syncmers(12345, 4, 2)


def test_generate_syncmers_type_error_k_not_int() -> None:
    """Test case 8: TypeError when k is not an integer."""
    with pytest.raises(TypeError, match="k must be int, got str"):
        generate_syncmers("ATGC", "4", 2)


def test_generate_syncmers_type_error_s_not_int() -> None:
    """Test case 9: TypeError when s is not an integer."""
    with pytest.raises(TypeError, match="s must be int, got float"):
        generate_syncmers("ATGC", 4, 2.0)


def test_generate_syncmers_type_error_method_not_string() -> None:
    """Test case 10: TypeError when method is not a string."""
    with pytest.raises(TypeError, match="method must be str, got int"):
        generate_syncmers("ATGC", 4, 2, method=1)


def test_generate_syncmers_value_error_empty() -> None:
    """Test case 11: ValueError when sequence is empty."""
    with pytest.raises(ValueError, match="seq cannot be empty"):
        generate_syncmers("", 4, 2)


def test_generate_syncmers_value_error_k_negative() -> None:
    """Test case 12: ValueError when k is negative."""
    with pytest.raises(ValueError, match="k must be positive"):
        generate_syncmers("ATGC", -1, 2)


def test_generate_syncmers_value_error_s_negative() -> None:
    """Test case 13: ValueError when s is negative."""
    with pytest.raises(ValueError, match="s must be positive"):
        generate_syncmers("ATGC", 4, -1)


def test_generate_syncmers_value_error_s_not_less_than_k() -> None:
    """Test case 14: ValueError when s >= k."""
    with pytest.raises(ValueError, match="s must be less than k"):
        generate_syncmers("ATGC", 4, 4)


def test_generate_syncmers_value_error_k_too_long() -> None:
    """Test case 15: ValueError when k is longer than sequence."""
    with pytest.raises(ValueError, match="k cannot be longer than sequence"):
        generate_syncmers("ATG", 5, 2)


def test_generate_syncmers_value_error_invalid_method() -> None:
    """Test case 16: ValueError when method is not 'open' or 'closed'."""
    with pytest.raises(ValueError, match="method must be 'open' or 'closed'"):
        generate_syncmers("ATGC", 4, 2, method="invalid")


def test_generate_syncmers_returns_list() -> None:
    """Test case 17: Verify function returns a list."""
    result = generate_syncmers("ATGCGATCG", 4, 2)
    assert isinstance(result, list)


def test_generate_syncmers_all_strings() -> None:
    """Test case 18: All syncmers are strings."""
    result = generate_syncmers("ATGCGATCG", 4, 2)
    assert all(isinstance(s, str) for s in result)


def test_generate_syncmers_correct_length() -> None:
    """Test case 19: All syncmers have correct length."""
    k = 5
    result = generate_syncmers("ATGCGATCGAA", k, 2)
    assert all(len(s) == k for s in result)


def test_generate_syncmers_large_sequence() -> None:
    """Test case 20: Handle large sequence."""
    large_seq = "ATGC" * 1000
    result = generate_syncmers(large_seq, 5, 3)
    assert isinstance(result, list)
    assert all(isinstance(s, str) for s in result)
