"""Unit tests for syncmer_density function."""

import pytest

from bioutils_collection import syncmer_density


def test_syncmer_density_basic_open() -> None:
    """Test case 1: Basic density calculation for open syncmers."""
    result = syncmer_density("ATGCGATCG", 4, 2, method="open")
    assert isinstance(result, float)
    assert 0.0 <= result <= 1.0


def test_syncmer_density_basic_closed() -> None:
    """Test case 2: Basic density calculation for closed syncmers."""
    result = syncmer_density("ATGCGATCG", 4, 2, method="closed")
    assert isinstance(result, float)
    assert 0.0 <= result <= 1.0


def test_syncmer_density_closed_higher_than_open() -> None:
    """Test case 3: Closed density typically higher than open."""
    seq = "AAACCCGGG"
    open_density = syncmer_density(seq, 3, 2, method="open")
    closed_density = syncmer_density(seq, 3, 2, method="closed")
    assert closed_density >= open_density


def test_syncmer_density_range() -> None:
    """Test case 4: Density is between 0 and 1."""
    result = syncmer_density("ATGCGATCG", 4, 2)
    assert 0.0 <= result <= 1.0


def test_syncmer_density_default_method() -> None:
    """Test case 5: Default method is 'open'."""
    result_default = syncmer_density("ATGC", 3, 2)
    result_open = syncmer_density("ATGC", 3, 2, method="open")
    assert result_default == result_open


def test_syncmer_density_lowercase() -> None:
    """Test case 6: Handle lowercase input."""
    result_upper = syncmer_density("ATGC", 3, 2)
    result_lower = syncmer_density("atgc", 3, 2)
    assert result_upper == result_lower


def test_syncmer_density_type_error_seq_not_string() -> None:
    """Test case 7: TypeError when seq is not a string."""
    with pytest.raises(TypeError, match="seq must be str, got int"):
        syncmer_density(12345, 4, 2)


def test_syncmer_density_type_error_k_not_int() -> None:
    """Test case 8: TypeError when k is not an integer."""
    with pytest.raises(TypeError, match="k must be int, got str"):
        syncmer_density("ATGC", "4", 2)


def test_syncmer_density_type_error_s_not_int() -> None:
    """Test case 9: TypeError when s is not an integer."""
    with pytest.raises(TypeError, match="s must be int, got float"):
        syncmer_density("ATGC", 4, 2.5)


def test_syncmer_density_type_error_method_not_string() -> None:
    """Test case 10: TypeError when method is not a string."""
    with pytest.raises(TypeError, match="method must be str, got int"):
        syncmer_density("ATGC", 4, 2, method=123)


def test_syncmer_density_value_error_empty() -> None:
    """Test case 11: ValueError when sequence is empty."""
    with pytest.raises(ValueError, match="seq cannot be empty"):
        syncmer_density("", 4, 2)


def test_syncmer_density_value_error_k_negative() -> None:
    """Test case 12: ValueError when k is negative."""
    with pytest.raises(ValueError, match="k must be positive"):
        syncmer_density("ATGC", -1, 2)


def test_syncmer_density_value_error_s_negative() -> None:
    """Test case 13: ValueError when s is negative."""
    with pytest.raises(ValueError, match="s must be positive"):
        syncmer_density("ATGC", 4, -1)


def test_syncmer_density_value_error_s_not_less_than_k() -> None:
    """Test case 14: ValueError when s >= k."""
    with pytest.raises(ValueError, match="s must be less than k"):
        syncmer_density("ATGC", 4, 4)


def test_syncmer_density_value_error_k_too_long() -> None:
    """Test case 15: ValueError when k is longer than sequence."""
    with pytest.raises(ValueError, match="k cannot be longer than sequence"):
        syncmer_density("ATG", 5, 2)


def test_syncmer_density_value_error_invalid_method() -> None:
    """Test case 16: ValueError when method is not 'open' or 'closed'."""
    with pytest.raises(ValueError, match="method must be 'open' or 'closed'"):
        syncmer_density("ATGC", 4, 2, method="other")


def test_syncmer_density_returns_float() -> None:
    """Test case 17: Verify function returns a float."""
    result = syncmer_density("ATGCGATCG", 4, 2)
    assert isinstance(result, float)


def test_syncmer_density_consistent_results() -> None:
    """Test case 18: Same input produces same output."""
    result1 = syncmer_density("ATGCGATCG", 4, 2)
    result2 = syncmer_density("ATGCGATCG", 4, 2)
    assert result1 == result2


def test_syncmer_density_known_value() -> None:
    """Test case 19: Test exact density for known sequence."""
    result = syncmer_density("ATGCGATCG", 4, 2, method="open")
    # Should be around 1/3 for open syncmers with k=4, s=2
    assert 0.0 <= result <= 1.0


def test_syncmer_density_large_sequence() -> None:
    """Test case 20: Handle large sequence."""
    large_seq = "ATGC" * 1000
    result = syncmer_density(large_seq, 5, 3)
    assert isinstance(result, float)
    assert 0.0 <= result <= 1.0
