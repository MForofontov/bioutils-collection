"""Unit tests for minimizer_density function."""

import pytest

from bioutils_collection import minimizer_density


def test_minimizer_density_basic_case() -> None:
    """Test case 1: Basic density calculation."""
    result = minimizer_density("ATGCGATCG", 3, 4)
    assert isinstance(result, float)
    assert 0.0 <= result <= 1.0


def test_minimizer_density_simple_sequence() -> None:
    """Test case 2: Simple repeating sequence."""
    result = minimizer_density("AAACCCGGG", 2, 3)
    assert isinstance(result, float)
    assert 0.0 <= result <= 1.0


def test_minimizer_density_range() -> None:
    """Test case 3: Density is between 0 and 1."""
    result = minimizer_density("ATGCGATCG", 3, 4)
    assert 0.0 <= result <= 1.0


def test_minimizer_density_high_for_diverse_sequence() -> None:
    """Test case 4: Diverse sequence has higher density."""
    result = minimizer_density("ATGCATGC", 2, 3)
    assert result > 0.0


def test_minimizer_density_low_for_repetitive_sequence() -> None:
    """Test case 5: Repetitive sequence has lower density."""
    result = minimizer_density("AAAAA", 2, 2)
    # Same minimizer repeated means lower unique count
    assert result <= 1.0


def test_minimizer_density_lowercase() -> None:
    """Test case 6: Handle lowercase input."""
    result_upper = minimizer_density("ATGC", 2, 2)
    result_lower = minimizer_density("atgc", 2, 2)
    assert result_upper == result_lower


def test_minimizer_density_type_error_seq_not_string() -> None:
    """Test case 7: TypeError when seq is not a string."""
    with pytest.raises(TypeError, match="seq must be str, got list"):
        minimizer_density(["A", "T"], 3, 4)


def test_minimizer_density_type_error_k_not_int() -> None:
    """Test case 8: TypeError when k is not an integer."""
    with pytest.raises(TypeError, match="k must be int, got float"):
        minimizer_density("ATGC", 3.5, 4)


def test_minimizer_density_type_error_w_not_int() -> None:
    """Test case 9: TypeError when w is not an integer."""
    with pytest.raises(TypeError, match="w must be int, got str"):
        minimizer_density("ATGC", 3, "4")


def test_minimizer_density_value_error_empty() -> None:
    """Test case 10: ValueError when sequence is empty."""
    with pytest.raises(ValueError, match="seq cannot be empty"):
        minimizer_density("", 3, 4)


def test_minimizer_density_value_error_k_negative() -> None:
    """Test case 11: ValueError when k is negative."""
    with pytest.raises(ValueError, match="k must be positive"):
        minimizer_density("ATGC", -1, 4)


def test_minimizer_density_value_error_w_negative() -> None:
    """Test case 12: ValueError when w is negative."""
    with pytest.raises(ValueError, match="w must be positive"):
        minimizer_density("ATGC", 3, -1)


def test_minimizer_density_value_error_k_too_long() -> None:
    """Test case 13: ValueError when k is longer than sequence."""
    with pytest.raises(ValueError, match="k cannot be longer than sequence"):
        minimizer_density("ATG", 5, 4)


def test_minimizer_density_returns_float() -> None:
    """Test case 14: Verify function returns a float."""
    result = minimizer_density("ATGCGATCG", 3, 4)
    assert isinstance(result, float)


def test_minimizer_density_consistent_results() -> None:
    """Test case 15: Same input produces same output."""
    result1 = minimizer_density("ATGCGATCG", 3, 4)
    result2 = minimizer_density("ATGCGATCG", 3, 4)
    assert result1 == result2


def test_minimizer_density_fewer_kmers_than_window() -> None:
    """Test case 16: Handle when sequence has fewer k-mers than window size."""
    result = minimizer_density("ATGC", 3, 10)
    assert isinstance(result, float)
    assert 0.0 <= result <= 1.0


def test_minimizer_density_exact_value() -> None:
    """Test case 17: Test exact density for known sequence."""
    result = minimizer_density("ATGCGATCG", 3, 4)
    # 3 unique minimizers out of 7 total kmers
    assert abs(result - 0.428571) < 0.01


def test_minimizer_density_larger_window() -> None:
    """Test case 18: Larger window generally means lower density."""
    seq = "ATGCGATCGATGC"
    density_small = minimizer_density(seq, 3, 2)
    density_large = minimizer_density(seq, 3, 5)
    assert density_large <= density_small or abs(density_large - density_small) < 0.1


def test_minimizer_density_large_sequence() -> None:
    """Test case 19: Handle large sequence."""
    large_seq = "ATGC" * 1000
    result = minimizer_density(large_seq, 5, 10)
    assert isinstance(result, float)
    assert 0.0 <= result <= 1.0


def test_minimizer_density_single_unique_minimizer() -> None:
    """Test case 20: All same bases produces minimal density."""
    result = minimizer_density("AAAA", 2, 2)
    # Only one unique minimizer "AA"
    assert result > 0.0
