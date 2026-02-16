"""Unit tests for generate_minimizers function."""

import pytest

from bioutils_collection import generate_minimizers


def test_generate_minimizers_basic_case() -> None:
    """Test case 1: Basic minimizer generation."""
    result = generate_minimizers("ATGCGATCG", 3, 4)
    assert isinstance(result, list)
    assert len(result) == 4
    assert all(len(m) == 3 for m in result)


def test_generate_minimizers_simple_sequence() -> None:
    """Test case 2: Simple repeating sequence."""
    result = generate_minimizers("AAACCCGGG", 2, 3)
    # Windows: [AA,AA,AC], [AA,AC,CC], [AC,CC,CC], [CC,CC,CG], [CC,CG,GG], [CG,GG,GG]
    assert len(result) == 6  # Number of windows


def test_generate_minimizers_identical_sequence() -> None:
    """Test case 3: All same bases."""
    result = generate_minimizers("AAAAA", 2, 2)
    # Windows: [AA,AA], [AA,AA], [AA,AA] = 3 windows
    assert result == ["AA", "AA", "AA"]


def test_generate_minimizers_small_window() -> None:
    """Test case 4: Window size of 2."""
    result = generate_minimizers("ATGC", 2, 2)
    # Windows: [AT,TG], [TG,GC] = 2 windows, minimums are AT and GC
    assert len(result) == 2
    assert result == ["AT", "GC"]


def test_generate_minimizers_exact_window_size() -> None:
    """Test case 5: Sequence length equals window requirement."""
    result = generate_minimizers("ATGCG", 2, 3)
    assert len(result) == 2


def test_generate_minimizers_lowercase_input() -> None:
    """Test case 6: Handle lowercase input."""
    result_upper = generate_minimizers("ATGC", 2, 2)
    result_lower = generate_minimizers("atgc", 2, 2)
    assert result_upper == result_lower


def test_generate_minimizers_type_error_not_string() -> None:
    """Test case 7: TypeError when seq is not a string."""
    with pytest.raises(TypeError, match="seq must be str, got int"):
        generate_minimizers(12345, 3, 4)


def test_generate_minimizers_type_error_k_not_int() -> None:
    """Test case 8: TypeError when k is not an integer."""
    with pytest.raises(TypeError, match="k must be int, got str"):
        generate_minimizers("ATGC", "3", 4)


def test_generate_minimizers_type_error_w_not_int() -> None:
    """Test case 9: TypeError when w is not an integer."""
    with pytest.raises(TypeError, match="w must be int, got float"):
        generate_minimizers("ATGC", 3, 4.0)


def test_generate_minimizers_value_error_empty() -> None:
    """Test case 10: ValueError when sequence is empty."""
    with pytest.raises(ValueError, match="seq cannot be empty"):
        generate_minimizers("", 3, 4)


def test_generate_minimizers_value_error_k_negative() -> None:
    """Test case 11: ValueError when k is negative."""
    with pytest.raises(ValueError, match="k must be positive"):
        generate_minimizers("ATGC", -1, 4)


def test_generate_minimizers_value_error_k_zero() -> None:
    """Test case 12: ValueError when k is zero."""
    with pytest.raises(ValueError, match="k must be positive"):
        generate_minimizers("ATGC", 0, 4)


def test_generate_minimizers_value_error_w_negative() -> None:
    """Test case 13: ValueError when w is negative."""
    with pytest.raises(ValueError, match="w must be positive"):
        generate_minimizers("ATGC", 3, -1)


def test_generate_minimizers_value_error_w_zero() -> None:
    """Test case 14: ValueError when w is zero."""
    with pytest.raises(ValueError, match="w must be positive"):
        generate_minimizers("ATGC", 3, 0)


def test_generate_minimizers_value_error_k_too_long() -> None:
    """Test case 15: ValueError when k is longer than sequence."""
    with pytest.raises(ValueError, match="k cannot be longer than sequence"):
        generate_minimizers("ATG", 5, 4)


def test_generate_minimizers_returns_list() -> None:
    """Test case 16: Verify function returns a list."""
    result = generate_minimizers("ATGCGATCG", 3, 4)
    assert isinstance(result, list)


def test_generate_minimizers_all_strings() -> None:
    """Test case 17: All minimizers are strings."""
    result = generate_minimizers("ATGCGATCG", 3, 4)
    assert all(isinstance(m, str) for m in result)


def test_generate_minimizers_correct_kmer_length() -> None:
    """Test case 18: All k-mers have correct length."""
    k = 4
    result = generate_minimizers("ATGCGATCGAA", k, 3)
    assert all(len(m) == k for m in result)


def test_generate_minimizers_fewer_kmers_than_window() -> None:
    """Test case 19: Handle when sequence has fewer k-mers than window size."""
    result = generate_minimizers("ATGC", 3, 10)
    assert len(result) == 1
    assert result[0] == "ATG"


def test_generate_minimizers_large_sequence() -> None:
    """Test case 20: Handle large sequence."""
    large_seq = "ATGC" * 1000
    result = generate_minimizers(large_seq, 5, 10)
    assert isinstance(result, list)
    assert all(isinstance(m, str) for m in result)
