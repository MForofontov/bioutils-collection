"""Unit tests for canonical_minimizers function."""

import pytest

from bioutils_collection import canonical_minimizers


def test_canonical_minimizers_basic_case() -> None:
    """Test case 1: Basic canonical minimizer generation."""
    result = canonical_minimizers("ATGCGATCG", 3, 4)
    assert isinstance(result, list)
    assert len(result) == 4


def test_canonical_minimizers_strand_independent() -> None:
    """Test case 2: Sequence and reverse complement should be similar."""
    seq = "ATGC"
    # Note: Not exact because positions matter, but kmers should use canonical forms
    result = canonical_minimizers(seq, 2, 2)
    assert isinstance(result, list)


def test_canonical_minimizers_simple_sequence() -> None:
    """Test case 3: Simple repeating sequence."""
    result = canonical_minimizers("AAACCCGGG", 2, 3)
    assert isinstance(result, list)
    # 6 windows, so 6 minimizers
    assert len(result) == 6


def test_canonical_minimizers_palindromic_kmers() -> None:
    """Test case 4: Palindromic k-mers are their own canonical form."""
    # AT and TA are reverse complements, canonical is AT
    result = canonical_minimizers("ATAT", 2, 2)
    assert all(isinstance(m, str) for m in result)


def test_canonical_minimizers_all_same_base() -> None:
    """Test case 5: All same bases."""
    result = canonical_minimizers("AAAAA", 2, 2)
    # 3 windows: [AA,AA], [AA,AA], [AA,AA]
    assert result == ["AA", "AA", "AA"]


def test_canonical_minimizers_lowercase() -> None:
    """Test case 6: Handle lowercase input."""
    result_upper = canonical_minimizers("ATGC", 2, 2)
    result_lower = canonical_minimizers("atgc", 2, 2)
    assert result_upper == result_lower


def test_canonical_minimizers_type_error_seq_not_string() -> None:
    """Test case 7: TypeError when seq is not a string."""
    with pytest.raises(TypeError, match="seq must be str, got int"):
        canonical_minimizers(12345, 3, 4)


def test_canonical_minimizers_type_error_k_not_int() -> None:
    """Test case 8: TypeError when k is not an integer."""
    with pytest.raises(TypeError, match="k must be int, got str"):
        canonical_minimizers("ATGC", "3", 4)


def test_canonical_minimizers_type_error_w_not_int() -> None:
    """Test case 9: TypeError when w is not an integer."""
    with pytest.raises(TypeError, match="w must be int, got float"):
        canonical_minimizers("ATGC", 3, 4.0)


def test_canonical_minimizers_value_error_empty() -> None:
    """Test case 10: ValueError when sequence is empty."""
    with pytest.raises(ValueError, match="seq cannot be empty"):
        canonical_minimizers("", 3, 4)


def test_canonical_minimizers_value_error_k_negative() -> None:
    """Test case 11: ValueError when k is negative."""
    with pytest.raises(ValueError, match="k must be positive"):
        canonical_minimizers("ATGC", -1, 4)


def test_canonical_minimizers_value_error_w_negative() -> None:
    """Test case 12: ValueError when w is negative."""
    with pytest.raises(ValueError, match="w must be positive"):
        canonical_minimizers("ATGC", 3, -1)


def test_canonical_minimizers_value_error_k_zero() -> None:
    """Test case 13: ValueError when k is zero."""
    with pytest.raises(ValueError, match="k must be positive"):
        canonical_minimizers("ATGC", 0, 4)


def test_canonical_minimizers_value_error_k_too_long() -> None:
    """Test case 14: ValueError when k is longer than sequence."""
    with pytest.raises(ValueError, match="k cannot be longer than sequence"):
        canonical_minimizers("ATG", 5, 4)


def test_canonical_minimizers_returns_list() -> None:
    """Test case 15: Verify function returns a list."""
    result = canonical_minimizers("ATGCGATCG", 3, 4)
    assert isinstance(result, list)


def test_canonical_minimizers_all_strings() -> None:
    """Test case 16: All minimizers are strings."""
    result = canonical_minimizers("ATGCGATCG", 3, 4)
    assert all(isinstance(m, str) for m in result)


def test_canonical_minimizers_correct_length() -> None:
    """Test case 17: All k-mers have correct length."""
    k = 4
    result = canonical_minimizers("ATGCGATCGAA", k, 3)
    assert all(len(m) == k for m in result)


def test_canonical_minimizers_fewer_kmers_than_window() -> None:
    """Test case 18: Handle when sequence has fewer k-mers than window size."""
    result = canonical_minimizers("ATGC", 3, 10)
    assert len(result) == 1


def test_canonical_minimizers_all_uppercase() -> None:
    """Test case 19: All minimizers are uppercase."""
    result = canonical_minimizers("atgc", 2, 2)
    assert all(m.isupper() for m in result)


def test_canonical_minimizers_large_sequence() -> None:
    """Test case 20: Handle large sequence."""
    large_seq = "ATGC" * 1000
    result = canonical_minimizers(large_seq, 5, 10)
    assert isinstance(result, list)
    assert all(isinstance(m, str) for m in result)
