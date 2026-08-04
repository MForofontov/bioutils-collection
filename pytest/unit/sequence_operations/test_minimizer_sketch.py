"""Unit tests for minimizer_sketch function."""

import pytest

from bioutils_collection import minimizer_sketch


def test_minimizer_sketch_basic_case() -> None:
    """Test case 1: Basic sketch generation."""
    result = minimizer_sketch("ATGCGATCG", 3, 4)
    assert isinstance(result, dict)
    assert all(isinstance(k, str) for k in result.keys())
    assert all(isinstance(v, list) for v in result.values())


def test_minimizer_sketch_has_positions() -> None:
    """Test case 2: Sketch contains positions."""
    result = minimizer_sketch("ATGCGATCG", 3, 4)
    for positions in result.values():
        assert all(isinstance(p, int) for p in positions)
        assert all(p >= 0 for p in positions)


def test_minimizer_sketch_simple_sequence() -> None:
    """Test case 3: Simple repeating sequence."""
    result = minimizer_sketch("AAACCCGGG", 2, 3)
    assert "AA" in result
    assert "CC" in result
    # CG is selected from the last windows, not GG
    assert isinstance(result, dict)


def test_minimizer_sketch_with_density() -> None:
    """Test case 4: Test with density parameter."""
    result = minimizer_sketch("ATGCGATCG", 3, 4, density=0.5)
    assert isinstance(result, dict)


def test_minimizer_sketch_full_density() -> None:
    """Test case 5: Test with full density (1.0)."""
    result = minimizer_sketch("ATGCGATCG", 3, 4, density=1.0)
    assert isinstance(result, dict)


def test_minimizer_sketch_low_density() -> None:
    """Test case 6: Test with low density."""
    result = minimizer_sketch("ATGCGATCG", 3, 4, density=0.3)
    full_result = minimizer_sketch("ATGCGATCG", 3, 4, density=1.0)
    # Low density should have fewer or equal minimizers
    assert len(result) <= len(full_result)


def test_minimizer_sketch_lowercase() -> None:
    """Test case 7: Handle lowercase input."""
    result_upper = minimizer_sketch("ATGC", 2, 2)
    result_lower = minimizer_sketch("atgc", 2, 2)
    assert result_upper == result_lower


def test_minimizer_sketch_returns_dict() -> None:
    """Test case 8: Verify function returns a dict."""
    result = minimizer_sketch("ATGCGATCG", 3, 4)
    assert isinstance(result, dict)


def test_minimizer_sketch_all_positions_valid() -> None:
    """Test case 9: All positions are within sequence bounds."""
    seq = "ATGCGATCG"
    result = minimizer_sketch(seq, 3, 4)
    for positions in result.values():
        assert all(0 <= p < len(seq) for p in positions)


def test_minimizer_sketch_fewer_kmers_than_window() -> None:
    """Test case 10: Handle when sequence has fewer k-mers than window size."""
    result = minimizer_sketch("ATGC", 3, 10)
    assert len(result) == 1
    assert "ATG" in result


def test_minimizer_sketch_type_error_seq_not_string() -> None:
    """Test case 11: TypeError when seq is not a string."""
    with pytest.raises(TypeError, match="seq must be str, got int"):
        minimizer_sketch(12345, 3, 4)


def test_minimizer_sketch_type_error_k_not_int() -> None:
    """Test case 12: TypeError when k is not an integer."""
    with pytest.raises(TypeError, match="k must be int, got str"):
        minimizer_sketch("ATGC", "3", 4)


def test_minimizer_sketch_type_error_w_not_int() -> None:
    """Test case 13: TypeError when w is not an integer."""
    with pytest.raises(TypeError, match="w must be int, got str"):
        minimizer_sketch("ATGC", 3, "4")


def test_minimizer_sketch_type_error_density_not_numeric() -> None:
    """Test case 14: TypeError when density is not numeric."""
    with pytest.raises(TypeError, match="density must be float, got str"):
        minimizer_sketch("ATGC", 3, 4, density="0.5")


def test_minimizer_sketch_value_error_empty() -> None:
    """Test case 15: ValueError when sequence is empty."""
    with pytest.raises(ValueError, match="seq cannot be empty"):
        minimizer_sketch("", 3, 4)


def test_minimizer_sketch_value_error_k_negative() -> None:
    """Test case 16: ValueError when k is negative."""
    with pytest.raises(ValueError, match="k must be positive"):
        minimizer_sketch("ATGC", -1, 4)


def test_minimizer_sketch_value_error_w_negative() -> None:
    """Test case 17: ValueError when w is negative."""
    with pytest.raises(ValueError, match="w must be positive"):
        minimizer_sketch("ATGC", 3, -1)


def test_minimizer_sketch_value_error_k_too_long() -> None:
    """Test case 18: ValueError when k is longer than sequence."""
    with pytest.raises(ValueError, match="k cannot be longer than sequence"):
        minimizer_sketch("ATG", 5, 4)


def test_minimizer_sketch_value_error_density_too_high() -> None:
    """Test case 19: ValueError when density > 1.0."""
    with pytest.raises(ValueError, match="density must be between 0.0 and 1.0"):
        minimizer_sketch("ATGC", 2, 2, density=1.5)


def test_minimizer_sketch_value_error_density_negative() -> None:
    """Test case 20: ValueError when density is negative."""
    with pytest.raises(ValueError, match="density must be between 0.0 and 1.0"):
        minimizer_sketch("ATGC", 2, 2, density=-0.1)
