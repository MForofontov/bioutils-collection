import pytest

from bioutils_collection.repeat_functions.palindromic_sequence_finder import (
    palindromic_sequence_finder,
)

pytestmark = [pytest.mark.unit, pytest.mark.repeat]


def test_palindromic_sequence_finder_basic() -> None:
    """
    Test case 1: Basic biological palindrome finding (EcoRI site).
    """
    # GAATTC is a biological palindrome (restriction enzyme recognition site)
    sequence = "GAATTC"
    result = palindromic_sequence_finder(sequence, min_length=4)
    assert isinstance(result, list)
    assert len(result) > 0
    # Should find GAATTC as a biological palindrome
    assert any(item[2] == "GAATTC" for item in result)


def test_palindromic_sequence_finder_no_palindromes() -> None:
    """
    Test case 2: No palindromes found.
    """
    sequence = "ATGC"  # Changed to avoid TGCA which is a palindrome
    result = palindromic_sequence_finder(sequence, min_length=4)
    assert result == []


def test_palindromic_sequence_finder_short_palindrome() -> None:
    """
    Test case 3: Find short palindromes with min_length=2.
    """
    sequence = "ATTA"
    result = palindromic_sequence_finder(sequence, min_length=2)
    assert len(result) > 0


def test_palindromic_sequence_finder_empty_sequence() -> None:
    """
    Test case 4: Empty sequence returns empty list.
    """
    result = palindromic_sequence_finder("", min_length=4)
    assert result == []


def test_palindromic_sequence_finder_entire_palindrome() -> None:
    """
    Test case 5: String palindrome (biological=False).
    """
    sequence = "ATTA"
    # ATTA is a string palindrome but not a biological one
    result = palindromic_sequence_finder(sequence, min_length=4, biological=False)
    assert len(result) >= 1  # ATTA is a string palindrome
    
    # With biological=True (default), ATTA is NOT a palindrome
    result_bio = palindromic_sequence_finder(sequence, min_length=4, biological=True)
    # ATTA reverse complement is TAAT, not ATTA
    assert len(result_bio) == 0


def test_palindromic_sequence_finder_biological_vs_string() -> None:
    """
    Test case 6: Difference between biological and string palindromes.
    """
    # GAATTC is both biological and string palindrome
    sequence = "GAATTC"
    result_bio = palindromic_sequence_finder(sequence, min_length=6, biological=True)
    result_str = palindromic_sequence_finder(sequence, min_length=6, biological=False)
    assert len(result_bio) >= 1  # Biological palindrome
    assert len(result_str) == 0  # Not a string palindrome


def test_palindromic_sequence_finder_invalid_min_length_error() -> None:
    """
    Test case 7: ValueError for invalid min_length.
    """
    with pytest.raises(ValueError, match="min_length must be > 0"):
        palindromic_sequence_finder("ATGC", min_length=0)
    with pytest.raises(ValueError, match="min_length must be > 0"):
        palindromic_sequence_finder("ATGC", min_length=-1)


def test_palindromic_sequence_finder_type_errors() -> None:
    """
    Test case 8: TypeError for invalid input types.
    """
    with pytest.raises(TypeError, match="sequence must be str"):
        palindromic_sequence_finder(12345, min_length=4)  # type: ignore[arg-type]
    
    with pytest.raises(TypeError, match="min_length must be int"):
        palindromic_sequence_finder("ATGC", min_length="4")  # type: ignore[arg-type]


def test_palindromic_sequence_finder_invalid_bases() -> None:
    """
    Test case 9: ValueError for invalid DNA bases when biological=True.
    """
    with pytest.raises(ValueError, match="Sequence contains invalid DNA bases"):
        palindromic_sequence_finder("ATGCX", min_length=4, biological=True)
    
    # Should work fine with biological=False
    result = palindromic_sequence_finder("XYZYX", min_length=5, biological=False)
    assert len(result) >= 1  # XYZYX is a string palindrome
