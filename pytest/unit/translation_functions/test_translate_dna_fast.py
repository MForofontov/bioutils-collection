"""Unit tests for translate_dna_fast function."""

try:
    from bioutils_collection.translation_functions import (
        translate_dna_fast,
        translate_dna_to_protein,
        NUMBA_AVAILABLE,
    )
    TRANSLATE_FAST_AVAILABLE = True
except ImportError:
    TRANSLATE_FAST_AVAILABLE = False
    translate_dna_fast = None  # type: ignore
    translate_dna_to_protein = None  # type: ignore
    NUMBA_AVAILABLE = False

import pytest

pytestmark = pytest.mark.skipif(
    not TRANSLATE_FAST_AVAILABLE, reason="translate_dna_fast not available"
)
pytestmark = [pytestmark, pytest.mark.unit, pytest.mark.translation]


def test_translate_dna_fast_basic() -> None:
    """
    Test case 1: Basic DNA to protein translation.
    """
    # Arrange
    seq = "ATGGCC"
    expected = "MA"

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert result == expected


def test_translate_dna_fast_start_codon() -> None:
    """
    Test case 2: Start codon ATG translates to M.
    """
    # Arrange
    seq = "ATG"
    expected = "M"

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert result == expected


def test_translate_dna_fast_stop_codon() -> None:
    """
    Test case 3: Stop codon translates to asterisk.
    """
    # Arrange
    seq = "ATGTAA"
    expected = "M*"

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert result == expected


def test_translate_dna_fast_lowercase() -> None:
    """
    Test case 4: Lowercase input sequence.
    """
    # Arrange
    seq = "atggcc"
    expected = "MA"

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert result == expected


def test_translate_dna_fast_mixed_case() -> None:
    """
    Test case 5: Mixed case input.
    """
    # Arrange
    seq = "AtGgCc"
    expected = "MA"

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert result == expected


def test_translate_dna_fast_empty_sequence() -> None:
    """
    Test case 6: Empty sequence returns empty string.
    """
    # Arrange
    seq = ""
    expected = ""

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert result == expected


def test_translate_dna_fast_longer_sequence() -> None:
    """
    Test case 7: Translation of longer sequence.
    """
    # Arrange
    seq = "ATGGCCAAATTTTGA"
    expected = "MAKF*"

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert result == expected


def test_translate_dna_fast_all_codons() -> None:
    """
    Test case 8: Sequence with various codons.
    """
    # Arrange
    seq = "GCTTGCGGCGGTTATTCATTAACACCCGGCGTAGCACTTCTTGG"  # 44 bases, need 45 for multiple of 3
    seq = seq + "A"  # Now 45 bases
    
    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert len(result) == 15
    assert isinstance(result, str)


def test_translate_dna_fast_multiple_stops() -> None:
    """
    Test case 9: Multiple stop codons.
    """
    # Arrange
    seq = "ATGTAATAGTGA"  # M, TAA(stop), TAG(stop), TGA(stop)
    expected = "M***"

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert result == expected


def test_translate_dna_fast_consistency_default() -> None:
    """
    Test case 10: Results match translate_dna_to_protein.
    """
    # Arrange
    sequences = [
        "ATGGCC",
        "ATGAAATTCTGA",
        "GCTTGCGGCGGTTATTCA",  # Fixed: 18 bases (multiple of 3)
    ]

    # Act & Assert
    for seq in sequences:
        result_fast = translate_dna_fast(seq)
        result_default = translate_dna_to_protein(seq)
        assert result_fast == result_default


def test_translate_dna_fast_vertebrate_mitochondrial() -> None:
    """
    Test case 11: Vertebrate mitochondrial genetic code (table 2).
    """
    # Arrange
    seq = "ATGAGA"  # AGA is stop in vertebrate mito

    # Act
    result = translate_dna_fast(seq, table=2)
    expected = translate_dna_to_protein(seq, table=2)

    # Assert
    assert result == expected


def test_translate_dna_fast_bacterial_table() -> None:
    """
    Test case 12: Bacterial genetic code (table 11).
    """
    # Arrange
    seq = "ATGGCCAAATTT"

    # Act
    result = translate_dna_fast(seq, table=11)
    expected = translate_dna_to_protein(seq, table=11)

    # Assert
    assert result == expected


def test_translate_dna_fast_custom_table() -> None:
    """
    Test case 13: Custom codon table.
    """
    # Arrange
    from bioutils_collection.translation_functions import STANDARD_CODE
    custom = STANDARD_CODE.copy()
    custom['ATG'] = 'X'
    seq = "ATGGCC"
    expected = "XA"

    # Act
    result = translate_dna_fast(seq, table=custom)

    # Assert
    assert result == expected


def test_translate_dna_fast_large_sequence() -> None:
    """
    Test case 14: Large sequence translation.
    """
    # Arrange
    seq = "ATGGCCAAA" * 10000 + "TGA"  # ~90KB

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert len(result) == 30001
    assert result.endswith("*")


def test_translate_dna_fast_all_bases() -> None:
    """
    Test case 15: Sequence with all four bases.
    """
    # Arrange
    seq = "ATGCGATCGATG"  # Fixed: 12 bases
    expected_len = 4

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert len(result) == expected_len
    assert result[0] == "M"


def test_translate_dna_fast_repeating_pattern() -> None:
    """
    Test case 16: Repeating codon pattern.
    """
    # Arrange
    seq = "GCCGCCGCCGCC"  # All Alanine
    expected = "AAAA"

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert result == expected


def test_translate_dna_fast_various_tables() -> None:
    """
    Test case 17: Various genetic code tables.
    """
    # Arrange
    seq = "ATGGCCAAA"
    tables = [1, 2, 3, 4, 5, 11]

    # Act & Assert
    for table in tables:
        result_fast = translate_dna_fast(seq, table=table)
        result_default = translate_dna_to_protein(seq, table=table)
        assert result_fast == result_default


def test_translate_dna_fast_invalid_length() -> None:
    """
    Test case 18: ValueError when sequence length not multiple of 3.
    """
    # Arrange
    seq = "ATGC"  # 4 bases

    # Act & Assert
    with pytest.raises(ValueError, match="Sequence length must be a multiple of 3"):
        translate_dna_fast(seq)


def test_translate_dna_fast_invalid_length_longer() -> None:
    """
    Test case 19: ValueError with longer invalid sequence.
    """
    # Arrange
    seq = "ATGGCCTA"  # 8 bases

    # Act & Assert
    with pytest.raises(ValueError, match="Sequence length must be a multiple of 3"):
        translate_dna_fast(seq)


def test_translate_dna_fast_invalid_bases() -> None:
    """
    Test case 20: Invalid bases produce 'X' (unknown amino acid).
    """
    # Arrange
    seq = "ATGNNN"  # N is invalid base

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert result == "MX"


def test_translate_dna_fast_mixed_invalid_bases() -> None:
    """
    Test case 21: Mixed valid and invalid bases.
    """
    # Arrange
    seq = "ATGXYZGCC"  # XYZ are invalid

    # Act
    result = translate_dna_fast(seq)

    # Assert
    assert result == "MXA"


def test_translate_dna_fast_invalid_table_name() -> None:
    """
    Test case 22: ValueError for invalid table name.
    """
    # Arrange
    seq = "ATGGCC"

    # Act & Assert
    with pytest.raises((ValueError, KeyError)):
        translate_dna_fast(seq, table="nonexistent_table")


def test_translate_dna_fast_non_string_sequence() -> None:
    """
    Test case 23: TypeError for non-string sequence.
    """
    # Arrange
    seq = 12345

    # Act & Assert
    with pytest.raises(TypeError, match="dna_sequence must be str"):
        translate_dna_fast(seq)
