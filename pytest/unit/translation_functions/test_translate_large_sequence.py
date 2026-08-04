"""Unit tests for translate_large_sequence function."""

from bioutils_collection.translation_functions import (
    translate_large_sequence,
    translate_dna_to_protein,
    translate_dna_fast,
)

import pytest

pytestmark = [pytest.mark.unit, pytest.mark.translation]


def test_translate_large_sequence_small() -> None:
    """
    Test case 1: Small sequence (smaller than chunk_size).
    """
    # Arrange
    seq = "ATGGCC"
    expected = "MA"

    # Act
    result = translate_large_sequence(seq, chunk_size=1000)

    # Assert
    assert result == expected


def test_translate_large_sequence_empty() -> None:
    """
    Test case 2: Empty sequence.
    """
    # Arrange
    seq = ""
    expected = ""

    # Act
    result = translate_large_sequence(seq)

    # Assert
    assert result == expected


def test_translate_large_sequence_single_chunk() -> None:
    """
    Test case 3: Sequence fits in single chunk.
    """
    # Arrange
    seq = "ATG" * 100  # 300 bases
    expected = "M" * 100

    # Act
    result = translate_large_sequence(seq, chunk_size=1000)

    # Assert
    assert result == expected


def test_translate_large_sequence_multiple_chunks() -> None:
    """
    Test case 4: Sequence split into multiple chunks.
    """
    # Arrange
    seq = "ATGGCC" * 1000  # 6000 bases
    chunk_size = 2000

    # Act
    result = translate_large_sequence(seq, chunk_size=chunk_size)
    expected = translate_dna_to_protein(seq)

    # Assert
    assert result == expected


def test_translate_large_sequence_chunk_boundary() -> None:
    """
    Test case 5: Chunk size adjusted to multiple of 3.
    """
    # Arrange
    seq = "ATGGCC" * 500  # 3000 bases
    chunk_size = 1001  # Not multiple of 3

    # Act
    result = translate_large_sequence(seq, chunk_size=chunk_size)
    expected = translate_dna_to_protein(seq)

    # Assert
    assert result == expected


def test_translate_large_sequence_with_stops() -> None:
    """
    Test case 6: Large sequence with stop codons.
    """
    # Arrange
    seq = "ATG" + "GCC" * 500 + "TAA" + "GCC" * 500 + "TGA"

    # Act
    result = translate_large_sequence(seq, chunk_size=2000)
    expected = translate_dna_to_protein(seq)

    # Assert
    assert result == expected


def test_translate_large_sequence_consistency() -> None:
    """
    Test case 7: Results match single-threaded translation.
    """
    # Arrange
    seq = "ATGGCCAAA" * 10000  # ~90KB

    # Act
    result_chunked = translate_large_sequence(seq, chunk_size=10000, n_processes=2)
    result_direct = translate_dna_fast(seq)

    # Assert
    assert result_chunked == result_direct


def test_translate_large_sequence_single_process() -> None:
    """
    Test case 8: With n_processes=1.
    """
    # Arrange
    seq = "ATGGCC" * 2000  # 12KB

    # Act
    result = translate_large_sequence(seq, chunk_size=5000, n_processes=1)
    expected = translate_dna_to_protein(seq)

    # Assert
    assert result == expected


def test_translate_large_sequence_custom_table() -> None:
    """
    Test case 9: Custom genetic code table.
    """
    # Arrange
    from bioutils_collection.translation_functions import STANDARD_CODE
    custom = STANDARD_CODE.copy()
    custom['ATG'] = 'X'
    seq = "ATGGCC" * 1000

    # Act
    result = translate_large_sequence(seq, table=custom, chunk_size=2000)

    # Assert
    assert result.startswith("X")
    assert "M" not in result[:100]  # ATG should be X, not M


def test_translate_large_sequence_vertebrate_mito() -> None:
    """
    Test case 10: Vertebrate mitochondrial genetic code.
    """
    # Arrange
    seq = "ATGAGA" * 500

    # Act
    result = translate_large_sequence(seq, table=2, chunk_size=2000)
    expected = translate_dna_to_protein(seq, table=2)

    # Assert
    assert result == expected


def test_translate_large_sequence_use_fast_false() -> None:
    """
    Test case 11: With use_fast=False.
    """
    # Arrange
    seq = "ATGGCC" * 1000

    # Act
    result = translate_large_sequence(seq, chunk_size=2000, use_fast=False)
    expected = translate_dna_to_protein(seq)

    # Assert
    assert result == expected


def test_translate_large_sequence_use_fast_true() -> None:
    """
    Test case 12: With use_fast=True (Numba is always available).
    """
    # Arrange
    seq = "ATGGCC" * 1000

    # Act
    result = translate_large_sequence(seq, chunk_size=2000, use_fast=True)
    expected = translate_dna_fast(seq)

    # Assert
    assert result == expected


def test_translate_large_sequence_very_large() -> None:
    """
    Test case 13: Very large sequence (10MB).
    """
    # Arrange
    seq = "ATGGCCAAA" * 1_111_111  # ~10MB

    # Act
    result = translate_large_sequence(seq, chunk_size=1_000_000, n_processes=2)

    # Assert
    assert len(result) == 3_333_333
    assert result[0] == "M"


def test_translate_large_sequence_exact_chunk_fit() -> None:
    """
    Test case 14: Sequence length is exact multiple of chunk_size.
    """
    # Arrange
    seq = "ATG" * 1000  # 3000 bases
    chunk_size = 1500  # Adjusted to 1500 (multiple of 3)

    # Act
    result = translate_large_sequence(seq, chunk_size=chunk_size)
    expected = "M" * 1000

    # Assert
    assert result == expected


def test_translate_large_sequence_lowercase() -> None:
    """
    Test case 15: Lowercase input sequence.
    """
    # Arrange
    seq = "atggcc" * 1000

    # Act
    result = translate_large_sequence(seq, chunk_size=2000)
    expected = translate_dna_to_protein(seq.upper())

    # Assert
    assert result == expected


def test_translate_large_sequence_bacterial_code() -> None:
    """
    Test case 16: Bacterial genetic code (table 11).
    """
    # Arrange
    seq = "ATGGCCAAA" * 1000

    # Act
    result = translate_large_sequence(seq, table=11, chunk_size=2000)
    expected = translate_dna_to_protein(seq, table=11)

    # Assert
    assert result == expected


def test_translate_large_sequence_invalid_length() -> None:
    """
    Test case 17: ValueError for length not multiple of 3.
    """
    # Arrange
    seq = "ATGC"  # 4 bases, not multiple of 3

    # Act & Assert
    with pytest.raises(ValueError, match="Sequence length must be a multiple of 3"):
        translate_large_sequence(seq)


def test_translate_large_sequence_invalid_chunk_size_one() -> None:
    """
    Test case 18a: ValueError for chunk_size 1 (rounds to 0).
    """
    seq = "ATG" * 100
    with pytest.raises(ValueError, match="chunk_size must be at least 3"):
        translate_large_sequence(seq, chunk_size=1)


def test_translate_large_sequence_invalid_chunk_size_two() -> None:
    """
    Test case 18b: ValueError for chunk_size 2 (rounds to 0).
    """
    seq = "ATG" * 100
    with pytest.raises(ValueError, match="chunk_size must be at least 3"):
        translate_large_sequence(seq, chunk_size=2)


def test_translate_large_sequence_invalid_chunk_size_zero() -> None:
    """
    Test case 18: ValueError for zero chunk_size.
    """
    # Arrange
    seq = "ATG" * 1000

    # Act & Assert
    with pytest.raises((ValueError, ZeroDivisionError)):
        translate_large_sequence(seq, chunk_size=0)


def test_translate_large_sequence_invalid_table() -> None:
    """
    Test case 19: ValueError for invalid table parameter.
    """
    # Arrange
    seq = "ATG" * 1000

    # Act & Assert
    with pytest.raises((ValueError, KeyError)):
        translate_large_sequence(seq, table="nonexistent")


def test_translate_large_sequence_non_string() -> None:
    """
    Test case 20: TypeError for non-string sequence.
    """
    # Arrange
    seq = 12345

    # Act & Assert
    with pytest.raises(TypeError, match="sequence must be str"):
        translate_large_sequence(seq)  # type: ignore[arg-type]


def test_translate_large_sequence_invalid_chunk_size_type() -> None:
    """
    Test case 21: TypeError for non-integer chunk_size.
    """
    # Arrange
    seq = "ATG" * 1000

    # Act & Assert
    with pytest.raises(TypeError, match="chunk_size must be int"):
        translate_large_sequence(seq, chunk_size="1000")  # type: ignore[arg-type]


def test_translate_large_sequence_invalid_chunk_size_negative() -> None:
    """
    Test case 22: ValueError for negative chunk_size.
    """
    # Arrange
    seq = "ATG" * 1000

    # Act & Assert
    with pytest.raises(ValueError, match="chunk_size must be positive"):
        translate_large_sequence(seq, chunk_size=-1000)


def test_translate_large_sequence_invalid_n_processes_negative() -> None:
    """
    Test case 23: ValueError for negative n_processes.
    """
    # Arrange
    seq = "ATG" * 1000

    # Act & Assert
    with pytest.raises(ValueError, match="n_processes must be positive"):
        translate_large_sequence(seq, n_processes=-1)
