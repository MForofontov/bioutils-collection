"""Unit tests for translate_batch function."""

try:
    from bioutils_collection.translation_functions import (
        translate_batch,
        translate_dna_to_protein,
        NUMBA_AVAILABLE,
    )
    TRANSLATE_BATCH_AVAILABLE = True
except ImportError:
    TRANSLATE_BATCH_AVAILABLE = False
    translate_batch = None  # type: ignore
    translate_dna_to_protein = None  # type: ignore
    NUMBA_AVAILABLE = False

import pytest

pytestmark = pytest.mark.skipif(
    not TRANSLATE_BATCH_AVAILABLE, reason="translate_batch not available"
)
pytestmark = [pytestmark, pytest.mark.unit, pytest.mark.translation]


def test_translate_batch_empty() -> None:
    """
    Test case 1: Empty sequence list.
    """
    # Arrange
    sequences = []
    expected = []

    # Act
    result = translate_batch(sequences)

    # Assert
    assert result == expected


def test_translate_batch_single() -> None:
    """
    Test case 2: Single sequence.
    """
    # Arrange
    sequences = ["ATGGCC"]
    expected = ["MA"]

    # Act
    result = translate_batch(sequences)

    # Assert
    assert result == expected


def test_translate_batch_multiple() -> None:
    """
    Test case 3: Multiple sequences.
    """
    # Arrange
    sequences = ["ATGGCC", "ATGAAA", "ATTTTT"]
    expected = ["MA", "MK", "IF"]

    # Act
    result = translate_batch(sequences)

    # Assert
    assert result == expected


def test_translate_batch_different_lengths() -> None:
    """
    Test case 4: Sequences of different lengths.
    """
    # Arrange
    sequences = ["ATG", "ATGGCC", "ATGGCCAAATTT"]
    expected = ["M", "MA", "MAKF"]

    # Act
    result = translate_batch(sequences)

    # Assert
    assert result == expected


def test_translate_batch_with_stops() -> None:
    """
    Test case 5: Sequences with stop codons.
    """
    # Arrange
    sequences = ["ATGTAA", "ATGTAGTGA"]
    expected = ["M*", "M**"]

    # Act
    result = translate_batch(sequences)

    # Assert
    assert result == expected


def test_translate_batch_consistency() -> None:
    """
    Test case 6: Batch results match individual translations.
    """
    # Arrange
    sequences = ["ATGGCC", "ATGAAATTC", "GCTTGCGGC", "TTTGGGCCC"]

    # Act
    batch_result = translate_batch(sequences)
    individual_results = [translate_dna_to_protein(seq) for seq in sequences]

    # Assert
    assert batch_result == individual_results


def test_translate_batch_lowercase() -> None:
    """
    Test case 7: Lowercase sequences.
    """
    # Arrange
    sequences = ["atggcc", "atgaaa"]
    expected = ["MA", "MK"]

    # Act
    result = translate_batch(sequences)

    # Assert
    assert result == expected


def test_translate_batch_use_fast_false() -> None:
    """
    Test case 8: With use_fast=False.
    """
    # Arrange
    sequences = ["ATGGCC", "ATGAAA"]
    expected = ["MA", "MK"]

    # Act
    result = translate_batch(sequences, use_fast=False)

    # Assert
    assert result == expected


def test_translate_batch_single_process() -> None:
    """
    Test case 9: With n_processes=1.
    """
    # Arrange
    sequences = ["ATGGCC", "ATGAAA", "ATTTTT"]
    expected = ["MA", "MK", "IF"]

    # Act
    result = translate_batch(sequences, n_processes=1)

    # Assert
    assert result == expected


def test_translate_batch_multiple_processes() -> None:
    """
    Test case 10: With multiple processes.
    """
    # Arrange
    sequences = ["ATGGCC"] * 20
    expected = ["MA"] * 20

    # Act
    result = translate_batch(sequences, n_processes=4)

    # Assert
    assert result == expected


def test_translate_batch_vertebrate_mito() -> None:
    """
    Test case 11: Vertebrate mitochondrial genetic code.
    """
    # Arrange
    sequences = ["ATGAGA", "ATGGCC"]

    # Act
    result = translate_batch(sequences, table=2)

    # Assert
    assert result[0] == translate_dna_to_protein("ATGAGA", table=2)
    assert result[1] == translate_dna_to_protein("ATGGCC", table=2)


def test_translate_batch_custom_table() -> None:
    """
    Test case 12: Custom codon table.
    """
    # Arrange
    from bioutils_collection.translation_functions import STANDARD_CODE
    custom = STANDARD_CODE.copy()
    custom['ATG'] = 'X'
    sequences = ["ATGGCC", "ATGAAA"]
    expected = ["XA", "XK"]

    # Act
    result = translate_batch(sequences, table=custom)

    # Assert
    assert result == expected


def test_translate_batch_large_count() -> None:
    """
    Test case 13: Large number of sequences.
    """
    # Arrange
    sequences = ["ATGGCCAAA"] * 100
    expected = ["MAK"] * 100

    # Act
    result = translate_batch(sequences, n_processes=4)

    # Assert
    assert len(result) == 100
    assert result == expected


def test_translate_batch_large_sequences() -> None:
    """
    Test case 14: Batch of large sequences.
    """
    # Arrange
    large_seq = "ATGGCCAAA" * 1000
    sequences = [large_seq] * 10

    # Act
    result = translate_batch(sequences, n_processes=2, use_fast=False)

    # Assert
    assert len(result) == 10
    assert all(len(r) == 3000 for r in result)


def test_translate_batch_mixed_content() -> None:
    """
    Test case 15: Sequences with different codon patterns.
    """
    # Arrange
    sequences = ["ATGATGATG", "GCCGCCGCC", "AAAAAACCC"]
    expected = ["MMM", "AAA", "KKP"]

    # Act
    result = translate_batch(sequences)

    # Assert
    assert result == expected


def test_translate_batch_small_no_multiprocessing() -> None:
    """
    Test case 16: Small batches don't use multiprocessing.
    """
    # Arrange
    sequences = ["ATGGCC"] * 5
    expected = ["MA"] * 5

    # Act
    result = translate_batch(sequences, n_processes=4)

    # Assert
    assert result == expected


def test_translate_batch_all_stops() -> None:
    """
    Test case 17: Sequences that are all stop codons.
    """
    # Arrange
    sequences = ["TAATAGTGA", "TGATGATAA"]
    expected = ["***", "***"]

    # Act
    result = translate_batch(sequences)

    # Assert
    assert result == expected


def test_translate_batch_bacterial_code() -> None:
    """
    Test case 18: Batch translation with bacterial code.
    """
    # Arrange
    sequences = ["ATGGCC", "AAATTT", "GGGTTT"]

    # Act
    result = translate_batch(sequences, table=11)
    expected = [translate_dna_to_protein(seq, table=11) for seq in sequences]

    # Assert
    assert result == expected


@pytest.mark.skipif(not NUMBA_AVAILABLE, reason="Numba not installed")
def test_translate_batch_fast_vs_default() -> None:
    """
    Test case 19: use_fast=True and False give same results.
    """
    # Arrange
    sequences = ["ATGGCC", "ATGAAATTCTGA", "GCTTGCGGCGGTTATTCA"]  # Fixed: all multiples of 3

    # Act
    result_default = translate_batch(sequences, use_fast=False, n_processes=1)
    result_fast = translate_batch(sequences, use_fast=True, n_processes=1)

    # Assert
    assert result_default == result_fast


def test_translate_batch_invalid_sequence_length() -> None:
    """
    Test case 20: ValueError for invalid sequence length in batch.
    """
    # Arrange
    sequences = ["ATGGCC", "ATGC", "GCTAGC"]  # Middle one is 4 bases

    # Act & Assert
    with pytest.raises(ValueError, match="Sequence length must be a multiple of 3"):
        translate_batch(sequences)


def test_translate_batch_invalid_table() -> None:
    """
    Test case 21: ValueError for invalid table parameter.
    """
    # Arrange
    sequences = ["ATGGCC", "GCTAGC"]

    # Act & Assert
    with pytest.raises((ValueError, KeyError)):
        translate_batch(sequences, table="invalid_table")


def test_translate_batch_non_string_sequence() -> None:
    """
    Test case 22: TypeError for non-string sequence.
    """
    # Arrange
    sequences = ["ATGGCC", None, "GCTAGC"]

    # Act & Assert
    with pytest.raises((AttributeError, TypeError)):
        translate_batch(sequences)


def test_translate_batch_non_list_input() -> None:
    """
    Test case 23: TypeError for non-list input.
    """
    # Arrange
    sequences = "ATGGCC"

    # Act & Assert
    with pytest.raises(TypeError, match="sequences must be list"):
        translate_batch(sequences)


def test_translate_batch_invalid_n_processes_negative() -> None:
    """
    Test case 24: ValueError for negative n_processes.
    """
    # Arrange
    sequences = ["ATGGCC"]

    # Act & Assert
    with pytest.raises(ValueError, match="n_processes must be positive"):
        translate_batch(sequences, n_processes=-1)


def test_translate_batch_invalid_n_processes_zero() -> None:
    """
    Test case 25: ValueError for zero n_processes.
    """
    # Arrange
    sequences = ["ATGGCC"]

    # Act & Assert
    with pytest.raises(ValueError, match="n_processes must be positive"):
        translate_batch(sequences, n_processes=0)
