import pytest
from bioutils_collection.sequence_operations.sequence_to_kmers_with_positions import (
    sequence_to_kmers_with_positions,
)

pytestmark = [pytest.mark.unit, pytest.mark.sequence_operations]


def test_sequence_to_kmers_with_positions_basic() -> None:
    """Test case 1: Basic k-mer extraction with positions."""
    result = sequence_to_kmers_with_positions("ATGCGA", 3)
    assert result == [("ATG", 0), ("TGC", 1), ("GCG", 2), ("CGA", 3)]


def test_sequence_to_kmers_with_positions_k2() -> None:
    """Test case 2: k=2 returns dinucleotides with correct positions."""
    result = sequence_to_kmers_with_positions("ATGC", 2)
    assert result == [("AT", 0), ("TG", 1), ("GC", 2)]


def test_sequence_to_kmers_with_positions_k_equals_length() -> None:
    """Test case 3: k equal to sequence length returns single kmer at position 0."""
    result = sequence_to_kmers_with_positions("AAAA", 4)
    assert result == [("AAAA", 0)]
    assert len(result) == 1


def test_sequence_to_kmers_with_positions_returns_list_of_tuples() -> None:
    """Test case 4: Return type is list of (str, int) tuples."""
    result = sequence_to_kmers_with_positions("ATGC", 2)
    assert isinstance(result, list)
    assert all(isinstance(item, tuple) for item in result)
    assert all(isinstance(kmer, str) and isinstance(pos, int) for kmer, pos in result)


def test_sequence_to_kmers_with_positions_count() -> None:
    """Test case 5: Number of k-mers is len(seq) - k + 1."""
    seq = "ATGCATGC"
    k = 3
    result = sequence_to_kmers_with_positions(seq, k)
    assert len(result) == len(seq) - k + 1


def test_sequence_to_kmers_with_positions_positions_are_sequential() -> None:
    """Test case 6: Positions increment by 1 each step."""
    result = sequence_to_kmers_with_positions("ATGCGATCG", 3)
    positions = [pos for _, pos in result]
    assert positions == list(range(len(positions)))


def test_sequence_to_kmers_with_positions_kmer_content() -> None:
    """Test case 7: Each kmer matches its slice in the original sequence."""
    seq = "ATGCGATCG"
    k = 4
    result = sequence_to_kmers_with_positions(seq, k)
    for kmer, pos in result:
        assert kmer == seq[pos : pos + k]


def test_sequence_to_kmers_with_positions_k1() -> None:
    """Test case 8: k=1 returns each character with its index."""
    result = sequence_to_kmers_with_positions("ATGC", 1)
    assert result == [("A", 0), ("T", 1), ("G", 2), ("C", 3)]


def test_sequence_to_kmers_with_positions_lowercase_preserved() -> None:
    """Test case 9: Lowercase input is preserved as-is."""
    result = sequence_to_kmers_with_positions("atgc", 2)
    assert result == [("at", 0), ("tg", 1), ("gc", 2)]


def test_sequence_to_kmers_with_positions_positions_zero_indexed() -> None:
    """Test case 10: First position is always 0."""
    result = sequence_to_kmers_with_positions("ATGCATGC", 3)
    assert result[0][1] == 0


def test_sequence_to_kmers_with_positions_last_position() -> None:
    """Test case 11: Last position is len(seq) - k."""
    seq = "ATGCATGC"
    k = 3
    result = sequence_to_kmers_with_positions(seq, k)
    assert result[-1][1] == len(seq) - k


def test_sequence_to_kmers_with_positions_repeated_sequence() -> None:
    """Test case 12: Repeated sequence produces repeated kmers with different positions."""
    result = sequence_to_kmers_with_positions("ATATAT", 2)
    kmers = [kmer for kmer, _ in result]
    positions = [pos for _, pos in result]
    assert kmers == ["AT", "TA", "AT", "TA", "AT"]
    assert positions == [0, 1, 2, 3, 4]


def test_sequence_to_kmers_with_positions_long_sequence() -> None:
    """Test case 13: Long sequence returns correct count and positions."""
    seq = "ATGC" * 100
    k = 5
    result = sequence_to_kmers_with_positions(seq, k)
    assert len(result) == len(seq) - k + 1
    assert result[0][1] == 0
    assert result[-1][1] == len(seq) - k


def test_sequence_to_kmers_with_positions_single_char_seq() -> None:
    """Test case 14: Single character sequence with k=1."""
    result = sequence_to_kmers_with_positions("A", 1)
    assert result == [("A", 0)]


def test_sequence_to_kmers_with_positions_type_error_not_string() -> None:
    """Test case 15: TypeError when seq is not a string."""
    with pytest.raises(TypeError, match="seq must be str, got int"):
        sequence_to_kmers_with_positions(12345, 3)  # type: ignore[arg-type]


def test_sequence_to_kmers_with_positions_type_error_none() -> None:
    """Test case 16: TypeError when seq is None."""
    with pytest.raises(TypeError, match="seq must be str, got NoneType"):
        sequence_to_kmers_with_positions(None, 3)  # type: ignore[arg-type]


def test_sequence_to_kmers_with_positions_type_error_k_not_int() -> None:
    """Test case 17: TypeError when k is not an integer."""
    with pytest.raises(TypeError, match="k must be int, got float"):
        sequence_to_kmers_with_positions("ATGC", 2.0)  # type: ignore[arg-type]


def test_sequence_to_kmers_with_positions_value_error_k_zero() -> None:
    """Test case 18: ValueError when k is zero."""
    with pytest.raises(ValueError, match="k must be positive"):
        sequence_to_kmers_with_positions("ATGC", 0)


def test_sequence_to_kmers_with_positions_value_error_k_negative() -> None:
    """Test case 19: ValueError when k is negative."""
    with pytest.raises(ValueError, match="k must be positive"):
        sequence_to_kmers_with_positions("ATGC", -1)


def test_sequence_to_kmers_with_positions_value_error_k_too_long() -> None:
    """Test case 20: ValueError when k exceeds sequence length."""
    with pytest.raises(ValueError, match="k cannot be longer than sequence"):
        sequence_to_kmers_with_positions("ATGC", 5)
