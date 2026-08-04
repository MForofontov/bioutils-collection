import pytest

from bioutils_collection.sequence_operations.find_orfs import find_orfs

pytestmark = [pytest.mark.unit, pytest.mark.sequence_operations]


def test_find_orfs_basic() -> None:
    """
    Test case 1: Basic ORF finding.
    """
    seq = "ATGAAATAGATGTAA"
    result = list(find_orfs(seq))
    assert isinstance(result, list)
    assert all(isinstance(item, tuple) for item in result)
    assert all(len(item) == 3 for item in result)


def test_find_orfs_no_orfs() -> None:
    """
    Test case 2: No ORFs found.
    """
    seq = "ATCATCATC"
    result = list(find_orfs(seq))
    assert result == []


def test_find_orfs_multiple_orfs() -> None:
    """
    Test case 3: Multiple ORFs in sequence.
    """
    seq = "ATGAAATAGATGTAA"
    result = list(find_orfs(seq))
    assert len(result) >= 1


def test_find_orfs_empty_sequence() -> None:
    """
    Test case 4: Empty sequence returns no ORFs.
    """
    result = list(find_orfs(""))
    assert result == []


def test_find_orfs_lowercase() -> None:
    """
    Test case 5: Lowercase input sequence.
    """
    seq = "atgaaatag"
    result = list(find_orfs(seq))
    assert len(result) >= 1


def test_find_orfs_invalid_type_error() -> None:
    """
    Test case 6: TypeError for invalid input type.
    """
    with pytest.raises(TypeError, match="seq must be str"):
        list(find_orfs(12345))  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="seq must be str"):
        list(find_orfs(None))  # type: ignore[arg-type]


def test_find_orfs_invalid_base_error() -> None:
    """
    Test case 7: ValueError for invalid DNA bases.
    """
    with pytest.raises(ValueError, match="Sequence contains invalid DNA bases"):
        list(find_orfs("ATGCX"))
    with pytest.raises(ValueError, match="Sequence contains invalid DNA bases"):
        list(find_orfs("ATGCU"))


def test_find_orfs_no_stop_codon() -> None:
    """
    Test case 8: ORF starting but no stop codon (now included in results).
    """
    # ATG followed by codons but no stop codon - should now be included
    seq = "ATGAAACCCTTT"
    result = list(find_orfs(seq))
    # Should return the ORF reaching end without stop codon
    assert len(result) == 1
    assert result[0][0] == 0  # Start position
    assert result[0][1] == 12  # End position (full sequence)
    assert result[0][2] == "ATGAAACCCTTT"  # ORF sequence


def test_find_orfs_partial_codon_truncated() -> None:
    """
    Test case 10: ORF truncated to last complete codon when sequence ends mid-codon.
    """
    assert list(find_orfs("ATGAA")) == [(0, 3, "ATG")]


def test_find_orfs_with_stop_codon() -> None:
    """
    Test case 9: ORF with proper stop codon.
    """
    seq = "ATGAAATAA"  # ATG + AAA + TAA (stop)
    result = list(find_orfs(seq))
    assert len(result) == 1
    assert result[0] == (0, 9, "ATGAAATAA")
