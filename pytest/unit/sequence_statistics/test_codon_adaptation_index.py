import pytest
import numpy
from bioutils_collection.sequence_statistics.codon_adaptation_index import (
    codon_adaptation_index,
)

pytestmark = [pytest.mark.unit, pytest.mark.sequence_statistics]

UNIFORM_ATG_WEIGHTS = {"ATG": 1.0}


def test_codon_adaptation_index_uniform() -> None:
    """Test case 1: Test CAI with uniform reference weights."""
    result = codon_adaptation_index("ATGATGATG", reference_weights=UNIFORM_ATG_WEIGHTS)
    assert 0.0 <= result <= 1.0
    assert result == 1.0


def test_codon_adaptation_index_with_weights() -> None:
    """Test case 2: Test CAI with custom reference weights."""
    weights = {"ATG": 1.0, "ATT": 0.5, "ATC": 0.8}
    result = codon_adaptation_index("ATGATT", reference_weights=weights)
    assert 0.0 <= result <= 1.0
    assert result < 1.0


def test_codon_adaptation_index_stop_codon() -> None:
    """Test case 3: Test CAI excludes stop codons."""
    result = codon_adaptation_index(
        "ATGATGTAA", reference_weights=UNIFORM_ATG_WEIGHTS
    )
    assert result > 0.0


def test_codon_adaptation_index_lowercase() -> None:
    """Test case 4: Test CAI with lowercase sequence."""
    result = codon_adaptation_index(
        "atgatgatg", reference_weights=UNIFORM_ATG_WEIGHTS
    )
    assert result == 1.0


def test_codon_adaptation_index_not_multiple_of_3() -> None:
    """Test case 5: Test ValueError for sequence not multiple of 3."""
    with pytest.raises(ValueError, match="Sequence length must be multiple of 3"):
        codon_adaptation_index("ATGATG A", reference_weights=UNIFORM_ATG_WEIGHTS)


def test_codon_adaptation_index_invalid_bases() -> None:
    """Test case 6: Test ValueError for invalid DNA bases."""
    with pytest.raises(ValueError, match="Invalid DNA bases found"):
        codon_adaptation_index("ATXATGATG", reference_weights=UNIFORM_ATG_WEIGHTS)


def test_codon_adaptation_index_empty() -> None:
    """Test case 7: Test ValueError for empty sequence."""
    with pytest.raises(ValueError, match="Sequence cannot be empty"):
        codon_adaptation_index("", reference_weights=UNIFORM_ATG_WEIGHTS)


def test_codon_adaptation_index_type_error() -> None:
    """Test case 8: Test TypeError for non-string input."""
    with pytest.raises(TypeError, match="seq must be a string"):
        codon_adaptation_index(123, reference_weights=UNIFORM_ATG_WEIGHTS)  # type: ignore[arg-type]


def test_codon_adaptation_index_weights_type_error() -> None:
    """Test case 9: Test TypeError for non-dict reference_weights."""
    with pytest.raises(TypeError, match="reference_weights must be a dict"):
        codon_adaptation_index("ATGATG", reference_weights="not_a_dict")  # type: ignore[arg-type]


def test_codon_adaptation_index_missing_reference_weights() -> None:
    """Test case 10: Test ValueError when reference_weights is None."""
    with pytest.raises(ValueError, match="reference_weights is required"):
        codon_adaptation_index("ATGATGATG")
