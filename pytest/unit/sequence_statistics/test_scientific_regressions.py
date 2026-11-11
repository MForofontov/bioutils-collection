"""Small analytical fixtures, independent of implementation branches."""

import math

import pytest
from bioutils_collection.data_validation import validate_protein_sequence
from bioutils_collection.sequence_statistics.calculate_isoelectric_point import (
    calculate_isoelectric_point,
)
from bioutils_collection.sequence_statistics.calculate_molecular_weight import (
    calculate_molecular_weight,
)
from bioutils_collection.sequence_statistics.codon_adaptation_index import (
    codon_adaptation_index,
)


def test_zero_cai_weight_is_not_discarded() -> None:
    assert codon_adaptation_index("ATGATT", {"ATG": 1.0, "ATT": 0.0}) == 0.0


@pytest.mark.parametrize("weight", [-0.1, 1.1, math.inf, math.nan])
def test_cai_rejects_invalid_weights(weight: float) -> None:
    with pytest.raises(ValueError, match="weight"):
        codon_adaptation_index("ATG", {"ATG": weight})


def test_terminal_charges_determine_neutral_peptide_pi() -> None:
    assert calculate_isoelectric_point("AG") == pytest.approx(
        (9.69 + 2.34) / 2, abs=0.0051
    )


@pytest.mark.parametrize(
    "seq", ["CCCC", "YYYY", "DEKR", "ACDEFGH", "KDEKDEKD", "KKKKK"]
)
def test_pi_solves_charge_balance(seq: str) -> None:
    ph = calculate_isoelectric_point(seq)
    # Independent Henderson-Hasselbalch charge evaluation at the rounded root.
    positive = 1 / (1 + 10 ** (ph - 9.69))
    positive += sum(
        seq.count(aa) / (1 + 10 ** (ph - pk))
        for aa, pk in [("K", 10.53), ("R", 12.48), ("H", 6.0)]
    )
    negative = 1 / (1 + 10 ** (2.34 - ph))
    negative += sum(
        seq.count(aa) / (1 + 10 ** (pk - ph))
        for aa, pk in [("D", 3.86), ("E", 4.25), ("C", 8.33), ("Y", 10.07)]
    )
    assert abs(positive - negative) < 0.02


def test_polymer_mass_accounts_for_water_and_termini() -> None:
    assert calculate_molecular_weight("GG", "protein") == pytest.approx(
        132.18, abs=0.01
    )
    assert calculate_molecular_weight("A", "dna") == pytest.approx(331.22, abs=0.02)
    assert calculate_molecular_weight("A", "rna") == pytest.approx(347.22, abs=0.02)
    assert validate_protein_sequence("GG")["molecular_weight"] == pytest.approx(
        132.18, abs=0.01
    )
