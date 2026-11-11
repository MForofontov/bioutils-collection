"""Cross-check translation backends and file-level configuration failures."""

from itertools import product
from pathlib import Path

import pytest
from bioutils_collection.translation_functions import (
    CODON_TABLES,
    translate_dna_fast,
    translate_dna_to_protein,
    translate_fasta,
)

pytestmark = [
    pytest.mark.unit,
    pytest.mark.translation,
    pytest.mark.filterwarnings(
        "error::numba.core.errors.NumbaPendingDeprecationWarning"
    ),
]


@pytest.mark.parametrize("table", sorted(CODON_TABLES))
def test_fast_translation_matches_all_codons(table: str) -> None:
    sequence = "".join("".join(codon) for codon in product("ACGT", repeat=3))
    assert translate_dna_fast(sequence, table) == translate_dna_to_protein(
        sequence, table
    )


def test_fast_translation_sparse_custom_table() -> None:
    table = {"ATG": "start", "TAA": "*"}
    assert translate_dna_fast("ATGGCCTAA", table) == "startX*"


@pytest.mark.parametrize("source_text", ["", ">first\nATGGCC\n"])
def test_fasta_invalid_table_preserves_output(tmp_path: Path, source_text: str) -> None:
    source = tmp_path / "input.fa"
    source.write_text(source_text)
    destination = tmp_path / "output.fa"
    destination.write_text("previous results\n")

    with pytest.raises(ValueError, match="Unknown table"):
        translate_fasta(source, destination, table="invalid")

    assert destination.read_text() == "previous results\n"


@pytest.mark.parametrize("use_fast", [False, True])
def test_fasta_records_and_errors(tmp_path: Path, use_fast: bool) -> None:
    source = tmp_path / "input.fa"
    source.write_text(">first description\nATG\nGCC\n>bad\nAT\n>last\nTAA\n")
    destination = tmp_path / "output.fa"

    result = translate_fasta(source, destination, use_fast=use_fast)

    assert list(result) == ["first", "bad", "last"]
    assert result["first"] == "MA"
    assert result["last"] == "*"
    assert result["bad"].startswith("ERROR:")
    assert destination.read_text() == (
        ">first\nMA\n>bad\n" + result["bad"] + "\n>last\n*\n"
    )
