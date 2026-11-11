"""Every supported NCBI table against an independently retrieved AAs row."""

import json
from pathlib import Path

from bioutils_collection.translation_functions.genetic_code_tables import (
    get_codon_table,
)
from bioutils_collection.translation_functions.translate_dna_to_protein import (
    translate_dna_to_protein,
)


def test_all_registered_ncbi_codes() -> None:
    fixture = json.loads(
        (Path(__file__).parents[2] / "fixtures/ncbi_genetic_codes.json").read_text()
    )
    codons = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
    for identifier, expected in fixture["tables"].items():
        table = get_codon_table(identifier)
        assert "".join(table[codon] for codon in codons) == expected, identifier
        assert translate_dna_to_protein("".join(codons), identifier) == expected
