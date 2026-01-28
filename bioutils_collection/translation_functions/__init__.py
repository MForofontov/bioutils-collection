"""
translation_functions: DNA/RNA translation utilities.
"""

from .rna_to_dna import rna_to_dna
from .transcribe_dna_to_rna import transcribe_dna_to_rna
from .translate_dna_to_protein import translate_dna_to_protein
from .translate_dna_fast import translate_dna_fast, NUMBA_AVAILABLE
from .translate_batch import translate_batch
from .translate_large_sequence import translate_large_sequence
from .translate_fasta import translate_fasta
from .translate_fasta_streaming import translate_fasta_streaming
from .genetic_code_tables import (
    get_codon_table,
    CODON_TABLES,
    STANDARD_CODE,
    VERTEBRATE_MITOCHONDRIAL,
    YEAST_MITOCHONDRIAL,
)

__all__ = [
    "rna_to_dna",
    "transcribe_dna_to_rna",
    "translate_dna_to_protein",
    "translate_dna_fast",
    "translate_batch",
    "translate_large_sequence",
    "translate_fasta",
    "translate_fasta_streaming",
    "get_codon_table",
    "CODON_TABLES",
    "STANDARD_CODE",
    "VERTEBRATE_MITOCHONDRIAL",
    "YEAST_MITOCHONDRIAL",
    "NUMBA_AVAILABLE",
]
