"""Generate reverse complement of FASTA sequences."""

from collections.abc import Iterator

from .fasta_parser import parse_fasta
from ..sequence_operations.reverse_complement import reverse_complement


def fasta_reverse_complement(fasta_str: str) -> Iterator[tuple[str, str]]:
    """
    Generate reverse complement for all DNA sequences in a FASTA string.

    Parameters
    ----------
    fasta_str : str
        FASTA formatted string.

    Yields
    ------
    Tuple[str, str]
        (header, reverse_complement_sequence) pairs.

    Examples
    --------
    >>> list(fasta_reverse_complement('>seq1\nATGC'))
    [('seq1', 'GCAT')]
    """
    for header, seq in parse_fasta(fasta_str):
        yield (header, reverse_complement(seq))


__all__ = ["fasta_reverse_complement"]
