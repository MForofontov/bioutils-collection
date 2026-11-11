"""Calculate statistics for FASTA sequences."""

from .fasta_parser import parse_fasta


def fasta_stats(fasta_str: str) -> tuple[int, float, float]:
    r"""
    Compute statistics for a FASTA string: (num_sequences, avg_length, avg_gc_content).

    Parameters
    ----------
    fasta_str : str
        FASTA formatted string.

    Returns
    -------
    Tuple[int, float, float]
        Number of sequences, average length, and mean per-record GC fraction
        (0 to 1). GC counting is case-insensitive; all sequence characters
        contribute to the denominator. Empty records contribute zero GC.

    Examples
    --------
    >>> fasta_stats('>seq1\nATGC\n>seq2\nGGTT')
    (2, 4.0, 0.5)
    """
    num = 0
    total_len = 0
    total_gc: float = 0.0
    for _, seq in parse_fasta(fasta_str):
        num += 1
        total_len += len(seq)
        gc_count = seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c")
        total_gc += gc_count / len(seq) if seq else 0
    avg_len = total_len / num if num else 0.0
    avg_gc = total_gc / num if num else 0.0
    return num, avg_len, avg_gc


__all__ = ["fasta_stats"]
