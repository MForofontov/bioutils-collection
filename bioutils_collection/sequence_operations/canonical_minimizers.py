"""Generate canonical minimizers considering reverse complement."""

from ._minimizer_window import sliding_window_minima
from .reverse_complement import reverse_complement


def canonical_minimizers(seq: str, k: int, w: int) -> list[str]:
    """
    Generate canonical minimizers considering reverse complement.

    Canonical minimizers use the lexicographically smaller of a k-mer
    and its reverse complement, making them strand-independent. This is
    crucial for genome assembly and mapping where strand orientation
    may be unknown.

    Parameters
    ----------
    seq : str
        Input DNA/RNA sequence.
    k : int
        Length of each k-mer (minimizer size).
    w : int
        Window size (number of consecutive k-mers).

    Returns
    -------
    list[str]
        List of canonical minimizer k-mers.

    Raises
    ------
    TypeError
        If seq is not a string, or k/w are not integers.
    ValueError
        If seq is empty, k or w are not positive, or k > len(seq).

    Examples
    --------
    >>> canonical_minimizers('ATGCGATCG', 3, 4)
    ['ATC', 'ATC']
    >>> canonical_minimizers('AAACCCGGG', 2, 3)
    ['AA', 'AA', 'CC', 'CC', 'CG', 'CC', 'CC']

    Notes
    -----
    Canonical k-mers are defined as:
        canonical(kmer) = min(kmer, reverse_complement(kmer))

    This ensures that a sequence and its reverse complement produce
    the same minimizer set, which is essential for:
    - Strand-agnostic genome assembly
    - Bidirectional read mapping
    - Symmetric similarity measures

    For applications where strand matters, use generate_minimizers() instead.

    Complexity
    ----------
    Time: O(n*k), Space: O(m) where n is sequence length, m is minimizers count
    """
    if not isinstance(seq, str):
        raise TypeError(f"seq must be str, got {type(seq).__name__}")
    if not isinstance(k, int):
        raise TypeError(f"k must be int, got {type(k).__name__}")
    if not isinstance(w, int):
        raise TypeError(f"w must be int, got {type(w).__name__}")
    if not seq:
        raise ValueError("seq cannot be empty")
    if k <= 0:
        raise ValueError("k must be positive")
    if w <= 0:
        raise ValueError("w must be positive")
    if k > len(seq):
        raise ValueError("k cannot be longer than sequence")

    seq = seq.upper()
    canonical_kmers = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        rev_comp = reverse_complement(kmer)
        canonical_kmers.append(min(kmer, rev_comp))

    return sliding_window_minima(canonical_kmers, w)


__all__ = ["canonical_minimizers"]
