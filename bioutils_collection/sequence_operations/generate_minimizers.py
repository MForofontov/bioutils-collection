"""Generate minimizers from sequences using sliding window approach."""

from ._minimizer_window import sliding_window_minima
from .sequence_to_kmers import sequence_to_kmers


def generate_minimizers(seq: str, k: int, w: int) -> list[str]:
    """
    Generate minimizers from a sequence using a sliding window approach.

    A minimizer is the lexicographically smallest k-mer in a window of size w.
    Minimizers reduce data redundancy while maintaining sequence information,
    making them useful for genome assembly, read mapping, and similarity search.

    Parameters
    ----------
    seq : str
        Input DNA/RNA sequence.
    k : int
        Length of each k-mer (minimizer size).
    w : int
        Window size (number of consecutive k-mers to consider).

    Returns
    -------
    list[str]
        List of minimizer k-mers in order of appearance.

    Raises
    ------
    TypeError
        If seq is not a string, or k/w are not integers.
    ValueError
        If seq is empty, k or w are not positive, or k > len(seq).

    Examples
    --------
    >>> generate_minimizers('ATGCGATCG', 3, 4)
    ['ATC', 'ATC']
    >>> generate_minimizers('AAACCCGGG', 2, 3)
    ['AA', 'AA', 'CC', 'CC', 'CG', 'GG', 'GG']
    >>> generate_minimizers('ATGC', 2, 2)
    ['AT', 'AT', 'GC']

    Notes
    -----
    Minimizers were introduced by Roberts et al. (2004) for reducing
    sequence redundancy while maintaining sensitivity for sequence comparison.
    The lexicographically smallest k-mer in each window is selected.

    For canonical minimizers that consider reverse complement, use
    canonical_minimizers() instead.

    References
    ----------
    Roberts, M., Hayes, W., Hunt, B.R., Mount, S.M., Yorke, J.A. (2004).
    Reducing storage requirements for biological sequence comparison.
    Bioinformatics 20(18), 3363-3369.

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
    kmers = sequence_to_kmers(seq, k)
    return sliding_window_minima(kmers, w)


__all__ = ["generate_minimizers"]
