"""Get positions of minimizers in sequences."""

from ._minimizer_window import sliding_window_minima
from .sequence_to_kmers_with_positions import sequence_to_kmers_with_positions


def minimizer_positions(seq: str, k: int, w: int) -> list[int]:
    """
    Get positions of all minimizers in a sequence.

    Returns the starting position of each minimizer selected from
    sliding windows across the sequence. Useful for identifying
    important subsequence locations.

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
    list[int]
        List of starting positions (0-indexed) of minimizers.

    Raises
    ------
    TypeError
        If seq is not a string, or k/w are not integers.
    ValueError
        If seq is empty, k or w are not positive, or k > len(seq).

    Examples
    --------
    >>> minimizer_positions('ATGCGATCG', 3, 4)
    [4, 4]
    >>> minimizer_positions('AAACCCGGG', 2, 3)
    [0, 1, 3, 4, 5, 6, 7]

    Notes
    -----
    Positions indicate where important k-mers (minimizers) occur,
    which can be used for:
    - Anchoring sequence alignments
    - Identifying conserved regions
    - Sparse dynamic programming
    - Seed-and-extend algorithms

    The same minimizer at different positions will have multiple entries.

    Complexity
    ----------
    Time: O(n*k), Space: O(m) where m is number of windows
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
    kmers = sequence_to_kmers_with_positions(seq, k)
    selected = sliding_window_minima(kmers, w, key=lambda item: item[0])
    return [pos for _, pos in selected]


__all__ = ["minimizer_positions"]
