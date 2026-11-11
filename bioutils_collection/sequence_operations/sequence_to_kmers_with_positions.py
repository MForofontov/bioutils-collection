"""Split sequences into k-mers with their starting positions."""


def sequence_to_kmers_with_positions(seq: str, k: int) -> list[tuple[str, int]]:
    """
    Split a sequence into k-mers of length k, each paired with its start position.

    Similar to sequence_to_kmers(), but returns (kmer, position) tuples
    instead of just k-mers. Useful when position information is needed
    for downstream processing such as minimizer computation.

    Parameters
    ----------
    seq : str
        Input sequence.
    k : int
        Length of each k-mer.

    Returns
    -------
    list[tuple[str, int]]
        List of (kmer, start_position) tuples in order of appearance.
        Positions are 0-indexed.

    Raises
    ------
    TypeError
        If seq is not a string or k is not an integer.
    ValueError
        If k is not positive or longer than sequence.

    Examples
    --------
    >>> sequence_to_kmers_with_positions('ATGCGA', 3)
    [('ATG', 0), ('TGC', 1), ('GCG', 2), ('CGA', 3)]
    >>> sequence_to_kmers_with_positions('ATGC', 2)
    [('AT', 0), ('TG', 1), ('GC', 2)]
    >>> sequence_to_kmers_with_positions('AAAA', 4)
    [('AAAA', 0)]

    Notes
    -----
    The position returned is the 0-indexed start position of the k-mer
    in the original sequence. For a sequence of length n and k-mer size k,
    positions range from 0 to n-k inclusive.

    For k-mers without position information, use sequence_to_kmers() instead.

    Complexity
    ----------
    Time: O(n*k), Space: O(n) where n is sequence length
    """
    if not isinstance(seq, str):
        raise TypeError(f"seq must be str, got {type(seq).__name__}")
    if not isinstance(k, int):
        raise TypeError(f"k must be int, got {type(k).__name__}")
    if k <= 0:
        raise ValueError("k must be positive")
    if k > len(seq):
        raise ValueError("k cannot be longer than sequence")
    return [(seq[i : i + k], i) for i in range(len(seq) - k + 1)]


__all__ = ["sequence_to_kmers_with_positions"]
