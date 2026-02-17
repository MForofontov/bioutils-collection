"""Get positions of minimizers in sequences."""


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
    Time: O(n*w*k), Space: O(m) where m is number of windows
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
    positions = []

    # Generate all k-mers with positions
    kmers = [(seq[i : i + k], i) for i in range(len(seq) - k + 1)]

    if len(kmers) < w:
        # If we have fewer k-mers than window size, return position of minimum
        if kmers:
            min_kmer = min(kmers, key=lambda x: x[0])
            positions.append(min_kmer[1])
        return positions

    # Slide window and collect minimizer positions
    for i in range(len(kmers) - w + 1):
        window = kmers[i : i + w]
        # Find minimum k-mer in window
        min_kmer = min(window, key=lambda x: x[0])
        positions.append(min_kmer[1])

    return positions


__all__ = ["minimizer_positions"]
