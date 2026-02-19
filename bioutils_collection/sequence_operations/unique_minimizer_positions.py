"""Get unique minimizers and all their positions in sequence."""

from .sequence_to_kmers_with_positions import sequence_to_kmers_with_positions


def unique_minimizer_positions(seq: str, k: int, w: int) -> dict[str, list[int]]:
    """
    Get unique minimizers and all their positions in sequence.

    Unlike minimizer_positions(), this returns a mapping of each unique
    minimizer to all positions where it was selected as a minimizer,
    removing duplicates.

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
    dict[str, list[int]]
        Dictionary mapping each unique minimizer to its positions.

    Raises
    ------
    TypeError
        If seq is not a string, or k/w are not integers.
    ValueError
        If seq is empty, k or w are not positive, or k > len(seq).

    Examples
    --------
    >>> unique_minimizer_positions('ATGCGATCG', 3, 4)
    {'ATC': [4]}
    >>> unique_minimizer_positions('AAACCCGGG', 2, 3)
    {'AA': [0, 1], 'CC': [3, 4], 'CG': [5], 'GG': [6, 7]}

    Notes
    -----
    This function groups minimizer positions by k-mer sequence,
    which is useful for:
    - Finding all occurrences of important k-mers
    - Building position indices
    - Analyzing minimizer density
    - Clustering similar regions

    Positions are sorted and deduplicated for each minimizer.

    Complexity
    ----------
    Time: O(n*w*k), Space: O(m) where m is unique minimizers
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
    position_map: dict[str, set[int]] = {}

    # Generate all k-mers with positions
    kmers = sequence_to_kmers_with_positions(seq, k)

    if len(kmers) < w:
        # If we have fewer k-mers than window size, add minimum
        if kmers:
            min_kmer = min(kmers, key=lambda x: x[0])
            position_map[min_kmer[0]] = {min_kmer[1]}
        return {k: sorted(v) for k, v in position_map.items()}

    # Slide window and collect minimizer positions
    for i in range(len(kmers) - w + 1):
        window = kmers[i : i + w]
        # Find minimum k-mer in window
        min_kmer = min(window, key=lambda x: x[0])
        if min_kmer[0] not in position_map:
            position_map[min_kmer[0]] = set()
        position_map[min_kmer[0]].add(min_kmer[1])

    # Convert sets to sorted lists
    return {kmer: sorted(positions) for kmer, positions in position_map.items()}


__all__ = ["unique_minimizer_positions"]
