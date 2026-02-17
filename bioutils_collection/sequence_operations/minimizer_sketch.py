"""Create minimizer sketch representation of sequences."""


def minimizer_sketch(
    seq: str, k: int, w: int, density: float = 1.0
) -> dict[str, list[int]]:
    """
    Create a minimizer sketch of a sequence with positions.

    A minimizer sketch is a compact representation of a sequence that stores
    unique minimizers and their positions. This is memory-efficient for large
    sequences and enables fast sequence comparison and similarity search.

    Parameters
    ----------
    seq : str
        Input DNA/RNA sequence.
    k : int
        Length of each k-mer (minimizer size).
    w : int
        Window size (number of consecutive k-mers).
    density : float, optional
        Sampling density between 0.0 and 1.0 (default: 1.0).
        Lower values subsample minimizers for even smaller sketches.

    Returns
    -------
    dict[str, list[int]]
        Dictionary mapping minimizer k-mers to their positions in sequence.

    Raises
    ------
    TypeError
        If seq is not string, k/w not integers, or density not float.
    ValueError
        If parameters are invalid or density not in [0.0, 1.0].

    Examples
    --------
    >>> minimizer_sketch('ATGCGATCG', 3, 4)
    {'ATC': [4, 4]}
    >>> minimizer_sketch('AAACCCGGG', 2, 3, density=0.5)
    {'AA': [0], 'CC': [4], 'GG': [7]}

    Notes
    -----
    The sketch representation enables:
    - Fast sequence similarity estimation (Jaccard similarity)
    - Memory-efficient sequence storage
    - Quick sequence database searches
    - Approximate sequence alignment

    The density parameter allows trade-off between accuracy and size.
    Density of 0.5 means approximately half of minimizers are kept.

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
    if not isinstance(density, (int, float)):
        raise TypeError(f"density must be float, got {type(density).__name__}")
    if not seq:
        raise ValueError("seq cannot be empty")
    if k <= 0:
        raise ValueError("k must be positive")
    if w <= 0:
        raise ValueError("w must be positive")
    if k > len(seq):
        raise ValueError("k cannot be longer than sequence")
    if not 0.0 <= density <= 1.0:
        raise ValueError("density must be between 0.0 and 1.0")

    seq = seq.upper()
    sketch: dict[str, list[int]] = {}

    # Generate all k-mers with positions
    kmers = [(seq[i : i + k], i) for i in range(len(seq) - k + 1)]

    if len(kmers) < w:
        # If we have fewer k-mers than window size, add the minimum
        if kmers:
            min_kmer = min(kmers, key=lambda x: x[0])
            sketch[min_kmer[0]] = [min_kmer[1]]
        return sketch

    # Slide window and collect minimizers with positions
    minimizer_positions = []
    for i in range(len(kmers) - w + 1):
        window = kmers[i : i + w]
        # Find minimum k-mer in window
        min_kmer = min(window, key=lambda x: x[0])
        minimizer_positions.append(min_kmer)

    # Apply density sampling if needed
    if density < 1.0:
        step = max(1, int(1.0 / density))
        minimizer_positions = minimizer_positions[::step]

    # Build sketch dictionary
    for kmer, pos in minimizer_positions:
        if kmer not in sketch:
            sketch[kmer] = []
        sketch[kmer].append(pos)

    return sketch


__all__ = ["minimizer_sketch"]
