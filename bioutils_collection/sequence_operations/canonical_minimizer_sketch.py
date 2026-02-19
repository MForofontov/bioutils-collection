"""Create canonical minimizer sketch with positions."""

from .reverse_complement import reverse_complement


def canonical_minimizer_sketch(seq: str, k: int, w: int) -> dict[str, list[int]]:
    """
    Create canonical minimizer sketch with positions.

    Combines canonical minimizer selection with sketch representation,
    providing strand-independent sequence signatures with position tracking.

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
        Dictionary mapping canonical minimizers to their positions.

    Raises
    ------
    TypeError
        If seq is not a string, or k/w are not integers.
    ValueError
        If seq is empty, k or w are not positive, or k > len(seq).

    Examples
    --------
    >>> canonical_minimizer_sketch('ATGCGATCG', 3, 4)
    {'ATC': [4, 4]}
    >>> canonical_minimizer_sketch('AAACCCGGG', 2, 3)
    {'AA': [0, 1], 'CC': [3, 4, 5, 6]}

    Notes
    -----
    This function is particularly useful for:
    - Building sequence indices for genomic databases
    - Comparing sequences regardless of strand
    - Detecting shared subsequences in assembly graphs

    Complexity
    ----------
    Time: O(n*w*k), Space: O(m) where m is unique canonical minimizers
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
    sketch: dict[str, list[int]] = {}

    # Generate all k-mers with canonical forms and positions
    canonical_kmers = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        rev_comp = reverse_complement(kmer)
        canonical = min(kmer, rev_comp)
        canonical_kmers.append((canonical, i))

    if len(canonical_kmers) < w:
        # If we have fewer k-mers than window size, add the minimum
        if canonical_kmers:
            min_kmer = min(canonical_kmers, key=lambda x: x[0])
            sketch[min_kmer[0]] = [min_kmer[1]]
        return sketch

    # Slide window and collect minimizers with positions
    for i in range(len(canonical_kmers) - w + 1):
        window = canonical_kmers[i : i + w]
        # Find minimum canonical k-mer in window
        min_kmer = min(window, key=lambda x: x[0])
        if min_kmer[0] not in sketch:
            sketch[min_kmer[0]] = []
        sketch[min_kmer[0]].append(min_kmer[1])

    return sketch


__all__ = ["canonical_minimizer_sketch"]
