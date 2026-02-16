"""Create minimizer sketch from k-mer generator (memory efficient)."""

from collections.abc import Iterator


def minimizer_sketch_from_generator(
    kmers_with_pos: Iterator[tuple[str, int]], w: int
) -> dict[str, list[int]]:
    """
    Create minimizer sketch from k-mer generator (memory efficient).

    Useful for streaming large sequences without loading entire sequence
    into memory at once.

    Parameters
    ----------
    kmers_with_pos : Iterator[tuple[str, int]]
        Iterator yielding (kmer, position) tuples.
    w : int
        Window size (number of consecutive k-mers).

    Returns
    -------
    dict[str, list[int]]
        Dictionary mapping minimizer k-mers to their positions.

    Raises
    ------
    TypeError
        If w is not an integer.
    ValueError
        If w is not positive.

    Examples
    --------
    >>> kmers = [('ATG', 0), ('TGC', 1), ('GCG', 2)]
    >>> minimizer_sketch_from_generator(iter(kmers), 2)
    {'ATG': [0], 'GCG': [2]}

    Notes
    -----
    This function is designed for streaming applications where
    the full sequence cannot fit in memory.

    Complexity
    ----------
    Time: O(n*w), Space: O(w + m) where m is unique minimizers
    """
    if not isinstance(w, int):
        raise TypeError(f"w must be int, got {type(w).__name__}")
    if w <= 0:
        raise ValueError("w must be positive")

    sketch: dict[str, list[int]] = {}
    window: list[tuple[str, int]] = []
    
    for kmer, pos in kmers_with_pos:
        window.append((kmer, pos))
        
        if len(window) == w:
            # Find minimizer in window
            min_kmer = min(window, key=lambda x: x[0])
            if min_kmer[0] not in sketch:
                sketch[min_kmer[0]] = []
            sketch[min_kmer[0]].append(min_kmer[1])
            
            # Slide window
            window.pop(0)
    
    # Process last window if exists
    if window and len(window) > 0:
        min_kmer = min(window, key=lambda x: x[0])
        if min_kmer[0] not in sketch:
            sketch[min_kmer[0]] = []
        sketch[min_kmer[0]].append(min_kmer[1])
    
    return sketch


__all__ = ["minimizer_sketch_from_generator"]
