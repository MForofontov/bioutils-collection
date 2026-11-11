"""Create minimizer sketch from k-mer generator (memory efficient)."""

from collections import deque
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
    the full sequence cannot fit in memory. Selects one leftmost minimum per
    complete window. If the entire input is shorter than w, selects one minimum;
    incomplete trailing windows after a complete window are skipped.

    Complexity
    ----------
    Time: O(n*k), Space: O(w*k + r) where k is k-mer length and r is output size
    """
    if not isinstance(w, int):
        raise TypeError(f"w must be int, got {type(w).__name__}")
    if w <= 0:
        raise ValueError("w must be positive")

    sketch: dict[str, list[int]] = {}
    candidates: deque[tuple[int, str, int]] = deque()
    count = 0
    for index, (kmer, pos) in enumerate(kmers_with_pos):
        count += 1
        while candidates and candidates[0][0] <= index - w:
            candidates.popleft()
        while candidates and candidates[-1][1] > kmer:
            candidates.pop()
        candidates.append((index, kmer, pos))
        if count >= w:
            _, minimum, position = candidates[0]
            sketch.setdefault(minimum, []).append(position)
    if 0 < count < w:
        _, minimum, position = candidates[0]
        sketch.setdefault(minimum, []).append(position)
    return sketch


__all__ = ["minimizer_sketch_from_generator"]
