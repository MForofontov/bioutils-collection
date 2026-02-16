"""Calculate minimizer density (ratio of minimizers to k-mers)."""

from .unique_minimizer_positions import unique_minimizer_positions


def minimizer_density(seq: str, k: int, w: int) -> float:
    """
    Calculate minimizer density (ratio of minimizers to k-mers).

    Minimizer density indicates how efficiently the minimizer scheme
    compresses the sequence. Lower density means better compression.

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
    float
        Density ratio between 0.0 and 1.0.

    Raises
    ------
    TypeError
        If seq is not a string, or k/w are not integers.
    ValueError
        If seq is empty, k or w are not positive, or k > len(seq).

    Examples
    --------
    >>> round(minimizer_density('ATGCGATCG', 3, 4), 2)
    0.29
    >>> round(minimizer_density('AAACCCGGG', 2, 3), 2)
    0.88

    Notes
    -----
    Density is calculated as:
        density = unique_minimizers / total_kmers

    Lower density indicates:
    - Better sequence compression
    - More repetitive minimizers
    - Higher data reduction

    Typical densities range from 0.1 to 0.5 for genomic sequences.

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

    total_kmers = len(seq) - k + 1
    if total_kmers == 0:
        return 0.0
    
    # Get unique minimizers
    unique_mins = unique_minimizer_positions(seq, k, w)
    unique_count = len(unique_mins)
    
    return unique_count / total_kmers


__all__ = ["minimizer_density"]
