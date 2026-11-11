"""Calculate syncmer density (ratio of syncmers to k-mers)."""

from ._syncmer_selection import iter_syncmer_positions, validate_syncmer_inputs


def syncmer_density(seq: str, k: int, s: int, method: str = "open") -> float:
    """
    Calculate syncmer density (ratio of syncmers to k-mers).

    Syncmer density indicates the proportion of k-mers selected as syncmers.
    This metric helps compare different syncmer parameters and methods.

    Parameters
    ----------
    seq : str
        Input DNA/RNA sequence.
    k : int
        Length of each k-mer (syncmer size).
    s : int
        Length of s-mer used for selection (s < k).
    method : str, optional
        Selection method: "open" or "closed" (default: "open").

    Returns
    -------
    float
        Density ratio between 0.0 and 1.0.

    Raises
    ------
    TypeError
        If seq is not string, k/s not integers, or method not string.
    ValueError
        If parameters are invalid or s >= k or method not recognized.

    Examples
    --------
    >>> round(syncmer_density('ATGCGATCG', 4, 2, method='open'), 2)
    0.33
    >>> round(syncmer_density('AAACCCGGG', 3, 2, method='closed'), 2)
    1.0

    Notes
    -----
    Density is calculated as:
        density = num_syncmers / total_kmers

    Typical open syncmer density: ~1/(k-s+1)
    Closed syncmer density is typically higher than open.

    Lower density means:
    - Better data reduction
    - Fewer selected k-mers
    - Faster downstream processing

    Complexity
    ----------
    Time: O(n*s), Space: O(n + (k-s)*s), including uppercase input
    """
    validate_syncmer_inputs(seq, k, s, method)
    count = sum(1 for _ in iter_syncmer_positions(seq.upper(), k, s, method))
    return count / (len(seq) - k + 1)


__all__ = ["syncmer_density"]
