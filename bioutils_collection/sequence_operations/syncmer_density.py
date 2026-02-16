"""Calculate syncmer density (ratio of syncmers to k-mers)."""

from .syncmer_positions import syncmer_positions


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
    0.71

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
    Time: O(n*k*s), Space: O(1)
    """
    if not isinstance(seq, str):
        raise TypeError(f"seq must be str, got {type(seq).__name__}")
    if not isinstance(k, int):
        raise TypeError(f"k must be int, got {type(k).__name__}")
    if not isinstance(s, int):
        raise TypeError(f"s must be int, got {type(s).__name__}")
    if not isinstance(method, str):
        raise TypeError(f"method must be str, got {type(method).__name__}")
    if not seq:
        raise ValueError("seq cannot be empty")
    if k <= 0:
        raise ValueError("k must be positive")
    if s <= 0:
        raise ValueError("s must be positive")
    if s >= k:
        raise ValueError("s must be less than k")
    if k > len(seq):
        raise ValueError("k cannot be longer than sequence")
    if method not in ("open", "closed"):
        raise ValueError("method must be 'open' or 'closed'")

    total_kmers = len(seq) - k + 1
    if total_kmers == 0:
        return 0.0
    
    num_syncmers = len(syncmer_positions(seq, k, s, method))
    
    return num_syncmers / total_kmers


__all__ = ["syncmer_density"]
