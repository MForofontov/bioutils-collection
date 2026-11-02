"""Get positions of all syncmers in a sequence."""

from ._syncmer_selection import iter_syncmer_positions, validate_syncmer_inputs


def syncmer_positions(seq: str, k: int, s: int, method: str = "open") -> list[int]:
    """
    Get positions of all syncmers in a sequence.

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
    list[int]
        List of starting positions (0-indexed) of syncmers.

    Raises
    ------
    TypeError
        If seq is not string, k/s not integers, or method not string.
    ValueError
        If parameters are invalid or s >= k or method not recognized.

    Examples
    --------
    >>> syncmer_positions('ATGCGATCG', 4, 2, method='open')
    [0, 5]
    >>> syncmer_positions('AAACCCGGG', 3, 2, method='closed')
    [0, 1, 2, 3, 4, 5, 6]

    Notes
    -----
    Syncmer positions are useful for:
    - Anchoring sequence alignments
    - Building sequence indices
    - Identifying conserved regions
    - Sparse sequence comparison

    Complexity
    ----------
    Time: O(n*s), Space: O(n + (k-s)*s + m) including uppercase input and output
    """
    validate_syncmer_inputs(seq, k, s, method)
    return list(iter_syncmer_positions(seq.upper(), k, s, method))


__all__ = ["syncmer_positions"]
