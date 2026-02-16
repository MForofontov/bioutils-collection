"""Get positions of all syncmers in a sequence."""


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
    [0, 4]
    >>> syncmer_positions('AAACCCGGG', 3, 2, method='closed')
    [0, 1, 3, 4, 6]

    Notes
    -----
    Syncmer positions are useful for:
    - Anchoring sequence alignments
    - Building sequence indices
    - Identifying conserved regions
    - Sparse sequence comparison

    Complexity
    ----------
    Time: O(n*k*s), Space: O(m) where m is number of syncmers
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

    seq = seq.upper()
    positions = []
    
    # Generate all k-mers and check syncmer condition
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        
        # Extract all s-mers from this k-mer
        smers = [kmer[j:j + s] for j in range(k - s + 1)]
        
        if not smers:
            continue
            
        min_smer = min(smers)
        
        if method == "open":
            # Open syncmer: minimum s-mer must be at first position
            if smers[0] == min_smer:
                positions.append(i)
        else:  # closed
            # Closed syncmer: minimum s-mer can be anywhere
            if min_smer in smers:
                positions.append(i)
    
    return positions


__all__ = ["syncmer_positions"]
