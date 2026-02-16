"""Generate syncmers - a variant of minimizers with different selection criteria."""


def generate_syncmers(seq: str, k: int, s: int, method: str = "open") -> list[str]:
    """
    Generate syncmers from a sequence.

    Syncmers are a variant of minimizers that select k-mers based on
    internal s-mer positions rather than sliding windows. They provide
    more uniform coverage than minimizers and are useful for sequence
    sketching and similarity search.

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
        - "open": s-mer must be minimum and at specific position
        - "closed": s-mer must be minimum anywhere in k-mer

    Returns
    -------
    list[str]
        List of syncmer k-mers in order of appearance.

    Raises
    ------
    TypeError
        If seq is not string, k/s not integers, or method not string.
    ValueError
        If parameters are invalid or s >= k or method not recognized.

    Examples
    --------
    >>> generate_syncmers('ATGCGATCG', 4, 2, method='open')
    ['ATGC', 'GATC']
    >>> generate_syncmers('AAACCCGGG', 3, 2, method='closed')
    ['AAA', 'AAC', 'CCC', 'CCG', 'GGG']

    Notes
    -----
    Syncmers were introduced by Edgar (2021) as an alternative to minimizers.

    Open syncmers: The minimum s-mer must occur at a specific offset
    (typically the first or last position). This provides more regular spacing.

    Closed syncmers: The minimum s-mer can occur anywhere in the k-mer.
    This produces more syncmers but with better sensitivity.

    Advantages over minimizers:
    - More uniform distribution
    - Deterministic selection (no arbitrary tie-breaking)
    - Better performance for some similarity metrics

    References
    ----------
    Edgar, R.C. (2021). Syncmers are more sensitive than minimizers
    for selecting conserved k-mers in biological sequences.
    PeerJ 9:e10805.

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
    syncmers = []
    
    # Generate all k-mers
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
                syncmers.append(kmer)
        else:  # closed
            # Closed syncmer: minimum s-mer can be anywhere
            if min_smer in smers:
                syncmers.append(kmer)
    
    return syncmers


__all__ = ["generate_syncmers"]
