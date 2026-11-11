"""Generate syncmers - a variant of minimizers with different selection criteria."""

from ._syncmer_selection import iter_syncmer_positions, validate_syncmer_inputs


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
        - "closed": leftmost minimum s-mer must be at the first or last position

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
    ['ATGC', 'ATCG']
    >>> generate_syncmers('AAACCCGGG', 3, 2, method='closed')
    ['AAA', 'AAC', 'ACC', 'CCC', 'CCG', 'CGG', 'GGG']

    Notes
    -----
    Syncmers were introduced by Edgar (2021) as an alternative to minimizers.

    Open syncmers use offset zero: the first s-mer must be the leftmost
    minimum. Selection uses lexicographic ordering of uppercase s-mers.

    Closed syncmers: The leftmost minimum s-mer must be first or last.
    Equal minima are resolved from left to right, including boundary ties.

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
    Time: O(n*s + m*k), Space: O(n + (k-s)*s + m*k) for generated k-mers
    """
    validate_syncmer_inputs(seq, k, s, method)
    seq = seq.upper()
    return [seq[i : i + k] for i in iter_syncmer_positions(seq, k, s, method)]


__all__ = ["generate_syncmers"]
