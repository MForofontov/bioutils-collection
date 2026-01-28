"""Compute evolutionary distance between sequences."""


def phylogenetic_distance(seq1: str, seq2: str) -> float:
    """
    Compute simple evolutionary distance (proportion of differences) between two sequences.

    Parameters
    ----------
    seq1 : str
        First sequence.
    seq2 : str
        Second sequence.

    Returns
    -------
    float
        Proportion of differing positions (0.0 to 1.0).

    Raises
    ------
    TypeError
        If seq1 or seq2 is not a string.
    ValueError
        If sequences are not the same length.
        If sequences are empty.

    Examples
    --------
    >>> phylogenetic_distance("ATGC", "ATGT")
    0.25
    
    >>> phylogenetic_distance("AAAA", "AAAA")
    0.0
    
    >>> phylogenetic_distance("ATGC", "CGTA")
    1.0
    
    Notes
    -----
    This is a simple p-distance metric. For more sophisticated
    evolutionary distance calculations, consider using specialized
    phylogenetic software.

    Complexity
    ----------
    Time: O(n), Space: O(1)
    """
    # Input validation
    if not isinstance(seq1, str):
        raise TypeError(f"seq1 must be str, got {type(seq1).__name__}")
    if not isinstance(seq2, str):
        raise TypeError(f"seq2 must be str, got {type(seq2).__name__}")
    
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be the same length")
    if len(seq1) == 0:
        raise ValueError("Sequences must not be empty")
    differences = sum(a != b for a, b in zip(seq1, seq2, strict=False))
    return differences / len(seq1)


__all__ = ["phylogenetic_distance"]
