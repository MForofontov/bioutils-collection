"""Calculate Levenshtein (edit) distance between sequences."""


def levenshtein_distance(seq1: str, seq2: str) -> int:
    """
    Calculate the Levenshtein distance (edit distance) between two sequences.

    The Levenshtein distance is the minimum number of single-character edits
    (insertions, deletions, or substitutions) required to change one sequence
    into the other.

    Parameters
    ----------
    seq1 : str
        First sequence.
    seq2 : str
        Second sequence.

    Returns
    -------
    int
        The Levenshtein distance between the two sequences.

    Raises
    ------
    TypeError
        If seq1 or seq2 is not a string.

    Examples
    --------
    >>> levenshtein_distance("ACGT", "ACGT")
    0
    >>> levenshtein_distance("ACGT", "ACT")
    1
    >>> levenshtein_distance("kitten", "sitting")
    3
    >>> levenshtein_distance("", "ABC")
    3
    >>> levenshtein_distance("ABC", "")
    3

    Notes
    -----
    The algorithm uses dynamic programming with O(n*m) time complexity.
    Space complexity is O(min(n,m)) using a two-row rolling array.

    References
    ----------
    Levenshtein, V.I. (1966).
    Binary codes capable of correcting deletions, insertions, and reversals.
    Soviet Physics Doklady 10(8):707-710.

    Complexity
    ----------
    Time: O(n*m), Space: O(min(n,m)) where n, m are sequence lengths
    """
    # Input validation
    if not isinstance(seq1, str):
        raise TypeError(f"seq1 must be a string, got {type(seq1).__name__}")
    if not isinstance(seq2, str):
        raise TypeError(f"seq2 must be a string, got {type(seq2).__name__}")

    # Handle empty strings
    if len(seq1) == 0:
        return len(seq2)
    if len(seq2) == 0:
        return len(seq1)

    # Use shorter sequence for columns to minimize memory
    if len(seq1) < len(seq2):
        seq1, seq2 = seq2, seq1

    n, m = len(seq1), len(seq2)
    previous = list(range(m + 1))
    current = [0] * (m + 1)

    for i in range(1, n + 1):
        current[0] = i
        for j in range(1, m + 1):
            cost = 0 if seq1[i - 1] == seq2[j - 1] else 1
            current[j] = min(
                previous[j] + 1,
                current[j - 1] + 1,
                previous[j - 1] + cost,
            )
        previous, current = current, previous

    return previous[m]


__all__ = ["levenshtein_distance"]
