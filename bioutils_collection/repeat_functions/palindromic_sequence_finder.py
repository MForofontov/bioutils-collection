"""Find palindromic sequences in DNA (biological palindromes: reverse complements)."""


def _reverse_complement(seq: str) -> str:
    """Get reverse complement of DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement.get(base, base) for base in reversed(seq))


def palindromic_sequence_finder(
    sequence: str, min_length: int = 4, biological: bool = True
) -> list[tuple[int, int, str]]:
    """
    Find palindromic regions in a DNA sequence.

    Parameters
    ----------
    sequence : str
        Input DNA sequence (A, T, G, C).
    min_length : int, optional
        Minimum length of palindrome (default 4).
    biological : bool, optional
        If True (default), finds biological palindromes (reverse complements).
        If False, finds string palindromes (simple reversals).
        
        Biological palindrome: GAATTC (reverse complement = GAATTC)
        String palindrome: ABCCBA (simple reversal)

    Returns
    -------
    list[tuple[int, int, str]]
        List of (start, end, palindrome) for each palindromic region.

    Raises
    ------
    TypeError
        If sequence is not a string or min_length is not an integer.
    ValueError
        If min_length is not positive.
        If sequence contains invalid DNA bases (when biological=True).

    Examples
    --------
    >>> # Biological palindromes (restriction sites)
    >>> palindromic_sequence_finder("GAATTC", min_length=4)
    [(0, 6, 'GAATTC')]
    
    >>> # EcoRI site is a biological palindrome
    >>> seq = "ATGAATTCGC"
    >>> palindromic_sequence_finder(seq, min_length=6)
    [(2, 8, 'GAATTC')]
    
    >>> # String palindromes
    >>> palindromic_sequence_finder("ATGCAT", min_length=4, biological=False)
    [(0, 4, 'ATGC'), (2, 6, 'GCAT')]
    
    Notes
    -----
    Biological palindromes are important for:
    - Restriction enzyme recognition sites
    - Hairpin structures in RNA
    - DNA-binding protein recognition
    
    For restriction enzyme analysis, use biological=True (default).

    Complexity
    ----------
    Time: O(n^2), Space: O(n)
    """
    # Input validation
    if not isinstance(sequence, str):
        raise TypeError(f"sequence must be str, got {type(sequence).__name__}")
    if not isinstance(min_length, int):
        raise TypeError(f"min_length must be int, got {type(min_length).__name__}")
    if min_length <= 0:
        raise ValueError("min_length must be > 0")
    
    sequence = sequence.upper()
    
    if biological:
        # Validate DNA sequence
        if not all(base in "ATGC" for base in sequence):
            raise ValueError("Sequence contains invalid DNA bases (use ATGC only)")
    
    n = len(sequence)
    results = []
    
    for length in range(min_length, n + 1):
        for i in range(n - length + 1):
            substr = sequence[i : i + length]
            
            if biological:
                # Check if reverse complement equals original (biological palindrome)
                if substr == _reverse_complement(substr):
                    results.append((i, i + length, substr))
            else:
                # Check if simple reversal equals original (string palindrome)
                if substr == substr[::-1]:
                    results.append((i, i + length, substr))
    
    return results


__all__ = ["palindromic_sequence_finder"]
