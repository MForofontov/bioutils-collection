"""Find open reading frames (ORFs) in sequences."""

from collections.abc import Iterator


def find_orfs(seq: str) -> Iterator[tuple[int, int, str]]:
    """
    Find open reading frames (ORFs) in a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence (A, T, C, G).

    Yields
    ------
    Tuple[int, int, str]
        (start, end, orf_sequence) for each ORF found.
        ORFs are yielded from 5' to 3' direction.
        Includes ORFs that reach the end without a stop codon.

    Raises
    ------
    TypeError
        If seq is not a string.
    ValueError
        If seq contains invalid characters.

    Examples
    --------
    >>> list(find_orfs('ATGAAATAGATGTAA'))
    [(0, 9, 'ATGAAATAG'), (9, 15, 'ATGTAA')]
    
    >>> # ORF without stop codon (reaches end)
    >>> list(find_orfs('ATGGCCAAA'))
    [(0, 9, 'ATGGCCAAA')]
    
    Notes
    -----
    - Start codon: ATG
    - Stop codons: TAA, TAG, TGA
    - ORFs without stop codons (reaching end) are included
    - Only searches in the 5' to 3' direction (use with reverse_complement for other strand)
    - Returns ORFs in forward reading frame only
    """
    if not isinstance(seq, str):
        raise TypeError(f"seq must be str, got {type(seq).__name__}")
    seq = seq.upper()
    if not all(base in "ATCG" for base in seq):
        raise ValueError("Sequence contains invalid DNA bases")
    
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    i = 0
    
    while i <= len(seq) - 3:  # Fixed: was len(seq) - 2
        if seq[i : i + 3] == start_codon:
            # Search for stop codon
            found_stop = False
            for j in range(i + 3, len(seq), 3):  # Fixed: removed - 2 to scan to end
                if j + 3 > len(seq):
                    # Incomplete codon at end - yield ORF without stop
                    yield (i, len(seq), seq[i:])
                    i = len(seq)
                    found_stop = True
                    break
                    
                codon = seq[j : j + 3]
                if codon in stop_codons:
                    yield (i, j + 3, seq[i : j + 3])
                    i = j + 3
                    found_stop = True
                    break
            
            if not found_stop:
                # ORF reaches end without stop codon
                yield (i, len(seq), seq[i:])
                i = len(seq)
        else:
            i += 1


__all__ = ["find_orfs"]
