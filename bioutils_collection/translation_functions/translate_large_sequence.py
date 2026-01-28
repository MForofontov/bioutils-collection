"""Parallel translation of a single large sequence by splitting into chunks."""

from typing import Optional
from multiprocessing import cpu_count
from .translate_batch import translate_batch


def translate_large_sequence(
    sequence: str,
    table: str | int | dict[str, str] = "standard",
    chunk_size: int = 1_000_000,
    n_processes: Optional[int] = None,
    use_fast: bool = True,
) -> str:
    """
    Translate a single large DNA sequence using parallel processing.
    
    Splits the sequence into chunks and processes them in parallel,
    then combines the results. Useful for very large sequences (>10MB).
    
    Parameters
    ----------
    sequence : str
        DNA sequence to translate (must be multiple of 3)
    table : str, int, or dict[str, str], optional
        Genetic code table (default: "standard")
    chunk_size : int, optional
        Size of each chunk in bases (default: 1,000,000)
        Will be adjusted to nearest multiple of 3
    n_processes : int, optional
        Number of parallel processes (default: CPU count)
    use_fast : bool, optional
        Use Numba JIT implementation (default: True)
        
    Returns
    -------
    str
        Translated protein sequence
        
    Examples
    --------
    >>> # Translate a 50MB sequence using parallel processing
    >>> large_seq = 'ATG' + 'GCC' * 16_666_666 + 'TAA'
    >>> protein = translate_large_sequence(large_seq, chunk_size=5_000_000)
    
    >>> # Use specific number of processes
    >>> protein = translate_large_sequence(large_seq, n_processes=8, use_fast=True)
    
    Notes
    -----
    - Recommended for sequences larger than 10MB
    - chunk_size will be adjusted to maintain codon boundaries (multiple of 3)
    - For smaller sequences, use translate_dna_to_protein() or translate_dna_fast()
    - Set use_fast=True (default) for best performance with large chunks
    """
    if not sequence:
        return ""
    
    seq_len = len(sequence)
    
    # Validate sequence length
    if seq_len % 3 != 0:
        raise ValueError("Sequence length must be a multiple of 3")
    
    # Adjust chunk_size to be multiple of 3
    chunk_size = (chunk_size // 3) * 3
    
    # If sequence is smaller than chunk_size, process directly
    if seq_len <= chunk_size:
        from .translate_dna_to_protein import translate_dna_to_protein
        from .translate_dna_fast import translate_dna_fast, NUMBA_AVAILABLE
        
        translate_func = translate_dna_fast if (use_fast and NUMBA_AVAILABLE) else translate_dna_to_protein
        return translate_func(sequence, table=table)
    
    # Split sequence into chunks
    chunks = []
    for i in range(0, seq_len, chunk_size):
        chunks.append(sequence[i:i + chunk_size])
    
    # Process chunks in parallel
    if n_processes is None:
        n_processes = min(cpu_count(), len(chunks))
    
    protein_chunks = translate_batch(
        chunks,
        table=table,
        n_processes=n_processes,
        use_fast=use_fast,
    )
    
    # Combine results
    return ''.join(protein_chunks)


__all__ = ["translate_large_sequence"]
