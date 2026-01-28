"""Batch translation with multiprocessing support."""

from typing import Optional
from multiprocessing import Pool, cpu_count
from functools import partial

from .translate_dna_to_protein import translate_dna_to_protein
from .translate_dna_fast import translate_dna_fast, NUMBA_AVAILABLE


def translate_batch(
    sequences: list[str],
    table: str | int | dict[str, str] = "standard",
    n_processes: Optional[int] = None,
    use_fast: bool = False,
) -> list[str]:
    """
    Translate multiple sequences in parallel using multiprocessing.
    
    Parameters
    ----------
    sequences : list[str]
        List of DNA sequences to translate
    table : str, int, or dict[str, str], optional
        Genetic code table
    n_processes : int, optional
        Number of parallel processes (default: CPU count)
    use_fast : bool, optional
        Use Numba JIT implementation for 2.7x speedup (default: False)
        Requires: pip install numba
        
    Returns
    -------
    list[str]
        List of translated protein sequences
        
    Examples
    --------
    >>> seqs = ['ATGGCC', 'ATGAGA', 'TTTGGG']
    >>> translate_batch(seqs)
    ['MA', 'MR', 'FG']
    
    >>> # Use multiprocessing for many sequences
    >>> big_batch = ['ATGGCC' * 100] * 1000
    >>> proteins = translate_batch(big_batch, n_processes=4)
    """
    if not sequences:
        return []
    
    if n_processes is None:
        n_processes = min(cpu_count(), len(sequences))
    
    # For small batches, don't bother with multiprocessing
    if len(sequences) < 10 or n_processes == 1:
        translate_func = translate_dna_fast if (use_fast and NUMBA_AVAILABLE) else translate_dna_to_protein
        return [translate_func(seq, table=table) for seq in sequences]
    
    # Multiprocessing for large batches
    translate_func = translate_dna_fast if (use_fast and NUMBA_AVAILABLE) else translate_dna_to_protein
    # Use partial to bind the table parameter (works with custom dicts too)
    worker = partial(translate_func, table=table)
    
    with Pool(n_processes) as pool:
        return pool.map(worker, sequences)


__all__ = ["translate_batch"]
