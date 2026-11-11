"""Batch translation with multiprocessing support."""

from functools import partial
from multiprocessing import Pool, cpu_count

from .translate_dna_fast import translate_dna_fast
from .translate_dna_to_protein import translate_dna_to_protein


def translate_batch(
    sequences: list[str],
    table: str | int | dict[str, str] = "standard",
    n_processes: int | None = None,
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

    Raises
    ------
    TypeError
        If sequences is not a list or contains non-string elements
    ValueError
        If any sequence length is not a multiple of 3
        If table name/number is unrecognized
        If n_processes is negative

    Examples
    --------
    >>> seqs = ['ATGGCC', 'ATGAGA', 'TTTGGG']
    >>> translate_batch(seqs)
    ['MA', 'MR', 'FG']

    >>> # Use multiprocessing for many sequences
    >>> big_batch = ['ATGGCC' * 100] * 1000
    >>> proteins = translate_batch(big_batch, n_processes=4)
    """
    # Input validation
    if not isinstance(sequences, list):
        raise TypeError(f"sequences must be list, got {type(sequences).__name__}")

    if not sequences:
        return []

    # Validate all sequences are strings
    for i, seq in enumerate(sequences):
        if not isinstance(seq, str):
            raise TypeError(f"sequences[{i}] must be str, got {type(seq).__name__}")

    if n_processes is not None:
        if not isinstance(n_processes, int):
            raise TypeError(
                f"n_processes must be int, got {type(n_processes).__name__}"
            )
        if n_processes < 1:
            raise ValueError(f"n_processes must be positive, got {n_processes}")

    if n_processes is None:
        n_processes = min(cpu_count(), len(sequences))

    # For small batches, don't bother with multiprocessing
    if len(sequences) < 10 or n_processes == 1:
        translate_func = translate_dna_fast if use_fast else translate_dna_to_protein
        return [translate_func(seq, table=table) for seq in sequences]

    # Multiprocessing for large batches
    translate_func = translate_dna_fast if use_fast else translate_dna_to_protein
    # Use partial to bind the table parameter (works with custom dicts too)
    worker = partial(translate_func, table=table)

    with Pool(n_processes) as pool:
        return pool.map(worker, sequences)


__all__ = ["translate_batch"]
