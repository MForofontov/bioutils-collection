"""FASTA file translation with in-memory result storage."""

from pathlib import Path

from .translate_dna_fast import translate_dna_fast
from .translate_dna_to_protein import translate_dna_to_protein


def translate_fasta(
    fasta_file: str | Path,
    output_file: str | Path | None = None,
    table: str | int | dict[str, str] = "standard",
    use_fast: bool = False,
) -> dict[str, str]:
    """
    Translate sequences from a FASTA file.

    Parameters
    ----------
    fasta_file : str or Path
        Input FASTA file path
    output_file : str or Path, optional
        Output FASTA file path (if None, only returns dict)
    table : str, int, or dict[str, str], optional
        Genetic code table
    use_fast : bool, optional
        Use Numba JIT implementation for 2.7x speedup (default: False)
        Requires: pip install numba

    Returns
    -------
    dict[str, str]
        Dictionary mapping sequence IDs to protein sequences

    Raises
    ------
    TypeError
        If fasta_file or output_file are not str or Path
    FileNotFoundError
        If input FASTA file does not exist
    ValueError
        If table name/number is unrecognized
    OSError
        If unable to read input file or write output file

    Examples
    --------
    >>> proteins = translate_fasta('sequences.fasta', 'proteins.fasta')
    >>> len(proteins)
    42

    Notes
    -----
    - Invalid sequences (not multiple of 3) will have ERROR message in results
    - For huge files that don't fit in memory, use translate_fasta_streaming()
    - Set use_fast=True for 2.7x speedup on large sequences
    """
    # Input validation
    if not isinstance(fasta_file, (str, Path)):
        raise TypeError(
            f"fasta_file must be str or Path, got {type(fasta_file).__name__}"
        )

    if output_file is not None and not isinstance(output_file, (str, Path)):
        raise TypeError(
            f"output_file must be str or Path, got {type(output_file).__name__}"
        )

    fasta_file = Path(fasta_file)

    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    translate_func = translate_dna_fast if use_fast else translate_dna_to_protein

    results = {}
    current_id = None
    current_seq: list[str] = []

    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save previous sequence
                if current_id is not None and current_seq:
                    seq = "".join(current_seq)
                    try:
                        results[current_id] = translate_func(seq, table=table)
                    except ValueError as e:
                        results[current_id] = f"ERROR: {e}"

                # Start new sequence
                current_id = line[1:].split()[0]  # Take first word after >
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id is not None and current_seq:
            seq = "".join(current_seq)
            try:
                results[current_id] = translate_func(seq, table=table)
            except ValueError as e:
                results[current_id] = f"ERROR: {e}"

    # Write output file if requested
    if output_file is not None:
        output_file = Path(output_file)
        with open(output_file, "w") as f:
            for seq_id, protein in results.items():
                f.write(f">{seq_id}\n")
                # Write in 60-character lines
                for i in range(0, len(protein), 60):
                    f.write(f"{protein[i : i + 60]}\n")

    return results


__all__ = ["translate_fasta"]
