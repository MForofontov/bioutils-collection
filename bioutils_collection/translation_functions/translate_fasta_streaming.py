"""Memory-efficient streaming FASTA translation for huge files."""

from pathlib import Path

from .translate_dna_fast import translate_dna_fast
from .translate_dna_to_protein import translate_dna_to_protein


def translate_fasta_streaming(
    fasta_file: str | Path,
    output_file: str | Path,
    table: str | int | dict[str, str] = "standard",
    use_fast: bool = False,
) -> int:
    """
    Translate FASTA file in streaming mode (memory-efficient for huge files).

    Parameters
    ----------
    fasta_file : str or Path
        Input FASTA file path
    output_file : str or Path
        Output FASTA file path
    table : str, int, or dict[str, str], optional
        Genetic code table
    use_fast : bool, optional
        Use Numba JIT implementation for 2.7x speedup (default: False)
        Requires: pip install numba

    Returns
    -------
    int
        Number of sequences translated

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
    >>> n = translate_fasta_streaming('huge.fasta', 'proteins.fasta')
    >>> print(f"Translated {n} sequences")

    Notes
    -----
    - Memory-efficient: processes one sequence at a time
    - Invalid sequences (not multiple of 3) will have ERROR in output
    - Recommended for files larger than available RAM
    - Set use_fast=True for 2.7x speedup on large sequences
    """
    # Input validation
    if not isinstance(fasta_file, (str, Path)):
        raise TypeError(
            f"fasta_file must be str or Path, got {type(fasta_file).__name__}"
        )

    if not isinstance(output_file, (str, Path)):
        raise TypeError(
            f"output_file must be str or Path, got {type(output_file).__name__}"
        )

    fasta_file = Path(fasta_file)
    output_file = Path(output_file)

    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    translate_func = translate_dna_fast if use_fast else translate_dna_to_protein

    count = 0
    current_id = None
    current_seq: list[str] = []

    with open(fasta_file) as fin, open(output_file, "w") as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Translate and write previous sequence
                if current_id is not None and current_seq:
                    seq = "".join(current_seq)
                    try:
                        protein = translate_func(seq, table=table)
                        fout.write(f">{current_id}\n")
                        for i in range(0, len(protein), 60):
                            fout.write(f"{protein[i : i + 60]}\n")
                        count += 1
                    except ValueError as e:
                        fout.write(f">{current_id}\n")
                        fout.write(f"ERROR: {e}\n")

                # Start new sequence
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Translate and write last sequence
        if current_id is not None and current_seq:
            seq = "".join(current_seq)
            try:
                protein = translate_func(seq, table=table)
                fout.write(f">{current_id}\n")
                for i in range(0, len(protein), 60):
                    fout.write(f"{protein[i : i + 60]}\n")
                count += 1
            except ValueError as e:
                fout.write(f">{current_id}\n")
                fout.write(f"ERROR: {e}\n")

    return count


__all__ = ["translate_fasta_streaming"]
