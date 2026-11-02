"""FASTA file translation with in-memory result storage."""

from pathlib import Path

from bioutils_collection.fasta_misc.fasta_parser import _parse_fasta_lines

from ._atomic_output import atomic_text_output
from .genetic_code_tables import get_codon_table
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
        Use the Numba JIT implementation (default: False)
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
        If table name/number is unrecognized, FASTA headers are malformed,
        or input and output refer to the same file
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
    - Empty sequence records are skipped; IDs use the first header word.
    - Duplicate IDs overwrite earlier entries in the dictionary API; streaming
      preserves each record.
    - Output is replaced atomically on success; failed processing preserves it.
    - Set use_fast=True to use the Numba backend on large sequences
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

    if (
        output_file is not None
        and Path(output_file).exists()
        and fasta_file.samefile(output_file)
    ):
        raise ValueError("Input and output must refer to different files")

    # Configuration errors must not be converted into per-record errors or
    # overwrite an existing output, including when the input file is empty.
    table = get_codon_table(table)
    translate_func = translate_dna_fast if use_fast else translate_dna_to_protein

    results = {}
    with open(fasta_file, encoding="utf-8") as source:
        for header, seq in _parse_fasta_lines(source):
            if not seq:
                continue
            seq_id = header.split()[0]
            try:
                results[seq_id] = translate_func(seq, table=table)
            except ValueError as e:
                results[seq_id] = f"ERROR: {e}"

    if output_file is not None:
        with atomic_text_output(Path(output_file)) as output:
            for seq_id, protein in results.items():
                output.write(f">{seq_id}\n")
                for i in range(0, len(protein), 60):
                    output.write(f"{protein[i : i + 60]}\n")
    return results


__all__ = ["translate_fasta"]
