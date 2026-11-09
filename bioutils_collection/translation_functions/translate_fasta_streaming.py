"""Memory-efficient streaming FASTA translation for huge files."""

from pathlib import Path

from bioutils_collection.fasta_misc.fasta_parser import _parse_fasta_lines

from ._atomic_output import atomic_text_output
from .genetic_code_tables import get_codon_table
from .translate_dna_fast import translate_dna_fast
from .translate_dna_to_protein import translate_dna_to_protein


def translate_fasta_streaming(
    fasta_file: str | Path,
    output_file: str | Path,
    table: str | int | dict[str, str] = "standard",
    use_fast: bool = False,
) -> int:
    r"""
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
        Use the Numba JIT implementation (default: False)
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
        If input and output refer to the same file
    OSError
        If unable to read input file or write output file

    Examples
    --------
    >>> from tempfile import TemporaryDirectory
    >>> with TemporaryDirectory() as directory:
    ...     source = Path(directory) / 'sequences.fasta'
    ...     output = Path(directory) / 'proteins.fasta'
    ...     _ = source.write_text('>example\nATGGCCTAA\n', encoding='utf-8')
    ...     n = translate_fasta_streaming(source, output)
    >>> n
    1

    Notes
    -----
    - Memory-efficient: processes one sequence at a time
    - Invalid sequences (not multiple of 3) will have ERROR in output
    - Recommended for files larger than available RAM
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

    if not isinstance(output_file, (str, Path)):
        raise TypeError(
            f"output_file must be str or Path, got {type(output_file).__name__}"
        )

    fasta_file = Path(fasta_file)
    output_file = Path(output_file)

    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    if output_file.exists() and fasta_file.samefile(output_file):
        raise ValueError("Input and output must refer to different files")

    # Reject invalid configuration before opening (and truncating) the output.
    table = get_codon_table(table)
    translate_func = translate_dna_fast if use_fast else translate_dna_to_protein

    count = 0
    with (
        open(fasta_file, encoding="utf-8") as source,
        atomic_text_output(output_file) as output,
    ):
        for header, seq in _parse_fasta_lines(source):
            if not seq:
                continue
            seq_id = header.split()[0]
            try:
                protein = translate_func(seq, table=table)
            except ValueError as e:
                output.write(f">{seq_id}\nERROR: {e}\n")
                continue
            output.write(f">{seq_id}\n")
            for i in range(0, len(protein), 60):
                output.write(f"{protein[i : i + 60]}\n")
            count += 1
    return count


__all__ = ["translate_fasta_streaming"]
