"""Parse FASTA format sequences."""

from collections.abc import Iterable, Iterator


def parse_fasta(fasta_str: str) -> Iterator[tuple[str, str]]:
    """
    Parse a FASTA formatted string and yield (header, sequence) tuples.

    Parameters
    ----------
    fasta_str : str
        FASTA formatted string.

    Yields
    ------
    Tuple[str, str]
        (header, sequence) pairs.

    Raises
    ------
    TypeError
        If fasta_str is not a string.
    ValueError
        If FASTA format is invalid.

    Examples
    --------
    >>> list(parse_fasta('>seq1\nATGC\n>seq2\nGGTT'))
    [('seq1', 'ATGC'), ('seq2', 'GGTT')]
    """
    if not isinstance(fasta_str, str):
        raise TypeError(f"fasta_str must be str, got {type(fasta_str).__name__}")
    yield from _parse_fasta_lines(fasta_str.splitlines())


def _parse_fasta_lines(lines: Iterable[str]) -> Iterator[tuple[str, str]]:
    """Parse text lines without retaining the entire file; reject missing IDs."""
    header = None
    seq: list[str] = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield (header, "".join(seq))
            header = line[1:].strip()
            if not header:
                raise ValueError("FASTA format error: empty header")
            seq = []
        else:
            if header is None:
                raise ValueError("FASTA format error: sequence before header")
            seq.append(line)
    if header is not None:
        yield (header, "".join(seq))


__all__ = ["parse_fasta"]
