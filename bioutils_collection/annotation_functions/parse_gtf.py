"""Parse GTF format files."""

from collections.abc import Iterator
from typing import Any


def parse_gtf(gtf_file: str) -> Iterator[dict[str, Any]]:
    """
    Parse a GTF (Gene Transfer Format) file and yield annotation records as dictionaries.

    Parameters
    ----------
    gtf_file : str
        Path to the GTF file.

    Returns
    -------
    Iterator[Dict[str, Any]]
        Iterator over parsed annotation records.

    Raises
    ------
    TypeError
        If gtf_file is not a string.
    ValueError
        If a non-comment record has malformed columns or coordinates.
    FileNotFoundError
        If the file does not exist.

    Examples
    --------
    >>> for record in parse_gtf('example.gtf'):
    ...     print(record)
    {'seqname': 'chr1', 'source': 'ENSEMBL', ...}

    Notes
    -----
    Coordinates are one-based and inclusive; the attribute field is retained as text.
    Malformed records raise ValueError rather than being silently discarded.

    Complexity
    ----------
    Time: O(n), Space: O(1)
    """
    if not isinstance(gtf_file, str):
        raise TypeError(f"gtf_file must be a string, got {type(gtf_file).__name__}")
    try:
        with open(gtf_file, encoding="utf-8") as f:
            for line_number, line in enumerate(f, 1):
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) != 9:
                    raise ValueError(
                        f"Invalid GTF line {line_number}: expected 9 columns"
                    )
                record = {
                    "seqname": parts[0],
                    "source": parts[1],
                    "feature": parts[2],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "score": parts[5],
                    "strand": parts[6],
                    "frame": parts[7],
                    "attribute": parts[8],
                }
                yield record
    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found: {gtf_file}") from e


__all__ = ["parse_gtf"]
