"""Parse BED format files."""

from collections.abc import Iterator
from typing import Any


def parse_bed(bed_file: str) -> Iterator[dict[str, Any]]:
    r"""
    Parse BED records, preserving zero-based, half-open start/end coordinates.

    Parameters
    ----------
    bed_file : str
        Path to BED formatted file.

    Yields
    ------
    Dict[str, Any]
        Annotation dictionary for each feature.

    Raises
    ------
    TypeError
        If bed_file is not a string.
    FileNotFoundError
        If the file does not exist.
    ValueError
        If BED format is invalid.

    Examples
    --------
    >>> from pathlib import Path
    >>> from tempfile import TemporaryDirectory
    >>> with TemporaryDirectory() as directory:
    ...     path = Path(directory) / 'example.bed'
    ...     _ = path.write_text('chr1\t0\t100\tfeature1\n', encoding='utf-8')
    ...     records = list(parse_bed(str(path)))
    >>> records
    [{'chrom': 'chr1', 'start': 0, 'end': 100, 'name': 'feature1'}]

    Notes
    -----
    Blank lines, comments and track/browser metadata are skipped. Tab-delimited
    records with integer coordinates are retained even when their chromosome
    is named track or browser.
    """
    if not isinstance(bed_file, str):
        raise TypeError("bed_file must be str")
    try:
        with open(bed_file) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if line.split(maxsplit=1)[0] in ("track", "browser"):
                    # Metadata keywords are also legal chromosome names.
                    try:
                        int(parts[1])
                        int(parts[2])
                    except (IndexError, ValueError):
                        continue
                if len(parts) < 3:
                    raise ValueError("Invalid BED line: " + line)
                ann = {
                    "chrom": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                }
                if len(parts) > 3:
                    ann["name"] = parts[3]
                yield ann
    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found: {bed_file}") from e


__all__ = ["parse_bed"]
