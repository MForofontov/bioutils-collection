"""Parse GFF format files."""

from collections.abc import Iterator
from typing import Any


def parse_gff(gff_file: str) -> Iterator[dict[str, Any]]:
    r"""
    Parse GFF records with one-based inclusive coordinates and encoded attributes.

    Parameters
    ----------
    gff_file : str
        Path to GFF formatted file.

    Yields
    ------
    Dict[str, Any]
        Annotation dictionary for each feature.

    Raises
    ------
    TypeError
        If gff_file is not a string.
    FileNotFoundError
        If the file does not exist.
    ValueError
        If GFF format is invalid.

    Examples
    --------
    >>> from pathlib import Path
    >>> from tempfile import TemporaryDirectory
    >>> with TemporaryDirectory() as directory:
    ...     path = Path(directory) / 'example.gff'
    ...     _ = path.write_text('chr1\tRefSeq\tgene\t1\t100\t.\t+\t.\tID=gene1\n', encoding='utf-8')
    ...     records = list(parse_gff(str(path)))
    >>> records[0]['attributes']
    {'ID': 'gene1'}
    """
    if not isinstance(gff_file, str):
        raise TypeError("gff_file must be str")
    try:
        with open(gff_file) as f:
            for line in f:
                if line.rstrip() == "##FASTA":
                    break
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) != 9:
                    raise ValueError("Invalid GFF line: " + line)
                attr_dict = {}
                for attr in parts[8].split(";"):
                    if "=" in attr:
                        k, v = attr.split("=", 1)
                        attr_dict[k] = v
                yield {
                    "seqid": parts[0],
                    "source": parts[1],
                    "type": parts[2],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "score": parts[5],
                    "strand": parts[6],
                    "phase": parts[7],
                    "attributes": attr_dict,
                }
    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found: {gff_file}") from e


__all__ = ["parse_gff"]
