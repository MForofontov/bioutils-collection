"""Convert annotations to GFF format."""

from collections.abc import Sequence
from typing import Any

from ._annotation_fields import annotation_fields


def annotation_to_gff(annotations: Sequence[dict[str, Any]]) -> list[str]:
    r"""
    Convert annotation records to GFF format strings.

    Parameters
    ----------
    annotations : Sequence[Dict[str, Any]]
        Existing writer records or native parser records. Coordinates remain
        one-based, inclusive. Existing writer field names take precedence.
        GFF attribute dictionaries retain their existing percent encoding.

    Returns
    -------
    list of str
        List of GFF format strings.

    Raises
    ------
    TypeError
        If input types are incorrect.
    KeyError
        If required keys are missing in any record.

    Examples
    --------
    >>> annotation_to_gff([
    ...     {'seqid': 'chr1', 'source': '.', 'type': 'gene', 'start': 1, 'end': 100, 'score': '.', 'strand': '+', 'phase': '.', 'attributes': 'ID=gene1'}
    ... ])
    ['chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1']
    """
    if not isinstance(annotations, (list, tuple)):
        raise TypeError(
            f"annotations must be a list or tuple, got {type(annotations).__name__}"
        )
    return [
        "\t".join(map(str, annotation_fields(record, "GFF"))) for record in annotations
    ]


__all__ = ["annotation_to_gff"]
