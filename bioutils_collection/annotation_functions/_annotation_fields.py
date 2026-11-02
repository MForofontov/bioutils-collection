"""Accept existing writer fields and native parser records without mutation."""

from typing import Any


def annotation_fields(record: dict[str, Any], format_name: str) -> list[Any]:
    if not isinstance(record, dict):
        raise TypeError("record must be a dict")
    fields = []
    for key, alias in (
        ("seqid", "seqname"),
        ("source", "source"),
        ("feature", "type"),
        ("start", "start"),
        ("end", "end"),
        ("score", "score"),
        ("strand", "strand"),
        ("frame", "phase"),
        ("attribute", "attributes"),
    ):
        if key in record:
            value = record[key]
        elif alias in record:
            value = record[alias]
        else:
            raise KeyError(f"Missing required {format_name} key: {key}")
        if key == "attribute" and isinstance(value, dict):
            if format_name == "GFF":
                # parse_gff retains GFF3 percent encoding, including list separators.
                value = ";".join(f"{k}={v}" for k, v in value.items()) or "."
            else:
                raise TypeError("GTF attribute must be a string")
        fields.append(value)
    return fields
