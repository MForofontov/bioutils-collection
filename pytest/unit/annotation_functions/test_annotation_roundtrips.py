"""Parser/writer interoperability with native coordinate conventions."""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from bioutils_collection.annotation_functions.annotation_statistics import (
    annotation_statistics,
)
from bioutils_collection.annotation_functions.annotation_to_bed import annotation_to_bed
from bioutils_collection.annotation_functions.annotation_to_gff import annotation_to_gff
from bioutils_collection.annotation_functions.annotation_to_gtf import annotation_to_gtf
from bioutils_collection.annotation_functions.filter_annotations import (
    filter_annotations,
)
from bioutils_collection.annotation_functions.parse_bed import parse_bed
from bioutils_collection.annotation_functions.parse_gff import parse_gff
from bioutils_collection.annotation_functions.parse_gtf import parse_gtf
from bioutils_collection.annotation_functions.parse_vcf import parse_vcf


def test_gff_roundtrip(tmp_path: Path) -> None:
    path = tmp_path / "input.gff"
    line = "chr1\t.\tgene\t1\t1\t.\t+\t.\tID=a;Note=a%3Bb;Parent=p1,p2"
    path.write_text(line + "\n", encoding="utf-8")
    records = list(parse_gff(str(path)))
    assert annotation_to_gff(records) == [line]
    assert annotation_to_bed(records) == ["chr1\t0\t1\tgene"]
    assert (
        list(filter_annotations(records, feature_type="gene", chrom="chr1")) == records
    )
    assert annotation_statistics(records) == {
        "feature_counts": {"gene": 1},
        "total_coverage": 1,
    }


def test_gtf_roundtrip(tmp_path: Path) -> None:
    path = tmp_path / "input.gtf"
    line = 'chr1\t.\texon\t1\t1\t.\t-\t.\tgene_id "g"; transcript_id "t";'
    path.write_text(line + "\n", encoding="utf-8")
    assert annotation_to_gtf(list(parse_gtf(str(path)))) == [line]


def test_gff_blank_lines_and_embedded_fasta(tmp_path: Path) -> None:
    path = tmp_path / "input.gff"
    path.write_text("\nchr1\t.\tgene\t1\t1\t.\t+\t.\t.\n##FASTA\n>chr1\nA\n")
    assert len(list(parse_gff(str(path)))) == 1


def test_bed_headers_and_native_coordinates(tmp_path: Path) -> None:
    path = tmp_path / "input.bed"
    path.write_text('track name="test"\nbrowser position chr1:1-1\n\nchr1\t0\t1\tx\n')
    assert list(parse_bed(str(path))) == [
        {"chrom": "chr1", "start": 0, "end": 1, "name": "x"}
    ]


@pytest.mark.parametrize("writer", [annotation_to_gff, annotation_to_gtf])
def test_writers_reject_nonrecord(writer: Callable[..., Any]) -> None:
    with pytest.raises(TypeError, match="record must be a dict"):
        writer(["bad"])


@pytest.mark.parametrize("parser", [parse_gtf, parse_vcf])
def test_malformed_records_are_not_silently_discarded(
    tmp_path: Path, parser: Callable[..., Any]
) -> None:
    path = tmp_path / "malformed"
    path.write_text("# comment\nmalformed_record\n")
    with pytest.raises(ValueError, match="line 2"):
        list(parser(str(path)))
