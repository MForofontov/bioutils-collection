"""Independent parsing and numeric boundary checks for the sync audit."""

from pathlib import Path

import pytest
from bioutils_collection import blast_score_ratio, parse_bed

pytestmark = pytest.mark.unit


@pytest.mark.parametrize("chrom", ["track", "browser"])
def test_bed_metadata_keywords_can_be_chromosome_names(
    tmp_path: Path, chrom: str
) -> None:
    path = tmp_path / "keywords.bed"
    path.write_text(
        f'track name="example"\nbrowser position chr1:1-10\n{chrom}\t0\t10\tfeature\n',
        encoding="utf-8",
    )
    assert list(parse_bed(str(path))) == [
        {"chrom": chrom, "start": 0, "end": 10, "name": "feature"}
    ]


@pytest.mark.parametrize("score", [float("nan"), float("inf"), -float("inf")])
@pytest.mark.parametrize("argument", ["reference_score", "target_score"])
def test_blast_ratio_rejects_nonfinite_scores(score: float, argument: str) -> None:
    scores = {"reference_score": 100.0, "target_score": 50.0, argument: score}
    with pytest.raises(ValueError, match=argument):
        blast_score_ratio(**scores)


@pytest.mark.parametrize("target, expected", [(0, 0.0), (50, 0.5), (150, 1.0)])
def test_blast_ratio_preserves_finite_values(target: int, expected: float) -> None:
    assert blast_score_ratio(100, target) == expected
