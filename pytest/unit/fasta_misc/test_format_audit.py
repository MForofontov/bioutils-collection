import pytest
from bioutils_collection import fastq_to_fasta, parse_fasta


@pytest.mark.parametrize("text", [">\nATG\n", ">   \nATG\n", "ATG\n>a\nATG\n"])
def test_invalid_fasta_headers_are_explicit_errors(text: str) -> None:
    with pytest.raises(ValueError, match="header"):
        list(parse_fasta(text))


@pytest.mark.parametrize(
    "text",
    ["@\nAC\n+\nII\n", "@a\nAC\n+b\nII\n", "@a\nAC\n+\nI\x7f\n", "@a\nA\n+\nI \n"],
)
def test_malformed_fastq_is_rejected(text: str) -> None:
    with pytest.raises(ValueError):
        fastq_to_fasta(text)


def test_empty_and_duplicate_fasta_records_preserved() -> None:
    assert list(parse_fasta(">a\n>a\nATG\n")) == [("a", ""), ("a", "ATG")]
