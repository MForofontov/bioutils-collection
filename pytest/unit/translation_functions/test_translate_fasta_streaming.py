"""Regression tests for streaming translation and preservation of input files."""

from pathlib import Path

import pytest
from bioutils_collection.translation_functions import translate_fasta_streaming

pytestmark = [pytest.mark.unit, pytest.mark.translation]


def test_streaming_rejects_same_file(tmp_path: Path) -> None:
    source = tmp_path / "input.fa"
    original = ">first\nATGGCC\n"
    source.write_text(original)

    with pytest.raises(ValueError, match="different files"):
        translate_fasta_streaming(source, source)

    assert source.read_text() == original


def test_streaming_rejects_hard_link(tmp_path: Path) -> None:
    source = tmp_path / "input.fa"
    original = ">first\nATGGCC\n"
    source.write_text(original)
    destination = tmp_path / "alias.fa"
    destination.hardlink_to(source)

    with pytest.raises(ValueError, match="different files"):
        translate_fasta_streaming(source, destination)

    assert source.read_text() == original


def test_streaming_invalid_table_preserves_output(tmp_path: Path) -> None:
    source = tmp_path / "input.fa"
    source.write_text(">first\nATGGCC\n")
    destination = tmp_path / "output.fa"
    destination.write_text("previous results\n")

    with pytest.raises(ValueError, match="Unknown table"):
        translate_fasta_streaming(source, destination, table="invalid")

    assert destination.read_text() == "previous results\n"


@pytest.mark.parametrize("use_fast", [False, True])
def test_streaming_translates_multiline_records(tmp_path: Path, use_fast: bool) -> None:
    source = tmp_path / "input.fa"
    source.write_text(">first description\nATG\nGCC\n\n>second\nTAA\n")
    destination = tmp_path / "output.fa"

    assert translate_fasta_streaming(source, destination, use_fast=use_fast) == 2
    assert destination.read_text() == ">first\nMA\n>second\n*\n"
