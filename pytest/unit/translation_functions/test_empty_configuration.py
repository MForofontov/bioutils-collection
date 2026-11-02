"""Configuration errors must be visible even when there are no records."""

from collections.abc import Callable
from typing import Any

import pytest
from bioutils_collection.translation_functions import (
    translate_batch,
    translate_dna_fast,
    translate_dna_to_protein,
    translate_large_sequence,
)

TRANSLATORS = [translate_dna_to_protein, translate_dna_fast, translate_large_sequence]


@pytest.mark.parametrize("translate", TRANSLATORS)
@pytest.mark.parametrize(
    "table,exception", [("invalid", ValueError), (999, ValueError), (None, TypeError)]
)
def test_empty_scalar_validates_table(
    translate: Callable[..., str], table: Any, exception: type[Exception]
) -> None:
    with pytest.raises(exception):
        translate("", table=table)


@pytest.mark.parametrize("fast", [False, True])
@pytest.mark.parametrize(
    "table,exception", [("invalid", ValueError), (999, ValueError), (None, TypeError)]
)
def test_empty_batch_validates_table(
    fast: bool, table: Any, exception: type[Exception]
) -> None:
    with pytest.raises(exception):
        translate_batch([], table=table, use_fast=fast)


@pytest.mark.parametrize(
    "translate,empty", [(translate_batch, []), (translate_large_sequence, "")]
)
@pytest.mark.parametrize(
    "processes,exception", [(0, ValueError), (-1, ValueError), (1.5, TypeError)]
)
def test_empty_validates_processes(
    translate: Callable[..., Any],
    empty: Any,
    processes: Any,
    exception: type[Exception],
) -> None:
    with pytest.raises(exception):
        translate(empty, n_processes=processes)


@pytest.mark.parametrize("chunk", [1, 2])
def test_empty_validates_adjusted_chunk(chunk: int) -> None:
    with pytest.raises(ValueError, match="at least 3"):
        translate_large_sequence("", chunk_size=chunk)


@pytest.mark.parametrize("table", ["standard", 2, {}, {"ATG": "start"}])
def test_valid_empty_inputs_remain_empty(table: Any) -> None:
    assert translate_dna_fast("", table) == ""
    assert translate_dna_to_protein("", table) == ""
    assert translate_large_sequence("", table, chunk_size=3) == ""
    assert translate_batch([], table, n_processes=1) == []
