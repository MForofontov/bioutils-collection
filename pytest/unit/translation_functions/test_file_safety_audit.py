"""Output transactions must preserve files on failure."""

import importlib
import os
import stat
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from bioutils_collection.translation_functions.translate_fasta import translate_fasta
from bioutils_collection.translation_functions.translate_fasta_streaming import (
    translate_fasta_streaming,
)


@pytest.mark.parametrize("translate", [translate_fasta, translate_fasta_streaming])
@pytest.mark.parametrize("alias", [False, True])
def test_input_output_alias(
    tmp_path: Path, translate: Callable[..., Any], alias: bool
) -> None:
    source = tmp_path / "input.fa"
    source.write_text(">a\nATG\n")
    output = source
    if alias:
        output = tmp_path / "alias.fa"
        os.link(source, output)
    with pytest.raises(ValueError, match="different files"):
        translate(source, output)
    assert source.read_text() == ">a\nATG\n"


@pytest.mark.parametrize("translate", [translate_fasta, translate_fasta_streaming])
@pytest.mark.parametrize("header", [">", ">   "])
def test_malformed_header_preserves_output(
    tmp_path: Path, translate: Callable[..., Any], header: str
) -> None:
    source, output = tmp_path / "input.fa", tmp_path / "output.fa"
    source.write_text(f">valid\nATG\n{header}\nATG\n")
    output.write_text("original")
    with pytest.raises(ValueError, match="header"):
        translate(source, output)
    assert output.read_text() == "original"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["input.fa", "output.fa"]


@pytest.mark.parametrize("translate", [translate_fasta, translate_fasta_streaming])
def test_interrupted_translation_preserves_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, translate: Callable[..., Any]
) -> None:
    source, output = tmp_path / "input.fa", tmp_path / "output.fa"
    source.write_text(">a\nATG\n>b\nATG\n")
    output.write_text("original")
    module = importlib.import_module(translate.__module__)
    calls = 0

    def interrupted(*args: object, **kwargs: object) -> str:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("interrupted")
        return "M"

    monkeypatch.setattr(module, "translate_dna_to_protein", interrupted)
    with pytest.raises(OSError, match="interrupted"):
        translate(source, output)
    assert output.read_text() == "original"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["input.fa", "output.fa"]


@pytest.mark.parametrize("failure", ["fsync", "replace"])
@pytest.mark.parametrize("translate", [translate_fasta, translate_fasta_streaming])
def test_commit_failure_preserves_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure: str,
    translate: Callable[..., Any],
) -> None:
    source, output = tmp_path / "input.fa", tmp_path / "output.fa"
    source.write_text(">a\nATG\n")
    output.write_text("original")
    module = importlib.import_module(
        "bioutils_collection.translation_functions._atomic_output"
    )

    def fail(*args: object) -> None:
        raise OSError("commit failed")

    monkeypatch.setattr(module.os, failure, fail)
    with pytest.raises(OSError, match="commit failed"):
        translate(source, output)
    assert output.read_text() == "original"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["input.fa", "output.fa"]


@pytest.mark.parametrize("translate", [translate_fasta, translate_fasta_streaming])
def test_symlink_alias_is_rejected(
    tmp_path: Path, translate: Callable[..., Any]
) -> None:
    source, alias = tmp_path / "input.fa", tmp_path / "alias.fa"
    source.write_text(">a\nATG\n")
    try:
        alias.symlink_to(source)
    except OSError as error:
        if os.name != "nt":
            raise
        pytest.skip(f"Symlink creation is unavailable: {error}")
    with pytest.raises(ValueError, match="different files"):
        translate(source, alias)
    assert source.read_text() == ">a\nATG\n"


@pytest.mark.parametrize("translate", [translate_fasta, translate_fasta_streaming])
def test_output_symlink_is_preserved(
    tmp_path: Path, translate: Callable[..., Any]
) -> None:
    source, target, alias = (
        tmp_path / "input.fa",
        tmp_path / "target.fa",
        tmp_path / "alias.fa",
    )
    source.write_text(">a\nATG\n")
    target.write_text("old")
    try:
        alias.symlink_to(target)
    except OSError as error:
        if os.name != "nt":
            raise
        pytest.skip(f"Symlink creation is unavailable: {error}")
    translate(source, alias)
    assert alias.is_symlink()
    assert target.read_text() == ">a\nM\n"


@pytest.mark.parametrize("translate", [translate_fasta, translate_fasta_streaming])
@pytest.mark.parametrize("failure", ["write", "close"])
@pytest.mark.parametrize("existing", [False, True])
def test_output_io_failure_is_transactional(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    translate: Callable[..., Any],
    failure: str,
    existing: bool,
) -> None:
    source, output = tmp_path / "input.fa", tmp_path / "output.fa"
    source.write_text(">a\nATG\n>b\nGCC\n")
    if existing:
        output.write_text("original")
    module = importlib.import_module(
        "bioutils_collection.translation_functions._atomic_output"
    )
    original_fdopen = module.os.fdopen

    class FailingOutput:
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            self.file = original_fdopen(*args, **kwargs)
            self.writes = 0

        def __enter__(self) -> Any:
            return self

        def __exit__(self, *args: Any) -> None:
            self.file.__exit__(*args)
            if failure == "close" and args[0] is None:
                raise OSError("simulated close failure")

        def __getattr__(self, name: str) -> Any:
            return getattr(self.file, name)

        def write(self, value: str) -> int:
            self.writes += 1
            if failure == "write" and self.writes == 2:
                raise OSError("simulated write failure")
            return self.file.write(value)

    monkeypatch.setattr(module.os, "fdopen", FailingOutput)
    with pytest.raises(OSError, match=f"simulated {failure} failure"):
        translate(source, output)
    assert output.read_text() == "original" if existing else not output.exists()
    assert {p.name for p in tmp_path.iterdir()} == (
        {"input.fa", "output.fa"} if existing else {"input.fa"}
    )


@pytest.mark.parametrize("translate", [translate_fasta, translate_fasta_streaming])
def test_cancellation_preserves_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, translate: Callable[..., Any]
) -> None:
    source, output = tmp_path / "input.fa", tmp_path / "output.fa"
    source.write_text(">a\nATG\n")
    output.write_text("original")
    module = importlib.import_module(translate.__module__)

    def cancel(*args: Any, **kwargs: Any) -> str:
        raise KeyboardInterrupt("cancelled")

    monkeypatch.setattr(module, "translate_dna_to_protein", cancel)
    with pytest.raises(KeyboardInterrupt, match="cancelled"):
        translate(source, output)
    assert output.read_text() == "original"
    assert {p.name for p in tmp_path.iterdir()} == {"input.fa", "output.fa"}


@pytest.mark.parametrize("translate", [translate_fasta, translate_fasta_streaming])
def test_output_permissions_preserved(
    tmp_path: Path, translate: Callable[..., Any]
) -> None:
    source, output = tmp_path / "input.fa", tmp_path / "output.fa"
    source.write_text(">a\nATG\n")
    output.write_text("old")
    output.chmod(0o640)
    expected_mode = stat.S_IMODE(output.stat().st_mode)
    translate(source, output)
    assert output.read_text() == ">a\nM\n"
    assert stat.S_IMODE(output.stat().st_mode) == expected_mode
