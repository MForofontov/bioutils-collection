"""Large-input translation equivalence checks. Timing gates live in explicit benchmarks."""

import pytest
from bioutils_collection.translation_functions import (
    translate_batch,
    translate_dna_fast,
    translate_dna_to_protein,
    translate_large_sequence,
)

pytestmark = [pytest.mark.unit, pytest.mark.slow, pytest.mark.translation]


def test_performance_small_sequence() -> None:
    """
    Test case 1: Performance comparison for small sequences (<1KB).
    """
    seq = "ATGGCCAAA" * 100
    result_default = translate_dna_to_protein(seq)
    _ = translate_dna_fast("ATGGCC" * 10)
    result_fast = translate_dna_fast(seq)
    assert result_default == result_fast


def test_performance_medium_sequence() -> None:
    """
    Test case 2: Performance comparison for medium sequences (~100KB).
    """
    seq = "ATGGCCAAA" * 11111
    result_default = translate_dna_to_protein(seq)
    _ = translate_dna_fast("ATGGCC" * 100)
    result_fast = translate_dna_fast(seq)
    assert result_default == result_fast


def test_performance_large_sequence() -> None:
    """
    Test case 3: Performance comparison for large sequences (~1MB).
    """
    seq = "ATGGCCAAA" * 111111
    result_default = translate_dna_to_protein(seq)
    _ = translate_dna_fast("ATGGCC" * 1000)
    result_fast = translate_dna_fast(seq)
    assert result_default == result_fast


def test_performance_batch_small_sequences() -> None:
    """
    Test case 4: Batch processing performance with many small sequences.
    """
    sequences = ["ATGGCCAAA" * 10] * 100
    result_default = translate_batch(sequences, use_fast=False, n_processes=1)
    _ = translate_dna_fast("ATGGCC" * 10)
    result_fast = translate_batch(sequences, use_fast=True, n_processes=1)
    assert result_default == result_fast


def test_performance_batch_large_sequences() -> None:
    """
    Test case 5: Batch processing with fewer large sequences.
    """
    sequences = ["ATGGCCAAA" * 10000] * 20
    result_default = translate_batch(sequences, use_fast=False, n_processes=1)
    _ = translate_dna_fast("ATGGCC" * 1000)
    result_fast = translate_batch(sequences, use_fast=True, n_processes=1)
    assert result_default == result_fast


def test_performance_batch_multiprocessing() -> None:
    """
    Test case 6: Multiprocessing performance comparison.
    """
    sequences = ["ATGGCCAAA" * 10000] * 40
    result_single = translate_batch(sequences, n_processes=1, use_fast=False)
    result_multi = translate_batch(sequences, n_processes=4, use_fast=False)
    assert result_single == result_multi


def test_performance_large_sequence_chunked() -> None:
    """
    Test case 7: Large sequence with chunked parallel processing.
    """
    seq = "ATGGCCAAA" * 1111111
    _ = translate_dna_fast("ATGGCC" * 1000)
    result_single = translate_dna_fast(seq)
    result_chunked = translate_large_sequence(
        seq, chunk_size=2000000, n_processes=4, use_fast=True
    )
    assert result_single == result_chunked


def test_performance_consistency_across_methods() -> None:
    """
    Test case 8: All methods produce identical results.
    """
    seq = "ATGGCCAAATTTGGG" * 1000
    result_default = translate_dna_to_protein(seq)
    result_batch = translate_batch([seq], n_processes=1, use_fast=False)[0]
    result_large = translate_large_sequence(seq, chunk_size=5000, use_fast=False)
    result_fast = translate_dna_fast(seq)
    result_batch_fast = translate_batch([seq], n_processes=1, use_fast=True)[0]
    result_large_fast = translate_large_sequence(seq, chunk_size=5000, use_fast=True)
    assert result_default == result_fast
    assert result_default == result_batch_fast
    assert result_default == result_large_fast
    assert result_default == result_batch
    assert result_default == result_large


def test_performance_various_sizes() -> None:
    """
    Test case 9: Performance scaling across various sequence sizes.
    """
    sizes = [100, 1000, 10000, 100000]
    for size in sizes:
        seq = "ATGGCCAAA" * (size // 9)
        result_default = translate_dna_to_protein(seq)
        result_fast = translate_dna_fast(seq)
        assert result_default == result_fast
        assert len(result_default) == len(seq) // 3
