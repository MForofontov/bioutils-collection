"""Performance comparison tests for translation functions.

These tests compare different translation implementations to validate
performance characteristics. Mark as slow tests.
"""

try:
    from bioutils_collection.translation_functions import (
        translate_dna_to_protein,
        translate_dna_fast,
        translate_batch,
        translate_large_sequence,
        NUMBA_AVAILABLE,
    )
    import time
    TRANSLATE_AVAILABLE = True
except ImportError:
    TRANSLATE_AVAILABLE = False
    translate_dna_to_protein = None  # type: ignore
    translate_dna_fast = None  # type: ignore
    translate_batch = None  # type: ignore
    translate_large_sequence = None  # type: ignore
    NUMBA_AVAILABLE = False

import pytest

pytestmark = pytest.mark.skipif(
    not TRANSLATE_AVAILABLE, reason="translation functions not available"
)
pytestmark = [pytestmark, pytest.mark.unit, pytest.mark.slow, pytest.mark.translation]


def test_performance_small_sequence() -> None:
    """
    Test case 1: Performance comparison for small sequences (<1KB).
    """
    # Arrange
    seq = "ATGGCCAAA" * 100  # ~900 bases
    iterations = 1000

    # Act - Default
    start = time.perf_counter()
    for _ in range(iterations):
        result_default = translate_dna_to_protein(seq)
    time_default = time.perf_counter() - start

    # Act - Fast
    if NUMBA_AVAILABLE:
        # Warm up JIT
        _ = translate_dna_fast("ATGGCC" * 10)
        
        start = time.perf_counter()
        for _ in range(iterations):
            result_fast = translate_dna_fast(seq)
        time_fast = time.perf_counter() - start
        
        # Assert
        assert result_default == result_fast
        # For small sequences, difference should be reasonable
        assert time_default > 0
        assert time_fast > 0


def test_performance_medium_sequence() -> None:
    """
    Test case 2: Performance comparison for medium sequences (~100KB).
    """
    # Arrange
    seq = "ATGGCCAAA" * 11_111  # ~100KB
    iterations = 100

    # Act - Default
    start = time.perf_counter()
    for _ in range(iterations):
        result_default = translate_dna_to_protein(seq)
    time_default = time.perf_counter() - start

    # Act - Fast
    if NUMBA_AVAILABLE:
        # Warm up JIT
        _ = translate_dna_fast("ATGGCC" * 100)
        
        start = time.perf_counter()
        for _ in range(iterations):
            result_fast = translate_dna_fast(seq)
        time_fast = time.perf_counter() - start
        
        # Assert
        assert result_default == result_fast
        # Fast should be faster for medium sequences
        speedup = time_default / time_fast
        assert speedup > 1.0, f"Expected speedup >1.0x, got {speedup:.2f}x"


def test_performance_large_sequence() -> None:
    """
    Test case 3: Performance comparison for large sequences (~1MB).
    """
    # Arrange
    seq = "ATGGCCAAA" * 111_111  # ~1MB
    iterations = 10

    # Act - Default
    start = time.perf_counter()
    for _ in range(iterations):
        result_default = translate_dna_to_protein(seq)
    time_default = time.perf_counter() - start

    # Act - Fast
    if NUMBA_AVAILABLE:
        # Warm up JIT
        _ = translate_dna_fast("ATGGCC" * 1000)
        
        start = time.perf_counter()
        for _ in range(iterations):
            result_fast = translate_dna_fast(seq)
        time_fast = time.perf_counter() - start
        
        # Assert
        assert result_default == result_fast
        # Fast should show significant speedup for large sequences
        speedup = time_default / time_fast
        assert speedup > 1.5, f"Expected speedup >1.5x for large sequences, got {speedup:.2f}x"


def test_performance_batch_small_sequences() -> None:
    """
    Test case 4: Batch processing performance with many small sequences.
    """
    # Arrange
    sequences = ["ATGGCCAAA" * 10] * 100  # 100 sequences of ~90 bases

    # Act - use_fast=False
    start = time.perf_counter()
    result_default = translate_batch(sequences, use_fast=False, n_processes=1)
    time_default = time.perf_counter() - start

    # Act - use_fast=True
    if NUMBA_AVAILABLE:
        # Warm up
        _ = translate_dna_fast("ATGGCC" * 10)
        
        start = time.perf_counter()
        result_fast = translate_batch(sequences, use_fast=True, n_processes=1)
        time_fast = time.perf_counter() - start
        
        # Assert
        assert result_default == result_fast
        assert time_default > 0
        assert time_fast > 0


def test_performance_batch_large_sequences() -> None:
    """
    Test case 5: Batch processing with fewer large sequences.
    """
    # Arrange
    sequences = ["ATGGCCAAA" * 10_000] * 20  # 20 sequences of ~90KB each

    # Act - use_fast=False
    start = time.perf_counter()
    result_default = translate_batch(sequences, use_fast=False, n_processes=1)
    time_default = time.perf_counter() - start

    # Act - use_fast=True
    if NUMBA_AVAILABLE:
        # Warm up
        _ = translate_dna_fast("ATGGCC" * 1000)
        
        start = time.perf_counter()
        result_fast = translate_batch(sequences, use_fast=True, n_processes=1)
        time_fast = time.perf_counter() - start
        
        # Assert
        assert result_default == result_fast
        # Fast should be faster for large sequences
        if time_default > 0 and time_fast > 0:
            speedup = time_default / time_fast
            assert speedup > 1.0, f"Expected speedup for large sequences, got {speedup:.2f}x"


def test_performance_batch_multiprocessing() -> None:
    """
    Test case 6: Multiprocessing performance comparison.
    """
    # Arrange
    sequences = ["ATGGCCAAA" * 10_000] * 40  # 40 sequences of ~90KB

    # Act - Single process
    start = time.perf_counter()
    result_single = translate_batch(sequences, n_processes=1, use_fast=False)
    time_single = time.perf_counter() - start

    # Act - Multiple processes
    start = time.perf_counter()
    result_multi = translate_batch(sequences, n_processes=4, use_fast=False)
    time_multi = time.perf_counter() - start

    # Assert
    assert result_single == result_multi
    assert time_single > 0
    assert time_multi > 0
    # Multiprocessing should show some benefit with many large sequences
    # (Though overhead may dominate with moderate sizes)


def test_performance_large_sequence_chunked() -> None:
    """
    Test case 7: Large sequence with chunked parallel processing.
    """
    # Arrange
    seq = "ATGGCCAAA" * 1_111_111  # ~10MB

    # Act - Single threaded fast
    if NUMBA_AVAILABLE:
        # Warm up
        _ = translate_dna_fast("ATGGCC" * 1000)
        
        start = time.perf_counter()
        result_single = translate_dna_fast(seq)
        time_single = time.perf_counter() - start

        # Act - Chunked parallel
        start = time.perf_counter()
        result_chunked = translate_large_sequence(
            seq, chunk_size=2_000_000, n_processes=4, use_fast=True
        )
        time_chunked = time.perf_counter() - start

        # Assert
        assert result_single == result_chunked
        assert time_single > 0
        assert time_chunked > 0
        # Chunked may or may not be faster depending on overhead


def test_performance_consistency_across_methods() -> None:
    """
    Test case 8: All methods produce identical results.
    """
    # Arrange
    seq = "ATGGCCAAATTTGGG" * 1000

    # Act
    result_default = translate_dna_to_protein(seq)
    result_batch = translate_batch([seq], n_processes=1, use_fast=False)[0]
    result_large = translate_large_sequence(seq, chunk_size=5000, use_fast=False)
    
    if NUMBA_AVAILABLE:
        result_fast = translate_dna_fast(seq)
        result_batch_fast = translate_batch([seq], n_processes=1, use_fast=True)[0]
        result_large_fast = translate_large_sequence(seq, chunk_size=5000, use_fast=True)
        
        # Assert - All methods produce same result
        assert result_default == result_fast
        assert result_default == result_batch_fast
        assert result_default == result_large_fast
    
    # Assert - Non-fast methods all match
    assert result_default == result_batch
    assert result_default == result_large


def test_performance_various_sizes() -> None:
    """
    Test case 9: Performance scaling across various sequence sizes.
    """
    # Arrange
    sizes = [100, 1_000, 10_000, 100_000]  # 100 bases to 100KB
    
    # Act & Assert
    for size in sizes:
        seq = "ATGGCCAAA" * (size // 9)
        
        # Should complete without errors
        result_default = translate_dna_to_protein(seq)
        
        if NUMBA_AVAILABLE:
            result_fast = translate_dna_fast(seq)
            assert result_default == result_fast
        
        # Verify result is reasonable
        assert len(result_default) == len(seq) // 3
