"""Calculate GC content using sliding windows."""


def gc_content_windows(
    seq: str, window_size: int = 100, step_size: int | None = None
) -> list[tuple[int, int, float]]:
    """
    Calculate GC content across a sequence using sliding windows.

    Parameters
    ----------
    seq : str
        DNA sequence to analyze.
    window_size : int, optional
        Size of sliding window in base pairs (by default 100).
    step_size : int | None, optional
        Step size for sliding window. If None, uses window_size (non-overlapping).
        Use smaller values for overlapping windows (by default None).

    Returns
    -------
    list[tuple[int, int, float]]
        List of tuples containing (start_pos, end_pos, gc_content_percentage).

    Raises
    ------
    TypeError
        If seq is not a string.
        If window_size or step_size is not an integer or None.
    ValueError
        If window_size is less than 1.
        If step_size is less than 1.
        If seq contains invalid DNA bases.
        If seq is shorter than window_size.

    Examples
    --------
    >>> gc_content_windows("ATGCATGCATGC", window_size=4, step_size=4)
    [(0, 4, 50.0), (4, 8, 50.0), (8, 12, 50.0)]
    >>> gc_content_windows("AAAAGGGCCCCTTTTT", window_size=4, step_size=2)
    [(0, 4, 0.0), (2, 6, 50.0), (4, 8, 100.0), ...]
    >>> gc_content_windows("ATGCATGC", window_size=4)
    [(0, 4, 50.0), (4, 8, 50.0)]

    Notes
    -----
    GC content is calculated as: (G + C) / window_size * 100
    Non-overlapping windows when step_size == window_size.
    Overlapping windows when step_size < window_size.
    Incomplete trailing windows are skipped.

    Complexity
    ----------
    Time: O(n), Space: O(n) including normalized input and results.
    """
    # Input validation
    if not isinstance(seq, str):
        raise TypeError(f"seq must be a string, got {type(seq).__name__}")
    if not isinstance(window_size, int):
        raise TypeError(
            f"window_size must be an integer, got {type(window_size).__name__}"
        )
    if step_size is not None and not isinstance(step_size, int):
        raise TypeError(
            f"step_size must be an integer or None, got {type(step_size).__name__}"
        )

    if window_size < 1:
        raise ValueError(f"window_size must be at least 1, got {window_size}")
    if step_size is not None and step_size < 1:
        raise ValueError(f"step_size must be at least 1, got {step_size}")

    if len(seq) < window_size:
        raise ValueError(
            f"Sequence length ({len(seq)}) must be at least window_size ({window_size})"
        )

    # Validate DNA sequence
    seq_upper = seq.upper()
    valid_bases = set("ATGCN")
    invalid_bases = set(seq_upper) - valid_bases
    if invalid_bases:
        raise ValueError(f"Invalid DNA bases found: {', '.join(sorted(invalid_bases))}")

    # Set step size to window size if not specified (non-overlapping)
    if step_size is None:
        step_size = window_size

    results = []

    gc_count = seq_upper.count("G", 0, window_size) + seq_upper.count(
        "C", 0, window_size
    )
    results.append((0, window_size, round(gc_count / window_size * 100.0, 2)))

    for start in range(step_size, len(seq) - window_size + 1, step_size):
        end = start + window_size
        if step_size < window_size:
            # Update only the outgoing and incoming bases of overlapping windows.
            previous_start = start - step_size
            previous_end = end - step_size
            gc_count -= seq_upper.count("G", previous_start, start)
            gc_count -= seq_upper.count("C", previous_start, start)
            gc_count += seq_upper.count("G", previous_end, end)
            gc_count += seq_upper.count("C", previous_end, end)
        else:
            gc_count = seq_upper.count("G", start, end) + seq_upper.count(
                "C", start, end
            )
        results.append((start, end, round(gc_count / window_size * 100.0, 2)))

    return results


__all__ = ["gc_content_windows"]
