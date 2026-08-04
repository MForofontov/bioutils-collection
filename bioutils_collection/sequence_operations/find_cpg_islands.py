"""Find CpG islands in DNA sequences."""


def find_cpg_islands(
    seq: str, window: int = 200, min_gc: float = 0.5, min_obs_exp: float = 0.6
) -> list[tuple[int, int]]:
    """
    Identify CpG islands in a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence (A, T, C, G).
    window : int, optional
        Sliding window size (default: 200).
    min_gc : float, optional
        Minimum GC content (default: 0.5).
    min_obs_exp : float, optional
        Minimum observed/expected CpG ratio (default: 0.6).

    Returns
    -------
    List[Tuple[int, int]]
        List of (start, end) positions of CpG islands.

    Raises
    ------
    TypeError
        If seq is not a string or parameters have wrong types.
    ValueError
        If seq contains invalid characters or parameters are invalid.

    Examples
    --------
    >>> find_cpg_islands('GCGCGCGCGC' * 20, window=50)
    [(0, 200)]

    Notes
    -----
    CpG islands are regions with high GC content and CpG dinucleotide frequency.
    Often found near gene promoters.

    References
    ----------
    Gardiner-Garden, M., Frommer, M. (1987).
    CpG islands in vertebrate genomes.
    Journal of Molecular Biology 196(2):261-282.

    Complexity
    ----------
    Time: O(n), Space: O(k) where k is number of islands
    """
    if not isinstance(seq, str):
        raise TypeError(f"seq must be str, got {type(seq).__name__}")
    if not isinstance(window, int):
        raise TypeError(f"window must be int, got {type(window).__name__}")
    if not isinstance(min_gc, (int, float)):
        raise TypeError(f"min_gc must be a number, got {type(min_gc).__name__}")
    if not isinstance(min_obs_exp, (int, float)):
        raise TypeError(
            f"min_obs_exp must be a number, got {type(min_obs_exp).__name__}"
        )

    seq = seq.upper()
    if not all(base in "ATCG" for base in seq):
        raise ValueError("Sequence contains invalid DNA bases")

    if window <= 0 or window > len(seq):
        raise ValueError("window must be positive and <= sequence length")
    if not 0 <= min_gc <= 1:
        raise ValueError("min_gc must be between 0 and 1")
    if min_obs_exp < 0:
        raise ValueError("min_obs_exp must be non-negative")

    gc_count = sum(1 for base in seq[:window] if base in "GC")
    c_count = seq[:window].count("C")
    g_count = seq[:window].count("G")
    cpg_count = sum(
        1 for i in range(window - 1) if seq[i : i + 2] == "CG"
    )

    def _window_qualifies() -> bool:
        gc_content = gc_count / window
        if gc_content < min_gc:
            return False
        if c_count > 0 and g_count > 0:
            expected_cpg = (c_count * g_count) / window
            obs_exp_ratio = cpg_count / expected_cpg if expected_cpg > 0 else 0
            return obs_exp_ratio >= min_obs_exp
        return False

    qualifying_starts: list[int] = []
    if _window_qualifies():
        qualifying_starts.append(0)

    for start in range(1, len(seq) - window + 1):
        left = start - 1
        right = start + window - 1

        left_base = seq[left]
        right_base = seq[right]
        if left_base in "GC":
            gc_count -= 1
        if right_base in "GC":
            gc_count += 1
        if left_base == "C":
            c_count -= 1
        elif left_base == "G":
            g_count -= 1
        if right_base == "C":
            c_count += 1
        elif right_base == "G":
            g_count += 1

        if seq[left : left + 2] == "CG":
            cpg_count -= 1
        if seq[right - 1 : right + 1] == "CG":
            cpg_count += 1

        if _window_qualifies():
            qualifying_starts.append(start)

    if not qualifying_starts:
        return []

    merged: list[tuple[int, int]] = []
    island_start = qualifying_starts[0]
    prev_start = qualifying_starts[0]
    for start in qualifying_starts[1:]:
        if start <= prev_start + 1:
            prev_start = start
        else:
            merged.append((island_start, prev_start + window))
            island_start = start
            prev_start = start
    merged.append((island_start, prev_start + window))
    return merged


__all__ = ["find_cpg_islands"]
