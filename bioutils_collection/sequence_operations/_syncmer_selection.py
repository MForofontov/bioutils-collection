"""Shared syncmer validation and streaming, leftmost minimum selection."""

from collections import deque
from collections.abc import Iterator


def validate_syncmer_inputs(seq: str, k: int, s: int, method: str) -> None:
    if not isinstance(seq, str):
        raise TypeError(f"seq must be str, got {type(seq).__name__}")
    if not isinstance(k, int):
        raise TypeError(f"k must be int, got {type(k).__name__}")
    if not isinstance(s, int):
        raise TypeError(f"s must be int, got {type(s).__name__}")
    if not isinstance(method, str):
        raise TypeError(f"method must be str, got {type(method).__name__}")
    if not seq:
        raise ValueError("seq cannot be empty")
    if k <= 0:
        raise ValueError("k must be positive")
    if s <= 0:
        raise ValueError("s must be positive")
    if s >= k:
        raise ValueError("s must be less than k")
    if k > len(seq):
        raise ValueError("k cannot be longer than sequence")
    if method not in ("open", "closed"):
        raise ValueError("method must be 'open' or 'closed'")


def iter_syncmer_positions(seq: str, k: int, s: int, method: str) -> Iterator[int]:
    """Scan uppercase input with a monotonic deque; retain equal keys for ties."""
    width = k - s + 1
    candidates: deque[tuple[int, str]] = deque()
    for i in range(len(seq) - s + 1):
        start = i - width + 1
        while candidates and candidates[0][0] < start:
            candidates.popleft()
        smer = seq[i : i + s]
        while candidates and candidates[-1][1] > smer:
            candidates.pop()
        candidates.append((i, smer))
        if start >= 0:
            minimum_position = candidates[0][0]
            if minimum_position == start or (
                method == "closed" and minimum_position == i
            ):
                yield start
