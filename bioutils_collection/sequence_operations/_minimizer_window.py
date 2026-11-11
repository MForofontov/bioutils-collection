"""Sliding-window minimizer selection utilities."""

from collections import deque
from collections.abc import Callable
from typing import TypeVar, cast

T = TypeVar("T")


def sliding_window_minima(
    items: list[T],
    window_size: int,
    key: Callable[[T], str] | None = None,
) -> list[T]:
    """
    Select the leftmost lexicographic minimum in each sliding window.

    Parameters
    ----------
    items : list[T]
        Items to scan (e.g. k-mers or (kmer, position) tuples).
    window_size : int
        Number of consecutive items per window.
    key : Callable[[T], str] | None
        Function extracting the comparable string from each item.
        Defaults to identity for str items.

    Returns
    -------
    list[T]
        One selected item per window (or a single minimum if len(items) < window_size).
    """
    if not items:
        return []

    if key is not None:
        selected_key: Callable[[T], str] = key
    else:

        def _default_key(item: T) -> str:
            return cast(str, item)

        selected_key = cast(Callable[[T], str], _default_key)

    n = len(items)
    if n < window_size:
        return [min(items, key=selected_key)]

    result: list[T] = []
    dq: deque[int] = deque()

    for i in range(n):
        while dq and dq[0] <= i - window_size:
            dq.popleft()
        while dq and selected_key(items[dq[-1]]) > selected_key(items[i]):
            dq.pop()
        dq.append(i)
        if i >= window_size - 1:
            result.append(items[dq[0]])

    return result


__all__ = ["sliding_window_minima"]
