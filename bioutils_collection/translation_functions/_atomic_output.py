"""Commit complete text output only after successful processing and close."""

import os
import stat
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO


@contextmanager
def atomic_text_output(path: Path) -> Iterator[TextIO]:
    # Follow existing symlinks as open(path, 'w') did, without replacing the link.
    target = path.resolve()
    fd, temporary = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as output:
            yield output
            output.flush()
            os.fsync(output.fileno())
        if target.exists():
            os.chmod(temporary, stat.S_IMODE(target.stat().st_mode))
        os.replace(temporary, target)
    finally:
        Path(temporary).unlink(missing_ok=True)
