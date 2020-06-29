import gzip
import sys
from typing import TextIO
from contextlib import contextmanager


def open_path(path: str, args="rt") -> TextIO:
    """Open a stream and return the handle"""
    if path == "stdin":
        handle = sys.stdin
    elif path == "stdout":
        handle = sys.stdout
    elif path == "stderr":
        handle = sys.stderr
    elif path.endswith(".gz"):
        handle = gzip.open(path, args)
    else:
        handle = open(path, args)
    return handle


@contextmanager
def open_stream(path: str, args="rt") -> TextIO:
    """Open a stream or return a context manager."""
    handle = open_path(path, args)
    try:
        yield handle
    finally:
        if path not in ("stdin", "stdout", "stderr"):
            handle.close()
