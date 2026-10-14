"""Accept paths pasted from a shell or entered in portable CSV tables."""
from __future__ import annotations

import os
from pathlib import Path, PureWindowsPath


def clean_path_text(value: str | os.PathLike) -> str:
    text = os.fspath(value).strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {'"', "'"}:
        text = text[1:-1]
    text = os.path.expanduser(os.path.expandvars(text))
    # Relative paths in portable tables may use either separator. Absolute
    # Windows paths keep their drive/UNC semantics for native Windows.
    if os.name != "nt" and not PureWindowsPath(text).drive and not text.startswith("\\"):
        text = text.replace("\\", "/")
    return text


def local_path(value: str | os.PathLike) -> Path:
    text = clean_path_text(value)
    if os.name != "nt" and (PureWindowsPath(text).drive or text.startswith("\\")):
        raise ValueError(f"Windows drive or network path needs Windows: {text}")
    if not text:
        raise ValueError("Enter a file or directory path")
    return Path(text)
