"""Waters SYNAPT-specific I/O: STS binary, _extern.inf, CDCReader.

These functions read proprietary Waters binary formats that no other
vendor uses. Generic text file readers (read_im_txt, read_ms_txt) are
in deconvovo/io.py since the text format is universal.
"""
from __future__ import annotations

import math
import struct
from pathlib import Path


def find_pusher_period(raw_path: Path) -> float | None:
    """Read pusher period (μs) from Waters _FUNC001.STS binary.

    Stat code 76 = Pusher Frequency in Hz. Period = floor(1e6 / freq).
    """
    sts = raw_path / "_FUNC001.STS"
    if not sts.exists():
        return None
    try:
        data = sts.read_bytes()
        if len(data) < 0x22:
            return None
        n_desc = struct.unpack_from('<H', data, 4)[0]
        data_start = 0x20 + n_desc * 32
        ds = data[data_start:]
        if len(ds) < n_desc + 16:
            return None
        pusher_doff = None
        for i in range(n_desc):
            off = 0x20 + i * 32
            code = struct.unpack_from('<H', data, off)[0]
            if code == 76:
                pusher_doff = struct.unpack_from('<H', data, off + 4)[0]
                break
        if pusher_doff is None:
            return None
        header_size = len(ds) % n_desc
        val = struct.unpack_from('<i', ds, header_size + pusher_doff)[0]
        if 1000 < val < 100000:
            return float(math.floor(1e6 / val))
        return None
    except Exception:
        return None


def read_extern_inf(raw_path: Path) -> dict[str, str]:
    """Parse Waters _extern.inf instrument parameters."""
    path = raw_path / "_extern.inf"
    params: dict[str, str] = {}
    if not path.exists():
        return params
    with open(path, errors="replace") as f:
        for line in f:
            parts = [p.strip() for p in line.strip().split("\t") if p.strip()]
            if len(parts) >= 2:
                params[parts[0]] = parts[-1]
    return params
