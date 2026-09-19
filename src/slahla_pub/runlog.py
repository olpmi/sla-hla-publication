"""Record outputs that could not be produced, and why.

A reproduction that quietly emits fewer files than it promises is worse than one
that fails: the reader has no way to tell a missing output from an output that
was never meant to exist. Every stage that cannot run because an external input
is absent records the fact here, and ``make reproduce`` prints the result in a
block that cannot be mistaken for success.
"""
from __future__ import annotations

import json

from .paths import REPORTS, ensure_outputs

PATH = REPORTS / "missing_inputs.json"


def clear() -> None:
    ensure_outputs()
    if PATH.exists():
        PATH.unlink()


def record(what: str, blocks: list[str], how: str) -> None:
    ensure_outputs()
    rec = read()
    if not any(r["input"] == what for r in rec):
        rec.append({"input": what, "blocked_outputs": blocks, "how_to_supply_it": how})
    PATH.write_text(json.dumps(rec, indent=2))


def read() -> list[dict]:
    return json.loads(PATH.read_text()) if PATH.exists() else []


def blocked() -> list[str]:
    out: list[str] = []
    for r in read():
        out.extend(r["blocked_outputs"])
    return sorted(set(out))


def banner() -> str:
    rec = read()
    if not rec:
        return ""
    lines = ["", "=" * 72, "MISSING EXTERNAL INPUTS - some published outputs were NOT produced", "=" * 72]
    for r in rec:
        lines.append(f"  input   : {r['input']}")
        lines.append(f"  blocks  : {', '.join(r['blocked_outputs'])}")
        for ln in r["how_to_supply_it"].splitlines():
            lines.append(f"  {ln}" if ln.startswith(" ") else f"  supply  : {ln}")
    lines.append("=" * 72)
    return "\n".join(lines)
