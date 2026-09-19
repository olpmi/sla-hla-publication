"""Path resolution.

Every path is configurable. Nothing here resolves relative to the current
working directory, to a user's home, or to a sibling checkout: the repository
root is found from this file's own location, and any of it can be overridden by
``configs/reproduce.yaml`` or by the ``SLAHLA_PUB_*`` environment variables.

The invariant the reproduction depends on: ``SOURCE_DATA`` and everything else
under ``data/`` is **read-only input**. Every generated file goes under
``OUTPUTS``. ``make verify`` asserts the inputs are byte-unchanged after a run.
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _cfg() -> dict:
    path = Path(os.environ.get("SLAHLA_PUB_CONFIG", ROOT / "configs" / "reproduce.yaml"))
    if not path.exists():
        return {}
    import yaml

    return yaml.safe_load(path.read_text()) or {}


_C = _cfg()


def _p(key: str, default: str) -> Path:
    env = os.environ.get("SLAHLA_PUB_" + key.upper())
    if env:
        return Path(env).expanduser().resolve()
    value = _C.get("paths", {}).get(key, default)
    p = Path(value).expanduser()
    return p if p.is_absolute() else (ROOT / p).resolve()


DATA = _p("data", "data")
SOURCE_DATA = _p("source_data", "data/source_data")
PANELS = _p("structural_panels", "data/structural_panels")
PREDICTIONS = _p("archived_predictions", "data/archived_predictions")
EPLET_AGREEMENT = _p("eplet_agreement", "data/eplet_agreement")
MODELS = _p("structure_models", "data/structures/models")
REFERENCE = _p("structure_reference", "data/structures/reference")
HAPLOTYPES = _p("haplotypes", "data/haplotypes")
ARTWORK = _p("artwork", "publication_artwork")

OUTPUTS = _p("outputs", "outputs")
FIGURES = OUTPUTS / "figures"
TABLES = OUTPUTS / "tables"
AUDIT = OUTPUTS / "audit"
REPORTS = OUTPUTS / "reports"


def ensure_outputs() -> None:
    for d in (OUTPUTS, FIGURES, TABLES, AUDIT, REPORTS):
        d.mkdir(parents=True, exist_ok=True)


def settings() -> dict:
    return _C.get("settings", {})
