"""The curated SLA haplotype table behind Table S1.

Methods: "Selection of SLA alleles for structural and position-level analyses was
guided by a curated haplotype table listing common SLA haplotypes across pig
breeds/herds and SLA-defined experimental miniature-pig lines."

Input
-----
The curated workbooks in ``data/haplotypes/``, laid out as::

    SLA Haplotype | Class I Haplotype | SLA-1 | SLA-3 | SLA-2 || Class II Haplotype | DRB1 | DQB1 | Local Names

A blank-ish row carrying only free text introduces a section -- ``Yucatan
Miniature Pigs``, ``NIH/MGH Mini Pigs`` -- recorded as the herd for the rows
beneath it. Allele cells hold the fields only, so the locus prefix is added here
(``01:01`` under ``SLA-1`` becomes ``SLA-1*01:01``). A cell may hold two alleles
(``09:01, 15:01``); both are kept. ``Null`` and ``NT`` mean no typed allele, and
an ``XX`` second field (``07:XX``) is ambiguous -- both are recorded but flagged,
never silently expanded.

Ported from the development package's ``slahla.haplotypes``; the parsing is
unchanged. The workbook filenames carry their haplotype count so the larger
table wins where the two disagree, which is what the published table used.
"""
from __future__ import annotations

import re

import pandas as pd

from .paths import DATA, HAPLOTYPES

CLASS_I_COLUMNS = {"SLA-1": "SLA-1", "SLA-3": "SLA-3", "SLA-2": "SLA-2"}
CLASS_II_COLUMNS = {"DRB1": "SLA-DRB1", "DQB1": "SLA-DQB1",
                    "DQA": "SLA-DQA", "DRA": "SLA-DRA"}
SECTION_WORDS = ("miniature", "mini pig", "minipig", "commercial", "yucatan", "nih", "mgh")
NO_ALLELE = {"null", "nt", "n/a", "na", "-", ""}
AMBIGUOUS = re.compile(r"XX", re.IGNORECASE)


def workbooks() -> list:
    """Candidate workbooks, most-haplotypes-first so the larger table wins ties."""
    found = sorted(HAPLOTYPES.glob("*.xlsx"))
    return sorted(found, key=lambda p: ("_14_" not in p.name, p.name))


def _clean(value) -> str:
    return "" if pd.isna(value) else str(value).strip()


def _compose_name(class1: str, class2: str) -> str:
    """Class I ``1a.0`` plus class II ``0.1`` is named ``Hp-1a.1``.

    The class I trailing ``.0`` and the class II leading ``0.`` are both dropped;
    composing them naively yields ``Hp-1a.0.0.1``, which then fails to match the
    same haplotype named explicitly in another workbook.
    """
    left = re.sub(r"\.0$", "", class1.strip())
    right = re.sub(r"^0\.", "", class2.strip())
    if left and right:
        return f"Hp-{left}.{right}"
    return f"Hp-{left or right}" if (left or right) else ""


def _split_alleles(cell: str) -> list[str]:
    parts = [p.strip() for p in re.split(r"[,;/]", cell) if p.strip()]
    return [p for p in parts if p.lower() not in NO_ALLELE]


def parse_workbook(path) -> pd.DataFrame:
    """Parse one curated haplotype workbook into tidy rows."""
    raw = pd.read_excel(path, header=None)
    header_row = None
    for i in range(min(6, len(raw))):
        joined = " ".join(_clean(v) for v in raw.iloc[i]).lower()
        if "haplotype" in joined and "sla-1" in joined:
            header_row = i
            break
    if header_row is None:
        raise ValueError(f"{path}: could not find the header row")
    headers = [_clean(v) for v in raw.iloc[header_row]]
    col_of = {h: j for j, h in enumerate(headers) if h}
    body = raw.iloc[header_row + 1:]

    def column(name: str):
        for h, j in col_of.items():
            if h.lower() == name.lower():
                return j
        return None

    hap_col, class1_col = column("SLA Haplotype"), column("Class I Haplotype")
    class2_col, local_col = column("Class II Haplotype"), column("Local Names")
    rows, herd = [], "commercial"
    for _, line in body.iterrows():
        cells = [_clean(v) for v in line]
        if not any(cells):
            continue
        text = " ".join(cells).lower()
        typed = sum(1 for c in cells if re.search(r"\d+:\d+", c))
        if typed == 0 and any(w in text for w in SECTION_WORDS):
            herd = " ".join(c for c in cells if c)
            continue
        if typed == 0:
            continue
        haplotype = (cells[hap_col] if hap_col is not None and hap_col < len(cells)
                     and cells[hap_col] else "")
        class1 = cells[class1_col] if class1_col is not None and class1_col < len(cells) else ""
        class2 = cells[class2_col] if class2_col is not None and class2_col < len(cells) else ""
        if not haplotype:
            haplotype = _compose_name(class1, class2)
        local = cells[local_col] if local_col is not None and local_col < len(cells) else ""
        for group, mapping in (("I", CLASS_I_COLUMNS), ("II", CLASS_II_COLUMNS)):
            for header, locus in mapping.items():
                j = column(header)
                if j is None or j >= len(cells):
                    continue
                cell = cells[j]
                if not cell or cell.lower() in NO_ALLELE:
                    continue
                for fields in _split_alleles(cell):
                    rows.append(dict(
                        haplotype=haplotype, locus=locus, allele=f"{locus}*{fields}",
                        mhc_class=group, herd=herd, class1_haplotype=class1,
                        class2_haplotype=class2, local_name=local,
                        ambiguous=bool(AMBIGUOUS.search(fields)),
                        source_ref=path.name))
    return pd.DataFrame(rows)


def load() -> pd.DataFrame:
    """All curated workbooks, concatenated, duplicates resolved to the larger table."""
    books = workbooks()
    if not books:
        raise FileNotFoundError(f"no curated haplotype workbook in {HAPLOTYPES}")
    frames = [parse_workbook(p) for p in books]
    out = pd.concat(frames, ignore_index=True)
    return out.drop_duplicates(subset=["haplotype", "locus", "allele"], keep="first")


#: Alleles whose reported haplotype cannot be regenerated, with the archived
#: selection carried in its place. See the file's own ``reason`` column.
CARRIED_SELECTION = DATA / "carried" / "tableS1_haplotype_selection.csv"


def carried_selection() -> pd.DataFrame:
    """The archived haplotype selections this repository cannot regenerate.

    Three Table S1 alleles sit in more than one curated haplotype, and the
    approved table reports one. No candidate rule reproduces the choice --
    first row, last row, first commercial row and last commercial row were all
    tested against the 34 multi-haplotype alleles and each fails on at least one.
    The rule is therefore unrecoverable, and rather than invent one the archived
    value is carried and labelled as carried.
    """
    if not CARRIED_SELECTION.exists():
        return pd.DataFrame(columns=["sla_allele", "archived_haplotype",
                                     "archived_haplotypes", "source", "reason"])
    return pd.read_csv(CARRIED_SELECTION)


def alleles_by_haplotype() -> pd.DataFrame:
    """One row per SLA allele, with every haplotype it appears in.

    ``haplotype`` and ``haplotypes`` are regenerated from the curated workbooks,
    except for the alleles listed in ``carried_selection()``, where the archived
    value is used. ``haplotype_selection_source`` says which applies to each row,
    and ``haplotypes_all_curated`` always carries the full regenerated set so no
    information is lost to the carry.
    """
    d = load()
    g = (d.groupby("allele")
           .agg(sla_locus=("locus", "first"),
                mhc_class=("mhc_class", "first"),
                haplotype=("haplotype", "first"),
                haplotypes=("haplotype", lambda s: ", ".join(sorted(set(s)))),
                herd=("herd", "first"),
                # local_name and herd describe the SELECTED haplotype row, not any
                # row the allele appears in. Taking the first non-empty value across
                # all of them imports a miniature-line label onto a commercial
                # haplotype: five Table S1 rows differ from the approved table that
                # way, all of them alleles shared between a commercial haplotype
                # (no local name) and a Yucatan or NIH/MGH one (X, Y or Z).
                local_name=("local_name", "first"),
                ambiguous=("ambiguous", "any"))
           .reset_index().rename(columns={"allele": "sla_allele"}))
    g["haplotypes_all_curated"] = g["haplotypes"]
    g["haplotype_selection_source"] = "regenerated from the curated workbooks"
    carried = carried_selection()
    if len(carried):
        idx = g.sla_allele.isin(carried.sla_allele)
        m = carried.set_index("sla_allele")
        g.loc[idx, "haplotype"] = g.loc[idx, "sla_allele"].map(m.archived_haplotype)
        g.loc[idx, "haplotypes"] = g.loc[idx, "sla_allele"].map(m.archived_haplotypes)
        g.loc[idx, "haplotype_selection_source"] = (
            "carried from the approved Table S1 (selection rule not recorded)")
    return g.sort_values(["sla_locus", "sla_allele"]).reset_index(drop=True)
