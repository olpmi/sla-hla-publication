"""Structural measurements: per-model confidence (Table S4a) and superpositions (S4b).

Two workflows, deliberately separate
------------------------------------
**Default (``make reproduce``)** assembles Figures 5/6 from the archived rendered
panels and formats Tables S4a/S4b from the archived measurements. Its starting
point is an archived upstream result and the documentation says so; it needs no
PyMOL and no structure prediction.

**Optional (``make structures``)** recomputes the measurements from the 24
deposited model PDBs and 3 experimental reference structures. It requires PyMOL
(``environments/structures.yml``) and is what a reader runs to check the numbers
rather than the formatting.

Why PyMOL and not a pure-Python aligner
---------------------------------------
Table S4b reports three *different* superpositions of the same pair --
``cealign`` (CE, C-alpha, no outlier rejection), ``super`` (sequence-independent,
all atoms, 5 cycles, 2.0 A cutoff) and ``align`` (Needleman-Wunsch, C-alpha, 5
cycles, 2.0 A cutoff) -- plus groove and non-groove domain RMSDs. Only PyMOL
produces that set: ``super`` and ``align`` have no equivalent in Biopython, and a
substituted CE implementation would give a different placement, different aligned
residues and different rounded values. The published procedure is retained.

The superposition code is a faithful port of ``recompute_structural_rmsd.py``
from the manuscript support package, which is the authoritative generator for
Table S4b. Only the paths were made configurable.

A note the published table carries: the nine RMSDs printed in the manuscript were
produced interactively and typed in, so they have no recoverable computational
provenance. The table therefore reports all three methods with their units and
denominators beside the printed value, and the ``delta_*_minus_reported`` columns
state the difference rather than hiding it.
"""
from __future__ import annotations

import hashlib
import json
import sys

import numpy as np
import pandas as pd

from .paths import (ARTWORK, DATA, MODELS, REFERENCE, REPORTS, SOURCE_DATA, TABLES,
                    ensure_outputs)
from .paths import PANELS as PANEL_DIR   # `PANELS` below is the panel definition list

R = "_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb"

#: The nine published panels plus the alternative allele investigated for 5F.
PANELS = [
    dict(final="Figure 5", panel="A", pn="class1_A", cls="class1", hla="A*02:01", sla="SLA-1*02:01",
         hf=f"HLA_A02_A_02_01{R}", sf=f"SLA_A02_SLA-1_02_01{R}", reported=2.8, match="exact"),
    dict(final="Figure 5", panel="B", pn="class1_B", cls="class1", hla="A*03:01", sla="SLA-2*01:02",
         hf=f"HLA_A03_A_03_01{R}", sf=f"SLA_A03_SLA-2_01_02{R}", reported=3.8, match="exact"),
    dict(final="Figure 5", panel="C", pn="class1_C", cls="class1", hla="B*08:01", sla="SLA-3*01:01",
         hf=f"HLA_B08_B_08_01{R}", sf=f"SLA_B08_SLA-3_01_01{R}", reported=2.6, match="exact"),
    dict(final="Figure 5", panel="D", pn="class1_D", cls="class1", hla="B*15:01", sla="SLA-1*01:01",
         hf=f"HLA_B15_B_15_01{R}", sf=f"SLA_B15_SLA-1_01_01{R}", reported=2.9, match="exact"),
    dict(final="Figure 5", panel="E", pn="class1_E", cls="class1", hla="B*57:01", sla="SLA-2*01:01",
         hf=f"HLA_B57_B_57_01{R}", sf=f"SLA_B57_SLA-2_01_01{R}", reported=2.3, match="exact"),
    dict(final="Figure 5", panel="F", pn="class1_F", cls="class1", hla="C*07:02", sla="SLA-3*06:01",
         hf=f"HLA_C07_C_07_02{R}", sf=f"SLA_C07_SLA-3_06_01{R}", reported=3.4,
         match="SUBSTITUTED - the legend names C*07:01", legend_allele="C*07:01"),
    dict(final="Figure 5", panel="F", pn="class1_F_alt", cls="class1", hla="C*07:01", sla="SLA-3*06:01",
         hf="C07_C0701_unrelaxed_CUSTOM_B.pdb", sf=f"SLA_C07_SLA-3_06_01{R}", reported=3.4,
         match="the allele the legend names, from the authors' CUSTOM_B model set"),
    dict(final="Figure 6", panel="A", pn="class2_A", cls="class2", hla="DRB1*01:02", sla="SLA-DRB1*01:01",
         hf=f"HLA_DRB101_DRB1_01_02{R}", sf=f"SLA_DRB101_SLA-DRB1_01_01{R}", reported=2.7, match="exact"),
    dict(final="Figure 6", panel="B", pn="class2_B", cls="class2", hla="DRB1*11:04", sla="SLA-DRB1*04:02",
         hf=f"HLA_DRB111_DRB1_11_04{R}", sf=f"SLA_DRB111_SLA-DRB1_04_02{R}", reported=3.3, match="exact"),
    dict(final="Figure 6", panel="C", pn="class2_C", cls="class2", hla="DQB1*02:02", sla="SLA-DQB1*02:02",
         hf=f"HLA_DQB102_DQB1_02_02{R}", sf=f"SLA_DQB102_SLA-DQB1_02_02{R}", reported=1.5, match="exact"),
]

DOMAINS = {"class1": dict(groove=(1, 182), rest=(183, 275)),
           "class2": dict(groove=(1, 95), rest=(96, 190))}

#: What each experimental reference structure is, and what the monomer models
#: therefore lack. Descriptive constants, not measurements.
REFERENCE_INFO = {
    "1HHK": ("A", "HLA-A*02:01 with beta-2-microglobulin and the HTLV-1 Tax nonamer",
             "beta-2-microglobulin (chain B) and peptide (chain C)"),
    "1AQD": ("B", "HLA-DR1 (DRA*01:01 / DRB1*01:01) with a bound peptide",
             "the DR alpha chain (chain A) and peptide (chain C)"),
    "1S9V": ("B", "HLA-DQ2 (DQA1*05:01 / DQB1*02:01) with an alpha-I gliadin peptide",
             "the DQ alpha chain (chain A) and peptide (chain C)"),
}


def allele_of(model_name: str) -> str | None:
    """Recover the allele a deposited model file is for, from its filename.

    ``HLA_A02_A_02_01_unrelaxed_...``      -> ``A*02:01``
    ``SLA_DRB101_SLA-DRB1_01_01_unrel...`` -> ``SLA-DRB1*01:01``

    Returns None for files that do not follow the ColabFold naming (the one
    ``*_CUSTOM_B.pdb`` model, which is the alternative allele investigated for
    Figure 5F and is not a Table S4a row).
    """
    if "_unrelaxed_rank_001" not in model_name:
        return None
    stem = model_name.split("_unrelaxed_rank_001")[0]
    parts = stem.split("_", 2)
    if len(parts) < 3:
        return None
    rest = parts[2].split("_")
    if len(rest) < 3:
        return None
    return f"{rest[0]}*{rest[1]}:{rest[2]}"


def reference_map() -> dict:
    return json.loads((DATA / "structures" / "structure_reference_map.json").read_text())


def sha256(p) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _resi_spec(refmap: dict, allele: str, lo: int, hi: int):
    """Model residue numbers whose canonical number falls in [lo, hi]."""
    e = refmap.get(allele)
    if not e:
        return None, 0
    keep = sorted(int(m) for m, c in e["pairs"].items() if lo <= c <= hi)
    if len(keep) < 20:
        return None, len(keep)
    runs, start, prev = [], keep[0], keep[0]
    for v in keep[1:]:
        if v != prev + 1:
            runs.append((start, prev)); start = v
        prev = v
    runs.append((start, prev))
    return "+".join(f"{a}-{b}" for a, b in runs), len(keep)


def plddt(path) -> np.ndarray:
    """Per-residue pLDDT from a model's C-alpha B-factor column.

    ColabFold writes pLDDT into the B-factor column, so this needs no structural
    library. Note that the ``*_CUSTOM_B.pdb`` models have that column overwritten
    and their confidence statistics are therefore not pLDDT.
    """
    vals = []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            vals.append(float(line[60:66]))
    return np.asarray(vals)


def plddt_by_residue(path) -> dict[int, float]:
    """Model residue number -> C-alpha pLDDT."""
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out[int(line[22:26])] = float(line[60:66])
    return out


def confidence() -> pd.DataFrame:
    """Table S4a's per-model confidence block, from the deposited model files.

    Whole-chain pLDDT needs nothing but the B-factor column. The groove/rest
    split and the ``unmapped`` residues need the reference map, which says which
    model residue carries which canonical position; residues with no canonical
    position are the ``unmapped`` set (mostly disordered termini).
    """
    refmap = reference_map()
    rows = []
    for p in sorted(MODELS.glob("*.pdb")):
        by_resi = plddt_by_residue(p)
        if not by_resi:
            continue
        v = np.array(list(by_resi.values()))
        allele = allele_of(p.name)
        rec = {
            "record_level": "individual model", "allele": allele, "model": p.name,
            "n_residues": len(v),
            "plddt_mean": round(float(v.mean()), 1),
            "plddt_median": round(float(np.median(v)), 1),
            "plddt_min": round(float(v.min()), 1),
            "pct_ge_70": round(float((v >= 70).mean() * 100), 1),
            "pct_ge_90": round(float((v >= 90).mean() * 100), 1),
            "predictor": ("ColabFold (AlphaFold2, monomer)" if "rank_001" in p.name
                          else "authors' CUSTOM_B model (B-factor column overwritten)"),
            "sha256": sha256(p),
        }
        entry = refmap.get(allele) if allele else None
        if entry:
            cls = entry["mhc_class"]
            canon = {int(m): c for m, c in entry["pairs"].items()}
            lo_g, hi_g = DOMAINS[cls]["groove"]
            groove = [by_resi[m] for m, c in canon.items() if m in by_resi and lo_g <= c <= hi_g]
            rest = [by_resi[m] for m, c in canon.items() if m in by_resi and c > hi_g]
            mapped = [by_resi[m] for m in canon if m in by_resi]
            unmapped = [b for m, b in by_resi.items() if m not in canon]
            ref = entry["reference"]
            chain, desc, omitted = REFERENCE_INFO[ref]
            rec.update({
                "system": "SLA" if str(allele).startswith("SLA") else "HLA",
                "mhc_class": cls,
                "groove_plddt_mean": round(float(np.mean(groove)), 1), "groove_n": len(groove),
                "rest_plddt_mean": round(float(np.mean(rest)), 1), "rest_n": len(rest),
                "unmapped_plddt_mean": (round(float(np.mean(unmapped)), 1) if unmapped else np.nan),
                "unmapped_n": len(unmapped),
                "mapped_plddt_mean": round(float(np.mean(mapped)), 1),
                "n_model": len(by_resi), "n_mapped": len(mapped), "n_unmapped": len(unmapped),
                "reference_pdb": ref, "reference_chain": chain,
                "reference_description": desc, "omitted_from_model": omitted,
                "coverage_model": round(len(mapped) / len(by_resi), 3),
            })
        rows.append(rec)
    return pd.DataFrame(rows)


#: The two sequence-identity columns of Table S4b. The data dictionary states
#: their method -- "MAFFT identity between the two modelled sequences" over
#: "positions at which both chains have a residue" -- so they are recomputed here
#: rather than carried, whenever MAFFT is available. Without MAFFT they fall back
#: to the approved values and the provenance record says so.
IDENTITY_S4B_COLUMNS = ("percent_identity_over_shared_positions", "n_positions_both_present")

#: Historical fields describing the SUBMITTED version of the manuscript. They are
#: not current publication values and are kept in a separate provenance table so
#: they cannot be read as such. See ``tables.table_s4b_provenance``.
HISTORICAL_S4B_COLUMNS = (
    "hla_allele_in_legend", "allele_match", "rmsd_reported_in_manuscript",
    "delta_cealign_minus_reported", "delta_super_minus_reported",
    "delta_align_minus_reported",
)

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def modelled_sequence(path) -> str:
    """The modelled chain A sequence, in residue order, from the C-alpha records."""
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == "A":
            out.append(AA3.get(line[17:20].strip().upper(), "X"))
    return "".join(out)


def _mafft() -> str | None:
    """Locate MAFFT: on PATH, or beside the running interpreter.

    A conda/micromamba environment puts its tools next to the python binary but
    does not add them to PATH unless the environment is activated, so an
    unactivated `env/bin/python -m slahla_pub.structural` would otherwise miss a
    MAFFT that is installed and usable.
    """
    import os
    import shutil
    found = shutil.which("mafft")
    if found:
        return found
    beside = os.path.join(os.path.dirname(sys.executable), "mafft")
    return beside if os.path.exists(beside) and os.access(beside, os.X_OK) else None


def pair_identity(hla_path, sla_path):
    """(percent identity, positions where both chains have a residue), via MAFFT.

    The dictionary defines these as MAFFT identity between the two modelled
    sequences over the positions at which both chains have a residue. Returns
    (None, None) if MAFFT is unavailable, so the caller can fall back and record
    that it did.
    """
    exe = _mafft()
    if exe is None:
        return None, None
    import subprocess
    import tempfile
    a, b = modelled_sequence(hla_path), modelled_sequence(sla_path)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        fh.write(f">hla\n{a}\n>sla\n{b}\n")
        path = fh.name
    res = subprocess.run([exe, "--quiet", "--auto", path], capture_output=True, text=True)
    if res.returncode != 0:
        return None, None
    seqs, name = {}, None
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].strip(); seqs[name] = ""
        elif name:
            seqs[name] += line.strip()
    both = [(x, y) for x, y in zip(seqs["hla"], seqs["sla"]) if x != "-" and y != "-"]
    if not both:
        return None, None
    same = sum(1 for x, y in both if x.upper() == y.upper())
    return round(100.0 * same / len(both), 2), len(both)


def superpositions() -> pd.DataFrame:
    """Table S4b, recomputed. Requires PyMOL.

    Faithful port of ``recompute_structural_rmsd.py`` from the manuscript support
    package; only the paths differ.
    """
    try:
        from pymol import cmd
    except ImportError as exc:                                   # pragma: no cover
        raise SystemExit(
            "PyMOL is required to recompute Table S4b. Install the optional "
            "environment:  micromamba env create -f environments/structures.yml\n"
            "The default reproduction workflow does not need this; it formats the "
            "archived measurements instead.") from exc
    refmap = reference_map()
    rows = []
    for p in PANELS:
        cmd.reinitialize()
        hp, sp = MODELS / p["hf"], MODELS / p["sf"]
        cmd.load(str(hp), "hla"); cmd.load(str(sp), "sla")
        cmd.create("h", "hla and chain A and polymer")
        cmd.create("s", "sla and chain A and polymer")
        n_hla_ca = cmd.count_atoms("h and name CA")
        n_sla_ca = cmd.count_atoms("s and name CA")
        n_hla_at = cmd.count_atoms("h")
        n_sla_at = cmd.count_atoms("s")
        hi = cmd.get_model("h and name CA").atom
        si = cmd.get_model("s and name CA").atom
        rec = dict(
            record_level="model pair",
            final_figure=p["final"], panel=p["panel"], panel_name=p["pn"], mhc_class=p["cls"],
            hla_allele_in_model=p["hla"], sla_allele_in_model=p["sla"],
            hla_allele_in_legend=p.get("legend_allele", p["hla"]),
            allele_match=p["match"], rmsd_reported_in_manuscript=p["reported"],
            hla_model_file=p["hf"], sla_model_file=p["sf"],
            hla_model_sha256=sha256(hp), sla_model_sha256=sha256(sp),
            hla_modelled_residue_range=f'{hi[0].resi}-{hi[-1].resi}',
            sla_modelled_residue_range=f'{si[0].resi}-{si[-1].resi}',
            n_hla_CA=n_hla_ca, n_sla_CA=n_sla_ca, n_hla_atoms=n_hla_at, n_sla_atoms=n_sla_at,
            predictor=("ColabFold (AlphaFold2, monomer, rank_001)" if "rank_001" in p["hf"]
                       else "authors' CUSTOM_B model (B-factor column overwritten)"),
        )
        ce = cmd.cealign("h", "s")
        rec.update(rmsd_cealign_CA=round(ce["RMSD"], 3), n_aligned_CA_pairs=ce["alignment_length"],
                   coverage_of_hla_chain_CA=round(ce["alignment_length"] / n_hla_ca, 4),
                   coverage_of_sla_chain_CA=round(ce["alignment_length"] / n_sla_ca, 4))
        sup = cmd.super("s", "h")
        rec.update(rmsd_super_all_atom=round(sup[0], 3), n_super_ATOMS=sup[1],
                   n_super_cycles=5, super_outlier_cutoff_A=2.0,
                   coverage_of_sla_chain_ATOMS=round(sup[1] / n_sla_at, 4))
        aln = cmd.align("s", "h")
        rec.update(rmsd_align_CA=round(aln[0], 3), n_align_ATOMS=aln[1],
                   align_raw_rmsd_before_refinement=round(aln[3], 3),
                   n_align_residues_aligned=aln[6])
        for dom, (lo, hi_) in DOMAINS[p["cls"]].items():
            hs, hn = _resi_spec(refmap, p["hla"], lo, hi_)
            ss, sn = _resi_spec(refmap, p["sla"], lo, hi_)
            if hs and ss:
                cmd.create("hd", f"h and resi {hs}")
                cmd.create("sd", f"s and resi {ss}")
                d = cmd.cealign("hd", "sd")
                rec[f"rmsd_{dom}_cealign_CA"] = round(d["RMSD"], 3)
                rec[f"n_{dom}_aligned_CA_pairs"] = d["alignment_length"]
                rec[f"n_{dom}_residues_selected_hla"] = hn
                rec[f"n_{dom}_residues_selected_sla"] = sn
                cmd.delete("hd"); cmd.delete("sd")
            else:
                rec[f"rmsd_{dom}_cealign_CA"] = None
        pid, npos = pair_identity(hp, sp)
        rec["percent_identity_over_shared_positions"] = pid
        rec["n_positions_both_present"] = npos
        rec["identity_method"] = ("MAFFT --auto over positions where both chains have a residue"
                                  if pid is not None else "not recomputed (MAFFT unavailable)")
        # What the final panel prints, and whether it agrees with the measurement.
        rec["rmsd_displayed_in_final_figure"] = round(rec["rmsd_cealign_CA"], 1)
        for k in ("cealign_CA", "super_all_atom", "align_CA"):
            short = k.split("_")[0]
            rec[f"delta_{short}_minus_reported"] = round(rec[f"rmsd_{k}"] - p["reported"], 3)
        rows.append(rec)
        print(f'{p["pn"]:15} reported {p["reported"]:.1f}  cealign {rec["rmsd_cealign_CA"]:.3f}'
              f'  super {rec["rmsd_super_all_atom"]:.3f}  align {rec["rmsd_align_CA"]:.3f}', flush=True)
    keys = list(rows[0])
    for r in rows:
        for k in keys:
            r.setdefault(k, None)
    return pd.DataFrame(rows)


def reference_superpositions() -> pd.DataFrame:
    """Each deposited model against its experimental reference chain. Requires PyMOL.

    The reference structures are complete complexes; the models are monomers, so
    the comparison is restricted to the reference's own MHC chain. ``coverage_*``
    state the denominators explicitly, because a coverage figure without its
    denominator is the easiest number here to misread.
    """
    try:
        from pymol import cmd
    except ImportError as exc:                                   # pragma: no cover
        raise SystemExit("PyMOL is required; see environments/structures.yml") from exc
    refmap = reference_map()
    rows = []
    for p in sorted(MODELS.glob("*.pdb")):
        allele = allele_of(p.name)
        entry = refmap.get(allele) if allele else None
        if not entry:
            continue
        cls = entry["mhc_class"]
        ref = entry["reference"]
        chain, _, _ = REFERENCE_INFO[ref]
        cmd.reinitialize()
        cmd.load(str(p), "mdl"); cmd.load(str(REFERENCE / f"{ref}.pdb"), "ref")
        cmd.create("m", "mdl and chain A and polymer")
        cmd.create("r", f"ref and chain {chain} and polymer")
        n_reference = cmd.count_atoms("r and name CA")
        ce = cmd.cealign("r", "m")
        rec = {"model": p.name, "rmsd_ce": round(ce["RMSD"], 3),
               "n_aligned": ce["alignment_length"], "n_reference": n_reference}
        for dom, (lo, hi) in DOMAINS[cls].items():
            # The model is selected through the reference map, because its own
            # residue numbering is arbitrary. The reference structure is already
            # in canonical numbering, so its domain is the plain range.
            spec, _ = _resi_spec(refmap, allele, lo, hi)
            if spec:
                cmd.create("md", f"m and resi {spec}")
                cmd.create("rd", f"r and resi {lo}-{hi}")
                d = cmd.cealign("rd", "md")
                rec[f"{dom}_rmsd"] = round(d["RMSD"], 3)
                rec[f"{dom}_n_aligned"] = d["alignment_length"]
                cmd.delete("md"); cmd.delete("rd")
        rows.append(rec)
    return pd.DataFrame(rows)


def panel_audit() -> pd.DataFrame:
    """The nine displayed panels: measured, displayed, and whether they agree.

    The number each panel prints is the whole-chain C-alpha CE RMSD from the
    render report, rounded to one decimal. This audit re-derives it from the
    recomputed superpositions and checks the three independent records agree:
    the render report the figure builder reads, the final structural manifest,
    and the recomputed measurement.

    ``class1_F_alt`` is deliberately absent: it is a historical comparison
    against the allele the *submitted* legend named, not a displayed panel.
    """
    rep = json.loads((PANEL_DIR / "render_report.json").read_text())
    man = pd.read_csv(SOURCE_DATA / "structural_panel_manifest_v6.csv")
    recomputed = TABLES / "structural_superpositions_recomputed.csv"
    sup = pd.read_csv(recomputed) if recomputed.exists() else None
    rows = []
    for name in [k for k in rep]:
        m = man[man.panel_name == name]
        if not len(m):
            continue
        m = m.iloc[0]
        ce_report = float(rep[name]["ce_rmsd"])
        rec = {
            "final_figure": m.final_figure, "panel": m.panel, "panel_name": name,
            "mhc_class": m.mhc_class,
            "hla_allele": m.hla_allele, "sla_allele": m.sla_allele,
            "ce_rmsd_render_report": round(ce_report, 4),
            "ce_rmsd_manifest_adopted": round(float(m.rmsd_adopted), 4),
            "rmsd_displayed_in_final_figure": round(ce_report, 1),
        }
        if sup is not None and (sup.panel_name == name).any():
            r = sup[sup.panel_name == name].iloc[0]
            rec["ce_rmsd_recomputed"] = round(float(r.rmsd_cealign_CA), 4)
            rec["n_aligned_CA_pairs"] = int(r.n_aligned_CA_pairs)
        else:
            rec["ce_rmsd_recomputed"] = None
            rec["n_aligned_CA_pairs"] = None
        vals = [rec["ce_rmsd_render_report"], rec["ce_rmsd_manifest_adopted"]]
        if rec["ce_rmsd_recomputed"] is not None:
            vals.append(rec["ce_rmsd_recomputed"])
        rec["max_abs_spread_between_records_A"] = round(max(vals) - min(vals), 4)
        rec["displayed_equals_round1_of_measurement"] = (
            rec["rmsd_displayed_in_final_figure"] == round(vals[-1], 1))
        rec["agreement"] = ("agree" if rec["max_abs_spread_between_records_A"] <= 0.001
                            and rec["displayed_equals_round1_of_measurement"] else "REVIEW")
        rec["record_level"] = "displayed panel"
        rows.append(rec)
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    ensure_outputs()
    dest = TABLES / "structural_model_confidence_recomputed.csv"
    no_pymol = "--no-pymol" in argv
    complete = dest.exists() and "rmsd_ce" in pd.read_csv(dest, nrows=0).columns
    if no_pymol and complete:
        # A previous `make structures` produced the full table, including the
        # reference superpositions. The pLDDT-only pass must not clobber it.
        print("confidence: keeping the complete table from an earlier `make structures`")
    else:
        conf = confidence()
        conf.to_csv(dest, index=False)
        print(f"confidence: {len(conf)} deposited models scored from their B-factor column"
              + (" (pLDDT only; the reference superpositions need PyMOL)" if no_pymol else ""))
    if not no_pymol:
        refsup = reference_superpositions()
        conf = conf.merge(refsup, on="model", how="left")
        conf["coverage_reference"] = (conf.n_mapped / conf.n_reference).round(3)
        conf.to_csv(dest, index=False)
        print(f"reference superpositions: {len(refsup)} models against "
              f"{refsup.n_reference.nunique()} experimental reference chains")
        sup = superpositions()
        sup.to_csv(TABLES / "structural_superpositions_recomputed.csv", index=False)
        print(f"superpositions: {len(sup)} model pairs")
    audit = panel_audit()
    audit.to_csv(REPORTS / "structural_panel_audit.csv", index=False)
    bad = audit[audit.agreement != "agree"]
    print(f"panel audit: {len(audit)} displayed panels, "
          f"{len(audit) - len(bad)} agree" + (f", {len(bad)} need review" if len(bad) else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
