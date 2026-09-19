"""The registry importer: the locus is in the allele, never in the eplet name."""
from __future__ import annotations

import pandas as pd
import pytest

from slahla_pub import eplet_mask


# ------------------------------------------------- the demonstrated failure
def test_eplet_name_column_alone_is_rejected(tmp_path):
    """`eplet,resi` with an eplet named `62GE` must not quietly yield nothing.

    This is the reported defect: the old importer picked the eplet-name column as
    its locus source, matched `62GE` against `^(A|B|C)`, found nothing and wrote
    an empty mask without complaint. A file with no `Allele` column cannot say
    which locus a position belongs to, so it is refused.
    """
    (tmp_path / "ABC.csv").write_text("eplet,resi\n62GE,[[62, 63]]\n1C,[[1]]\n")
    with pytest.raises(SystemExit) as e:
        eplet_mask.from_registry(tmp_path)
    assert "Allele" in str(e.value)


def test_allele_column_after_the_eplet_column_still_resolves_the_locus(tmp_path):
    """Column order must not matter: the locus comes from `Allele` by name."""
    (tmp_path / "ABC.csv").write_text(
        "Name,Evidence,Exposition,Allele,resi\n"
        "62GE,A1,High,A*01:01,[[62]]\n"
        "1C,,High,C*01:02,[[1]]\n")
    d = eplet_mask.read_export(tmp_path)
    assert list(d.Allele) == ["A*01:01", "C*01:02"]
    assert [eplet_mask.locus_of(a) for a in d.Allele] == ["A", "C"]


def test_eplet_names_are_never_read_as_loci():
    """`62GE`, `1C` and `193PI` are eplet names; none names a locus."""
    for name in ("62GE", "1C", "193PI", "144KR+151H"):
        assert eplet_mask.locus_of(name) == name        # no '*' -> unchanged
        assert eplet_mask.locus_of(name) not in eplet_mask.CLASS_LOCI["class1"]


# --------------------------------------------------------- the real schema
def test_real_export_schema_is_the_supported_one():
    assert eplet_mask.REQUIRED_COLUMNS == ("Name", "Evidence", "Exposition", "Allele", "resi")


@pytest.mark.parametrize("allele,locus", [
    ("C*01:02", "C"), ("A*01:01", "A"), ("B*07:02", "B"),
    ("DQA1*01:01", "DQA1"), ("DQB1*02:01", "DQB1"),
    ("DRB1*01:01", "DRB1"), ("DRB3*01:01", "DRB3"), ("DPB1*01:01", "DPB1"),
])
def test_locus_is_parsed_from_the_allele(allele, locus):
    assert eplet_mask.locus_of(allele) == locus


def test_class_loci_are_the_analysis_loci():
    """Class II is the beta chains the cohort has comparators for; DP and the
    alpha chains are excluded, and including DP would add five positions."""
    assert eplet_mask.CLASS_LOCI["class1"] == ("A", "B", "C")
    assert eplet_mask.CLASS_LOCI["class2"] == ("DRB1", "DRB3", "DRB4", "DRB5", "DQB1")
    assert "DPB1" not in eplet_mask.CLASS_LOCI["class2"]
    assert "DQA1" not in eplet_mask.CLASS_LOCI["class2"]


def test_the_published_filter_is_a_disjunction():
    assert eplet_mask.EVIDENCE_KEEP == frozenset({"A1", "A2", "B"})
    assert eplet_mask.EXPOSITION_KEEP == "High"


@pytest.mark.parametrize("cell,expected", [
    ("[[1]]", [1]),
    ("[[62, 63]]", [62, 63]),
    ("[[9], [11]]", [9, 11]),
    ("[['62'], ['63']]", [62, 63]),
    ("", []),
    (None, []),
])
def test_nested_resi_is_flattened(cell, expected):
    assert eplet_mask.flatten_resi(cell) == expected


# ------------------------------------------------------- refuses bad input
def test_empty_directory_is_refused(tmp_path):
    with pytest.raises(SystemExit) as e:
        eplet_mask.from_registry(tmp_path)
    assert "no .csv" in str(e.value)


def test_table_with_no_rows_is_refused(tmp_path):
    (tmp_path / "ABC.csv").write_text("Name,Evidence,Exposition,Allele,resi\n")
    with pytest.raises(SystemExit):
        eplet_mask.from_registry(tmp_path)


def test_export_without_the_target_loci_is_refused(tmp_path):
    """A DP-only export cannot produce the published mask; say so, do not write one."""
    (tmp_path / "DP.csv").write_text(
        "Name,Evidence,Exposition,Allele,resi\n84GD,A1,High,DPB1*01:01,[[84]]\n")
    with pytest.raises(SystemExit) as e:
        eplet_mask.from_registry(tmp_path)
    assert "class1 loci" in str(e.value) or "no allele at" in str(e.value)


def test_nothing_passing_the_filter_is_refused(tmp_path):
    (tmp_path / "ABC.csv").write_text(
        "Name,Evidence,Exposition,Allele,resi\n62GE,C,Low,A*01:01,[[62]]\n")
    with pytest.raises(SystemExit) as e:
        eplet_mask.from_registry(tmp_path)
    assert "filter" in str(e.value)


# ------------------------------------------- end-to-end against the real export
@pytest.mark.skipif(not eplet_mask.available(),
                    reason="no locally supplied mask to compare against")
def test_generated_mask_has_the_published_membership():
    """Counts are not enough: the position sets themselves must match."""
    m = eplet_mask.load()
    assert len(m["class1"]) == 61 and len(m["class2"]) == 52
    agr_positions = set(pd.read_parquet(
        eplet_mask.EPLET_AGREEMENT / "identity_eplet_agreement.parquet").position)
    for cls, s in m.items():
        assert s <= agr_positions, f"{cls} names positions the analysis never scores"
