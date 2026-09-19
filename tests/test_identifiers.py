"""Identifier parsing. A bug in any of these silently mislabels a published row."""
from __future__ import annotations

import pandas as pd
import pytest

from slahla_pub import haplotypes as H
from slahla_pub.ensemble import _top_level, label_space, probs_path, run_tag
from slahla_pub.structural import allele_of
from slahla_pub.paths import SOURCE_DATA


@pytest.mark.parametrize("label,locus", [
    ("A01", "A"), ("B40", "B"), ("C08", "C"), ("C_Other", "C"),
])
def test_class1_label_collapses_to_its_locus(label, locus):
    labels, _ = label_space("class1")
    assert _top_level(label, labels) == locus


@pytest.mark.parametrize("label,locus", [
    ("DQB102", "DQB1"), ("DRB101", "DRB1"), ("DRB3_Other", "DRB3"),
    ("DRB4_Other", "DRB4"), ("DPB1_Other", "DPB1"),
])
def test_class2_label_collapses_to_its_locus(label, locus):
    """DRB1/DRB3/DRB4/DRB5 share a prefix; a shortest-match rule would fold them
    together and silently merge four loci into one."""
    labels, _ = label_space("class2")
    assert _top_level(label, labels) == locus


def test_label_space_matches_the_archived_prediction_columns():
    import pandas as pd
    for which in ("class1", "class2"):
        labels, _ = label_space(which)
        d = pd.read_parquet(probs_path(which, run_tag(42)))
        assert not set(labels) - set(d.columns)


@pytest.mark.parametrize("filename,allele", [
    ("HLA_A02_A_02_01_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb", "A*02:01"),
    ("SLA_A02_SLA-1_02_01_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb", "SLA-1*02:01"),
    ("HLA_DRB101_DRB1_01_02_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb", "DRB1*01:02"),
    ("SLA_DQB102_SLA-DQB1_02_02_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb",
     "SLA-DQB1*02:02"),
    ("C07_C0701_unrelaxed_CUSTOM_B.pdb", None),
])
def test_model_filename_to_allele(filename, allele):
    assert allele_of(filename) == allele


def test_haplotype_name_composition():
    """Class I ``1a.0`` plus class II ``0.1`` is ``Hp-1a.1``, not ``Hp-1a.0.0.1``."""
    assert H._compose_name("1a.0", "0.1") == "Hp-1a.1"
    assert H._compose_name("22.0", "0.15b") == "Hp-22.15b"


def test_haplotype_cell_with_two_alleles_keeps_both():
    assert H._split_alleles("09:01, 15:01") == ["09:01", "15:01"]
    assert H._split_alleles("Null") == []
    assert H._split_alleles("NT") == []


def test_haplotype_workbook_covers_every_published_table_s1_allele():
    from slahla_pub.paths import ARTWORK
    s1 = pd.read_csv(ARTWORK / "tables" /
                     "TableS1_haplotype_nearest_and_most_distant_HLA.csv")
    parsed = set(H.alleles_by_haplotype().sla_allele)
    assert not set(s1.sla_allele) - parsed
