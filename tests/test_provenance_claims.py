"""Reproduction claims must separate recomputed values from carried ones."""
from __future__ import annotations

import pandas as pd

from slahla_pub import output_provenance, s1_audit


def test_every_published_column_is_classified():
    d = output_provenance.classify()
    assert set(d.provenance) <= {"recomputed", "carried", "external_not_executed", "key"}
    for t in ["S1", "S2", "S3", "S4a", "S4b", "S5", "S6"]:
        assert (d.table == t).any(), f"no provenance rows for {t}"


def test_identity_columns_of_table_s1_are_not_claimed_as_recomputed():
    """They need DIAMOND over IPD-IMGT/HLA and CIWD, which the clean run does not have."""
    d = output_provenance.classify()
    s1 = d[d.table == "S1"].set_index("column")
    for c in s1_audit.CARRIED_COLUMNS:
        assert s1.loc[c, "provenance"] == "external_not_executed"


def test_registry_dependent_columns_of_table_s6_are_not_claimed_as_recomputed():
    d = output_provenance.classify()
    s6 = d[d.table == "S6"].set_index("column")
    for c in ["eplet_names", "eplet_loci_defining", "eplet_name_source"]:
        assert s6.loc[c, "provenance"] == "external_not_executed"


def test_table_s1_differences_are_enumerated_with_a_cause():
    d = s1_audit.compare()
    assert set(d.columns) >= {"sla_allele", "column", "approved_value",
                              "regenerated_value", "cause", "affects_a_manuscript_result"}
    assert (d.cause.str.len() > 0).all()
    # no difference may be left unexplained
    assert not (d.classification == "unclassified").any()


def test_no_table_s1_difference_affects_a_manuscript_result():
    d = s1_audit.compare()
    assert not len(d) or (d.affects_a_manuscript_result == "none").all()


def test_local_name_follows_the_selected_haplotype():
    """The five previously-differing local_name cells are resolved."""
    d = s1_audit.compare()
    assert not (d.column == "local_name").any()
