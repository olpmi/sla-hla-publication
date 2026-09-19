# Reproduction entry points for the SLA/HLA manuscript.
#
#   make reproduce   regenerate every published figure and table, then verify
#   make verify      compare what was regenerated against the approved materials
#   make test        the retained test set
#   make structures  OPTIONAL: recompute Tables S4a/S4b from the model files (PyMOL)
#
# Everything runs on CPU. Nothing downloads, trains, or predicts structures.

PYTHON ?= python
export PYTHONPATH := src

MODULES := ensemble eplet_stats tables figures_final
OUT     := outputs
SNAPSHOT := $(OUT)/reports/source_data_before.json

.PHONY: help reproduce classifier eplets figures tables verify checks audit mask workbook test structures clean

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
	  | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

$(SNAPSHOT):
	@mkdir -p $(OUT)/reports
	@$(PYTHON) -c "import json,hashlib,pathlib; r=pathlib.Path('data/source_data'); \
print(json.dumps({str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() \
for p in sorted(r.rglob('*')) if p.is_file()}))" > $@

classifier: ## recompute the five-seed ensemble from the archived predictions
	$(PYTHON) -m slahla_pub.ensemble

eplets: ## recompute the eplet-position statistics
	$(PYTHON) -m slahla_pub.eplet_stats

tables: ## assemble publication Tables S1-S6
	$(PYTHON) -m slahla_pub.tables

figures: ## regenerate Figures 1-6 and S1-S10
	$(PYTHON) -m slahla_pub.figures_final

structures: ## OPTIONAL: recompute S4a/S4b from the deposited models (needs PyMOL)
	$(PYTHON) -m slahla_pub.structural

audit: ## the nine-panel structural audit, S1 cell differences and column provenance
	$(PYTHON) -m slahla_pub.structural --no-pymol
	$(PYTHON) -m slahla_pub.s1_audit
	$(PYTHON) -m slahla_pub.output_provenance

mask: ## report whether the locally supplied eplet position mask is present
	$(PYTHON) -m slahla_pub.eplet_mask --check

checks: ## the manuscript's scientific endpoints
	$(PYTHON) -m slahla_pub.checks

reproduce: $(SNAPSHOT) ## the documented default workflow
	@$(PYTHON) -c "from slahla_pub import runlog; runlog.clear()"
	-@$(PYTHON) -m slahla_pub.eplet_mask --check
	$(PYTHON) -m slahla_pub.ensemble
	$(PYTHON) -m slahla_pub.eplet_stats
	$(PYTHON) -m slahla_pub.structural --no-pymol
	$(PYTHON) -m slahla_pub.tables
	$(PYTHON) -m slahla_pub.figures_final
	$(PYTHON) -m slahla_pub.s1_audit
	$(PYTHON) -m slahla_pub.output_provenance
	-$(PYTHON) scripts/build_corrected_workbook.py
	$(MAKE) verify

verify: ## compare outputs with the approved materials; writes the report
	$(PYTHON) -m slahla_pub.verify --inputs-before=$(SNAPSHOT)
	$(PYTHON) -m slahla_pub.report

workbook: ## corrected copy of the supplementary workbook (originals untouched)
	$(PYTHON) scripts/build_corrected_workbook.py

test: ## the retained test set
	$(PYTHON) -m pytest -q

clean: ## remove regenerated outputs only (never touches data/ or publication_artwork/)
	rm -rf $(OUT)
