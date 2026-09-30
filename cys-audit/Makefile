# make test       run the unit tests
# make demo       simulate a dataset with a planted cleavage artefact and audit it
# make reproduce  rebuild the public demonstration (needs the PRIDE files; see data/README.md)
PY ?= python
export PYTHONPATH := src
export PYTHONDONTWRITEBYTECODE := 1
DATA ?= data/PXD063463

test:
	$(PY) -m unittest discover -s tests -t tests

demo:
	$(PY) -m cys_audit simulate --artefact cleavage --strength 1 --out-dir examples/synthetic_cleavage
	$(PY) -m cys_audit audit --input examples/synthetic_cleavage/sites.tsv --fasta examples/synthetic_cleavage/proteome.fasta \
	  --protease trypsin --background proteome --output examples/synthetic_cleavage/audit_report

reproduce:
	$(PY) examples/convert_maxquant_abe.py --data-dir $(DATA) --enzyme Trypsin --hydroxylamine P --out examples/pxd063463_trypsin.tsv
	$(PY) -m cys_audit audit --input examples/pxd063463_trypsin.tsv --fasta $(DATA)/Mouse_2025_UP000000589_10090.fasta \
	  --expand-background --protease trypsin --modification S-palmitoylation --identity-readout chemistry_inferred \
	  --background proteome --dataset PXD063463_trypsin --output examples/pxd063463_trypsin_report

.PHONY: test demo reproduce
