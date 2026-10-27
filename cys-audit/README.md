# Cys-Audit

Audit a cysteine-modification site table for five measurement-driven artefacts before believing its
site preferences.

| test | artefact | statistic | needs |
|---|---|---|---|
| `cleavage_geometry` | protease cleavage geometry | log2 OR of a cleavage-competent residue within 1-3 and 6-12 residues, positive vs background | `--fasta` |
| `multi_cysteine` | multi-cysteine peptides | log2 OR of multi-Cys-only support, positive vs observed background | `n_cys_peptide`, `detected` |
| `protein_abundance` | protein abundance | rank AUC of log abundance, positive vs background | `abundance` |
| `positive_detected_overlap` | positive set = detected set | share of detected cysteines that are positives | `detected` |
| `search_space` | isobaric alternatives absent from the search | resolving mass of the closest alternative | `--precursor-ppm`, `--search-mods` |

Each test returns **PASS / WARNING / FAIL / UNDECIDABLE** from its interval and a fixed equivalence margin
(see `src/cys_audit/status.py`); a p-value never sets a status. Intervals are protein-clustered bootstrap,
5000 replicates, fixed seed. An optional claim retest (`--claim-feature`) returns
**PASS / ATTENUATED / UNDECIDABLE / FAIL** for one binary site feature after stratifying by every
detectability covariate the input supports, or (v0.2.0) on propensity quintiles of covariates you name with
`--claim-covariates COL1,COL2,...` (ridge logistic regression, C = 1). Without named covariates the control uses only
the built-in covariates, which can be far too coarse: in the authors' Phase 5 cross-check, a cleavage-band-only
control reproduced none of the five control-decided non-PASS verdicts of 18 published site-level claims.

## Install and run

```bash
pip install -e .            # or: export PYTHONPATH=src
make test
cys-audit schema
cys-audit audit --input sites.tsv --fasta proteome.fasta --protease trypsin \
    --modification persulfidation --background observed \
    --precursor-ppm 4.5 --search-mods Sulfide,Oxidation --output audit_report
```

Output: `audit_report/{audit.json, audit.html, summary.csv, figures/*.svg, manifest.json}`. Files contain no
timestamps or absolute paths; the same input and seed give byte-identical outputs, and `manifest.json` records
the sha256 of every input and output.

## What a status does and does not mean

A status describes one dataset on one axis. PASS means an artefact of meaningful size is excluded on that axis
with the given data; it does not certify a biological claim. FAIL means an artefact signature of meaningful
size is established; it does not show the claim is false, only that the claim needs a control on that axis.
The margins are conventions fixed in `constants.py` before any real dataset was audited.

## Versions

- 0.2.2: `search_space` carries one non-cysteine pair, phosphorylation against sulfation (delta 0.009515 Da, derived from element masses); no rule, threshold or statistic changes, and v0.2.1 is archived byte for byte in `archive/v0.2.1/`.
- 0.2.1: the size-matched random control is replaced by a permuted-stratum control when stratification
  retains more than 95% of the sample, where a size-matched draw cannot be a null; criteria registered in
  `protocols/cys_audit_random_control_preregistration_2026-09-22.json`. Archived in `archive/v0.2.0/`.
- 0.2.0: `--claim-covariates` (propensity-quintile claim control); HTML shows '-' for undefined counts.
- 0.1.0: first release; archived byte-for-byte in `archive/v0.1.0/`.

## Limits

- `search_space` knows only the isobaric pairs in its built-in table (+S / +2O, +3O / +S+O); anything else is
  reported UNDECIDABLE ("not covered"), never PASS.
- The claim retest controls only for covariates present in the input; an unmeasured detectability factor is not
  controlled.
- Validation so far: unit tests, synthetic calibration, one exact reproduction of a published-analysis number,
  and one public demonstration dataset (see the Phase 4 report in the analysis tree).
