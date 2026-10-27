# Minimal computational reproduction code

This directory contains code and dependency pins only. No experimental or derived data, predictions, fitted models, figures, manuscripts, or saved analysis results are distributed in this update. Inputs must be supplied separately. Nothing is automatically downloaded or published.

The five modules in `code/` preserve the original independent reconstruction source byte for byte. `run.py` is a new filesystem adapter: it accepts separately supplied prepared site membership, full protein sequences and fixed group folds, then calls the unchanged feature construction, model fitting and paired-unit bootstrap. It does not ingest MS spectra, reconstruct source-table labels, or infer chemical negatives.

## Install

Use a fresh environment with Python **3.10** for the strongest numerical comparability with the archived execution. No GPU is required. From this directory:

```sh
python3.10 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

Run with the explicit environment interpreter. The runner uses two OpenMP threads and one BLAS thread. It never overwrites an existing output directory.

## Separately supplied inputs

| Input | Required content |
| --- | --- |
| `--sites` | Tab-separated file with exactly `acc, site, unit, is_pos, in_A, in_B` in that order. Positions are one-based cysteines. Membership fields are 0/1; positives cannot also be references. Each accession-position pair is unique. |
| `--fasta` | Full sequences; the first token after each `>` is exactly the corresponding `acc`. Use the same sequence version as the prepared sites. |
| `--folds` | Tab-separated file with the same six identity/membership columns, in the same row order, plus integer `new_fold` values 0–4. Each group `unit` must belong to one fold. Other columns, including any historical predictions, are ignored. |

`in_A` and `in_B` are supplied reference memberships, not model-derived truth labels. Their experimental definitions must be documented by the data provider. Observed/unassigned is not equivalent to chemically unmodified. Fixed folds must be supplied: regenerating GroupKFold ties across hardware can change the split. There is no bundled synthetic or experimental dataset.

For author-held reproduction, the prepared `sites.tsv`, `reference.fasta` and archived `folds.tsv.gz` are in the local complete analysis backup. Those files are deliberately not uploaded to this repository. This code-only repository cannot reproduce the article's numerical values without those separately acquired, version-matched inputs.

## Execute

```sh
.venv/bin/python run.py \
  --sites /data/sites.tsv \
  --fasta /data/reference.fasta \
  --folds /data/fixed_folds.tsv.gz \
  --output /data/new_reproduction
```

Defaults run both `strict_mc_v1` and `probe_compatible_v1`, with the VIS10 feature set and 2,000 paired-unit bootstrap draws. Other prespecified feature sets can be requested with `--feature-sets VIS10 DIG25 CPL15`. Each definition/feature-set combination fits ten models: two training backgrounds × five folds.

All newly generated features, models, fold records, prediction arrays, held-out predictions, estimates and uncertainty intervals are written only to the supplied output directory. An execution identity records the code, consumed inputs, environment and analysis options. Bootstrap draws condition on fitted held-out predictions; they do not refit the models.

The runner performs analytical boundary and weighted-AUC checks before fitting. `tests/test_runner.py` generates temporary synthetic cases in memory/on a temporary filesystem to check the adapter; these are software tests, not biological validation or article results. Run `python -m unittest discover -s tests` with the scientific dependencies installed.

## Scope

This is an explicitly versioned adapter, not a claim of byte-identical rerunning of every historical analysis. The scientific implementation is unchanged, but the adapter's execution identity is new. Python version, feature arithmetic and fixed folds matter for numerical comparison. AV1 source-workbook processing, plotting and manuscript production are outside this minimum entry. The repository's older release folders predate this update and retain their historical scope.
