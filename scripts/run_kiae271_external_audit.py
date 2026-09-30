"""Zero-shot external audit of the deployed v2 rankers on the collaborator's own
persulfidome (Plant Physiol 2024, doi 10.1093/plphys/kiae271, Dataset S1).

WHY THIS SET IS WORTH A RUN
---------------------------
Every external audit so far used human PXD044043 and rice PXD072089 - another species and another
workflow.  This set is same species, same tissue type, and the
same laboratory that will run the wet-lab validation, and it was confirmed on
2026-09-15 NOT to be the source of the 15,215 training labels: of 167 unique
(accession, position) pairs only 40 appear in the frozen cohort (18 positives,
22 negatives) and the proteins of the remaining 127 are absent from the cohort
entirely, SlWRKY6 (A0A3Q7F586:396) and WRKY7 (A0A3Q7GA72:230/240) included.

LABEL SEMANTICS - state it, do not overclaim it
-----------------------------------------------
kiae271's Methods say the persulfidome "was identified according to the method
of Aroca et al. (2017)" and report no modification masses, no search engine and
no data-availability accession.  Dataset S1's MaxQuant column `S(C)` implies a
sulfur-addition modification on cysteine, and Aroca's own deposition declares
Sulfide (+31.972071) with carbamidomethyl absent, so this is PROBABLY a direct
+S adduct label.  It is not proven, so this report must call it
"the collaborator's reported persulfidation sites", never "a verified direct
+S label".  Aroca's direct channels also failed decoy QC in this project's own
census, which is a further reason not to treat it as a gold standard.

PRE-DECLARED DESIGN (fixed before any score was computed)
---------------------------------------------------------
Scoring: `scripts/score_v2_fasta.py` deployment path, unchanged, both feature
sets (chem = the predeclared wet-lab ranker, full = the benchmark model).  The
models are used zero-shot; nothing is refitted on this set.

Sites: Dataset S1 rows are parsed with `Proteins` and `Positions within
proteins` treated as PARALLEL arrays (MaxQuant semantics), never as a cartesian
product.  A pair is kept only when the accession is in `inputs/sly_proteome.tsv`
and the residue at that position is a cysteine.

Two strata, reported separately and never pooled:
* `clean_external` - proteins with NO row in `inputs/primary_site_folds.csv`,
  and additionally no cysteine whose 31-aa window is identical to a training
  positive's window.  This is the stratum the conclusion may rest on.
* `cohort_overlap` - proteins that do appear in the frozen cohort.  Reported for
  completeness only; agreement here is consistency with training data.

Proteins with a single cysteine are EXCLUDED from the endpoint: the ranking has
no freedom there and a hit is automatic (the project's standing rule).

Primary endpoint: per protein, is a kiae271-reported site ranked 1 (Top-1) or
within the top 2 (Top-2) by the model, against that protein's own random
baseline (Top-1: k/n; Top-2: 1 - C(n-k,2)/C(n,2), with k reported sites among n
cysteines).  Statistic = mean over proteins of (model hit - random expectation).
Intervals: paired bootstrap over proteins, 5000 replicates, seed 20260915.

Reading rule: lower bound above 0 means the ranker beats random on this set;
an interval containing 0 means undetermined at this sample size; an upper bound
below 0 means worse than random.  A null result is a legitimate outcome.

Boundaries to carry into any report: the label is an enrichment-workflow site
list of ~119 rows, so absence of a site is not evidence of absence of
persulfidation; the negatives here are "not reported by kiae271", which is a
detection-limited class, so no ROC-AUC or AP is computed - only top-of-ranking
coverage, exactly as in `reports/V2_RETROSPECTIVE_EXTERNAL_AUDIT.md`; and
homology to training proteins is screened only by identical 31-aa windows, not
by alignment, because the deployment venv has no Biopython.
"""
from __future__ import annotations

import csv
import json
import pathlib
import re
import subprocess
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
XLSX = pathlib.Path.home() / "Downloads/kiae271_supplementary_data/DSs.xlsx"
PROTEOME = ROOT / "inputs/sly_proteome.tsv"
FOLDS = ROOT / "inputs/primary_site_folds.csv"
WORK = ROOT / "external/kiae271_audit"
OUT_SITES = ROOT / "results/kiae271_external_sites.csv"
OUT_SUMMARY = ROOT / "results/kiae271_external_summary.csv"
OUT_AUDIT = ROOT / "results/kiae271_external_audit.json"
SEED = 20260915
REPLICATES = 5000
WINDOW = 15  # +-15 gives the project's 31-aa window


def sha256(path: pathlib.Path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def window(sequence: str, position: int) -> str:
    start, end = position - 1 - WINDOW, position + WINDOW
    left = "-" * max(0, -start) + sequence[max(0, start):position - 1]
    right = sequence[position:end] + "-" * max(0, end - len(sequence))
    return left + sequence[position - 1] + right


def main() -> None:
    started = time.time()
    import openpyxl

    proteome = {}
    with PROTEOME.open(encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            proteome[row["Entry"]] = row["Sequence"]

    folds = list(csv.DictReader(FOLDS.open(encoding="utf-8-sig")))
    cohort_proteins = {r["accession"] for r in folds}
    cohort_label = {(r["accession"], int(r["position"])): int(r["label"]) for r in folds}
    positive_windows = {
        window(proteome[r["accession"]], int(r["position"]))
        for r in folds if int(r["label"]) == 1 and r["accession"] in proteome
    }

    workbook = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    rows = [r for r in list(workbook["Supplementary Dataset S1"].iter_rows(values_only=True))[2:] if r[0]]
    pairs, dropped = set(), {"parallel_length_mismatch": 0, "accession_not_in_proteome": 0,
                             "position_not_cysteine": 0}
    for row in rows:
        accessions = re.findall(r"(?:tr|sp)\|([A-Z0-9]+)\|", str(row[0] or ""))
        positions = [p.strip() for p in str(row[1] or "").split(";")]
        if len(accessions) != len(positions):
            dropped["parallel_length_mismatch"] += 1
            continue
        for accession, position in zip(accessions, positions):
            if not position.isdigit():
                continue
            if accession not in proteome:
                dropped["accession_not_in_proteome"] += 1
                continue
            index = int(position)
            sequence = proteome[accession]
            if index > len(sequence) or sequence[index - 1] != "C":
                dropped["position_not_cysteine"] += 1
                continue
            pairs.add((accession, index))

    reported = {}
    for accession, position in sorted(pairs):
        reported.setdefault(accession, set()).add(position)

    strata = {}
    for accession, sites in reported.items():
        sequence = proteome[accession]
        shares_window = any(
            window(sequence, i + 1) in positive_windows
            for i, residue in enumerate(sequence) if residue == "C"
        )
        if accession in cohort_proteins or shares_window:
            strata[accession] = "cohort_overlap"
        else:
            strata[accession] = "clean_external"

    WORK.mkdir(parents=True, exist_ok=True)
    fasta = WORK / "kiae271_proteins.fasta"
    with fasta.open("w", encoding="utf-8") as handle:
        for accession in sorted(reported):
            handle.write(f">{accession}\n{proteome[accession]}\n")

    scores = {}
    for feature_set in ("chem", "full"):
        output = WORK / f"scores_{feature_set}.csv"
        subprocess.run([sys.executable, str(ROOT / "scripts/score_v2_fasta.py"),
                        "--fasta", str(fasta), "--feature-set", feature_set,
                        "--output", str(output)], check=True, cwd=ROOT / "scripts")
        table = {}
        with output.open(encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                table.setdefault(row["protein"], []).append(row)
        scores[feature_set] = table

    site_rows, summary_rows = [], []
    per_protein = {"chem": {}, "full": {}}
    for feature_set, table in scores.items():
        for accession, rows_for_protein in table.items():
            ranked = sorted(rows_for_protein, key=lambda r: int(r["rank_in_protein"]))
            n = len(ranked)
            sites = reported[accession]
            k = len(sites)
            hit1 = int(int(ranked[0]["position"]) in sites)
            hit2 = int(any(int(r["position"]) in sites for r in ranked[:2]))
            rand1 = k / n
            rand2 = 1.0 - ((n - k) * (n - k - 1)) / (n * (n - 1)) if n >= 2 and n - k >= 2 else 1.0
            per_protein[feature_set][accession] = {
                "stratum": strata[accession], "n_cys": n, "n_reported": k,
                "hit_top1": hit1, "hit_top2": hit2, "random_top1": rand1, "random_top2": rand2,
            }
            if feature_set == "chem":
                for row in ranked:
                    position = int(row["position"])
                    key = (accession, position)
                    site_rows.append({
                        "protein": accession, "position": position,
                        "stratum": strata[accession],
                        "reported_by_kiae271": int(position in sites),
                        "chem_rank": int(row["rank_in_protein"]), "n_cys": n,
                        "in_frozen_cohort": int(key in cohort_label),
                        "frozen_label": cohort_label.get(key, ""),
                    })

    rng = np.random.default_rng(SEED)
    for feature_set in ("chem", "full"):
        for stratum in ("clean_external", "cohort_overlap", "all"):
            selected = [v for v in per_protein[feature_set].values()
                        if (stratum == "all" or v["stratum"] == stratum) and v["n_cys"] >= 2]
            if not selected:
                continue
            for metric in ("top1", "top2"):
                model = np.array([v[f"hit_{metric}"] for v in selected], dtype=float)
                random = np.array([v[f"random_{metric}"] for v in selected], dtype=float)
                delta = model - random
                draws = rng.integers(0, len(delta), size=(REPLICATES, len(delta)))
                boot = delta[draws].mean(axis=1)
                low, high = np.percentile(boot, [2.5, 97.5])
                verdict = ("beats random" if low > 0 else
                           "worse than random" if high < 0 else "undetermined at this n")
                summary_rows.append({
                    "feature_set": feature_set, "stratum": stratum, "metric": metric,
                    "n_proteins": len(selected),
                    "model": round(float(model.mean()), 4),
                    "random": round(float(random.mean()), 4),
                    "delta": round(float(delta.mean()), 4),
                    "ci_low": round(float(low), 4), "ci_high": round(float(high), 4),
                    "verdict": verdict,
                })

    for rows_out, path in ((site_rows, OUT_SITES), (summary_rows, OUT_SUMMARY)):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows_out[0]))
            writer.writeheader()
            writer.writerows(rows_out)

    single_cys = sorted(a for a, v in per_protein["chem"].items() if v["n_cys"] == 1)
    # only sites kiae271 actually reports can contradict our label; without the
    # first condition this counts every cysteine of these proteins that our
    # cohort calls negative, which is not a contradiction at all
    contradictions = sorted(
        (r["protein"], r["position"]) for r in site_rows
        if r["reported_by_kiae271"] and r["in_frozen_cohort"] and r["frozen_label"] == 0)
    OUT_AUDIT.write_text(json.dumps({
        "question": "do the deployed v2 rankers put the collaborator's reported "
                    "persulfidation sites at the top of each protein, zero-shot",
        "script_sha256": sha256(pathlib.Path(__file__)),
        "source_xlsx_sha256": sha256(XLSX),
        "label_semantics": "kiae271 Dataset S1 reported sites; Methods cite Aroca et al. 2017 and "
                           "state no modification masses, so a direct +S label is probable but unproven",
        "rows_parsed": len(rows),
        "unique_pairs_kept": len(pairs),
        "dropped": dropped,
        "proteins_scored": len(reported),
        "proteins_by_stratum": {s: sum(1 for v in strata.values() if v == s)
                                for s in set(strata.values())},
        "single_cys_proteins_excluded_from_endpoint": single_cys,
        "sites_in_cohort_labelled_negative": contradictions,
        "n_sites_in_cohort_labelled_negative": len(contradictions),
        "seed": SEED, "bootstrap_replicates": REPLICATES,
        "runtime_seconds": round(time.time() - started, 1),
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"{'set':<6}{'stratum':<16}{'metric':<7}{'n':>5}{'model':>8}{'random':>8}"
          f"{'delta':>9}  95% CI            verdict")
    for row in summary_rows:
        print(f"{row['feature_set']:<6}{row['stratum']:<16}{row['metric']:<7}{row['n_proteins']:>5}"
              f"{row['model']:>8.3f}{row['random']:>8.3f}{row['delta']:>+9.4f}  "
              f"[{row['ci_low']:+.4f}, {row['ci_high']:+.4f}]  {row['verdict']}")
    print(f"\nreported sites kept {len(pairs)} | dropped {dropped}")
    print(f"single-Cys proteins excluded: {len(single_cys)}")
    print(f"kiae271 sites that are NEGATIVES in the frozen cohort: {len(contradictions)}")


if __name__ == "__main__":
    main()
