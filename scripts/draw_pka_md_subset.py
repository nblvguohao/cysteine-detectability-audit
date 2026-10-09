"""Draw the pre-declared protein subset for the pKa round's ensemble stage.

The rule is quoted from `protocols/pka_round_preregistration_2026-09-17.json`, fixed before any pKa
value was computed:

    every protein carrying at least one positive site in supplementary table 1, capped at 120
    proteins by a seeded random draw (seed 20260915) if more qualify, with all cysteines of each
    selected protein carried along - positives and in-protein background alike.

**Why the draw is seeded and made before instrument A is read.** Selecting proteins by the static
calculation's verdict would turn the ensemble stage into a confirmation of the static stage. The
draw here depends only on (a) which proteins carry a table-1 positive and (b) whether an AlphaFold
model exists, neither of which is a pKa value.

Inputs
  results/doulias2010_sno_sites_alltables.csv   parsed supplementary tables 1-4
  external/alphafold_mouse/pdb/                 AlphaFold DB v6 models, one per accession

Exclusions applied HERE (structural availability only, not chemistry):
  * an accession with no model is dropped and listed in the audit - P11352 (GPX1_MOUSE) is the
    known case, absent from AlphaFold DB because its sequence carries selenocysteine;
  * nothing else is dropped. Metal-coordinated cysteines (supplementary table 3) and cysteines in
    predicted disulfides are excluded from the pKa STATISTICS, not from the simulations, so that
    the exclusion can be applied identically to both instruments afterwards.

Writes results/pka_md_subset.csv (one row per selected protein) and
results/pka_md_subset_audit.json, plus results/pka_md_subset_accessions.txt for the job queue.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import random
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
MODELS = os.path.join(ROOT, "external", "alphafold_mouse", "pdb")
SCRIPT = os.path.abspath(__file__)
SITES = os.path.join(RESULTS, "doulias2010_sno_sites_alltables.csv")
TABLE1 = "st01_wt_mouse_liver"
SEED = 20260915
CAP = 120


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    started = time.time()
    with open(SITES, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))

    positives = [r for r in rows if r["table"] == TABLE1 and r["position"].strip()]
    by_protein = {}
    for r in positives:
        by_protein.setdefault(r["accession_normalised"], []).append(int(r["position"]))

    have_model, missing_model = [], []
    for acc in sorted(by_protein):
        path = os.path.join(MODELS, f"AF-{acc}-F1-model_v6.pdb")
        (have_model if os.path.exists(path) else missing_model).append(acc)

    rng = random.Random(SEED)
    selected = sorted(have_model) if len(have_model) <= CAP else sorted(rng.sample(sorted(have_model), CAP))

    out_rows = []
    for acc in selected:
        path = os.path.join(MODELS, f"AF-{acc}-F1-model_v6.pdb")
        residues = set()
        n_cys = 0
        with open(path) as fh:
            for line in fh:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    residues.add(int(line[22:26]))
                    if line[17:20].strip() == "CYS":
                        n_cys += 1
        sites = sorted(set(by_protein[acc]))
        out_rows.append({
            "accession": acc,
            "n_positive_sites_table1": len(sites),
            "positive_positions": ";".join(str(s) for s in sites),
            "n_residues_in_model": len(residues),
            "n_cysteines_in_model": n_cys,
            "positions_beyond_model_length": ";".join(
                str(s) for s in sites if s > max(residues, default=0)),
            "model": f"external/alphafold_mouse/pdb/AF-{acc}-F1-model_v6.pdb",
            "model_sha256": sha256_of(path),
        })

    cols = list(out_rows[0].keys())
    path = os.path.join(RESULTS, "pka_md_subset.csv")
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)
    os.replace(tmp, path)

    listing = os.path.join(RESULTS, "pka_md_subset_accessions.txt")
    tmp = listing + ".tmp"
    open(tmp, "w").write("\n".join(r["accession"] for r in out_rows) + "\n")
    os.replace(tmp, listing)

    mismatched = [r for r in out_rows if r["positions_beyond_model_length"]]
    audit = {
        "script": "scripts/draw_pka_md_subset.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "registration": "protocols/pka_round_preregistration_2026-09-17.json",
        "rule": "proteins with at least one table-1 positive site and an AlphaFold model, capped at "
                f"{CAP} by a seeded draw (seed {SEED})",
        "drawn_before_any_pka_was_computed": True,
        "n_proteins_with_a_table1_positive": len(by_protein),
        "n_with_model": len(have_model),
        "n_without_model": len(missing_model),
        "accessions_without_model": missing_model,
        "cap_applied": len(have_model) > CAP,
        "n_selected": len(out_rows),
        "n_positive_sites_in_selection": sum(r["n_positive_sites_table1"] for r in out_rows),
        "n_cysteines_in_selected_models": sum(r["n_cysteines_in_model"] for r in out_rows),
        "proteins_with_a_site_beyond_the_model_length": [
            {"accession": r["accession"], "positions": r["positions_beyond_model_length"],
             "n_residues_in_model": r["n_residues_in_model"]} for r in mismatched],
        "inputs": {"results/doulias2010_sno_sites_alltables.csv": sha256_of(SITES)},
        "versions": {"python": sys.version},
        "elapsed_seconds": round(time.time() - started, 2),
    }
    path = os.path.join(RESULTS, "pka_md_subset_audit.json")
    tmp = path + ".tmp"
    json.dump(audit, open(tmp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    os.replace(tmp, path)

    print(f"proteins with a table-1 positive: {len(by_protein)}; with a model: {len(have_model)}; "
          f"selected: {len(out_rows)}")
    print(f"positive sites in the selection: {audit['n_positive_sites_in_selection']}; "
          f"cysteines in those models: {audit['n_cysteines_in_selected_models']}")
    if missing_model:
        print("no model:", ", ".join(missing_model))
    if mismatched:
        print("sites beyond the model length:", [(r['accession'], r['positions_beyond_model_length'])
                                                 for r in mismatched])


if __name__ == "__main__":
    main()
