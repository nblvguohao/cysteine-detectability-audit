"""Instrument A of the pKa round: PROPKA3 on the static AlphaFold model, one value per cysteine.

Declared in `protocols/pka_round_preregistration_2026-09-17.json` before any value was computed.
This is the instrument that REPRODUCES the authors' design (PROPKA 2.0 on crystal structures), not
the one that improves on it; its weakness - empirical, single static conformer, unreliable for
cysteine - is declared there and must travel with every number it produces.

**What this script does and does not do.** It computes pKa and writes one row per cysteine with the
annotations a later round needs. It runs NO comparison, applies NO control and issues NO verdict;
the claim is re-tested only once instrument B (the ensemble) is in, so that a verdict cannot be
quietly taken from whichever instrument was available first.

Definitions fixed before the run:

  * **Cohort** - every cysteine of every mouse protein that carries at least one site in
    supplementary table 1 and has an AlphaFold DB v6 model. Positives are the cysteines listed in
    table 1; the background is the remaining cysteines of the SAME proteins, which is the authors'
    own comparison.
  * **Numbering** - AlphaFold models are numbered by UniProt sequence position, so a table-1
    position maps to the model directly. A position beyond the model length is recorded as
    `position_not_in_model` and is dropped from the positive set rather than silently matched.
  * **Metal coordination** - supplementary table 3 lists cysteines coordinated with metals. Those
    are flagged here and excluded from every pKa statistic later, because PROPKA returns extreme
    values for them. `is_metal_coordinated = 0` means "not listed in table 3", NOT "measured and
    found not to coordinate".
  * **Disulfide** - a cysteine whose SG is within 2.5 A of another cysteine's SG in the model.
    Flagged here and excluded from the statistics later.
  * **Value source** - the pKa is read from PROPKA's SUMMARY block, matched by regular expression
    rather than by column position, because a four-digit residue number is printed without a space
    after the residue name. The determinant table above it
    flags coupled or buried residues with an asterisk, which is not a float; reading that table cost
    pass 1 of this round 122 values and created a spurious cluster at 0.0 (see the audit's
    `corrections`).
  * **Model confidence** - the mean pLDDT over the residue's atoms, read from the model's B-factor
    column. Recorded, not filtered on: filtering by confidence conditions on order, and order is
    correlated with the accessibility this round is about.

Writes results/pka_static_propka.csv and results/pka_static_propka_audit.json.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
MODELS = os.path.join(ROOT, "external", "alphafold_mouse", "pdb")
SCRIPT = os.path.abspath(__file__)
SITES = os.path.join(RESULTS, "doulias2010_sno_sites_alltables.csv")
TABLE1 = "st01_wt_mouse_liver"
TABLE3 = "st03_metal_coordinated"
DISULFIDE_ANGSTROM = 2.5


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def model_atoms(path):
    """{residue number: {'resname': str, 'sg': (x,y,z) or None, 'plddt': mean b-factor}}"""
    out = {}
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        resi = int(line[22:26])
        entry = out.setdefault(resi, {"resname": line[17:20].strip(), "sg": None, "b": []})
        entry["b"].append(float(line[60:66]))
        if line[12:16].strip() == "SG":
            entry["sg"] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    for entry in out.values():
        entry["plddt"] = sum(entry["b"]) / len(entry["b"])
        del entry["b"]
    return out


def disulfide_partners(atoms):
    cys = {r: a["sg"] for r, a in atoms.items() if a["resname"] == "CYS" and a["sg"]}
    partners = {}
    items = sorted(cys.items())
    for i, (ri, si) in enumerate(items):
        for rj, sj in items[i + 1:]:
            d = math.dist(si, sj)
            if d <= DISULFIDE_ANGSTROM:
                partners.setdefault(ri, []).append((rj, round(d, 3)))
                partners.setdefault(rj, []).append((ri, round(d, 3)))
    return partners


def propka_cys(path, workdir):
    """{residue number: pKa} for cysteines, from a PROPKA3 run on a copy of the model."""
    local = os.path.join(workdir, os.path.basename(path))
    with open(path, "rb") as src, open(local, "wb") as dst:
        dst.write(src.read())
    run = subprocess.run([sys.executable, "-m", "propka", local],
                         cwd=workdir, capture_output=True, text=True, timeout=1800)
    pka_file = local[:-4] + ".pka" if local.endswith(".pdb") else local + ".pka"
    if run.returncode != 0 or not os.path.exists(pka_file):
        return None, (run.stderr or run.stdout)[-300:]
    # Read the SUMMARY block, not the determinant table. The determinant line marks coupled or
    # buried residues with an asterisk ("CYS 721 A   9.38*  100 % ..."), which does not parse as a
    # float; pass 1 of this round skipped those residues and then picked up 0.00 from the residue's
    # CONTINUATION line, producing 122 missing values and a spurious cluster at 0.0
    # (results/pka_static_propka_pass1.csv, kept for the record). The SUMMARY block carries one
    # clean line per group: "CYS 721 A     9.38       9.00".
    values = {}
    in_summary = False
    for line in open(pka_file):
        if line.startswith("SUMMARY OF THIS PREDICTION"):
            in_summary = True
            continue
        if in_summary:
            if line.startswith("-") or line.strip().startswith("Free energy"):
                break
            # A four-digit residue number is printed WITHOUT a space after the residue name
            # ("CYS1010 A"), so splitting on whitespace and testing parts[0] == "CYS" silently
            # drops every cysteine past position 999. That cost pass 2 of this round 122 values,
            # all of them in large proteins. Same class of error as this project's recorded
            # \bC\d{1,4}\b trap: never assume the column layout, match it.
            m = re.match(r"\s*CYS\s*(\d+)\s+\S+\s+([-\d.]+)", line)
            if m:
                try:
                    values[int(m.group(1))] = float(m.group(2))
                except ValueError:
                    continue
    return values, None


def main():
    started = time.time()
    with open(SITES, encoding="utf-8-sig", newline="") as fh:
        site_rows = list(csv.DictReader(fh))

    positives, metal = {}, {}
    for r in site_rows:
        if not r["position"].strip():
            continue
        acc, pos = r["accession_normalised"], int(r["position"])
        if r["table"] == TABLE1:
            positives.setdefault(acc, set()).add(pos)
        elif r["table"] == TABLE3:
            metal.setdefault(acc, set()).add(pos)

    rows_out, failures, not_in_model = [], [], []
    workdir = tempfile.mkdtemp(prefix="propka_")
    for acc in sorted(positives):
        model = os.path.join(MODELS, f"AF-{acc}-F1-model_v6.pdb")
        if not os.path.exists(model):
            failures.append({"accession": acc, "reason": "no AlphaFold model"})
            continue
        atoms = model_atoms(model)
        bonds = disulfide_partners(atoms)
        values, error = propka_cys(model, workdir)
        if values is None:
            failures.append({"accession": acc, "reason": f"propka failed: {error}"})
            continue
        for pos in sorted(positives[acc]):
            if pos not in atoms:
                not_in_model.append({"accession": acc, "position": pos,
                                     "model_length": max(atoms) if atoms else 0})
        model_cys = [r for r in sorted(atoms) if atoms[r]["resname"] == "CYS"]
        absent = [r for r in model_cys if r not in values]
        if absent:
            failures.append({"accession": acc, "reason": "cysteines absent from the PROPKA summary",
                             "positions": absent})
        for resi in sorted(atoms):
            if atoms[resi]["resname"] != "CYS":
                continue
            rows_out.append({
                "accession": acc,
                "position": resi,
                "is_positive": int(resi in positives[acc]),
                "is_metal_coordinated": int(resi in metal.get(acc, ())),
                "in_disulfide": int(resi in bonds),
                "disulfide_partner": ";".join(str(p) for p, _ in bonds.get(resi, [])),
                "pka_static": round(values[resi], 3) if resi in values else "",
                "plddt": round(atoms[resi]["plddt"], 2),
                "model_length": max(atoms),
                "model": f"external/alphafold_mouse/pdb/AF-{acc}-F1-model_v6.pdb",
            })

    cols = list(rows_out[0].keys())
    path = os.path.join(RESULTS, "pka_static_propka.csv")
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        for row in rows_out:
            writer.writerow(row)
    os.replace(tmp, path)

    with_value = [r for r in rows_out if r["pka_static"] != ""]
    usable = [r for r in with_value if not r["is_metal_coordinated"] and not r["in_disulfide"]]
    audit = {
        "script": "scripts/run_pka_static_propka.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "registration": "protocols/pka_round_preregistration_2026-09-17.json",
        "instrument": "A: PROPKA3 on the static AlphaFold DB v6 model",
        "no_comparison_no_control_no_verdict_in_this_round": True,
        "disulfide_criterion_angstrom": DISULFIDE_ANGSTROM,
        "n_proteins_attempted": len(positives),
        "n_proteins_failed": len(failures),
        "failures": failures,
        "n_cysteines": len(rows_out),
        "n_cysteines_with_a_value": len(with_value),
        "n_positive_cysteines": sum(r["is_positive"] for r in rows_out),
        "n_metal_coordinated": sum(r["is_metal_coordinated"] for r in rows_out),
        "n_in_disulfide": sum(r["in_disulfide"] for r in rows_out),
        "n_usable_after_exclusions": len(usable),
        "n_usable_positives": sum(r["is_positive"] for r in usable),
        "table1_positions_not_in_the_model": not_in_model,
        "inputs": {"results/doulias2010_sno_sites_alltables.csv": sha256_of(SITES)},
        "versions": {"python": sys.version},
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }
    path = os.path.join(RESULTS, "pka_static_propka_audit.json")
    tmp = path + ".tmp"
    json.dump(audit, open(tmp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    os.replace(tmp, path)

    print(f"proteins {len(positives)} attempted, {len(failures)} failed; "
          f"cysteines {len(rows_out)} ({len(with_value)} with a value)")
    print(f"positives {audit['n_positive_cysteines']}; metal {audit['n_metal_coordinated']}; "
          f"disulfide {audit['n_in_disulfide']}; usable {len(usable)} "
          f"({audit['n_usable_positives']} positive)")
    if not_in_model:
        print("table-1 positions beyond the model:", not_in_model[:5])
    print(f"done in {audit['elapsed_minutes']} min")


if __name__ == "__main__":
    main()
