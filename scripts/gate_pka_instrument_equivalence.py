"""Equivalence gate: is today's PROPKA the same instrument that produced instrument A?

**Why this gate exists.** `results/pka_static_propka.csv` was computed on 2026-09-17 with PROPKA
invoked as `sys.executable -m propka`. Its audit records the Python version and the script's own
sha256 - but NOT the version of PROPKA itself, and PROPKA is the instrument. When instrument B (the
ensemble) came to be computed, PROPKA was no longer present in any environment on this machine, so
the build that produced instrument A is not recoverable by inspection. Installing "a" PROPKA and
using it for instrument B would confound the comparison this round exists to make - static versus
ensemble - with a silent change of instrument between the two arms.

This script therefore does what the pathology project on the shared cluster does before it lets a
new preprocessing chain touch anything: it recomputes a frozen reference with the new build and
requires element-wise agreement before the new build is allowed downstream.

**The gate, declared before it was run.**

  * Input: the same 192 mouse accessions, the same AlphaFold DB v6 models, and the same code path -
    `propka_cys` and `model_atoms` are IMPORTED from `scripts/run_pka_static_propka.py`, not
    reimplemented, so that any difference measured here is a difference of PROPKA and not of our
    parsing.
  * Reference: `results/pka_static_propka.csv` as stored, read as text. Its `pka_static` column holds
    values rounded to three decimals; the gate compares the newly computed value rounded the same
    way, so that the comparison is between what was STORED and what would be stored today.
  * PASS: every cysteine present in the reference is present now, with an identical rounded value.
    Zero missing, zero extra, zero differing.
  * FAIL: anything else. The gate then reports the number of differing cells, the maximum absolute
    difference, and the ten largest, and instrument B is NOT computed until the discrepancy is
    resolved by amendment. A gate that fails is reported, not loosened - the tolerance is exact
    equality and was fixed before the run precisely so that it cannot be widened afterwards.

**What this gate does not do.** It does not re-issue instrument A's values, does not overwrite any
product, and does not compare the two classes of cysteine. It is an instrument check, not an
analysis: it never reads the `is_positive` column for anything other than carrying it through.

Run with the dedicated PROPKA environment, because that is the interpreter whose `-m propka` the
imported code will invoke:

    /path/to/propka_env/bin/python scripts/gate_pka_instrument_equivalence.py

Writes results/pka_instrument_equivalence_gate.csv and
results/pka_instrument_equivalence_gate_audit.json. Overwrites nothing else.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
SCRIPT = os.path.abspath(__file__)
sys.path.insert(0, os.path.dirname(SCRIPT))

from run_pka_static_propka import MODELS, model_atoms, propka_cys  # noqa: E402

REFERENCE = os.path.join(RESULTS, "pka_static_propka.csv")


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def propka_version():
    try:
        import propka
        return getattr(propka, "__version__", "unknown")
    except Exception as exc:  # pragma: no cover - reported, not raised
        return f"import failed: {exc}"


def main():
    started = time.time()
    with open(REFERENCE, encoding="utf-8-sig", newline="") as fh:
        reference = list(csv.DictReader(fh))

    wanted = {}
    for row in reference:
        wanted.setdefault(row["accession"], {})[int(row["position"])] = row["pka_static"]

    workdir = tempfile.mkdtemp(prefix="pka_gate_")
    rows, failures = [], []
    n_match = n_differ = n_missing = n_extra = 0
    largest = []

    for accession in sorted(wanted):
        model = os.path.join(MODELS, f"AF-{accession}-F1-model_v6.pdb")
        if not os.path.exists(model):
            failures.append({"accession": accession, "reason": "model absent now"})
            continue
        atoms = model_atoms(model)
        values, error = propka_cys(model, workdir)
        if values is None:
            failures.append({"accession": accession, "reason": f"propka failed: {error}"})
            continue
        model_cys = [r for r in sorted(atoms) if atoms[r]["resname"] == "CYS"]
        for resi in model_cys:
            stored = wanted[accession].get(resi)
            now = round(values[resi], 3) if resi in values else None
            now_text = "" if now is None else f"{now:g}"
            if stored is None:
                n_extra += 1
                state = "extra_now"
            elif stored == "" and now is None:
                n_match += 1
                state = "match_both_blank"
            elif stored == "" or now is None:
                n_missing += 1
                state = "missing_on_one_side"
            else:
                same = abs(float(stored) - now) < 1e-9
                n_match += int(same)
                n_differ += int(not same)
                state = "match" if same else "differ"
                if not same:
                    largest.append((abs(float(stored) - now), accession, resi, stored, now_text))
            rows.append({"accession": accession, "position": resi,
                         "stored_pka_static": stored if stored is not None else "",
                         "recomputed_pka_static": now_text, "state": state})
        for resi in wanted[accession]:
            if resi not in atoms:
                n_missing += 1
                rows.append({"accession": accession, "position": resi,
                             "stored_pka_static": wanted[accession][resi],
                             "recomputed_pka_static": "", "state": "absent_from_model_now"})

    largest.sort(reverse=True)
    passed = (n_differ == 0 and n_missing == 0 and n_extra == 0 and not failures)

    path = os.path.join(RESULTS, "pka_instrument_equivalence_gate.csv")
    with open(path + ".tmp", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["accession", "position", "stored_pka_static",
                                                "recomputed_pka_static", "state"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(path + ".tmp", path)

    audit = {
        "script": "scripts/gate_pka_instrument_equivalence.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "what": "element-wise equivalence gate between the stored instrument A values and the same "
                "computation repeated with the PROPKA build available today",
        "why": "instrument A's audit recorded the Python version but not the PROPKA version, and "
               "PROPKA was absent from every environment on this machine when instrument B came to "
               "be computed. Without this gate, 'static versus ensemble' would be confounded with "
               "'one PROPKA build versus another'.",
        "tolerance": "exact equality of the stored three-decimal value; fixed before the run",
        "interpreter": {"python": sys.version, "propka": propka_version(),
                        "environment": sys.prefix},
        "reference": {"path": "results/pka_static_propka.csv", "sha256": sha256_of(REFERENCE),
                      "rows": len(reference)},
        "counts": {"compared": len(rows), "match": n_match, "differ": n_differ,
                   "missing_on_one_side": n_missing, "extra_now": n_extra},
        "per_protein_failures": failures,
        "ten_largest_differences": [
            {"accession": a, "position": p, "stored": s, "recomputed": n, "abs_diff": round(d, 6)}
            for d, a, p, s, n in largest[:10]],
        "gate": "PASS" if passed else "FAIL",
        "consequence_if_fail": "instrument B is not computed until an amendment resolves the "
                               "discrepancy; the tolerance is not widened",
        "what_this_gate_did_not_do": [
            "it did not overwrite results/pka_static_propka.csv or its audit",
            "it did not compare the positive and background classes",
            "it did not issue any verdict",
        ],
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "corrections": [],
    }
    path = os.path.join(RESULTS, "pka_instrument_equivalence_gate_audit.json")
    with open(path + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)
    os.replace(path + ".tmp", path)

    print(f"propka {audit['interpreter']['propka']}  python {sys.version.split()[0]}")
    print(f"compared {len(rows)}: match {n_match}, differ {n_differ}, "
          f"missing {n_missing}, extra {n_extra}, protein failures {len(failures)}")
    if largest:
        print("largest differences:", largest[:5])
    print("GATE:", audit["gate"], f"({audit['elapsed_minutes']} min)")


if __name__ == "__main__":
    main()
