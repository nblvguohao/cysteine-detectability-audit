"""Instrument B of the pKa round: PROPKA3 averaged over an explicit-solvent MD ensemble.

Declared in `protocols/pka_round_preregistration_2026-09-17.json` and its amendments A1 and A2,
all written before any ensemble pKa existed. This script MEASURES instrument B. It runs no
comparison, applies no control, computes no class-wise statistic and issues no verdict; the round's
verdict is issued by a separate script whose sha256 is recorded before it is pointed at this
script's output, per the blind-prediction lock.

**What instrument B is, and what it is not.** It is the SAME PROPKA3 calculation as instrument A,
run on snapshots of a 1 ns explicit-solvent trajectory of the same protein and averaged. Its role is
to ask whether instrument A's reading is an artefact of using one static conformer. A 1 ns
equilibrium ensemble is NOT converged sampling and is NOT a free-energy calculation; the
pre-registration declares this and it travels with every number below.

Rules fixed before the run (A2):

  * **Instrument** - PROPKA 3.5.1 in `/Users/lyuguohao/opt/propka_env`, the build that reproduced
    all 1,761 stored instrument A values exactly (`results/pka_instrument_equivalence_gate.csv`,
    gate PASS). PROPKA is not installed on the compute host, so it is run here for both instruments.
  * **Completeness** - a protein enters instrument B only if its `md.log` records
    `Statistics over 500001 steps` and its `n_snapshots.txt` reads 11. Anything else is named in the
    audit and excluded; a shorter ensemble is never silently included. This criterion exists because
    the queue's own completeness test was "md.xtc is non-empty", which an interrupted run also
    satisfies - two proteins were caught that way and re-run.
  * **Hydrogens** - PROPKA 3.5.1 refuses the all-atom snapshot outright ("Unexpected number (9) of
    atoms in residue ASP 293"). Every snapshot is therefore reduced to heavy atoms by the PDB
    element column before PROPKA is called, and nothing else about it is altered. This also makes
    instrument B's input format identical to instrument A's, whose AlphaFold models carry no
    hydrogens, so that what differs between the instruments is the conformation.
  * **Aggregation** - the ensemble value is the ARITHMETIC MEAN of the per-snapshot pKa over the
    snapshots in which PROPKA returned a value for that cysteine, as the registration states
    ("pKa computed on each snapshot and averaged"). The standard deviation and the number of
    contributing snapshots are stored beside it. A cysteine missing from a snapshot's SUMMARY block
    is recorded as missing for that snapshot and is never imputed.
  * **Value source** - the SUMMARY block, matched by the same regular expression as instrument A;
    `propka_cys` is IMPORTED from `scripts/run_pka_static_propka.py` rather than reimplemented, so
    the two instruments cannot drift apart in their parsing.
  * **Raw snapshots stay on the compute host.** They are ~0.5 GB and are not products of this tree.
    What is stored here is the per-snapshot value together with the sha256 of the heavy-atom PDB
    that produced it, so any number can be traced back to a specific file.

Writes results/pka_ensemble_per_snapshot.csv (one row per cysteine per snapshot),
results/pka_ensemble_propka.csv (one row per cysteine) and results/pka_ensemble_propka_audit.json.
Overwrites nothing else.

Run with the pinned environment:
    /Users/lyuguohao/opt/propka_env/bin/python scripts/run_pka_ensemble_propka.py
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
SCRIPT = os.path.abspath(__file__)
sys.path.insert(0, os.path.dirname(SCRIPT))

from run_pka_static_propka import propka_cys  # noqa: E402

HOST = "a100-jump"
REMOTE = "/data/lgh/prod/jobs/persulf_pka"
SUBSET = os.path.join(RESULTS, "pka_md_subset_accessions.txt")
REQUIRED_STEPS = "500001"
REQUIRED_SNAPSHOTS = 11


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def remote(cmd, binary=False):
    run = subprocess.run(["ssh", "-o", "ConnectTimeout=30", HOST, cmd],
                         capture_output=True, timeout=1800)
    if run.returncode != 0:
        raise RuntimeError(run.stderr.decode("utf-8", "replace")[-300:])
    return run.stdout if binary else run.stdout.decode("utf-8", "replace")


def strip_hydrogens(src, dst):
    """Keep every non-ATOM line and every ATOM whose element column is not H."""
    kept = dropped = 0
    with open(src, encoding="utf-8", errors="replace") as fh, \
            open(dst, "w", encoding="utf-8") as out:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                element = line[76:78].strip().upper()
                if element == "H":
                    dropped += 1
                    continue
                kept += 1
            out.write(line)
    return kept, dropped


def main():
    started = time.time()
    subset = [a.strip() for a in open(SUBSET) if a.strip()]

    # completeness, measured on the compute host from md.log and n_snapshots.txt
    listing = remote(
        'B=%s; for d in $B/ensembles/*/; do a=$(basename $d); '
        's=$(grep -a "Statistics over" "$d/md.log" 2>/dev/null | tail -1 | awk "{print \\$3}"); '
        'n=$(cat "$d/n_snapshots.txt" 2>/dev/null); echo "$a|$s|$n"; done' % REMOTE)
    state = {}
    for line in listing.splitlines():
        if line.count("|") == 2:
            acc, steps, snaps = line.split("|")
            state[acc] = (steps.strip(), snaps.strip())

    complete, excluded = [], []
    for acc in subset:
        steps, snaps = state.get(acc, ("", ""))
        if steps == REQUIRED_STEPS and snaps == str(REQUIRED_SNAPSHOTS):
            complete.append(acc)
        else:
            excluded.append({"accession": acc, "statistics_over_steps": steps or "(none)",
                             "n_snapshots": snaps or "(none)",
                             "reason": "did not meet the declared completeness criterion"})

    per_snapshot, failures = [], []
    workroot = tempfile.mkdtemp(prefix="pka_ens_")
    for index, acc in enumerate(complete, 1):
        workdir = os.path.join(workroot, acc)
        os.makedirs(workdir, exist_ok=True)
        try:
            blob = remote(f"cd {REMOTE}/ensembles/{acc} && tar cf - snap*.pdb | gzip -1",
                          binary=True)
            tgz = os.path.join(workdir, "s.tgz")
            open(tgz, "wb").write(blob)
            subprocess.run(["tar", "xzf", tgz, "-C", workdir], check=True, capture_output=True)
        except Exception as exc:
            failures.append({"accession": acc, "stage": "fetch", "error": str(exc)[-200:]})
            shutil.rmtree(workdir, ignore_errors=True)
            continue
        for snap in range(REQUIRED_SNAPSHOTS):
            src = os.path.join(workdir, f"snap{snap}.pdb")
            if not os.path.exists(src):
                failures.append({"accession": acc, "stage": "snapshot", "snapshot": snap,
                                 "error": "file absent after extraction"})
                continue
            dst = os.path.join(workdir, f"heavy{snap}.pdb")
            kept, dropped = strip_hydrogens(src, dst)
            digest = sha256_of(dst)
            # propka_cys copies its input into the workdir it is given. If that is the directory the
            # input already sits in, the copy opens the same path for writing, truncates it, and
            # PROPKA then reads an empty file ("does not seem to contain any molecular
            # conformations"). Instrument A never hit this because its inputs lived outside the
            # temporary workdir. Give PROPKA a directory of its own.
            rundir = os.path.join(workdir, f"propka{snap}")
            os.makedirs(rundir, exist_ok=True)
            values, error = propka_cys(dst, rundir)
            if values is None:
                failures.append({"accession": acc, "stage": "propka", "snapshot": snap,
                                 "error": (error or "")[-200:]})
                continue
            for position, value in sorted(values.items()):
                per_snapshot.append({"accession": acc, "position": position, "snapshot": snap,
                                     "pka": round(value, 3), "heavy_atoms": kept,
                                     "hydrogens_dropped": dropped,
                                     "heavy_pdb_sha256": digest})
        shutil.rmtree(workdir, ignore_errors=True)
        if index % 20 == 0:
            print(f"  {index}/{len(complete)} proteins, {len(per_snapshot)} values, "
                  f"{round((time.time() - started) / 60, 1)} min")
    shutil.rmtree(workroot, ignore_errors=True)

    grouped = {}
    for row in per_snapshot:
        grouped.setdefault((row["accession"], row["position"]), []).append(row["pka"])

    aggregate = []
    for (acc, position), values in sorted(grouped.items()):
        aggregate.append({
            "accession": acc, "position": position,
            "pka_ensemble_mean": round(statistics.fmean(values), 3),
            "pka_ensemble_sd": round(statistics.pstdev(values), 3) if len(values) > 1 else "",
            "pka_ensemble_min": round(min(values), 3),
            "pka_ensemble_max": round(max(values), 3),
            "n_snapshots_with_a_value": len(values),
            "n_snapshots_expected": REQUIRED_SNAPSHOTS,
        })

    if not aggregate:
        raise SystemExit(
            "no ensemble value was produced. complete proteins: %d, failures: %d. First failures: %s"
            % (len(complete), len(failures), failures[:3]))

    for path, data in (
            (os.path.join(RESULTS, "pka_ensemble_per_snapshot.csv"), per_snapshot),
            (os.path.join(RESULTS, "pka_ensemble_propka.csv"), aggregate)):
        cols = list(data[0].keys())
        with open(path + ".tmp", "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=cols)
            writer.writeheader()
            for row in data:
                writer.writerow(row)
        os.replace(path + ".tmp", path)

    partial = [r for r in aggregate if r["n_snapshots_with_a_value"] != REQUIRED_SNAPSHOTS]
    audit = {
        "script": "scripts/run_pka_ensemble_propka.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "registration": "protocols/pka_round_preregistration_2026-09-17.json, amendments A1 and A2",
        "instrument": "B: PROPKA3 on heavy atoms of 11 snapshots of a 1 ns explicit-solvent "
                      "ensemble, arithmetic mean over snapshots",
        "no_comparison_no_control_no_verdict_in_this_round": True,
        "interpreter": {"python": sys.version, "environment": sys.prefix,
                        "propka": __import__("propka").__version__
                        if hasattr(__import__("propka"), "__version__") else "unknown"},
        "instrument_equivalence_gate": "results/pka_instrument_equivalence_gate_audit.json (PASS, "
                                       "1,761 of 1,761 stored instrument A values reproduced "
                                       "exactly with this build)",
        "completeness_criterion": f"md.log records 'Statistics over {REQUIRED_STEPS} steps' and "
                                  f"n_snapshots.txt reads {REQUIRED_SNAPSHOTS}",
        "n_proteins_in_the_registered_draw": len(subset),
        "n_proteins_complete_and_measured": len(complete),
        "n_proteins_excluded_as_incomplete": len(excluded),
        "excluded": excluded,
        "n_cysteines_with_an_ensemble_value": len(aggregate),
        "n_per_snapshot_values": len(per_snapshot),
        "cysteines_missing_from_at_least_one_snapshot": len(partial),
        "cysteines_missing_from_at_least_one_snapshot_detail": partial[:50],
        "failures": failures,
        "declared_limits": [
            "a 1 ns equilibrium ensemble is not converged sampling and is not a free-energy "
            "calculation; it is used only to ask whether instrument A's ranking moves",
            "PROPKA is an empirical calculation and is acknowledged in the registration to be "
            "unreliable for cysteine in absolute terms; both instruments share that weakness",
            "raw snapshots remain on the compute host; each value carries the sha256 of the "
            "heavy-atom PDB that produced it",
        ],
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "corrections": [],
    }
    path = os.path.join(RESULTS, "pka_ensemble_propka_audit.json")
    with open(path + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)
    os.replace(path + ".tmp", path)

    print(f"proteins: {len(complete)} measured of {len(subset)} drawn "
          f"({len(excluded)} excluded as incomplete)")
    print(f"cysteines with an ensemble value: {len(aggregate)}; "
          f"per-snapshot values: {len(per_snapshot)}; failures: {len(failures)}")
    print(f"cysteines missing from at least one snapshot: {len(partial)}")
    print(f"done in {audit['elapsed_minutes']} min")


if __name__ == "__main__":
    main()
