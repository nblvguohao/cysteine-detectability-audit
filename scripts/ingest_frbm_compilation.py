"""Normalise the compiled human SSH / SNO site lists into one site-level table.

Source: supplementary files mmc1-mmc4 of ScienceDirect article
`1-s2.0-S0891584924006361` (Free Radic Biol Med, 2024), obtained 2026-09-15
through the reader's institutional subscription. The files stay under
`external/intake/frbm2024_compilation/` and are not redistributed in any
deliverable; only derived counts and model scores leave this directory.

Why this table is needed. The QTRP test
(`scripts/run_qtrp_persulfidation_test.py`) builds its negative class from
cysteines that one run detected as free thiol only. A compilation of
persulfidation sites from the wider literature is the natural way to ask how
many of those negatives are known positives elsewhere, which is label noise in
the negative class rather than in the positive class.

Two things come out of the four workbooks:

* `mmc1` holds the two master lists, `SSH-AC_site` and `SNO-AC_site`;
* `mmc2`, `mmc3`, `mmc4` hold the three-way partition with per-site annotation,
  `SSH and SNO`, `SNO only`, `SSH only`, carrying `pPSE value` (a predicted
  solvent-exposure measure), the UniProt feature context and the gene symbol.

The partition is checked against the master lists rather than trusted: `mmc2`
must equal `SSH and SNO`, `mmc4` must equal `SSH minus SNO`, and any
disagreement is recorded in the audit instead of being silently repaired.

The SNO list is kept because it gives a control the persulfidation work has
never had: another cysteine modification, compiled the same way, on the same
proteome. A structural descriptor block that discriminates persulfidation
should not discriminate S-nitrosylation equally well, and if it does, the block
is reading "modifiable cysteine", not persulfidation.
"""
from __future__ import annotations

import collections
import json
import platform
import sys
import time
from pathlib import Path

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json

INTAKE = ROOT / "external/intake/frbm2024_compilation"
PREFIX = "1-s2.0-S0891584924006361"
QTRP_SITES = RESULTS / "qtrp_sites_normalised.csv"

MASTER = {"SSH": "SSH-AC_site", "SNO": "SNO-AC_site"}
PARTITION = {
    "mmc2": ("SSH_and_SNO", "AC_site"),
    "mmc3": ("SNO_only", "AC_site"),
    "mmc4": ("SSH_only", "AC_site"),
}


def sheet_rows(path: Path, sheet: str):
    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        worksheet = workbook[sheet]
        rows = worksheet.iter_rows(values_only=True)
        header = [str(c).strip() if c is not None else "" for c in next(rows)]
        for values in rows:
            yield dict(zip(header, values))
    finally:
        workbook.close()


def site_key(row):
    accession = row.get("AC")
    position = row.get("site")
    if accession in (None, "") or position in (None, ""):
        return None
    try:
        return str(accession).strip(), int(position)
    except (TypeError, ValueError):
        return None


def main():
    started = time.time()
    master = {}
    for label, sheet in MASTER.items():
        keys = set()
        for row in sheet_rows(INTAKE / f"{PREFIX}-mmc1.xlsx", sheet):
            key = site_key(row)
            if key is not None:
                keys.add(key)
        master[label] = keys
        print(f"{sheet}: {len(keys)} sites, {len({a for a, _ in keys})} proteins", flush=True)

    annotation = {}
    partition_keys = {}
    for stem, (label, sheet) in PARTITION.items():
        keys = set()
        for row in sheet_rows(INTAKE / f"{PREFIX}-{stem}.xlsx", sheet):
            key = site_key(row)
            if key is None:
                continue
            keys.add(key)
            ppse = row.get("pPSE value")
            annotation[key] = {
                "gene": (row.get("Gene") or "").strip(),
                "uniprot_annotation": (row.get("Uniprot annotation") or "").strip(),
                "ppse": "" if ppse in (None, "", "NA") else str(ppse).strip(),
                "partition_sheet": label,
            }
        partition_keys[label] = keys
        print(f"{stem} {label}: {len(keys)} sites", flush=True)

    checks = {
        "mmc2_equals_SSH_and_SNO": sorted(
            partition_keys["SSH_and_SNO"] ^ (master["SSH"] & master["SNO"]))[:20],
        "mmc4_equals_SSH_minus_SNO": sorted(
            partition_keys["SSH_only"] ^ (master["SSH"] - master["SNO"]))[:20],
        "mmc3_equals_SNO_minus_SSH": sorted(
            partition_keys["SNO_only"] ^ (master["SNO"] - master["SSH"]))[:20],
    }
    consistency = {name: len(diff) for name, diff in checks.items()}

    rows = []
    for key in sorted(master["SSH"] | master["SNO"]):
        accession, position = key
        extra = annotation.get(key, {})
        rows.append({
            "accession": accession,
            "site": position,
            "compiled_ssh": int(key in master["SSH"]),
            "compiled_sno": int(key in master["SNO"]),
            "gene": extra.get("gene", ""),
            "uniprot_annotation": extra.get("uniprot_annotation", ""),
            "ppse": extra.get("ppse", ""),
            "partition_sheet": extra.get("partition_sheet", ""),
        })
    write_csv(RESULTS / "frbm_compiled_sites.csv", rows)

    # How the compilation lines up with the run-matched QTRP label, which is the
    # reason this ingestion exists.
    overlap = {}
    if QTRP_SITES.exists():
        qtrp = read_csv(QTRP_SITES)
        classes = collections.defaultdict(set)
        per_arm = collections.defaultdict(lambda: collections.defaultdict(set))
        for row in qtrp:
            if (row.get("species") or "") != "human":
                continue
            label = row["class_from_adduct"] or row["class_from_sheet"]
            if label not in ("SH", "SSH"):
                continue
            key = (row["accession"], int(row["site"]))
            classes[label].add(key)
            per_arm[f"{row['source_table']}_pH{row['ph']}"][label].add(key)
        qtrp_ssh = classes["SSH"]
        qtrp_sh_only = classes["SH"] - qtrp_ssh
        overlap["pooled_human"] = {
            "qtrp_ssh": len(qtrp_ssh),
            "qtrp_ssh_in_compilation": len(qtrp_ssh & master["SSH"]),
            "qtrp_sh_only": len(qtrp_sh_only),
            "qtrp_sh_only_flagged_ssh_by_compilation": len(qtrp_sh_only & master["SSH"]),
            "qtrp_sh_only_flagged_sno_only": len(
                qtrp_sh_only & (master["SNO"] - master["SSH"])),
            "compiled_ssh_not_detected_at_all_by_qtrp": len(
                master["SSH"] - qtrp_ssh - qtrp_sh_only),
        }
        for arm, label_sets in sorted(per_arm.items()):
            positives = label_sets["SSH"]
            negatives = label_sets["SH"] - positives
            # Of the negatives flagged persulfidated by the compilation, how many
            # are explained by another QTRP arm rather than by an outside study.
            flagged = negatives & master["SSH"]
            overlap[arm] = {
                "positives": len(positives),
                "negatives": len(negatives),
                "negatives_flagged_ssh_by_compilation": len(flagged),
                "negatives_flagged_share": round(len(flagged) / len(negatives), 4) if negatives else None,
                "flagged_and_ssh_in_another_qtrp_arm": len(flagged & qtrp_ssh),
                "positives_in_compilation": len(positives & master["SSH"]),
            }

    write_json(RESULTS / "frbm_compilation_ingest_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "source": (
            "Free Radic Biol Med 2024 compiled human S-persulfidation and "
            "S-nitrosylation site lists, supplementary mmc1-mmc4 of "
            "ScienceDirect 1-s2.0-S0891584924006361, institutional subscription"
        ),
        "redistribution": "intake files stay under external/intake and are never copied into deliverables",
        "master_lists": {k: len(v) for k, v in master.items()},
        "master_proteins": {k: len({a for a, _ in v}) for k, v in master.items()},
        "ssh_and_sno": len(master["SSH"] & master["SNO"]),
        "partition_sizes": {k: len(v) for k, v in partition_keys.items()},
        "partition_consistency_mismatches": consistency,
        "partition_mismatch_examples": {k: v for k, v in checks.items() if v},
        "qtrp_overlap": overlap,
        "caveat": (
            "the compilation carries no per-site provenance, so a flagged site "
            "cannot be attributed to a specific study; the QTRP cross-tabulation "
            "below separates the part explained by another QTRP arm"
        ),
        "input_hashes": {
            f"external/intake/frbm2024_compilation/{p.name}": sha256(p)
            for p in sorted(INTAKE.glob("*.xlsx"))
        },
        "versions": {"python": sys.version},
    })
    print(json.dumps(overlap, indent=2))


if __name__ == "__main__":
    main()
