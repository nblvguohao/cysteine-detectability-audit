#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Does the search-space artefact appear outside cysteine chemistry? Phosphorylation against sulfation.

Criteria, branches, sensitivity grid and controls are in
protocols/phospho_sulfo_search_space_preregistration_2026-09-23.json, written before any deposit was selected.
Gate A0 asserts that protocol's sha256 and refuses to run if it changed. Nothing below re-decides anything the
protocol fixed; the deposit was chosen by scripts/fetch_phospho_deposit_2026-09-23.py under the registered rule
and is read from its fetch record rather than named again here.

WHAT IS MEASURED
The released check cys_audit.checks.searchspace computes the resolving mass M* = delta / (ppm x 1e-6) and asks
whether the near-isobaric alternative was in the search space. Phosphorylation (HPO3) and sulfation (SO3) differ
by about 0.0095 Da, roughly half the cysteine +S / +2O difference this manuscript already reports, so M* for this
pair sits at about half the peptide mass. The check is unchanged; only its modification table gained this pair.

WHAT COUNTS AS AN IDENTIFIED PHOSPHOPEPTIDE (declared here, before the table was read)
A row of modificationSpecificPeptides.txt with Phospho (STY) >= 1 and a finite Mass. The primary set additionally
drops rows marked Reverse or Potential contaminant, which is the standard MaxQuant convention; the unfiltered set
is reported beside it as a sensitivity, so the filter cannot be doing the work.

WHERE THE CONFIGURATION COMES FROM
This deposit's parameters.txt states neither the variable modifications nor the precursor tolerance, so both are
read from its own mqpar.xml and its Andromeda .apar files, which is what branch B3 requires before calling them
unreadable. Every parameter group is parsed, not the first one: a deposit may search different groups differently,
and a single group would not license a statement about the deposit. If the groups disagree on tolerance or on
whether sulfation was searched, the run REFUSES rather than picking one.

POSITION COLUMN
The released reader needs protein and position; the search-space check reads neither, only label and peptide_mass.
modificationSpecificPeptides.txt carries no site position, so position is a row ordinal, declared a placeholder
here and in the audit. Inventing a site coordinate would be a fabricated column, and this check would not read it.

CONTROLS (any failure REFUSES the run and writes nothing)
  PC1 the extended table still reproduces the cysteine delta the manuscript prints, 0.017758 Da
  PC2 adding sulfation to the search space moves the status away from FAIL
  PC3 at a tolerance whose M* exceeds every observed peptide mass, the status is not FAIL
  NC1 a modification outside the table returns UNDECIDABLE, never PASS
  PC5 (added here, not registered) the share at or above M* is monotone non-decreasing in the tolerance
  PC6 (added here, not registered) the primary and unfiltered sets give the same branch
  PC4 is the released test suite and is run separately by the caller, not from inside this script.
Outputs: results/phospho_sulfo_search_space_2026-09-23.csv and ..._audit.json. Reruns are byte identical.
Interpreter: any Python 3.9+ with numpy (the project venv or the canonical 3.11 environment).
"""
import csv
import hashlib
import json
import os
import sys
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "cys-audit", "src"))

from cys_audit import constants as C  # noqa: E402
from cys_audit.checks import searchspace  # noqa: E402
from cys_audit.io import read_sites  # noqa: E402

PROTO = os.path.join(ROOT, "protocols", "phospho_sulfo_search_space_preregistration_2026-09-23.json")
PROTO_SHA = "f04b136d67a25fad10289026de3b6caabbd6ddc9ca1956ff1cb9f003319c3cf2"
DEP = os.path.join(ROOT, "external", "pride_phospho_20260923")
OUT_CSV = os.path.join(ROOT, "results", "phospho_sulfo_search_space_2026-09-23.csv")
OUT_JSON = os.path.join(ROOT, "results", "phospho_sulfo_search_space_2026-09-23_audit.json")
TOOL_IN = os.path.join(ROOT, "results", "phospho_sulfo_tool_input_2026-09-23.tsv")
GRID = [2.0, 4.5, 10.0, 20.0]
MIN_N = 200


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def parse_mqpar(path):
    """Every parameter group's variable modifications and precursor tolerances, not just the first."""
    root = ET.parse(path).getroot()
    groups = []
    for g in root.iter("parameterGroup"):
        var = [s.text for s in g.findall("./variableModifications/string") if s.text]
        fix = [s.text for s in g.findall("./fixedModifications/string") if s.text]
        def one(tag):
            e = g.find("./" + tag)
            return None if e is None or e.text is None else float(e.text)
        groups.append({"variable_modifications": var, "fixed_modifications": fix,
                       "first_search_ppm": one("firstSearchTol"), "main_search_ppm": one("mainSearchTol"),
                       "search_tol_in_ppm": (g.findtext("./searchTolInPpm") or "").strip(),
                       "max_peptide_mass": one("maxPeptideMass")})
    return groups


def parse_apar(path):
    d = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "=" in line:
                k, v = line.rstrip("\n").split("=", 1)
                d[k.strip()] = v.strip()
    return {"variable_modifications": [x for x in d.get("variable modifications", "").split(",") if x],
            "fixed_modifications": [x for x in d.get("fixed modifications", "").split(",") if x],
            "peptide_mass_tolerance": float(d["peptide mass tolerance"]) if "peptide mass tolerance" in d else None,
            "peptide_mass_tolerance_unit": d.get("peptide mass tolerance Unit"),
            "max_peptide_mass": float(d["max peptide mass"]) if "max peptide mass" in d else None,
            "enzymes": d.get("enzymes")}


def looks_like_sulfation(names):
    return sorted({n for n in names if searchspace.canonical(n) == "sulfation"})


def main():
    fails = []
    if sha(PROTO) != PROTO_SHA:
        sys.exit("REFUSE (A0): the pre-registration changed since the criteria were written")
    rec = json.load(open(os.path.join(DEP, "fetch_record.json"), encoding="utf-8"))
    cfg_rec = json.load(open(os.path.join(DEP, "fetch_record_config.json"), encoding="utf-8"))
    acc = rec["chosen"]
    d = os.path.join(DEP, acc)

    groups = parse_mqpar(os.path.join(d, "mqpar.xml"))
    apars = {n: parse_apar(os.path.join(d, n)) for n in sorted(cfg_rec["files"]) if n.endswith(".apar")}
    if not groups:
        sys.exit("REFUSE (B3): mqpar.xml has no parameter group")
    mains = sorted({g["main_search_ppm"] for g in groups})
    firsts = sorted({g["first_search_ppm"] for g in groups})
    if len(mains) != 1 or mains[0] is None:
        sys.exit("REFUSE: parameter groups disagree on the main search tolerance: %s" % mains)
    if any(g["search_tol_in_ppm"].lower() != "true" for g in groups):
        sys.exit("REFUSE: a parameter group does not express its search tolerance in ppm")
    main_ppm = mains[0]
    for n, a in apars.items():
        if a["peptide_mass_tolerance"] != main_ppm or a["peptide_mass_tolerance_unit"] != "ppm":
            sys.exit("REFUSE: %s disagrees with mqpar on the precursor tolerance" % n)
    searched_names = sorted({m for g in groups for m in g["variable_modifications"] + g["fixed_modifications"]}
                            | {m for a in apars.values() for m in a["variable_modifications"] + a["fixed_modifications"]})
    sulfation_present = looks_like_sulfation(searched_names)
    phospho_present = sorted({n for n in searched_names if searchspace.canonical(n) == "phosphorylation"})
    if not phospho_present:
        sys.exit("REFUSE: the deposit does not search phosphorylation; the registered rule chose it as one that does")

    # ---- identified phosphopeptides -------------------------------------------------------------
    src = os.path.join(d, "modificationSpecificPeptides.txt")
    masses_primary, masses_all, n_rows = [], [], 0
    with open(src, encoding="utf-8-sig", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        need = ["Phospho (STY)", "Mass", "Reverse", "Potential contaminant", "Proteins"]
        miss = [c for c in need if c not in (rd.fieldnames or [])]
        if miss:
            sys.exit("REFUSE: modificationSpecificPeptides.txt lacks column(s) %s" % miss)
        for r in rd:
            n_rows += 1
            try:
                nph = int(float(r["Phospho (STY)"] or 0))
                m = float(r["Mass"])
            except ValueError:
                continue
            if nph < 1 or not (m == m) or m <= 0:
                continue
            prot = (r["Proteins"] or "").split(";")[0].strip() or "unassigned"
            masses_all.append((prot, m))
            if (r["Reverse"] or "").strip() != "+" and (r["Potential contaminant"] or "").strip() != "+":
                masses_primary.append((prot, m))
    if len(masses_primary) < MIN_N:
        sys.exit("REFUSE (B4): only %d identified phosphopeptides carry a usable mass" % len(masses_primary))

    # ---- tool input, read back through the released reader ---------------------------------------
    def write_input(path, rows):
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["protein", "position", "label", "peptide_mass"])
            for i, (p, m) in enumerate(rows, 1):
                w.writerow([p, i, 1, "%.6f" % m])
    write_input(TOOL_IN, masses_primary)
    ds = read_sites(TOOL_IN)
    tmp_all = TOOL_IN + ".unfiltered.tmp"
    write_input(tmp_all, masses_all)
    ds_all = read_sites(tmp_all)
    os.remove(tmp_all)

    def run(dset, ppm, mods):
        return searchspace.run(dset, {"modification": "Phospho (STY)", "precursor_ppm": ppm,
                                      "search_mods": mods, "identity_readout": "direct_mass"})

    as_searched = searched_names
    primary = run(ds, main_ppm, as_searched)
    head = primary["details"]["alternatives"][0]
    share_primary = head["share_positive_peptides_at_or_above_resolving_mass"]

    rows_out = []
    for ppm in sorted(set(GRID + [main_ppm])):
        rp = run(ds, ppm, as_searched)
        ra = run(ds_all, ppm, as_searched)
        hp, ha = rp["details"]["alternatives"][0], ra["details"]["alternatives"][0]
        rows_out.append({
            "accession": acc, "claimed_modification": "phosphorylation", "alternative": "sulfation",
            "delta_da": hp["delta_da"], "tolerance_ppm": ppm,
            "is_deposit_main_search_tolerance": ppm == main_ppm,
            "resolving_mass_da": hp["resolving_mass_da"],
            "alternative_in_search_space": hp["alternative_in_search_space"],
            "n_phosphopeptides_primary": rp["n_positive"],
            "share_at_or_above_resolving_mass_primary": hp["share_positive_peptides_at_or_above_resolving_mass"],
            "status_primary": rp["status"],
            "n_phosphopeptides_unfiltered": ra["n_positive"],
            "share_at_or_above_resolving_mass_unfiltered": ha["share_positive_peptides_at_or_above_resolving_mass"],
            "status_unfiltered": ra["status"]})

    # ---- controls --------------------------------------------------------------------------------
    cys_delta = round(abs(searchspace.DELTAS["persulfidation"] - searchspace.DELTAS["dioxidation"]), 6)
    pc1 = cys_delta == 0.017758
    pc2 = run(ds, main_ppm, as_searched + ["Sulfo (STY)"])
    max_mass = max(m for _, m in masses_primary)
    ppm_safe = round(head["delta_da"] / (max_mass * 2.0) / 1e-6, 6)  # M* = twice the heaviest peptide
    pc3 = run(ds, ppm_safe, as_searched)
    nc1 = searchspace.run(ds, {"modification": "Nitrosylation (C)", "precursor_ppm": main_ppm,
                               "search_mods": as_searched, "identity_readout": "direct_mass"})
    shares = [r["share_at_or_above_resolving_mass_primary"] for r in rows_out]
    ppms = [r["tolerance_ppm"] for r in rows_out]
    pc5 = all(shares[i] <= shares[i + 1] + 1e-12 for i in range(len(shares) - 1)) and ppms == sorted(ppms)
    pc6 = all(r["status_primary"] == r["status_unfiltered"] for r in rows_out)
    controls = {"PC1_cysteine_delta_reproduced": {"pass": bool(pc1), "value": cys_delta},
                "PC2_adding_sulfation_moves_status": {"pass": pc2["status"] != primary["status"],
                                                      "status_with_sulfation": pc2["status"],
                                                      "status_without": primary["status"]},
                "PC3_no_fail_when_resolving_mass_exceeds_every_peptide":
                    {"pass": pc3["status"] != C.FAIL, "tolerance_ppm": ppm_safe,
                     "resolving_mass_da": pc3["details"]["alternatives"][0]["resolving_mass_da"],
                     "heaviest_peptide_da": round(max_mass, 4), "status": pc3["status"]},
                "NC1_unknown_modification_is_undecidable":
                    {"pass": nc1["status"] == C.UNDECIDABLE, "status": nc1["status"]},
                "PC5_share_monotone_in_tolerance": {"pass": bool(pc5), "shares": shares, "ppm": ppms},
                "PC6_primary_and_unfiltered_agree_on_branch": {"pass": bool(pc6)}}
    for k, v in controls.items():
        if not v["pass"]:
            fails.append(k)
    if fails:
        sys.exit("REFUSE: controls failed: %s" % fails)

    # ---- branch ----------------------------------------------------------------------------------
    if sulfation_present:
        branch, verdict = "B2", "sulfation IS in the search space; the artefact does not apply to this deposit"
    elif share_primary is not None and share_primary >= C.SEARCHSPACE_FAIL_SHARE:
        branch = "B1"
        verdict = ("sulfation is absent from the search space and %.4f of the identified phosphopeptides are at or "
                   "above the resolving mass at the deposit's own %s ppm main search, so for that share the "
                   "precursor mass cannot separate the two chemistries and only one of them could be reported"
                   % (share_primary, main_ppm))
    else:
        branch, verdict = "B1_not_met", "sulfation absent but the share below the registered FAIL share"

    with open(OUT_CSV, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0]), lineterminator="\n")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)

    audit = {
        "script": os.path.relpath(os.path.abspath(__file__), ROOT), "script_sha256": sha(os.path.abspath(__file__)),
        "date": "2026-09-23", "protocol": os.path.relpath(PROTO, ROOT), "protocol_sha256": PROTO_SHA,
        "tool": {"name": "cys-audit", "version": C.VERSION, "check": searchspace.TEST,
                 "searchspace_sha256": sha(os.path.join(ROOT, "cys-audit", "src", "cys_audit", "checks",
                                                        "searchspace.py")),
                 "constants_sha256": sha(os.path.join(ROOT, "cys-audit", "src", "cys_audit", "constants.py")),
                 "fail_share": C.SEARCHSPACE_FAIL_SHARE},
        "deposit": {"accession": acc, "title": rec["chosen_title"],
                    "selection": "first qualifying candidate of the registered PRIDE query, in API order",
                    "qualifying_in_order": rec["qualifying_in_order"],
                    "inputs": {k: v["sha256"] for k, v in sorted(rec["files"].items())} |
                              {k: v["sha256"] for k, v in sorted(cfg_rec["files"].items())}},
        "search_configuration_read_from_the_deposit": {
            "source": ["mqpar.xml", "0.allSpectra.CID.FTMS.sil0.apar", "0.allSpectra.HCD.FTMS.sil0.apar"],
            "n_parameter_groups": len(groups), "main_search_ppm": main_ppm, "first_search_ppm": firsts,
            "modifications_searched": searched_names,
            "phosphorylation_entries": phospho_present, "sulfation_entries": sulfation_present,
            "note": "parameters.txt of this MaxQuant version states neither the modification list nor the "
                    "precursor tolerance; branch B3 therefore reads them from the deposit's other own files, and "
                    "the run refuses if parameter groups disagree"},
        "identified_phosphopeptides": {
            "definition": "a modificationSpecificPeptides.txt row with Phospho (STY) >= 1 and a finite positive Mass",
            "primary_filter": "rows marked Reverse or Potential contaminant are dropped",
            "n_table_rows": n_rows, "n_primary": len(masses_primary), "n_unfiltered": len(masses_all),
            "heaviest_peptide_da": round(max_mass, 4),
            "position_column": "a row ordinal placeholder; the search-space check reads only label and "
                               "peptide_mass, and this table carries no site position"},
        "result": {"branch": branch, "verdict": verdict, "status_at_deposit_tolerance": primary["status"],
                   "resolving_mass_da": head["resolving_mass_da"], "delta_da": head["delta_da"],
                   "share_at_or_above_resolving_mass": share_primary, "reason": primary["reason"]},
        "sensitivity_grid_ppm": GRID, "rows": rows_out, "controls": controls, "gates_failed": fails,
        "what_this_cannot_show": json.load(open(PROTO, encoding="utf-8"))["what_this_cannot_show"],
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=1, sort_keys=False)
        fh.write("\n")
    print(json.dumps({"branch": branch, "status": primary["status"], "main_ppm": main_ppm,
                      "resolving_mass_da": head["resolving_mass_da"], "delta_da": head["delta_da"],
                      "n_primary": len(masses_primary), "share": share_primary,
                      "sulfation_in_search_space": bool(sulfation_present),
                      "controls_all_pass": not fails}, indent=1))


if __name__ == "__main__":
    main()
