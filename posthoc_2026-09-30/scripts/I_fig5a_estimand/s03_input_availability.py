"""Step 3 - can a chemistry-aware denominator be fitted for ONE Figure 5a cohort from local files?

POST HOC revision analysis (2026-09-30), item I_fig5a_estimand. Not registered, not pre-specified.

The decision rule given to this item: run a minimal chemistry+detectability demonstration for
one cohort only if its per-site inputs AND protein sequences are available locally (no network).
This script only checks availability and writes the evidence; it fits nothing.

What such a demonstration needs, per cohort:
  (a) the per-site role table (accession, site, positive / observed_unmodified) - the census
      input results/ptm_census_sites.csv, or the normalised source tables it was built from;
  (b) the full reference sequence of every protein carrying a site (for the 25 digest features,
      for enumerating NEG_A, and for any sequence-window chemistry feature);
  (c) optionally structural features for the same cysteines.

Run:  python -B s03_input_availability.py
"""
from __future__ import annotations

import csv
import json
import os
import pathlib

import common as C

DESKTOP_ROOT = pathlib.Path(r"/path/to/local")


def status(path: pathlib.Path, expected_sha: str | None = None) -> dict:
    d = {"path": str(path), "exists": path.exists()}
    if path.exists() and path.is_file():
        d["bytes"] = path.stat().st_size
        d["sha256"] = C.sha256(path)
        if expected_sha:
            d["expected_sha256"] = expected_sha
            d["sha256_matches_expected"] = d["sha256"] == expected_sha
    elif expected_sha:
        d["expected_sha256"] = expected_sha
    return d


def find_sequence_files(root: pathlib.Path) -> list:
    hits = []
    pats = (".fasta", ".fasta.gz", ".fa", ".fa.gz", ".faa", ".faa.gz")
    for dirpath, dirnames, filenames in os.walk(root):
        # skip virtual environments and git internals
        dirnames[:] = [d for d in dirnames if d not in (".git", ".venv_tf", "node_modules", "__pycache__")]
        for fn in filenames:
            low = fn.lower()
            if low.endswith(pats) or ("proteome" in low and (low.endswith(".tsv") or low.endswith(".tsv.gz"))):
                p = pathlib.Path(dirpath) / fn
                hits.append({"path": str(p), "bytes": p.stat().st_size})
    return sorted(hits, key=lambda h: h["path"])


def main():
    census_audit = json.loads(C.CENSUS_AUDIT.read_text(encoding="utf-8"))
    sites_audit = json.loads(C.CENSUS_SITES_AUDIT.read_text(encoding="utf-8"))
    inp = sites_audit["inputs"]
    need = {
        "census_site_table results/ptm_census_sites.csv": status(C.REPO / "results/ptm_census_sites.csv",
                                                                  census_audit["sites_table_sha256"]),
        "human proteome external/proteomes/hsa.fasta.gz": status(C.REPO / "external/proteomes/hsa.fasta.gz",
                                                                 inp["hsa_proteome"]["sha256"]),
        "Arabidopsis proteome external/proteomes/ath.fasta.gz": status(C.REPO / "external/proteomes/ath.fasta.gz",
                                                                       inp["ath_proteome"]["sha256"]),
        "Arabidopsis UniProt-TAIR map external/proteomes/ath_uniprot_tair_map.tsv.gz":
            status(C.REPO / "external/proteomes/ath_uniprot_tair_map.tsv.gz"),
        "external supplementary tables external/intake/ptm_census_supp/": {
            "path": str(C.REPO / "external/intake/ptm_census_supp"),
            "exists": (C.REPO / "external/intake/ptm_census_supp").exists()},
        "QTRP normalised sites results/qtrp_sites_normalised.csv": status(C.REPO / "results/qtrp_sites_normalised.csv",
                                                                         inp["qtrp"]["sha256"]),
        "qPerS-SID normalised sites results/qpers_sid_sites_normalised.csv":
            status(C.REPO / "results/qpers_sid_sites_normalised.csv", inp["qpers"]["sha256"]),
        "human homology components results/v3_human_homology_components.csv":
            status(C.REPO / "results/v3_human_homology_components.csv"),
        "structural features results/structural_site_features.csv": status(C.REPO / "results/structural_site_features.csv"),
    }
    # does the structural table cover any Fig. 5a cohort?
    cohorts_in_structural = {}
    with open(C.REPO / "results/structural_site_features.csv", newline="", encoding="utf-8-sig") as h:
        for r in csv.DictReader(h):
            for t in r["cohorts"].split(";"):
                cohorts_in_structural[t.strip()] = cohorts_in_structural.get(t.strip(), 0) + 1
    # coverage of QTRP / qPerS-SID (accession, site) pairs by the structural table
    struct_keys = set()
    with open(C.REPO / "results/structural_site_features.csv", newline="", encoding="utf-8-sig") as h:
        for r in csv.DictReader(h):
            struct_keys.add((r["accession"], r["position"]))
    cover = {}
    for name, path in (("qtrp", C.REPO / "results/qtrp_sites_normalised.csv"),
                       ("qpers_sid", C.REPO / "results/qpers_sid_sites_normalised.csv")):
        keys = set()
        with open(path, newline="", encoding="utf-8-sig") as h:
            for r in csv.DictReader(h):
                keys.add((r["accession"], r["site"]))
        cover[name] = {"unique_accession_site_pairs": len(keys),
                       "covered_by_structural_table": len(keys & struct_keys),
                       "fraction": len(keys & struct_keys) / len(keys) if keys else None,
                       "columns_carry_full_protein_sequence": False,
                       "columns": list(csv.DictReader(open(path, newline="", encoding="utf-8-sig")).fieldnames)}
    seq_files = find_sequence_files(DESKTOP_ROOT)
    decision = {
        "demonstration_run": False,
        "reason": ("The census site table (results/ptm_census_sites.csv) and the human and Arabidopsis reference "
                   "proteomes it was built against (external/proteomes/*) are not present in the local repository copy. "
                   "The QTRP and qPerS-SID normalised tables are present and hash-identical to the census inputs, "
                   "but they carry accession, site and the identified peptide only, not the full protein sequence "
                   "needed for the 25 digest features, for enumerating NEG_A and for sequence-window chemistry "
                   "features. The only local proteomes are mouse and rice, and the only human sequences are "
                   "the pLMSNOSite benchmark, none of which covers a Figure 5a cohort; the structural feature table "
                   "covers other cohorts. Network access is not permitted for this item. Per the item's rule, no "
                   "demonstration was fitted."),
    }
    out = {"label": "POST HOC revision analysis 2026-09-30, item I_fig5a_estimand; input availability for task 3",
           "required_inputs": need, "structural_table_cohorts": cohorts_in_structural,
           "structural_coverage_of_normalised_tables": cover,
           "sequence_files_found_under_project_tree": seq_files, "decision": decision}
    (C.OUT / "t3_input_availability.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    C.record_inputs("s03_input_availability", [C.CENSUS_AUDIT, C.CENSUS_SITES_AUDIT,
                                               C.REPO / "results/structural_site_features.csv",
                                               C.REPO / "results/qtrp_sites_normalised.csv",
                                               C.REPO / "results/qpers_sid_sites_normalised.csv"])
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk in ("exists", "sha256_matches_expected")}
                      for k, v in need.items()}, indent=1, ensure_ascii=False))
    print(json.dumps(cohorts_in_structural))
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "columns"} for k, v in cover.items()}, indent=1))
    for s in seq_files:
        print(s)


if __name__ == "__main__":
    main()
