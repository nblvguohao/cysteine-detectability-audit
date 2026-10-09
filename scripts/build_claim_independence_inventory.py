"""How independent are the re-tested claims from one another? Count it, do not assume it.

**Why this round exists.** The manuscript quotes a tally over re-tested claims - 35 claims, 28 with a
measured baseline, of which 12 survive, 11 undecidable, and so on. Every one of those five categories
is written with 28 as its denominator, and a reader will take 28 as 28 independent tests. That has
never been checked. Claims are extracted from papers, and one paper can contribute several claims
that are all re-tested against the SAME site table, in the same species, through the same
enrichment chemistry. If they do, the tally's denominator is inflated in exactly the way this
manuscript accuses the field of inflating its own.

The practice is borrowed from a rare-disease project on the shared cluster, which publishes a
per-disease table of `patient_units / publication_units / family_safe_units / provenance_units` and
then reports the collapsed number rather than the raw one. This script is the same idea applied to
our claims.

**Definitions, fixed before the run.**

  * **L1, the claim** - one row of results/claim_verdict_tally_2026-09-17.csv. 35 of them. This is
    the unit the manuscript currently counts in.
  * **L2, the paper** - the DOI recorded for the claim in results/ptm_site_preference_claims.csv.
    Two claims sharing a DOI were extracted from the same article.
  * **L3, the site table** - `sites_table_locator` from the same file, normalised only by stripping
    whitespace and lower-casing. This is the string the extractor wrote to say WHERE the positives
    came from. It is the closest thing on disk to "the same underlying measurement", and it is
    coarser than a hash of the table: two claims may share a locator and still use different
    columns of it, and a claim whose locator is blank is counted as its own unit rather than pooled
    with anything. Both limitations are reported, not adjusted away.
  * **has_measured_baseline** - taken verbatim from the tally, not recomputed. The two
    administrative verdicts (`out_of_instrument_scope`, `blocked_material_unreachable`) are the ones
    with no baseline.

**What this script does NOT do.** It does not re-test anything, does not recompute a verdict, does
not change a number, and does not propose a corrected tally. It reports how many independent units
the existing tally rests on, so the manuscript can say so in its own sentence. Whether to re-weight
anything is a decision for the corresponding author, not for this script.

**Reading rule, fixed before the run.** A paper contributing three or more re-tested claims with a
measured baseline is flagged `concentration_hotspot`. Three is the smallest count at which a single
article can carry more than a tenth of the 28-claim denominator; the threshold is stated here rather
than chosen after seeing the distribution.

Writes results/claim_independence_inventory.csv (one row per claim),
results/claim_independence_by_source.csv (one row per paper) and
results/claim_independence_audit.json. Overwrites nothing else.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
import time
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
SCRIPT = os.path.abspath(__file__)

TALLY = os.path.join(RESULTS, "claim_verdict_tally_2026-09-17.csv")
CLAIMS = os.path.join(RESULTS, "ptm_site_preference_claims.csv")
VERDICT_TABLES = [
    ("phase2", "phase2_claim_retest.csv"),
    ("phase2b", "phase2b_claim_retest_combined.csv"),
    ("phase2d", "phase2d_claim_retest.csv"),
    ("phase2e", "phase2e_claim_retest.csv"),
]
HOTSPOT_MIN_CLAIMS = 3
NO_BASELINE = {"out_of_instrument_scope", "blocked_material_unreachable"}


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    started = time.time()
    tally = read_csv(TALLY)
    claims = {r["claim_id"]: r for r in read_csv(CLAIMS)}

    # size fields, taken from whichever round actually re-tested the claim
    sizes = {}
    for round_name, filename in VERDICT_TABLES:
        for row in read_csv(os.path.join(RESULTS, filename)):
            if row.get("caliber") not in (None, "", "primary"):
                continue
            key = (row["claim_id"], round_name)
            sizes[key] = {
                "unit": row.get("unit", ""),
                "n_observations": row.get("n_observations", ""),
                "n_positive": row.get("n_positive", ""),
                "n_proteins": row.get("n_proteins", ""),
            }

    rows, missing_claim_metadata = [], []
    for row in tally:
        cid = row["claim_id"]
        meta = claims.get(cid)
        if meta is None:
            missing_claim_metadata.append(cid)
            meta = {}
        locator = (meta.get("sites_table_locator") or "").strip()
        size = sizes.get((cid, row["round"]), {})
        rows.append({
            "claim_id": cid,
            "family": meta.get("family", ""),
            "paper_short": meta.get("paper_short", ""),
            "doi": (meta.get("doi") or "").strip(),
            "pmid": (meta.get("pmid") or "").strip(),
            "year": meta.get("year", ""),
            "round": row["round"],
            "verdict": row["verdict"],
            "has_measured_baseline": row["has_measured_baseline"],
            "unit": size.get("unit", ""),
            "n_observations": size.get("n_observations", ""),
            "n_positive": size.get("n_positive", ""),
            "n_proteins": size.get("n_proteins", ""),
            "site_table_locator": locator,
            "site_table_key": locator.lower() if locator else f"__blank__{cid}",
            "negative_set_class": meta.get("negative_set_class", ""),
            "species": meta.get("species", ""),
        })

    with_baseline = [r for r in rows if r["has_measured_baseline"] == "1"]

    def group(records, key):
        out = defaultdict(list)
        for r in records:
            out[r[key] or f"__blank__{r['claim_id']}"].append(r)
        return out

    by_doi_all = group(rows, "doi")
    by_doi_base = group(with_baseline, "doi")
    by_table_all = group(rows, "site_table_key")
    by_table_base = group(with_baseline, "site_table_key")

    source_rows = []
    for doi, members in sorted(by_doi_all.items(), key=lambda kv: -len(kv[1])):
        base = by_doi_base.get(doi, [])
        source_rows.append({
            "doi": doi,
            "paper_short": members[0]["paper_short"],
            "family": members[0]["family"],
            "year": members[0]["year"],
            "n_claims_retested": len(members),
            "n_claims_with_measured_baseline": len(base),
            "claim_ids": ";".join(sorted(m["claim_id"] for m in members)),
            "verdicts_with_baseline": ";".join(
                f"{k}={v}" for k, v in sorted(Counter(m["verdict"] for m in base).items())),
            "distinct_site_tables": len({m["site_table_key"] for m in members}),
            "concentration_hotspot": int(len(base) >= HOTSPOT_MIN_CLAIMS),
        })

    def share(groups):
        """largest group's share of the claims it contains"""
        if not groups:
            return 0.0, 0
        biggest = max(len(v) for v in groups.values())
        total = sum(len(v) for v in groups.values())
        return (biggest / total if total else 0.0), biggest

    share_doi_base, biggest_doi_base = share(by_doi_base)
    share_tab_base, biggest_tab_base = share(by_table_base)

    blank_locator_claims = [r["claim_id"] for r in rows if not r["site_table_locator"]]
    hotspots = [s for s in source_rows if s["concentration_hotspot"]]

    for path, data, cols in (
        (os.path.join(RESULTS, "claim_independence_inventory.csv"), rows, list(rows[0].keys())),
        (os.path.join(RESULTS, "claim_independence_by_source.csv"), source_rows,
         list(source_rows[0].keys())),
    ):
        with open(path + ".tmp", "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=cols)
            writer.writeheader()
            for item in data:
                writer.writerow(item)
        os.replace(path + ".tmp", path)

    audit = {
        "script": "scripts/build_claim_independence_inventory.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "what_this_counts": "how many independent papers and site tables the manuscript's claim "
                            "tally rests on; it re-tests nothing and changes no verdict",
        "levels": {"L1_claim": len(rows),
                   "L1_claim_with_measured_baseline": len(with_baseline),
                   "L2_paper": len(by_doi_all),
                   "L2_paper_with_measured_baseline": len(by_doi_base),
                   "L3_site_table": len(by_table_all),
                   "L3_site_table_with_measured_baseline": len(by_table_base)},
        "concentration": {
            "largest_paper_claims_with_baseline": biggest_doi_base,
            "largest_paper_share_of_the_28": round(share_doi_base, 4),
            "largest_site_table_claims_with_baseline": biggest_tab_base,
            "largest_site_table_share_of_the_28": round(share_tab_base, 4),
            "hotspot_threshold_claims": HOTSPOT_MIN_CLAIMS,
            "hotspots": [{"doi": h["doi"], "paper_short": h["paper_short"],
                          "n_with_baseline": h["n_claims_with_measured_baseline"],
                          "claim_ids": h["claim_ids"],
                          "verdicts": h["verdicts_with_baseline"]} for h in hotspots],
        },
        "verdict_tally_reproduced_from_the_inventory": dict(
            sorted(Counter(r["verdict"] for r in with_baseline).items())),
        "verdict_tally_all_claims": dict(sorted(Counter(r["verdict"] for r in rows).items())),
        "declared_limits": [
            "L3 is the extractor's locator string, not a hash of the table. Two claims sharing a "
            "locator may still use different columns of it; two claims with different locator "
            "wording may still rest on the same measurement. The count is therefore a bound on "
            "independence in one direction only and is reported as such.",
            f"{len(blank_locator_claims)} claims have no locator recorded and are each counted as "
            "their own site table rather than pooled, which can only make the independence count "
            "look BETTER than it is: " + (";".join(blank_locator_claims) or "none"),
            "papers are grouped by DOI as recorded; no attempt is made to detect two papers "
            "re-publishing one dataset, which would further reduce independence",
            "this inventory does not re-weight, correct or replace the manuscript's tally",
        ],
        "claims_in_the_tally_without_metadata": missing_claim_metadata,
        "inputs": {os.path.relpath(TALLY, ROOT): sha256_of(TALLY),
                   os.path.relpath(CLAIMS, ROOT): sha256_of(CLAIMS)},
        "versions": {"python": sys.version},
        "elapsed_seconds": round(time.time() - started, 2),
        "corrections": [],
    }
    path = os.path.join(RESULTS, "claim_independence_audit.json")
    with open(path + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)
    os.replace(path + ".tmp", path)

    lv = audit["levels"]
    print(f"claims {lv['L1_claim']} ({lv['L1_claim_with_measured_baseline']} with a measured "
          f"baseline) rest on {lv['L2_paper']} papers ({lv['L2_paper_with_measured_baseline']} "
          f"with a baseline) and {lv['L3_site_table']} site tables "
          f"({lv['L3_site_table_with_measured_baseline']} with a baseline)")
    c = audit["concentration"]
    print(f"largest single paper: {c['largest_paper_claims_with_baseline']} of "
          f"{lv['L1_claim_with_measured_baseline']} = {c['largest_paper_share_of_the_28']:.1%}; "
          f"largest single site table: {c['largest_site_table_claims_with_baseline']}")
    for h in c["hotspots"]:
        print(f"  hotspot {h['n_with_baseline']:2d}  {h['paper_short'][:44]:46s} {h['verdicts']}")
    print("tally reproduced:", audit["verdict_tally_reproduced_from_the_inventory"])


if __name__ == "__main__":
    main()
