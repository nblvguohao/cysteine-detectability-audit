"""SFE-006 on the authors' own cohort: reproduction, not transfer.

`reports/PHASE2_CLAIMS_UNDER_DETECTABILITY_CONTROL.md` reported SFE-006 - "sulfenylation sites
are flanked by an enrichment of the positively charged residues K and R", Bui 2016's most
pronounced feature - as the one claim that VANISHES under detectability control. That run was a
**transfer test** on an Arabidopsis sulfenylation cohort, because the authors' own positive set
was unreachable; `results/phase2b_claim_retest_combined.csv` records it as
`is_transfer=True, baseline_reproduction=not_reproduced`, and
`reports/SISTER_JOURNAL_GAP_REGISTER.md` names that as gap G1 - the single decisive gap, because
no claim changed under control was a reproduction on the authors' own data.

The material arrived: Yang et al. 2014 (Nat Commun 5:4776) Supplementary Dataset 2, which is the
positive source behind Bui 2016's dataset (SOHSite = Yang 2014 plus RedoxDB). This round rebuilds
the claim on the authors' own species, experiment and negative rule.

**Predeclared construction, fixed before the run.**

* Positives: the `Total` sheet of Supplementary Dataset 2 (1106 rows), mapped to UniProt by
  PEPTIDE rather than by the sheet's `ID` column, which carries gene symbols and numeric ids
  rather than accessions. For each row the modified cysteine's offset inside the peptide is taken
  from the position of `C#` in `Modified Peptide Sequence`; the peptide is then located in the
  human reference proteome. Where a peptide occurs in more than one protein, the candidate whose
  resulting position equals the sheet's stated `Modified Site` is taken; if none matches and the
  peptide is not unique, the row is DROPPED as ambiguous and counted.
* Negatives: exactly the authors' rule - every other cysteine in those same proteins
  ("the rest of the cysteine residues in these proteins"). These are NOT same-run observed
  negatives, so the cohort is a NEG_A construction by the census taxonomy and its
  label-semantics tier is T3.
* Attribute: the authors' named offsets, unchanged from the transfer run -
  K or R at any of -10, -8, -7, -6, -4, -2, +4, +5, +6, +7, +8.
* Covariates, matching, size-matched random control, verdict rule and interval method are the
  phase-2 instrument as registered; nothing about the estimand is re-tuned here.

**Reading rule, fixed before the run.** Three outcomes, and the verdict function decides:
  (1) the effect vanishes or attenuates on the authors' own cohort while the size-matched random
      control holds -> the transfer result is CONFIRMED as a reproduction, and G1 closes;
  (2) the effect survives on the authors' own cohort -> the transfer result does NOT generalise
      to the authors' data, the vanishing is a property of the Arabidopsis cohort, and the paper
      must say so; the flagship claim is then that a published claim survives control while a
      same-content claim in another species does not;
  (3) both the controlled and the random-control intervals cross zero -> power statement, not
      evidence either way.
Whatever comes out is what gets reported. The one thing this script may not do is pick the
cohort after seeing the answer.

Writes results/phase2c_sfe006_reproduction{,_audit}.{csv,json} and a site table.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import platform
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from phase2_claim_cohorts import (
    SFE006_KR_OFFSETS, flank_residue_flag, site_feature_matrix,
)
from run_phase2_claims_under_detectability_control import (
    REPLICATES, SEED, run_contingency_claim, sha256_of,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
INTAKE = os.path.join(ROOT, "external", "intake", "manual_fetch_2026-09-16", "yang2014")
DATASET2 = os.path.join(INTAKE, "41467_2014_BFncomms5776_MOESM785_ESM.xlsx")
PROTEOME = os.path.join(ROOT, "external", "proteomes", "hsa.fasta.gz")
SCRIPT = os.path.abspath(__file__)
SHEET = "Total"


def read_proteome():
    sequences, name, buf = {}, None, []
    with gzip.open(PROTEOME, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(buf).upper()
                parts = line[1:].split()[0].split("|")
                name = parts[1] if len(parts) > 2 else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        sequences[name] = "".join(buf).upper()
    return sequences


def build_index(sequences):
    """One concatenated string plus an offset table, so peptide lookup is a C-level find."""
    parts, offsets, accs, pos = [], [], [], 0
    for acc, seq in sequences.items():
        parts.append(seq)
        offsets.append(pos)
        accs.append(acc)
        pos += len(seq) + 1
        parts.append("*")
    blob = "".join(parts)
    return blob, np.asarray(offsets), accs


def locate(blob, offsets, accs, sequences, peptide):
    out, start = [], blob.find(peptide)
    while start != -1:
        i = int(np.searchsorted(offsets, start, side="right") - 1)
        acc = accs[i]
        local = start - int(offsets[i])
        if local + len(peptide) <= len(sequences[acc]):
            out.append((acc, local))
        start = blob.find(peptide, start + 1)
    return out


def cys_offset(modified_peptide):
    """Index of the modified cysteine inside the peptide, from the `C#` marker."""
    plain = []
    for ch in modified_peptide:
        if ch == "#":
            if plain and plain[-1] == "C":
                return len(plain) - 1
            return None
        if ch.isalpha():
            plain.append(ch)
    return None


def load_rows():
    import openpyxl
    wb = openpyxl.load_workbook(DATASET2, read_only=True, data_only=True)
    ws = wb[SHEET]
    header = None
    rows = []
    for row in ws.iter_rows(values_only=True):
        if header is None:
            header = [str(c).strip() if c is not None else "" for c in row]
            continue
        rows.append(dict(zip(header, row)))
    wb.close()
    return rows


def build_cohort():
    sequences = read_proteome()
    blob, offsets, accs = build_index(sequences)
    rows = load_rows()
    stat = {"rows": len(rows), "no_peptide": 0, "no_marker": 0, "unmapped": 0,
            "ambiguous_dropped": 0, "mapped_by_stated_site": 0, "mapped_unique": 0,
            "stated_site_mismatch_but_unique": 0}
    sites, detail = set(), []
    for r in rows:
        pep = (r.get("Peptide Sequence") or "").strip().upper()
        mod = (r.get("Modified Peptide Sequence") or "").strip().upper()
        stated = (r.get("Modified Site") or "").strip()
        m = re.fullmatch(r"C(\d{1,5})", stated)
        stated_pos = int(m.group(1)) if m else None
        if not pep or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", pep):
            stat["no_peptide"] += 1
            continue
        off = cys_offset(mod)
        if off is None:
            stat["no_marker"] += 1
            continue
        cand = locate(blob, offsets, accs, sequences, pep)
        if not cand:
            stat["unmapped"] += 1
            continue
        hits = [(a, local + off + 1) for a, local in cand]
        pick = None
        if stated_pos is not None:
            exact = [h for h in hits if h[1] == stated_pos]
            if exact:
                pick, how = exact[0], "stated_site"
                stat["mapped_by_stated_site"] += 1
        if pick is None:
            if len({h[0] for h in hits}) == 1:
                pick, how = hits[0], "unique_protein"
                stat["mapped_unique"] += 1
                if stated_pos is not None and pick[1] != stated_pos:
                    stat["stated_site_mismatch_but_unique"] += 1
            else:
                stat["ambiguous_dropped"] += 1
                continue
        if sequences[pick[0]][pick[1] - 1] != "C":
            stat["ambiguous_dropped"] += 1
            continue
        sites.add(pick)
        detail.append({"accession": pick[0], "position": pick[1], "stated_site": stated,
                       "peptide": pep, "mapped_how": how, "gene_or_id": str(r.get("ID") or "")})

    proteins = sorted({a for a, _ in sites})
    keys, y = [], []
    for acc in proteins:
        seq = sequences[acc]
        for i, ch in enumerate(seq):
            if ch != "C":
                continue
            keys.append((acc, i + 1))
            y.append(1 if (acc, i + 1) in sites else 0)
    y = np.asarray(y, dtype=int)
    attribute = flank_residue_flag(sequences, keys, SFE006_KR_OFFSETS, {"K", "R"})
    cohort = {
        "claim_id": "SFE-006", "unit": "site", "keys": keys, "y": y,
        "attribute": attribute,
        "groups": np.asarray([a for a, _ in keys]),
        "covariates": site_feature_matrix(sequences, keys, "VIS10"),
        "covariate_set": "VIS10", "observed_depth": None,
        "claim_direction": "positive_preference",
        "attribute_label": "侧翼 −10/−8~−6/−4/−2/+4~+8 上出现 K 或 R（作者点名偏移，与转移检验一致）",
        "author_statistic": {"kind": "counts_only", "value": None,
                             "n_positive": 1443, "n_negative": 10521,
                             "text": "1443 阳性 / 10521 阴性 / 987 蛋白（Bui 2016 = Yang 2014 + RedoxDB）"},
        "reproduction_notes": [
            "阳性取自 Yang 2014 Supplementary Dataset 2 的 Total 表，按肽段映射到人源参考蛋白组",
            "阴性严格按作者规则：同一批蛋白里其余全部 Cys（NEG_A 构造，标签语义 T3）",
            "属性偏移与转移检验完全一致，未重新调参",
            "与 Bui 2016 的 1443 阳性不同：本轮只含 Yang 2014 成分，不含 RedoxDB 并集部分",
        ],
        "_stat": stat, "_detail": detail, "_n_proteins": len(proteins),
    }
    return cohort, sequences


def main():
    started = time.time()
    cohort, sequences = build_cohort()
    stat, detail = cohort.pop("_stat"), cohort.pop("_detail")
    n_proteins = cohort.pop("_n_proteins")
    print("mapping:", json.dumps(stat), flush=True)
    print(f"positives={int(cohort['y'].sum())} observations={len(cohort['y'])} "
          f"proteins={n_proteins}", flush=True)

    site_csv = os.path.join(RESULTS, "phase2c_yang2014_sulfenyl_sites.csv")
    with open(site_csv, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(detail[0]))
        w.writeheader()
        w.writerows(detail)

    record = run_contingency_claim(cohort)
    out = {k: v for k, v in record.items()}
    json.dump(out, open(os.path.join(RESULTS, "phase2c_sfe006_reproduction.json"), "w"),
              ensure_ascii=False, indent=2, default=str)

    audit = {
        "script": "scripts/run_phase2c_sfe006_reproduction.py", "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": platform.node(),
        "round": "SFE-006 reproduction on the authors' own cohort (gap G1)",
        "is_transfer": False, "baseline_reproduction_target": "Bui 2016 / Yang 2014 human sulfenylome",
        "mapping": stat, "n_proteins": n_proteins,
        "n_observations": int(len(cohort["y"])), "n_positive": int(cohort["y"].sum()),
        "attribute_offsets": SFE006_KR_OFFSETS,
        "negative_rule": "every other cysteine in the same proteins (the authors' rule; NEG_A / tier T3)",
        "instrument": {"replicates": REPLICATES, "seed": SEED,
                       "source": "run_phase2_claims_under_detectability_control.run_contingency_claim"},
        "inputs": {
            "external/intake/manual_fetch_2026-09-16/yang2014/41467_2014_BFncomms5776_MOESM785_ESM.xlsx": sha256_of(DATASET2),
            "external/proteomes/hsa.fasta.gz": sha256_of(PROTEOME),
        },
        "outputs": {
            "results/phase2c_yang2014_sulfenyl_sites.csv": sha256_of(site_csv),
        },
        "declared_difference_from_phase2": (
            "phase 2 ran SFE-006 as a transfer test on fps2020_ath_sulfenyl; this run uses the authors' own "
            "species, experiment and negative rule. The positive set is the Yang 2014 component only, not "
            "Bui 2016's union with RedoxDB, so it is a reproduction of the claim on the authors' data rather "
            "than a reconstruction of their exact 1443-site table."),
        "versions": {"python": sys.version, "numpy": np.__version__},
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }
    json.dump(audit, open(os.path.join(RESULTS, "phase2c_sfe006_reproduction_audit.json"), "w"),
              ensure_ascii=False, indent=2)
    for k in ("baseline_log2_or", "baseline_interval", "propensity_auc", "stratified_log2_or",
              "stratified_interval", "matched_log2_or", "matched_interval",
              "random_control_log2_or", "random_control_interval", "verdict", "primary_reason"):
        print(f"  {k}: {record.get(k)}", flush=True)
    print(f"done in {audit['elapsed_minutes']} min", flush=True)


if __name__ == "__main__":
    main()
