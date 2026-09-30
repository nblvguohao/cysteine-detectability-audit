"""Step 2 (POST HOC): reproduce the stored Artifact 3 numbers and build the per-group table.

Reproduction follows scripts/recompute_artefact3_osa_abundance_2026-09-19.py and
scripts/artefact3_effect_sizes_and_matching_2026-09-19.py of the internal repository (read-only).
Their functions are COPIED VERBATIM below (the original modules execute on import and read paths
that do not exist in this workspace). The only input that cannot be the same is the FASTA: the
original read external/proteomes/osa.fasta.gz (sha256 9df83a98...), which is not on disk. Two
substitutes are run and reported side by side, never tuned to the stored numbers:
  SEQ_REF   UP000059680 (UniProt 2026_03) only, original "first accession found" rule
  SEQ_FULL  UP000059680 + sequences fetched in s01 (active UniProtKB or UniParc), same rule

Outputs:
  results/F_rice_artifact3/s02_reproduction.csv       stored vs reproduced, per quantity
  results/F_rice_artifact3/s02_groups.csv             per-group table used by later steps
  results/F_rice_artifact3/s02_reproduction.json
"""
import random

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa

# ------------------------------------------------------------------ load (as the original) -----
peps = read_tsv(PEP)
pos_accs = {(x.get("Leading razor protein") or "").strip() for x in peps}
pos_accs.discard("")
prot = read_tsv(PRO)
IBAQ = [c for c in prot[0] if c.endswith("PG.IBAQ")]
QTY = [c for c in prot[0] if c.endswith("PG.Quantity")]


def group_value(row, cols, first_segment):          # verbatim
    vals = []
    for c in cols:
        raw = (row.get(c) or "").strip()
        if not raw or raw == "NaN":
            continue
        tok = raw.split(";")[0] if first_segment else raw
        try:
            v = float(tok)
        except ValueError:
            continue
        if v > 0:
            vals.append(v)
    return float(np.median(vals)) if vals else None


def n_runs(row, cols, first_segment):
    n = 0
    for c in cols:
        raw = (row.get(c) or "").strip()
        if not raw or raw == "NaN":
            continue
        tok = raw.split(";")[0] if first_segment else raw
        try:
            if float(tok) > 0:
                n += 1
        except ValueError:
            pass
    return n


seq_ref = read_fasta_any(FASTA_REF)
seq_full, seq_src = load_all_sequences()


def build_rows(seqs):
    rows, n_no_abund = [], 0
    for p in prot:
        accs = [a.strip() for a in (p.get("PG.ProteinAccessions") or "").split(";") if a.strip()]
        if not accs:
            continue
        ib = group_value(p, IBAQ, True)
        qt = group_value(p, QTY, False)
        if ib is None:
            n_no_abund += 1
            continue
        s, s_acc = None, None
        for a in accs:
            if a in seqs:
                s, s_acc = seqs[a], a
                break
        rows.append({"group_id": accs[0], "accessions": ";".join(accs), "n_accessions": len(accs),
                     "persulfidated": int(any(a in pos_accs for a in accs)),
                     "ibaq_median": ib, "quantity_median": qt,
                     "n_runs_ibaq": n_runs(p, IBAQ, True),
                     "seq_accession": s_acc,
                     "length": len(s) if s else None,
                     "cys_count": s.count("C") if s else None})
    return rows, n_no_abund


def ratio_of_medians(a, b):                          # verbatim
    return float(np.median(a) / np.median(b)) if len(a) and len(b) and np.median(b) > 0 else float("nan")


def boot_ratio(values, labels, n=5000, seed=SEED_ORIGINAL):   # verbatim (seed name changed)
    rng = np.random.default_rng(seed)
    v, l = np.asarray(values, float), np.asarray(labels, int)
    idx = np.arange(len(v)); out = []
    for _ in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        vs, ls = v[s], l[s]
        if ls.sum() == 0 or (1 - ls).sum() == 0: continue
        mb = np.median(vs[ls == 0])
        if mb > 0: out.append(np.median(vs[ls == 1]) / mb)
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))) if out else (float("nan"),) * 2


def mwu_z(x, y):                                     # verbatim
    x, y = np.asarray(x, float), np.asarray(y, float)
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2: return float("nan")
    allv = np.concatenate([x, y])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(len(allv), float)
    sorted_v = allv[order]
    i = 0
    while i < len(sorted_v):
        j = i
        while j + 1 < len(sorted_v) and sorted_v[j + 1] == sorted_v[i]: j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    R1 = ranks[:nx].sum()
    U = R1 - nx * (nx + 1) / 2.0
    mu = nx * ny / 2.0
    _, counts = np.unique(allv, return_counts=True)
    tie = (counts ** 3 - counts).sum()
    N = nx + ny
    sd = np.sqrt(nx * ny / 12.0 * ((N + 1) - tie / (N * (N - 1))))
    return float((U - mu) / sd) if sd > 0 else float("nan")


def stratified_z(rows_, key, nstrata=10):            # verbatim
    vals = np.array([r["ibaq_median"] for r in rows_], float)
    edges = np.percentile(vals, np.linspace(0, 100, nstrata + 1))
    edges[0] -= 1e-9; edges[-1] += 1e-9
    zs, used, dropped = [], 0, 0
    for k in range(nstrata):
        sel = [r for r in rows_ if edges[k] < r["ibaq_median"] <= edges[k + 1] and r[key] is not None]
        a = [r[key] for r in sel if r["persulfidated"]]
        b = [r[key] for r in sel if not r["persulfidated"]]
        if len(a) < 10 or len(b) < 10: dropped += 1; continue
        z = mwu_z(a, b)
        if not np.isnan(z): zs.append(z); used += 1
    if not zs: return float("nan"), used, dropped
    return float(np.sum(zs) / np.sqrt(len(zs))), used, dropped


def rank_biserial(x, y):                             # verbatim (effect-size script)
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2 or len(y) < 2: return float("nan")
    allv = np.concatenate([x, y]); order = allv.argsort(kind="mergesort")
    ranks = np.empty(len(allv), float); sv = allv[order]; i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]: j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    U = ranks[:len(x)].sum() - len(x) * (len(x) + 1) / 2.0
    return float(2.0 * U / (len(x) * len(y)) - 1.0)


def by_decile(have, key):                            # verbatim, 'have' passed explicitly
    vals = np.array([r["ibaq_median"] for r in have], float)
    edges = np.percentile(vals, np.linspace(0, 100, 11)); edges[0] -= 1e-9; edges[-1] += 1e-9
    rs, ns = [], []
    for k in range(10):
        sel = [r for r in have if edges[k] < r["ibaq_median"] <= edges[k + 1]]
        a = [r[key] for r in sel if r["persulfidated"]]
        b = [r[key] for r in sel if not r["persulfidated"]]
        if len(a) < 10 or len(b) < 10: continue
        rs.append(rank_biserial(a, b)); ns.append(len(a) + len(b))
    return (float(np.mean(rs)) if rs else float("nan")), len(rs)


def matched_11(have, key):                           # verbatim, 'have' passed explicitly
    pos = sorted([r for r in have if r["persulfidated"]], key=lambda r: r["ibaq_median"])
    neg = sorted([r for r in have if not r["persulfidated"]], key=lambda r: r["ibaq_median"])
    negv = np.log10([r["ibaq_median"] for r in neg])
    used = np.zeros(len(neg), bool)
    A, B = [], []
    for p in pos:
        t = np.log10(p["ibaq_median"])
        d = np.abs(negv - t); d[used] = np.inf
        j = int(d.argmin())
        if not np.isfinite(d[j]): break
        used[j] = True
        A.append(p[key]); B.append(neg[j][key])
    return A, B


def seq_stats(rows):
    have = [r for r in rows if r["length"] is not None]
    out = {"n_groups_with_sequence": len(have),
           "n_pos_with_sequence": sum(r["persulfidated"] for r in have),
           "n_neg_with_sequence": sum(1 - r["persulfidated"] for r in have)}
    for key, lab in (("length", "length"), ("cys_count", "cys")):
        a = [r[key] for r in have if r["persulfidated"]]
        b = [r[key] for r in have if not r["persulfidated"]]
        out[f"r_before_{lab}"] = rank_biserial(a, b)
        out[f"z_before_{lab}"] = mwu_z(a, b)
        out[f"r_decile_{lab}"], out[f"deciles_used_{lab}"] = by_decile(have, key)
        out[f"z_strat_{lab}"] = stratified_z(have, key)[0]
        A, B = matched_11(have, key)
        out[f"r_1to1_{lab}"] = rank_biserial(A, B)
        out[f"z_1to1_{lab}"] = mwu_z(A, B)
    return out


# ------------------------------------------------------------------ run -----
rows_ref, n_no_abund = build_rows(seq_ref)
rows_full, _ = build_rows(seq_full)
pos = [r for r in rows_ref if r["persulfidated"]]
neg = [r for r in rows_ref if not r["persulfidated"]]
vals = [r["ibaq_median"] for r in rows_ref]
labs = [r["persulfidated"] for r in rows_ref]
ratio = ratio_of_medians([r["ibaq_median"] for r in pos], [r["ibaq_median"] for r in neg])
ci = boot_ratio(vals, labs)
q = [r for r in rows_ref if r["quantity_median"]]
q_ratio = ratio_of_medians([r["quantity_median"] for r in q if r["persulfidated"]],
                           [r["quantity_median"] for r in q if not r["persulfidated"]])
rng = random.Random(SEED_ORIGINAL)
sh = [r["persulfidated"] for r in rows_ref]; rng.shuffle(sh)
sh_ratio = ratio_of_medians([r["ibaq_median"] for r, s in zip(rows_ref, sh) if s],
                            [r["ibaq_median"] for r, s in zip(rows_ref, sh) if not s])
sh_ci = boot_ratio([r["ibaq_median"] for r in rows_ref], sh)
st_ref = seq_stats(rows_ref)
st_full = seq_stats(rows_full)

stored = {r["quantity"]: r for r in csv.DictReader(open(STORED_A3, encoding="utf-8"))}
eff = {r["quantity"]: r for r in csv.DictReader(open(REPO_EFFECT_SIZES, encoding="utf-8"))}


def f(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else x


table = [
    ("n_protein_groups_in_quantification", stored["n_protein_groups_in_quantification"]["value"], len(prot), None),
    ("n_groups_dropped_no_usable_iBAQ", stored["n_groups_dropped_no_usable_iBAQ"]["value"], n_no_abund, None),
    ("n_groups_analysed", stored["n_groups_analysed"]["value"], len(rows_ref), None),
    ("n_persulfidated", stored["n_persulfidated"]["value"], len(pos), None),
    ("n_unmodified", stored["n_unmodified"]["value"], len(neg), None),
    ("n_unique_leading_razor_proteins_in_peptide_table", stored["n_unique_leading_razor_proteins_in_peptide_table"]["value"], len(pos_accs), None),
    ("median_iBAQ_ratio", stored["median_iBAQ_ratio_persulfidated_over_unmodified"]["value"], ratio, None),
    ("median_iBAQ_ratio_CI_low", stored["median_iBAQ_ratio_persulfidated_over_unmodified"]["interval_95"].strip("[]").split(",")[0], ci[0], None),
    ("median_iBAQ_ratio_CI_high", stored["median_iBAQ_ratio_persulfidated_over_unmodified"]["interval_95"].strip("[]").split(",")[1], ci[1], None),
    ("median_Quantity_ratio_raw_intensity", stored["sensitivity_median_Quantity_ratio"]["value"], q_ratio, None),
    ("shuffled_label_ratio", stored["negative_control_shuffled_label_ratio"]["value"], sh_ratio, None),
    ("shuffled_label_ratio_CI_low", stored["negative_control_shuffled_label_ratio"]["interval_95"].strip("[]").split(",")[0], sh_ci[0], None),
    ("shuffled_label_ratio_CI_high", stored["negative_control_shuffled_label_ratio"]["interval_95"].strip("[]").split(",")[1], sh_ci[1], None),
    ("n_groups_with_sequence", stored["n_groups_with_sequence"]["value"], st_ref["n_groups_with_sequence"], st_full["n_groups_with_sequence"]),
    ("n_pos_with_sequence", 1062, st_ref["n_pos_with_sequence"], st_full["n_pos_with_sequence"]),
    ("n_neg_with_sequence", 5301, st_ref["n_neg_with_sequence"], st_full["n_neg_with_sequence"]),
    ("length_r_before", eff["protein_length_r_before_matching"]["value"], st_ref["r_before_length"], st_full["r_before_length"]),
    ("length_z_before", eff["protein_length_z_before_matching"]["value"], st_ref["z_before_length"], st_full["z_before_length"]),
    ("length_r_decile_mean", eff["protein_length_r_decile_weighted_after"]["value"], st_ref["r_decile_length"], st_full["r_decile_length"]),
    ("length_z_stratified", eff["protein_length_z_stratified_after"]["value"], st_ref["z_strat_length"], st_full["z_strat_length"]),
    ("length_r_1to1", eff["protein_length_r_1to1_matched"]["value"], st_ref["r_1to1_length"], st_full["r_1to1_length"]),
    ("length_z_1to1", eff["protein_length_z_1to1_matched"]["value"], st_ref["z_1to1_length"], st_full["z_1to1_length"]),
    ("cys_r_before", eff["cysteine_count_r_before_matching"]["value"], st_ref["r_before_cys"], st_full["r_before_cys"]),
    ("cys_z_before", eff["cysteine_count_z_before_matching"]["value"], st_ref["z_before_cys"], st_full["z_before_cys"]),
    ("cys_r_decile_mean", eff["cysteine_count_r_decile_weighted_after"]["value"], st_ref["r_decile_cys"], st_full["r_decile_cys"]),
    ("cys_z_stratified", eff["cysteine_count_z_stratified_after"]["value"], st_ref["z_strat_cys"], st_full["z_strat_cys"]),
    ("cys_r_1to1", eff["cysteine_count_r_1to1_matched"]["value"], st_ref["r_1to1_cys"], st_full["r_1to1_cys"]),
    ("cys_z_1to1", eff["cysteine_count_z_1to1_matched"]["value"], st_ref["z_1to1_cys"], st_full["z_1to1_cys"]),
]
out = []
for name, st, rep, rep_full in table:
    stf = float(st)
    repf = float(rep)
    # stored values carry 4 decimals (or are integers): agreement = equal after rounding to 4 dp
    same = abs(round(repf, 4) - stf) < 1e-9
    out.append({"quantity": name, "stored": st, "reproduced_SEQ_REF": rep,
                "reproduced_SEQ_FULL": rep_full if rep_full is not None else "",
                "agrees_at_stored_precision_SEQ_REF": bool(same)})
pd.DataFrame(out).to_csv(RES / "s02_reproduction.csv", index=False)

# ------------------------------------------------------------------ per-group table -----
g = pd.DataFrame(rows_full).rename(columns={"length": "length_first_found", "cys_count": "cys_first_found",
                                            "seq_accession": "seq_accession_first_found"})
g["log10_ibaq"] = np.log10(g["ibaq_median"])
g["log10_quantity"] = np.log10(g["quantity_median"].astype(float))
g["lead_seq_available"] = g["group_id"].map(lambda a: a in seq_full)
g["lead_seq_source"] = g["group_id"].map(lambda a: seq_src.get(a, ""))
g["lead_in_reference_proteome"] = g["group_id"].map(lambda a: a in seq_ref)
g["length_ref_rule"] = [r["length"] for r in rows_ref]
g["cys_ref_rule"] = [r["cys_count"] for r in rows_ref]
g.to_csv(RES / "s02_groups.csv", index=False)

write_json(RES / "s02_reproduction.json", {
    "analysis_label": "POST HOC revision analysis 2026-09-30 (not registered)",
    "stats_SEQ_REF": {k: f(v) for k, v in st_ref.items()},
    "stats_SEQ_FULL": {k: f(v) for k, v in st_full.items()},
    "iBAQ_ratio": ratio, "iBAQ_ratio_ci": ci, "quantity_ratio": q_ratio,
    "shuffled_ratio": sh_ratio, "shuffled_ci": sh_ci,
    "n_groups_lead_seq_available": int(g["lead_seq_available"].sum()),
    "lead_seq_source_counts": g["lead_seq_source"].value_counts().to_dict(),
})
print(pd.DataFrame(out).to_string())
