"""The pKa round's verdict, and the settlement of the blind-prediction lock.

Everything this script does was fixed before it was pointed at instrument B's output:
`protocols/pka_round_preregistration_2026-09-17.json` and its amendments A1 (the minimum meaningful
difference, the A-versus-B comparison scope, and nine numbered predictions), A2 (the pinned PROPKA
build, the heavy-atom snapshot convention, the completeness criterion) and A3 (the covariate set,
the matching rule with its demotion precheck, and the power arithmetic). This script's own sha256 is
recorded in the lock BEFORE it is run, so that a prediction cannot be scored against a rule written
after the answer was seen.

THE CLAIM UNDER TEST is a NULL claim. SNO-003 / SNO-007 (Doulias 2010, PNAS) report that the pKa of
S-nitrosylated cysteines does not differ from that of unmodified cysteines - 10.0 +/- 2.10 (n = 142)
against 9.88 +/- 2.20 (n = 559), "not significantly different", measured with PROPKA 2.0 on crystal
structures. A null claim cannot be confirmed by an interval that covers zero in a design too weak to
see anything, so the reading rule is:

    the claim SURVIVES only if the interval covers the null AND the 95% cluster-bootstrap half-width
    of the baseline contrast is smaller than the minimum meaningful difference. Otherwise the same
    null-covering interval reads UNDECIDABLE. MMD primary 1.0 pKa unit (a tenfold change in
    thiolate fraction at fixed pH), secondaries 0.5 and 0.12 reported in the same table.

WHAT IS COMPUTED, per instrument:

  baseline      mean pKa of positives minus mean pKa of the in-protein background, unmatched
  stratified    the same difference within quintiles of the VIS10 propensity, combined
                inverse-variance
  matched       1:1 nearest neighbour on the propensity logit, caliper 0.2 x SD(logit), with a
                MANDATORY same-size random control (20 draws, seeds 20260915+i). The random control
                is what separates "the control moved the effect" from "the sample got smaller".
  precheck      if the share of unmatched positives exceeds ORPHAN_LIMIT = 0.05, the matched caliber
                is demoted to secondary and the stratified caliber becomes primary.

  Intervals are cluster bootstrap percentile intervals, 5,000 replicates, seed 20260915, clustering
  on PROTEIN (per the registration - this cohort is 120 randomly drawn mouse proteins, not a
  homology-grouped benchmark), implemented with multinomial cluster weights.

  Exclusions, from instrument A's own flags and applied identically to both instruments:
  metal-coordinated cysteines (supplementary table 3) and cysteines in a predicted disulfide.

SCOPE (A1). Instrument A's claim verdict uses A's full cohort. Instrument B's uses the proteins B
covers. The A-versus-B AGREEMENT test uses only cysteines both instruments value, and A's difference
is recomputed on that common subset for that purpose alone and labelled as such. Neither number
stands in for the other.

OUTCOME BRANCHES (registration): 1 both agree and the null survives; 2 both agree and the null
breaks; 3 the instruments disagree - no verdict is issued on the biology and the disagreement enters
the manuscript's third axis; 4 fewer than 30 usable positives - underpowered, no verdict.

SMOKE MODE. `--smoke` permutes `is_positive` WITHIN each protein with a fixed seed and writes to
*_smoke filenames marked non-authoritative. It exists so the machinery can be exercised end to end
without reading the real contrast, the way a formal run is preceded by a single-step smoke test. A
smoke run is not a verdict and its audit says so.

Writes results/pka_round_verdict.csv, results/pka_round_verdict_strata.csv,
results/pka_blind_prediction_scorecard.csv and results/pka_round_verdict_audit.json.
Overwrites nothing else.

Run in the project venv (numpy/sklearn/scipy; PROPKA is not needed here):
    PYTHONPATH=scripts /Users/lyuguohao/opt/persulfidation_v2/venv/bin/python \
        scripts/run_pka_round_verdict.py [--smoke]
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import sys
import time

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
SCRIPT = os.path.abspath(__file__)
sys.path.insert(0, os.path.dirname(SCRIPT))

from phase2_claim_cohorts import VIS10_INDEX  # noqa: E402
from run_cross_protease_detectability_probe import (  # noqa: E402
    PROTEASES, boundaries_for, features,
)

TRYPSIN = PROTEASES["Trypsin"]
STATIC = os.path.join(RESULTS, "pka_static_propka.csv")
ENSEMBLE = os.path.join(RESULTS, "pka_ensemble_propka.csv")
LOCK = os.path.join(ROOT, "protocols", "pka_blind_prediction_lock_2026-09-17.json")
PROTEOME = os.path.join(ROOT, "external", "proteomes", "phase2_mmu_uniprot.tsv.gz")

SEED = 20260915
N_BOOT = 5000
CHUNK = 500
N_RANDOM_DRAWS = 20
N_STRATA = 5
CALIPER_SD = 0.2
ORPHAN_LIMIT = 0.05
CLASS_FLOOR = 30
MMD_PRIMARY = 1.0
MMD_CALIBERS = (1.0, 0.5, 0.12)
AGREEMENT_SPEARMAN_FLOOR = 0.60


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def read_sequences(accessions):
    wanted, out = set(accessions), {}
    with gzip.open(PROTEOME, "rt", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["Entry"] in wanted:
                out[row["Entry"]] = row["Sequence"]
    return out


def cluster_weights(cluster_of_row, n_clusters, rng, n_rep):
    counts = rng.multinomial(n_clusters, np.full(n_clusters, 1.0 / n_clusters), size=n_rep)
    return counts[:, cluster_of_row].astype(float)


def weighted_difference(values, pos_mask, weights):
    pos_w = weights * pos_mask
    neg_w = weights * (~pos_mask)
    ps, ns = pos_w.sum(axis=1), neg_w.sum(axis=1)
    out = np.full(weights.shape[0], np.nan)
    ok = (ps > 0) & (ns > 0)
    out[ok] = ((pos_w @ values)[ok] / ps[ok]) - ((neg_w @ values)[ok] / ns[ok])
    return out


def bootstrap(fn, cluster_of_row, n_clusters):
    rng = np.random.default_rng(SEED)
    done, out = 0, []
    while done < N_BOOT:
        size = min(CHUNK, N_BOOT - done)
        out.append(fn(cluster_weights(cluster_of_row, n_clusters, rng, size)))
        done += size
    return np.concatenate(out)


def interval(samples):
    finite = samples[np.isfinite(samples)]
    if finite.size == 0:
        return float("nan"), float("nan")
    return float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))


def null_claim_verdict(low, high, half_width, mmd):
    """The registration's reading rule for a claim that asserts NO difference."""
    if not np.isfinite(low) or not np.isfinite(high):
        return "undecidable", "interval not estimable"
    covers = low <= 0 <= high
    powered = np.isfinite(half_width) and half_width < mmd
    if covers and powered:
        return "survives", (f"interval covers zero and the design resolves {mmd} pKa units "
                            f"(half-width {half_width:.4f})")
    if covers and not powered:
        return "undecidable", (f"interval covers zero but the design cannot resolve {mmd} pKa units "
                               f"(half-width {half_width:.4f})")
    return "null_broken_by_control" if False else "baseline_contradicts_claim", (
        "the interval excludes zero, so the claimed absence of a difference does not hold in this "
        "reconstruction")


def analyse(tag, values, pos_mask, protein, design, rows, strata_rows, notes):
    """One instrument, one cohort: baseline, power, stratified, matched, random control."""
    uniq = {p: i for i, p in enumerate(sorted(set(protein)))}
    cluster = np.array([uniq[p] for p in protein])
    n_pos, n_neg = int(pos_mask.sum()), int((~pos_mask).sum())

    if n_pos < CLASS_FLOOR:
        rows.append({"instrument": tag, "caliber": "baseline", "n_positive": n_pos,
                     "n_background": n_neg, "n_proteins": len(uniq), "point_estimate": "",
                     "ci_low": "", "ci_high": "", "half_width": "", "random_control": "",
                     "verdict": "out_of_instrument_scope",
                     "reason": f"fewer than {CLASS_FLOOR} usable positives (branch 4)"})
        notes[tag] = {"underpowered": True, "n_positive": n_pos}
        return

    point = weighted_difference(values, pos_mask, np.ones((1, values.size)))[0]
    draws = bootstrap(lambda w: weighted_difference(values, pos_mask, w), cluster, len(uniq))
    low, high = interval(draws)
    half = (high - low) / 2.0

    for mmd in MMD_CALIBERS:
        verdict, reason = null_claim_verdict(low, high, half, mmd)
        rows.append({"instrument": tag, "caliber": f"baseline | MMD={mmd}", "n_positive": n_pos,
                     "n_background": n_neg, "n_proteins": len(uniq),
                     "point_estimate": round(point, 4), "ci_low": round(low, 4),
                     "ci_high": round(high, 4), "half_width": round(half, 4),
                     "random_control": "",
                     "verdict": verdict + ("" if mmd == MMD_PRIMARY else " (secondary caliber)"),
                     "reason": reason})

    notes[tag] = {"baseline_point": round(point, 4), "baseline_ci": [round(low, 4), round(high, 4)],
                  "baseline_half_width": round(half, 4), "n_positive": n_pos,
                  "n_background": n_neg, "n_proteins": len(uniq),
                  "primary_verdict": null_claim_verdict(low, high, half, MMD_PRIMARY)[0]}

    if design is None:
        rows.append({"instrument": tag, "caliber": "matched", "n_positive": "", "n_background": "",
                     "n_proteins": "", "point_estimate": "", "ci_low": "", "ci_high": "",
                     "half_width": "", "random_control": "", "verdict": "out_of_instrument_scope",
                     "reason": "no sequence coverage for the VIS10 covariates on this cohort"})
        return

    mean, sd = design.mean(axis=0), design.std(axis=0)
    sd[sd == 0] = 1.0
    z = (design - mean) / sd
    model = LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs")
    model.fit(z, pos_mask.astype(int))
    score = model.predict_proba(z)[:, 1]
    logit = np.log(np.clip(score, 1e-9, 1 - 1e-9) / np.clip(1 - score, 1e-9, 1 - 1e-9))

    cuts = np.quantile(score, np.linspace(0, 1, N_STRATA + 1))
    cuts[0], cuts[-1] = -np.inf, np.inf
    stratum = np.clip(np.searchsorted(cuts, score, side="right") - 1, 0, N_STRATA - 1)
    pieces, ivw = [], []
    for k in range(N_STRATA):
        sel = stratum == k
        kp, kn = int((sel & pos_mask).sum()), int((sel & ~pos_mask).sum())
        strata_rows.append({"instrument": tag, "stratum": k, "n_positive": kp, "n_background": kn,
                            "included": kp >= 1 and kn >= 1,
                            "note": "" if kp >= 1 and kn >= 1 else "dropped: a class is empty"})
        if kp < 1 or kn < 1:
            continue
        sub = weighted_difference(values[sel], pos_mask[sel], np.ones((1, int(sel.sum()))))[0]
        if np.isfinite(sub):
            pieces.append(sub)
            ivw.append(1.0 / (1.0 / kp + 1.0 / kn))
    combined = float(np.average(pieces, weights=ivw)) if pieces else float("nan")

    def stratified_from_weights(w):
        out = np.empty(w.shape[0])
        for i in range(w.shape[0]):
            acc, accw = [], []
            for k in range(N_STRATA):
                sel = stratum == k
                if not ((sel & pos_mask).any() and (sel & ~pos_mask).any()):
                    continue
                sub = weighted_difference(values[sel], pos_mask[sel], w[i:i + 1, sel])[0]
                if np.isfinite(sub):
                    acc.append(sub)
                    accw.append(float(w[i, sel].sum()))
            out[i] = np.average(acc, weights=accw) if acc else np.nan
        return out

    s_low, s_high = interval(bootstrap(stratified_from_weights, cluster, len(uniq)))
    verdict, reason = null_claim_verdict(s_low, s_high, (s_high - s_low) / 2.0, MMD_PRIMARY)
    rows.append({"instrument": tag, "caliber": "stratified (VIS10 quintiles)", "n_positive": n_pos,
                 "n_background": n_neg, "n_proteins": len(uniq),
                 "point_estimate": round(combined, 4), "ci_low": round(s_low, 4),
                 "ci_high": round(s_high, 4), "half_width": round((s_high - s_low) / 2.0, 4),
                 "random_control": "", "verdict": verdict,
                 "reason": "primary caliber whenever matching is blocked by the orphan limit. "
                           + reason})

    caliper = CALIPER_SD * float(np.std(logit))
    pos_idx = np.flatnonzero(pos_mask)
    neg_idx = np.flatnonzero(~pos_mask)
    order = neg_idx[np.argsort(logit[neg_idx], kind="mergesort")]
    neg_vals = logit[order]
    used = np.zeros(order.size, dtype=bool)
    pairs = []
    for i in pos_idx[np.argsort(logit[pos_idx], kind="mergesort")]:
        j = np.searchsorted(neg_vals, logit[i])
        best, best_d = -1, np.inf
        for c in range(max(0, j - 50), min(order.size, j + 50)):
            if used[c]:
                continue
            d = abs(neg_vals[c] - logit[i])
            if d < best_d:
                best, best_d = c, d
        if best >= 0 and best_d <= caliper:
            used[best] = True
            pairs.append((i, order[best]))

    orphan = 1 - len(pairs) / max(1, pos_idx.size)
    blocked = orphan > ORPHAN_LIMIT
    notes[tag]["matching"] = {"caliper_logit": round(caliper, 6), "n_pairs": len(pairs),
                              "orphan_fraction": round(orphan, 4), "orphan_limit": ORPHAN_LIMIT,
                              "blocking": bool(blocked)}
    if not pairs:
        return

    m_pos = np.array([p for p, _ in pairs])
    m_neg = np.array([n for _, n in pairs])
    idx = np.concatenate([m_pos, m_neg])
    m_mask = np.concatenate([np.ones(m_pos.size, bool), np.zeros(m_neg.size, bool)])
    raw = cluster[idx]
    remap = {c: i for i, c in enumerate(sorted(set(raw.tolist())))}
    m_cluster = np.array([remap[c] for c in raw])
    m_values = values[idx]

    m_point = weighted_difference(m_values, m_mask, np.ones((1, m_values.size)))[0]
    m_low, m_high = interval(bootstrap(
        lambda w: weighted_difference(m_values, m_mask, w), m_cluster, len(remap)))

    rng_rc = np.random.default_rng(SEED)
    rc_points, rc_draws = [], []
    per_pool = N_BOOT // N_RANDOM_DRAWS
    for draw in range(N_RANDOM_DRAWS):
        picked = np.random.default_rng(SEED + draw).choice(neg_idx, size=m_pos.size, replace=False)
        pool = np.concatenate([m_pos, picked])
        rc_points.append(weighted_difference(values[pool], m_mask,
                                             np.ones((1, pool.size)))[0])
        rraw = cluster[pool]
        rmap = {c: i for i, c in enumerate(sorted(set(rraw.tolist())))}
        rcl = np.array([rmap[c] for c in rraw])
        rc_draws.append(weighted_difference(values[pool], m_mask,
                                            cluster_weights(rcl, len(rmap), rng_rc, per_pool)))
    rc_low, rc_high = interval(np.concatenate(rc_draws))
    rc_mean = float(np.mean(rc_points))

    verdict, reason = null_claim_verdict(m_low, m_high, (m_high - m_low) / 2.0, MMD_PRIMARY)
    rows.append({"instrument": tag,
                 "caliber": ("matched 1:1 (SECONDARY - blocked by the orphan limit)"
                             if blocked else "matched 1:1 (primary)"),
                 "n_positive": int(m_pos.size), "n_background": int(m_neg.size),
                 "n_proteins": len(remap), "point_estimate": round(m_point, 4),
                 "ci_low": round(m_low, 4), "ci_high": round(m_high, 4),
                 "half_width": round((m_high - m_low) / 2.0, 4),
                 "random_control": f"{round(rc_mean, 4)} [{round(rc_low, 4)}, {round(rc_high, 4)}]",
                 "verdict": verdict,
                 "reason": (f"orphan fraction {round(orphan, 4)} against the limit {ORPHAN_LIMIT}: "
                            + ("BLOCKING, this row is secondary. " if blocked else "")
                            + reason + ". Same-size random control: matched positives kept, an "
                            f"equal number of background cysteines drawn at random "
                            f"{N_RANDOM_DRAWS} times, seeds {SEED}+i")})


def main():
    started = time.time()
    smoke = "--smoke" in sys.argv
    suffix = "_smoke" if smoke else ""

    static = read_csv(STATIC)
    ensemble = {(r["accession"], int(r["position"])): r for r in read_csv(ENSEMBLE)}

    usable = [r for r in static
              if r["pka_static"] != "" and r["is_metal_coordinated"] == "0"
              and r["in_disulfide"] == "0"]

    if smoke:
        rng = np.random.default_rng(20260915)
        by_protein = {}
        for r in usable:
            by_protein.setdefault(r["accession"], []).append(r)
        for group in by_protein.values():
            labels = [r["is_positive"] for r in group]
            rng.shuffle(labels)
            for r, lab in zip(group, labels):
                r["is_positive"] = lab

    sequences = read_sequences({r["accession"] for r in usable})
    bounds_cache, vis_cache = {}, {}

    def vis10(accession, position):
        key = (accession, position)
        if key in vis_cache:
            return vis_cache[key]
        seq = sequences.get(accession)
        if not seq or position > len(seq) or seq[position - 1] != "C":
            vis_cache[key] = None
            return None
        if accession not in bounds_cache:
            bounds_cache[accession] = boundaries_for(seq, TRYPSIN)
        row = features(seq, position - 1, bounds_cache[accession], TRYPSIN)
        vis_cache[key] = [row[i] for i in VIS10_INDEX]
        return vis_cache[key]

    def cohort(records, value_key, value_source):
        values, mask, protein, design, uncovered = [], [], [], [], 0
        for r in records:
            if value_source == "static":
                value = float(r["pka_static"])
            else:
                key = (r["accession"], int(r["position"]))
                if key not in ensemble:
                    continue
                value = float(ensemble[key][value_key])
            block = vis10(r["accession"], int(r["position"]))
            if block is None:
                uncovered += 1
                continue
            values.append(value)
            mask.append(r["is_positive"] == "1")
            protein.append(r["accession"])
            design.append(block)
        return (np.array(values), np.array(mask, dtype=bool), protein,
                np.array(design, dtype=float) if design else None, uncovered)

    b_accessions = {a for a, _ in ensemble}
    rows, strata_rows, notes, coverage = [], [], {}, {}

    for tag, records, source in (
            ("A_static_full_cohort", usable, "static"),
            ("B_ensemble", [r for r in usable if r["accession"] in b_accessions],
             "pka_ensemble_mean"),
            ("A_static_on_B_subset", [r for r in usable if r["accession"] in b_accessions],
             "static")):
        values, mask, protein, design, uncovered = cohort(records, "pka_ensemble_mean", source)
        coverage[tag] = {"rows_considered": len(records), "rows_used": int(values.size),
                         "dropped_no_vis10_coverage": uncovered}
        analyse(tag, values, mask, protein, design, rows, strata_rows, notes)

    # A versus B agreement, on cysteines both instruments value
    common = [r for r in usable
              if (r["accession"], int(r["position"])) in ensemble]
    a_vals = np.array([float(r["pka_static"]) for r in common])
    b_vals = np.array([float(ensemble[(r["accession"], int(r["position"]))]["pka_ensemble_mean"])
                       for r in common])
    if smoke:
        # P6 and P7 do not depend on the labels, so a smoke run would otherwise reveal them before
        # the formal run. Break the A-to-B pairing so the agreement block is meaningless in smoke
        # mode, and say so in the audit.
        b_vals = np.random.default_rng(20260915).permutation(b_vals)
    rho = float(spearmanr(a_vals, b_vals).statistic) if a_vals.size > 2 else float("nan")
    shift = b_vals - a_vals
    agreement = {"n_common_cysteines": int(a_vals.size),
                 "spearman_rho": round(rho, 4),
                 "median_shift_B_minus_A": round(float(np.median(shift)), 4),
                 "mean_shift_B_minus_A": round(float(np.mean(shift)), 4),
                 "sd_shift": round(float(np.std(shift, ddof=1)), 4) if a_vals.size > 1 else None}

    verdict_a = notes.get("A_static_full_cohort", {}).get("primary_verdict")
    verdict_b = notes.get("B_ensemble", {}).get("primary_verdict")
    if verdict_a == "out_of_instrument_scope" or verdict_b == "out_of_instrument_scope":
        branch = "4_too_few_sites_survive_the_exclusions"
    elif verdict_a != verdict_b:
        branch = "3_the_instruments_disagree"
    elif verdict_a == "survives":
        branch = "1_both_instruments_agree_and_the_null_survives"
    elif verdict_a in ("baseline_contradicts_claim", "null_broken_by_control"):
        branch = "2_both_instruments_agree_and_the_null_breaks"
    else:
        branch = "3_the_instruments_disagree" if verdict_a != verdict_b else \
                 "1_both_instruments_agree_and_the_null_survives"
        branch = ("undecidable_on_both_instruments" if verdict_a == "undecidable" else branch)

    # score the blind predictions
    lock = json.load(open(LOCK, encoding="utf-8"))
    a_note = notes.get("A_static_full_cohort", {})
    b_note = notes.get("B_ensemble", {})
    a_sub = notes.get("A_static_on_B_subset", {})
    unshuffled = 0  # this round shuffles nothing; kept for P4's wording
    score = []

    def add(pid, statement, observed, correct):
        score.append({"prediction_id": pid, "statement": statement, "observed": observed,
                      "correct": "" if correct is None else str(bool(correct))})

    p1 = a_note.get("baseline_point")
    add("P1", "instrument A baseline difference in [-0.60, +0.20]",
        f"{p1}", None if p1 is None else (-0.60 <= p1 <= 0.20))
    add("P2", "that difference is negative", f"{p1}", None if p1 is None else (p1 < 0))
    both_cover = None
    if a_note.get("baseline_ci") and b_note.get("baseline_ci"):
        both_cover = (a_note["baseline_ci"][0] <= 0 <= a_note["baseline_ci"][1]
                      and b_note["baseline_ci"][0] <= 0 <= b_note["baseline_ci"][1])
    add("P3", "both instruments' 95% intervals cover zero",
        f"A {a_note.get('baseline_ci')}, B {b_note.get('baseline_ci')}", both_cover)
    p4 = None
    if verdict_a and verdict_b:
        at_primary = verdict_a == "survives" and verdict_b == "survives"
        rows_012 = [r for r in rows if r["caliber"].endswith("MMD=0.12")]
        at_012 = all(r["verdict"].startswith("undecidable") for r in rows_012) and rows_012
        p4 = bool(at_primary and at_012)
    add("P4", "verdict `survives` on both at MMD 1.0 and `undecidable` on both at MMD 0.12",
        f"A {verdict_a}, B {verdict_b}", p4)
    p5 = None
    if b_note.get("baseline_point") is not None and a_sub.get("baseline_point") is not None:
        same_sign = np.sign(b_note["baseline_point"]) == np.sign(a_sub["baseline_point"])
        gap = abs(b_note["baseline_point"] - a_sub["baseline_point"])
        p5 = bool(same_sign and gap < 0.5)
    add("P5", "B and A-on-the-same-subset share a sign and differ by less than 0.5",
        f"B {b_note.get('baseline_point')}, A|subset {a_sub.get('baseline_point')}", p5)
    ms = agreement["median_shift_B_minus_A"]
    add("P6", "median per-cysteine shift (B minus A) in [-1.5, 0.0]", f"{ms}",
        -1.5 <= ms <= 0.0)
    add("P7", f"Spearman rho between A and B exceeds {AGREEMENT_SPEARMAN_FLOOR}",
        f"{agreement['spearman_rho']}", agreement["spearman_rho"] > AGREEMENT_SPEARMAN_FLOOR)
    add("P8", "probability assigned to each registered outcome", f"outcome taken: {branch}", None)
    matched_rows = [r for r in rows if r["caliber"].startswith("matched")]
    p9 = None
    if matched_rows:
        checks = []
        for r in matched_rows:
            covers = (r["ci_low"] != "" and float(r["ci_low"]) <= 0 <= float(r["ci_high"]))
            small = r["point_estimate"] != "" and abs(float(r["point_estimate"])) < 0.5
            checks.append(covers and small)
        p9 = bool(all(checks))
    add("P9", "matched |difference| below 0.5 with an interval covering zero",
        "; ".join(f"{r['instrument']}: {r['point_estimate']} [{r['ci_low']}, {r['ci_high']}]"
                  for r in matched_rows), p9)

    for path, data in ((os.path.join(RESULTS, f"pka_round_verdict{suffix}.csv"), rows),
                       (os.path.join(RESULTS, f"pka_round_verdict_strata{suffix}.csv"),
                        strata_rows),
                       (os.path.join(RESULTS, f"pka_blind_prediction_scorecard{suffix}.csv"),
                        score)):
        cols = list(data[0].keys())
        with open(path + ".tmp", "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=cols)
            writer.writeheader()
            for row in data:
                writer.writerow(row)
        os.replace(path + ".tmp", path)

    audit = {
        "script": "scripts/run_pka_round_verdict.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "execution_scope": "SMOKE_NON_AUTHORITATIVE" if smoke else "FORMAL",
        "is_a_verdict": not smoke,
        "smoke_note": ("labels were permuted within protein with seed 20260915; no number here is a "
                       "result and none may be quoted" if smoke else ""),
        "registration": "protocols/pka_round_preregistration_2026-09-17.json, amendments A1-A3",
        "blind_prediction_lock": {"path": os.path.relpath(LOCK, ROOT), "sha256": sha256_of(LOCK)},
        "interpreter": {"python": sys.version, "numpy": np.__version__},
        "seed": SEED, "bootstrap_replicates": N_BOOT,
        "mmd": {"primary": MMD_PRIMARY, "calibers": list(MMD_CALIBERS),
                "powered_rule": "baseline 95% bootstrap half-width < MMD"},
        "inputs": {"results/pka_static_propka.csv": sha256_of(STATIC),
                   "results/pka_ensemble_propka.csv": sha256_of(ENSEMBLE)},
        "coverage": coverage,
        "per_instrument": notes,
        "agreement_A_versus_B": agreement,
        "outcome_branch": branch,
        "blind_prediction_scorecard": score,
        "n_correct": sum(1 for s in score if s["correct"] == "True"),
        "n_wrong": sum(1 for s in score if s["correct"] == "False"),
        "n_not_scorable": sum(1 for s in score if s["correct"] == ""),
        "declared_limits": [
            "both instruments are PROPKA, which the registration declares unreliable for cysteine "
            "in absolute terms; they share that weakness and agreement between them is not evidence "
            "of accuracy",
            "a 1 ns ensemble is not converged sampling; instrument B tests whether A's reading moves "
            "under conformational change, nothing more",
            "the matched caliber's propensity model is fitted once and held fixed across bootstrap "
            "replicates, so the intervals cover the contrast given the fit, not the fit",
            "clustering is on protein, so paralogy within the draw is not absorbed",
        ],
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "corrections": [],
    }
    path = os.path.join(RESULTS, f"pka_round_verdict_audit{suffix}.json")
    with open(path + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)
    os.replace(path + ".tmp", path)

    print(("SMOKE (non-authoritative) " if smoke else "") + f"outcome: {branch}")
    for tag, note in notes.items():
        print(f"  {tag}: baseline {note.get('baseline_point')} "
              f"{note.get('baseline_ci')} half-width {note.get('baseline_half_width')} "
              f"-> {note.get('primary_verdict')}")
    print(f"  agreement: rho {agreement['spearman_rho']}, median shift "
          f"{agreement['median_shift_B_minus_A']}, n {agreement['n_common_cysteines']}")
    print(f"  blind predictions: {audit['n_correct']} right, {audit['n_wrong']} wrong, "
          f"{audit['n_not_scorable']} not scorable")
    print(f"done in {audit['elapsed_minutes']} min")


if __name__ == "__main__":
    main()
