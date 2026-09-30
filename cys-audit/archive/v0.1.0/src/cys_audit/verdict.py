"""Verdict engine: dataset summary and the optional claim retest.

claim_verdict() is copied verbatim (logic, branch order, reason strings) from
scripts/ptm_detectability_diagnostics.py in the analysis tree, where it was pre-declared for the
published-claim retests and regression-tested on all 13 return paths. The Phase 4 validation script
checks that this copy returns the same verdict and reason as the original on every branch; a vendored
copy is used only so that the tool installs without the analysis tree.

Claim retest (optional, --claim-feature COLUMN): the claim is that a binary site feature (the column,
1 = feature present) is enriched among positives (claim_direction positive) or depleted (negative) or
absent (null). Effects are log2 odds ratios, positives vs background.
  baseline    uncontrolled log2 OR, clustered bootstrap 95%
  controlled  Mantel-Haenszel log2 OR within detectability strata built from every covariate the input
              supports: cleavage-band flags (needs --fasta), abundance tertile (needs abundance),
              multi-Cys-only flag (needs n_cys_peptide; undetected rows form their own level).
              Strata holding only one class are dropped; if they hold more than ORPHANED_LIMIT of the
              positives the controlled estimate is not trusted (precheck BLOCKING -> undecidable).
  random      size-matched random control: uncontrolled log2 OR on a random subset with as many
              positives and background rows as the informative strata retain (fixed seed); it
              separates "the control removed the effect" from "the control removed the power".
Verdict map: survives -> PASS; attenuated -> ATTENUATED; undecidable / out_of_instrument_scope ->
UNDECIDABLE; vanishes / reverses / baseline_contradicts_claim -> FAIL.
"""
from __future__ import annotations

import numpy as np

from . import stats
from .constants import (ATTENUATION_FLOOR, BOOTSTRAP_REPS, CI_LEVEL, CLAIM_ATTENUATED, CLAIM_FAIL,
                        CLAIM_NOT_EVALUATED, CLAIM_PASS, CLAIM_UNDECIDABLE, CLEAVAGE_BANDS, MIN_BACKGROUND,
                        MIN_POSITIVES, ORPHANED_LIMIT, SEVERITY)
from .io import background_mask


def claim_verdict(baseline, baseline_ci, controlled_ci, random_control_ci,
                  controlled_point, claim_direction="positive_preference",
                  baseline_reproduced=True, precheck_blocking=False):
    """The rule predeclared in phase 2, with the phase-2b direction guard. Pure function."""
    def excludes_zero(ci):
        return ci is not None and (ci[0] > 0 or ci[1] < 0)
    if precheck_blocking:
        return {"verdict": "undecidable", "reason": "a precheck is BLOCKING on this caliber"}
    if baseline is None or baseline_ci is None:
        return {"verdict": "out_of_instrument_scope", "reason": "no baseline estimate"}
    if claim_direction == "null_no_preference":
        if excludes_zero(baseline_ci):
            return {"verdict": "undecidable", "reason": "the author's null does not reproduce here"}
        if excludes_zero(controlled_ci):
            return {"verdict": "reverses", "reason": "null broken by the control"}
        return {"verdict": "survives", "reason": "null survives; quote with the null resolution"}
    if excludes_zero(baseline_ci) and claim_direction == "positive_preference" and baseline < 0:
        return {"verdict": "baseline_contradicts_claim",
                "reason": "baseline effect is opposite to the claim; a surviving effect is not the claim"}
    if not baseline_reproduced and not excludes_zero(baseline_ci):
        return {"verdict": "undecidable", "reason": "baseline not reproducible and not established here"}
    if not excludes_zero(baseline_ci):
        return {"verdict": "undecidable", "reason": "baseline interval crosses zero in our hands"}
    if excludes_zero(controlled_ci):
        if controlled_point is not None and np.sign(controlled_point) != np.sign(baseline):
            return {"verdict": "reverses", "reason": "sign flipped with the interval clear of zero"}
        if controlled_point is not None and abs(controlled_point) < ATTENUATION_FLOOR * abs(baseline):
            return {"verdict": "attenuated", "reason": "sign kept, interval clear of zero, under half the baseline"}
        return {"verdict": "survives", "reason": "sign kept, interval clear of zero, at least half the baseline"}
    if excludes_zero(random_control_ci):
        return {"verdict": "vanishes", "reason": "controlled interval crosses zero while the size-matched random control does not"}
    return {"verdict": "undecidable", "reason": "both the controlled and the size-matched random control cross zero: power-limited"}


VERDICT_TO_CLAIM = {"survives": CLAIM_PASS, "attenuated": CLAIM_ATTENUATED, "undecidable": CLAIM_UNDECIDABLE,
                    "out_of_instrument_scope": CLAIM_UNDECIDABLE, "vanishes": CLAIM_FAIL, "reverses": CLAIM_FAIL,
                    "baseline_contradicts_claim": CLAIM_FAIL}


def dataset_summary(tests):
    """Worst status across tests plus counts; describes the dataset, not a claim."""
    counts = {}
    for t in tests.values():
        counts[t["status"]] = counts.get(t["status"], 0) + 1
    worst = max(tests.values(), key=lambda t: SEVERITY[t["status"]])["status"] if tests else None
    return {"worst_test_status": worst, "status_counts": counts}


def _strata(ds, cfg):
    parts, used = [], []
    if ds.fasta is not None:
        from .checks.cleavage import band_flags
        flags, has_seq = band_flags(ds, cfg["protease"])
        for b in CLEAVAGE_BANDS:
            parts.append(np.where(has_seq, flags[b].astype(int), 2))
            used.append(f"cleavage band {b}")
    if ds.abundance is not None:
        a = ds.abundance
        fin = np.isfinite(a)
        lvl = np.full(len(ds), 3)
        if fin.sum() >= 3:
            q1, q2 = np.quantile(a[fin], [1 / 3, 2 / 3])
            lvl[fin] = np.where(a[fin] <= q1, 0, np.where(a[fin] <= q2, 1, 2))
        parts.append(lvl)
        used.append("abundance tertile (missing = own level)")
    if ds.n_cys_peptide is not None:
        n = ds.n_cys_peptide
        parts.append(np.where(np.isfinite(n), (n >= 2).astype(int), 2))
        used.append("multi-Cys-only (no peptide = own level)")
    if not parts:
        return None, used
    key = np.zeros(len(ds), dtype=np.int64)
    for p in parts:
        key = key * 4 + p
    return key, used


def claim_retest(ds, cfg):
    feat = cfg.get("claim_feature")
    if not feat:
        return {"status": CLAIM_NOT_EVALUATED, "reason": "no --claim-feature given; only dataset-level tests were run"}
    if feat not in ds.extra:
        return {"status": CLAIM_UNDECIDABLE, "reason": f"claim feature column {feat!r} not found or not numeric"}
    direction = cfg.get("claim_direction", "positive")
    rule_dir = {"positive": "positive_preference", "negative": "positive_preference",
                "null": "null_no_preference"}[direction]
    sign = -1.0 if direction == "negative" else 1.0
    bg = background_mask(ds, cfg["background"])
    if bg is None:
        return {"status": CLAIM_UNDECIDABLE, "reason": "background 'observed' needs a 'detected' column"}
    x = ds.extra[feat]
    ok = np.isfinite(x) & np.isin(x, (0.0, 1.0))
    pos = (ds.label == 1) & ok
    bgm = bg & ok
    if pos.sum() < MIN_POSITIVES or bgm.sum() < MIN_BACKGROUND:
        return {"status": CLAIM_UNDECIDABLE, "reason": "too few positives or background rows carry the feature"}
    keep = pos | bgm
    k_idx = np.flatnonzero(keep)
    codes, labels = stats.cluster_index([ds.cluster[i] for i in k_idx])
    K = len(labels)
    reps = cfg.get("reps", BOOTSTRAP_REPS)
    seed = cfg["seed"]
    flag = x[k_idx] == 1
    ispos = pos[k_idx]
    base_counts = stats.two_by_two_counts(codes, K, ispos, flag)
    b_pt, b_ci, _, _ = stats.boot_log2_or(base_counts, reps, seed + 501, CI_LEVEL)

    key, used = _strata(ds, cfg)
    details = {"feature": feat, "direction": direction, "background": cfg["background"], "strata_from": used}
    if key is None:
        v = claim_verdict(sign * b_pt, sorted((sign * b_ci[0], sign * b_ci[1])), None, None, None, rule_dir,
                          precheck_blocking=True)
        details.update(baseline={"estimate": b_pt, "ci": list(b_ci)})
        return {"status": CLAIM_UNDECIDABLE, "verdict": v["verdict"],
                "reason": "no detectability covariate available to build a control", "details": details}
    sk = key[k_idx]
    levels = np.unique(sk)
    s_index = np.searchsorted(levels, sk)
    S = len(levels)
    # per-stratum class counts decide informativeness
    npos_s = np.bincount(s_index, weights=ispos.astype(float), minlength=S)
    nbg_s = np.bincount(s_index, weights=(~ispos).astype(float), minlength=S)
    informative = (npos_s > 0) & (nbg_s > 0)
    orphaned = float(npos_s[~informative].sum() / npos_s.sum())
    blocking = orphaned > ORPHANED_LIMIT
    inf_rows = informative[s_index]
    # controlled: MH over informative strata; bootstrap over clusters with per-cluster-per-stratum counts
    cell = np.zeros((K, S * 4), dtype=np.int64)
    for j, sel in enumerate((ispos & flag, ispos & ~flag, ~ispos & flag, ~ispos & ~flag)):
        m = sel & inf_rows
        np.add.at(cell, (codes[m], s_index[m] * 4 + j), 1)
    c_pt = stats.mh_log2_or(cell.sum(axis=0).reshape(S, 4)[informative])
    vals = []
    for w in stats.multiplicities(K, reps, seed + 502):
        tab = (w @ cell).reshape(w.shape[0], S, 4)
        for r in range(w.shape[0]):
            t = tab[r][informative]
            t = t[t.sum(axis=1) > 0]
            has_both = (t[:, 0] + t[:, 1] > 0) & (t[:, 2] + t[:, 3] > 0)
            vals.append(stats.mh_log2_or(t[has_both]) if has_both.any() else np.nan)
    c_ci = stats.percentile_interval(vals, CI_LEVEL)
    # size-matched random control
    rng = np.random.default_rng(seed + 503)
    n_pos_keep = int(npos_s[informative].sum())
    n_bg_keep = int(nbg_s[informative].sum())
    pos_rows = np.flatnonzero(ispos)
    bg_rows = np.flatnonzero(~ispos)
    pick = np.concatenate([rng.choice(pos_rows, n_pos_keep, replace=False), rng.choice(bg_rows, n_bg_keep, replace=False)])
    sub = np.zeros(len(k_idx), dtype=bool)
    sub[pick] = True
    r_counts = stats.two_by_two_counts(codes, K, ispos, flag, mask=sub)
    r_pt, r_ci, _, _ = stats.boot_log2_or(r_counts, reps, seed + 504, CI_LEVEL)

    def flip(ci):
        return sorted((sign * ci[0], sign * ci[1]))
    v = claim_verdict(sign * b_pt, flip(b_ci), flip(c_ci), flip(r_ci), sign * c_pt, rule_dir,
                      precheck_blocking=blocking)
    details.update(baseline={"estimate": float(b_pt), "ci": [float(b_ci[0]), float(b_ci[1])]},
                   controlled={"estimate": float(c_pt), "ci": [float(c_ci[0]), float(c_ci[1])],
                               "method": "Mantel-Haenszel within informative strata"},
                   size_matched_random={"estimate": float(r_pt), "ci": [float(r_ci[0]), float(r_ci[1])]},
                   strata=int(S), informative_strata=int(informative.sum()), orphaned_positive_fraction=orphaned,
                   precheck_blocking=bool(blocking), sign_convention="effects multiplied by -1 for a negative claim"
                   if sign < 0 else "effects as estimated")
    return {"status": VERDICT_TO_CLAIM[v["verdict"]], "verdict": v["verdict"], "reason": v["reason"], "details": details}
