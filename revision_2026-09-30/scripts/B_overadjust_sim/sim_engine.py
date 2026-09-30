# -*- coding: utf-8 -*-
"""B_overadjust_sim simulation engine (POST HOC revision analysis; not registered).

One replicate = one simulated site table on real mouse sequences with real detection
(PXD063463 global arm), analysed exactly as the published claim re-test does. Every statistical
step of the re-test is IMPORTED from the read-only repository script
`repo/scripts/run_phase2_claims_under_detectability_control.py` (module alias `rt`):

  rt.log_odds_ratio, rt.weighted_table        Haldane-Anscombe (0.5) log2 odds ratio
  rt.make_weighted_table                      weighted 2x2 closure for the cluster bootstrap
  rt.propensity_score                         LogisticRegression(C=1, lbfgs, max_iter=1000) on
                                              standardised VIS10, clipped score and logit
  rt.nearest_neighbour_match                  1:1 nearest neighbour WITHOUT replacement on the
                                              propensity logit, caliper rt.CALIPER_SD x SD(logit)
                                              = 0.2 SD, positives processed from high to low logit
  rt.quantile_strata, rt.make_mh_log_odds_ratio  secondary specification: propensity quintiles,
                                              Mantel-Haenszel log2 OR
  rt.bootstrap_interval                       protein-clustered bootstrap (multinomial cluster
                                              weights), percentile 2.5/97.5, seed rt.SEED
  rt.random_control_effect                    same-size random control: the matched positives plus
                                              as many negatives drawn at random from the whole
                                              background, 20 draws (seeds rt.SEED+i); interval
                                              integrates draw choice and cluster resampling
  rt.classify, rt.DOWNGRADE                   verdict rule and the specification-disagreement
                                              downgrade (survives->attenuated, vanishes->undecidable)
  audit_cleaning_and_grouping.cleaning_collinearity   the pre-check that can BLOCK matching
                                              (called with groups=None, as the pipeline does
                                              for > 2000 groups; blocking never depends on groups)

Documented departures, none of which changes an estimator of the pipeline:
  * random_control_point() is a line-for-line copy of the POINT part of rt.random_control_effect
    (its 20 seeded draws and the mean of their log2 ORs); the original function is called
    unchanged whenever an interval is needed.
  * bootstrap replicates for intervals are N_BOOT (default 1000) instead of 5000.
  * the phase-2b direction guard (`baseline_contradicts_claim`) is re-stated from
    run_phase2b_claims_backfill.run_one (lines 351-379) because it is inline code there.

EXPLORATORY estimators added here (post hoc, NOT part of the published instrument; point
estimates only), used to locate the mechanism of any attenuation:
  pair_conditional     conditional (McNemar-type) log2 OR over the pipeline's own matched pairs,
                       log2((n10 + 0.5) / (n01 + 0.5)); removes the bias of analysing a sample
                       matched on an exposure-correlated factor with an unmatched (pooled) 2x2
  detmatch_crude       same matching algorithm, but on the logit of a DETECTION model
                       (rt.propensity_score applied to detected ~ VIS10 over the whole universe)
                       instead of a model of the claim label; pooled 2x2 of the matched sample
  detmatch_pair        pair-conditional log2 OR over the detection-model matched pairs
  neg_b_oracle         positives vs detected non-positives of the same proteins (an
                       observed-negative background, available in the simulation by construction)

The released tool's claim retest (cys_audit v0.2.2, verdict.claim_retest with
--claim-covariates = VIS10, background 'proteome') is run on the same table as a cross-check.

REVISION AFTER VERIFICATION (post hoc, 2026-09-30), backward compatible (every stored table is
reproduced bit for bit; 07_rerun_check.py):
  * two further selection models for positives, read from results/.../universe_extra.csv.gz
    (01d_extra_universe.py): 'steep_vis' (selection on the logit of the VIS10 detection model, a
    score VIS10 represents exactly) and 'steep_nocomp' (selection on the step-1 peptide detection
    model refitted without the peptide basic/acidic fractions);
  * one exploratory Artifact-1 analogue attribute, a1c_KRH_p5p8 (K/R/H at +5..+8), analysed only when
    a task asks for it (analyse = 'extra_only' or 'with_extra'); attribute code 7;
  * new scenario codes (16-23) for the new groups; codes 0-15 are unchanged.

REVISION ROUND 2 (post hoc, 2026-09-30), backward compatible (07_rerun_check.py):
  * PROXY-CHEMISTRY model (the verifier's round-2 alternative): genuine chemistry planted on a CONTINUOUS
    per-cysteine feature instead of on the claimed attribute, P(positive | detected) = expit(a + b ln2 F),
    F standardised over the universe, b = log2 odds ratio per SD of F. Main feature F_KRcount20_z, the local
    K/R count within +/-20 residues (repository column pep_cleavage_sites_within_20, which VIS10 omits);
    sensitivity features F_KR5_z (K/R within +/-5) and F_AKRV10_z (A/K/R/V within +/-10). Read from
    results/.../universe_krcount.csv.gz (01f_krcount_universe.py). draw_positives() is unchanged: it already
    accepts a continuous planted vector.
  * two HYPOTHETICAL count-defined claims analysed like the published attributes (binary): x1_KRcount20_hi
    (K/R count within +/-20 above the universe median) and x2_KR5_ge3 (>= 3 of +/-5 in {K,R}).
  * attribute/feature codes 8-12 and scenario codes 24-30 for the new groups; every earlier code is unchanged,
    and analysing further attributes never consumes the replicate's random stream, so re-generated tables are
    identical to the stored ones.
  * a task may name the analysed attributes explicitly (task['attrs']).

REVISION ROUND 3 (post hoc, 2026-09-30), backward compatible (07_rerun_check.py):
  * two further HYPOTHETICAL count-defined attributes used only as discriminators, x3_AKRV10_hi (A/K/R/V count
    within +/-10 above its universe median) and x4_KR5_hi (K/R count within +/-5 above its median), read from
    results/.../universe_discrim.csv.gz (01g_discrim_universe.py); attribute codes 13-14. Loading them draws nothing.
  * MIXTURE chemistry: a task may carry a second planted term (task['attribute2'], task['effect2'], log2 OR per unit
    of that attribute or feature), added to the logit among detected cysteines. Only tasks that carry it add
    (attribute2 code, round(100 * effect2)) to their seed; every other task keeps its seed and its table.
  * scenario codes 31-43 for the round-3 groups.

REVISION ROUND 4 (post hoc, 2026-09-30), backward compatible (07_rerun_check.py):
  * a DETECTABILITY-ONLY selection model on tryptic-peptide length, re-implemented from the round-4 verifier's
    specification (verify_r4_02_engine.py, mode 'lensel7'; verify_r4_08_lensel_calibrated.py): among detected
    cysteines, P(positive) = expit(a + gamma * s) with s = -ln(max(pep_len, 7)), pep_len being the VIS10 column
    'pep_len' (length of the cysteine's tryptic peptide), so gamma > 0 favours cysteines on short peptides. No
    chemistry term. Detection mode 'lensel7'; it enters draw_positives() through the same gamma * score path as the
    steep models.
  * only 'lensel7' tasks add (LENSEL_SEED_TAG, round(1000 * gamma)) to their seed, so cells at different gamma are
    independent tables; every other task keeps its seed and its table.
  * scenario codes 44-47 for the round-4 groups.
"""
from __future__ import annotations

import os
import sys
import warnings

sys.dont_write_bytecode = True
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import brentq  # noqa: E402
from scipy.special import expit  # noqa: E402

import run_phase2_claims_under_detectability_control as rt  # noqa: E402  (repo, read-only)
from audit_cleaning_and_grouping import cleaning_collinearity  # noqa: E402
from cys_audit.io import Dataset  # noqa: E402  (released tool v0.2.2)
from cys_audit.verdict import claim_retest as tool_claim_retest  # noqa: E402
from cys_audit import constants as tool_constants  # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)

VIS10_NAMES = ["pep_len", "pep_mass", "pep_gravy", "pep_detectable_length", "pep_detectable_mass",
               "pep_detectable_both", "pep_log_len", "pep_mc1_detectable", "pep_mc2_detectable",
               "pep_detectable_any_missed_cleavage"]
SEQ_ATTRIBUTES = ["a1_SFE006_KR", "a1b_SNO021_K", "a2_SNO016_KRH", "a2b_PERS009_AKRV",
                  "a3_SNO006_DE3", "a3b_SFE002_E"]
LATENT = "a4_LATENT30"
LATENT_PREVALENCE = 0.30
ALL_ATTRIBUTES = SEQ_ATTRIBUTES + [LATENT]
EXTRA_ATTRIBUTES = ["a1c_KRH_p5p8"]          # revision: exploratory Artifact-1 analogue
ATTR_CODES = {a: i for i, a in enumerate(ALL_ATTRIBUTES + EXTRA_ATTRIBUTES)}
# revision round 2: hypothetical count-defined claims (binary) and proxy-chemistry features (continuous)
COUNT_ATTRIBUTES = ["x1_KRcount20_hi", "x2_KR5_ge3"]
CHEM_FEATURES = ["F_KRcount20_z", "F_KR5_z", "F_AKRV10_z"]
ATTR_CODES.update({"x1_KRcount20_hi": 8, "x2_KR5_ge3": 9, "F_KRcount20_z": 10, "F_KR5_z": 11, "F_AKRV10_z": 12})
# revision round 3: discriminator attributes (count above its universe median)
DISCRIM_ATTRIBUTES = ["x3_AKRV10_hi", "x4_KR5_hi"]
ATTR_CODES.update({"x3_AKRV10_hi": 13, "x4_KR5_hi": 14})
STEEP_MODES = {"steep": "steep_z", "steep_vis": "steep_z_vis", "steep_nocomp": "steep_z_nocomp",
               "lensel7": "lensel7_s"}           # revision round 4: selection on tryptic-peptide length
LENSEL_SEED_TAG = 15                              # revision round 4: seed tag of 'lensel7' tasks only
N_BOOT = 1000

_U = None


def universe():
    """Load the universe table once per process; fit the detection-model logit once."""
    global _U
    if _U is None:
        df = pd.read_csv(cb.UNIVERSE, dtype={"protein": str})
        _U = {
            "protein": df["protein"].to_numpy(),
            "position": df["position"].to_numpy(),
            "detected": df["detected"].to_numpy().astype(int),
            "vis10": df[VIS10_NAMES].to_numpy(dtype=float),
            "attr": {a: df[a].to_numpy().astype(int) for a in SEQ_ATTRIBUTES},
        }
        z = df["m3_best_logit"].to_numpy(dtype=float)
        det = _U["detected"] == 1
        zs = np.full(len(z), np.nan)
        zs[det] = (z[det] - z[det].mean()) / z[det].std()
        _U["steep_z"] = zs
        # revision round 4: selection score of the length-selection model (no random draw; unstandardised, as in
        # the round-4 verifier's parametrisation, so gamma is directly comparable)
        pep_len = df["pep_len"].to_numpy(dtype=float)
        assert np.all(np.isfinite(pep_len[det])) and np.all(pep_len[det] > 0)
        _U["lensel7_s"] = -np.log(np.maximum(pep_len, 7.0))
        uniq, codes = np.unique(_U["protein"], return_inverse=True)
        _U["prot_codes"] = codes
        _U["n_prot"] = len(uniq)
        _U["det_logit"] = rt.propensity_score(_U["vis10"], _U["detected"])[1]
        td_path = os.path.join(cb.RESULTS, "universe_theor_detectable.csv.gz")
        if os.path.exists(td_path):
            td = pd.read_csv(td_path, dtype={"protein": str})
            assert (td["protein"].to_numpy() == _U["protein"]).all()
            assert (td["position"].to_numpy() == _U["position"]).all()
            _U["theor"] = td["theor_detectable_ms"].to_numpy().astype(int)
        ex_path = os.path.join(cb.RESULTS, "universe_extra.csv.gz")
        _U["extra_attr"] = {}
        if os.path.exists(ex_path):
            ex = pd.read_csv(ex_path, dtype={"protein": str})
            assert (ex["protein"].to_numpy() == _U["protein"]).all()
            assert (ex["position"].to_numpy() == _U["position"]).all()
            for a in EXTRA_ATTRIBUTES:
                _U["extra_attr"][a] = ex[a].to_numpy().astype(int)
            for col, key in (("vis_det_logit", "steep_z_vis"), ("m3nc_best_logit", "steep_z_nocomp")):
                v = ex[col].to_numpy(dtype=float)
                zz = np.full(len(v), np.nan)
                zz[det] = (v[det] - v[det].mean()) / v[det].std()
                _U[key] = zz
        # revision round 2: K/R-count columns (01f_krcount_universe.py)
        kc_path = os.path.join(cb.RESULTS, "universe_krcount.csv.gz")
        _U["count_attr"] = {}
        _U["chem_feat"] = {}
        if os.path.exists(kc_path):
            kc = pd.read_csv(kc_path, dtype={"protein": str})
            assert (kc["protein"].to_numpy() == _U["protein"]).all()
            assert (kc["position"].to_numpy() == _U["position"]).all()
            for a in COUNT_ATTRIBUTES:
                _U["count_attr"][a] = kc[a].to_numpy().astype(int)
            for f in CHEM_FEATURES:
                _U["chem_feat"][f] = kc[f].to_numpy(dtype=float)
        # revision round 3: discriminator attributes (01g_discrim_universe.py)
        dc_path = os.path.join(cb.RESULTS, "universe_discrim.csv.gz")
        if os.path.exists(dc_path):
            dc = pd.read_csv(dc_path, dtype={"protein": str})
            assert (dc["protein"].to_numpy() == _U["protein"]).all()
            assert (dc["position"].to_numpy() == _U["position"]).all()
            for a in DISCRIM_ATTRIBUTES:
                _U["count_attr"][a] = dc[a].to_numpy().astype(int)
    return _U


def seed_for(*parts):
    return np.random.SeedSequence([cb.MASTER_SEED] + [int(p) for p in parts])


def latent_attribute(rng, n):
    """Sequence-independent binary attribute on a random 30% of cysteines (exact count)."""
    a = np.zeros(n, dtype=int)
    a[rng.choice(n, size=int(round(LATENT_PREVALENCE * n)), replace=False)] = 1
    return a


def intercept_for(attr_det, b, rate):
    """Intercept a with mean(expit(a + b*A)) == rate over the detected cysteines."""
    return brentq(lambda a: float(np.mean(expit(a + b * attr_det))) - rate, -30.0, 30.0)


def draw_positives(detected, scenario, attr_all, effect_log2, rate, rng, steep_z=None, gamma=0.0, extra_all=None):
    """Positives among detected cysteines (the only cysteines that can be reported).
    S0: P(positive | detected) = rate. S1: P = expit(a + b*A), b = effect_log2 * ln 2, a set so
    that the mean over detected cysteines equals `rate`. Equivalent to modification M ~ expit(a+bA)
    on every cysteine with positive = M AND detected, because M depends on A only.
    'steep' sensitivity mode: P = expit(a + b*A + gamma*z), z = standardised peptide-level
    detectability of the cysteine (m3_best_logit), so positives are a more detectability-selected
    subset of the detected cysteines; the chemical truth is still b (conditional on z).
    Revision round 3: `extra_all` (log-odds scale, one value per cysteine) is a second planted term (mixture
    chemistry); it is None for every earlier task, which therefore follows the unchanged code path."""
    det_idx = np.flatnonzero(detected == 1)
    b = 0.0 if (scenario.startswith("S0") or effect_log2 == 0.0) else effect_log2 * np.log(2.0)
    lin = np.zeros(len(det_idx))
    if b != 0.0:
        lin = lin + b * attr_all[det_idx]
    if extra_all is not None:
        lin = lin + extra_all[det_idx]
    if gamma:
        lin = lin + gamma * steep_z[det_idx]
    if b == 0.0 and not gamma and extra_all is None:
        p = np.full(len(det_idx), rate)
    else:
        a0 = brentq(lambda a: float(np.mean(expit(a + lin))) - rate, -40.0, 40.0)
        p = expit(a0 + lin)
    y_all = np.zeros(len(detected), dtype=int)
    y_all[det_idx[rng.random(len(det_idx)) < p]] = 1
    return y_all


def cohort_index(prot_codes, y_all):
    """NEG_A cohort: every cysteine of every protein carrying >= 1 positive."""
    pos_prot = np.unique(prot_codes[y_all == 1])
    return np.flatnonzero(np.isin(prot_codes, pos_prot))


def random_control_point(y, attribute, positive_index, n_negatives,
                         seed=rt.SEED, draws=rt.RANDOM_DRAWS):
    """Line-for-line copy of the point part of rt.random_control_effect (no interval)."""
    negatives = np.flatnonzero(y == 0)
    if n_negatives > len(negatives):
        n_negatives = len(negatives)
    point = []
    for i in range(draws):
        rng = np.random.default_rng(seed + i)
        picked = rng.choice(negatives, size=n_negatives, replace=False)
        index = np.concatenate([positive_index, picked])
        n11, n10, n01, n00 = rt.weighted_table(y[index], attribute[index], np.ones(len(index)))
        point.append(rt.log_odds_ratio(n11, n10, n01, n00))
    return float(np.mean(point))


def pair_conditional(a, pos_idx, neg_idx):
    """Conditional log2 OR for 1:1 matched pairs, Haldane 0.5 on the discordant counts."""
    if len(pos_idx) == 0:
        return np.nan
    ap, an = a[pos_idx], a[neg_idx]
    n10 = float(((ap == 1) & (an == 0)).sum())
    n01 = float(((ap == 0) & (an == 1)).sum())
    return float(np.log2((n10 + 0.5) / (n01 + 0.5)))


def within_protein_z(y, a, prot):
    """Artifact-1 statistic: observed attribute-carrying positives against the within-protein
    hypergeometric expectation, normal approximation (Experimental Procedures, binary features)."""
    uniq, inv = np.unique(prot, return_inverse=True)
    n = np.bincount(inv).astype(float)
    m = np.bincount(inv, weights=y.astype(float))
    k = np.bincount(inv, weights=a.astype(float))
    o = float(((y == 1) & (a == 1)).sum())
    e = float(np.sum(m * k / n))
    ok = n > 1
    var = float(np.sum((m * (k / n) * (1 - k / n) * (n - m) / np.where(ok, n - 1, 1.0))[ok]))
    return (o - e) / np.sqrt(var) if var > 0 else np.nan


def artifact1_estimators(y, a, prot, theor):
    """Within-protein z and protein-stratified Mantel-Haenszel log2 OR (rt.make_mh_log_odds_ratio
    with strata = protein), against all cysteines and against positives + theoretically detectable
    cysteines only (the Artifact-1 background restriction)."""
    codes = np.unique(prot, return_inverse=True)[1]
    keep = (y == 1) | (theor == 1)
    out = {"a1_wp_z_all": within_protein_z(y, a, prot),
           "a1_wp_z_restricted": within_protein_z(y[keep], a[keep], prot[keep]),
           "a1_wp_mh_all": rt.make_mh_log_odds_ratio(y, a, codes)(np.ones(len(y))),
           "a1_wp_mh_restricted": rt.make_mh_log_odds_ratio(y[keep], a[keep], codes[keep])(np.ones(int(keep.sum()))),
           "a1_restricted_crude": rt.log_odds_ratio(*rt.weighted_table(y[keep], a[keep], np.ones(int(keep.sum())))),
           "a1_restricted_n_neg": int(((y == 0) & keep).sum())}
    return out


def direction_guard(verdict, baseline, baseline_ci):
    """phase-2b supplement 3 (run_phase2b_claims_backfill.run_one): a positive_preference claim whose
    baseline is negative with an interval excluding zero is recorded as baseline_contradicts_claim."""
    if baseline_ci is None or not np.all(np.isfinite(baseline_ci)):
        return verdict
    crosses = baseline_ci[0] <= 0 <= baseline_ci[1]
    if (not crosses) and baseline < 0:
        return "baseline_contradicts_claim"
    return verdict


def analyse(y, attrs, cov, groups, detected, det_logit, want_ci, n_boot=N_BOOT, run_tool=False,
            tool_proteins=None, tool_positions=None, theor=None):
    """Run the re-test on one table for every attribute in `attrs` (dict name -> 0/1 vector).

    The propensity model, matching and strata depend only on y and VIS10, so they are fitted once
    and shared by all attributes, exactly as the pipeline shares them between a claim's primary and
    secondary attributes."""
    out = {}
    n = len(y)
    ones = np.ones(n)
    score, logit = rt.propensity_score(cov, y)
    strata = rt.quantile_strata(score)
    caliper = rt.CALIPER_SD * float(np.std(logit))
    pos_idx, neg_idx = rt.nearest_neighbour_match(logit, y, caliper)
    matched_index = np.concatenate([pos_idx, neg_idx])
    keep = np.zeros(n, dtype=bool)
    keep[matched_index] = True
    col = cleaning_collinearity(y, keep, covariate=score, strata=strata, groups=None, label="sim")
    blocked = bool(col["blocking_reasons"])
    # exploratory: matching on a detection model instead of a model of the claim label
    d_caliper = rt.CALIPER_SD * float(np.std(det_logit))
    dpos_idx, dneg_idx = rt.nearest_neighbour_match(det_logit, y, d_caliper)
    d_index = np.concatenate([dpos_idx, dneg_idx])
    common = {
        "n_obs": int(n), "n_pos": int(y.sum()), "n_proteins": int(len(np.unique(groups))),
        "pos_rate": float(y.mean()), "n_pairs": int(len(pos_idx)),
        "pos_retention": float(len(pos_idx) / max(1, int(y.sum()))),
        "prop_auc_before": col["covariate"]["auc_before"], "prop_auc_after": col["covariate"]["auc_after"],
        "precheck_blocked": blocked,
        "matched_neg_detected_frac": float(detected[neg_idx].mean()) if len(neg_idx) else np.nan,
        "all_neg_detected_frac": float(detected[y == 0].mean()),
        "detmatch_n_pairs": int(len(dpos_idx)),
    }
    m_y = y[matched_index]
    m_groups = groups[matched_index]
    det_mask = detected == 1
    d_y = y[d_index]
    for name, a in attrs.items():
        rec = dict(common)
        rec["attr_prev_pos"] = float(a[y == 1].mean())
        rec["attr_prev_neg"] = float(a[y == 0].mean())
        rec["attr_prev_matched_neg"] = float(a[neg_idx].mean()) if len(neg_idx) else np.nan
        table = rt.make_weighted_table(y, a)
        baseline = rt.log_odds_ratio(*table(ones))
        m_a = a[matched_index]
        m_table = rt.make_weighted_table(m_y, m_a)
        matched = rt.log_odds_ratio(*m_table(np.ones(len(m_y)))) if len(pos_idx) >= 10 else np.nan
        mh_fn = rt.make_mh_log_odds_ratio(y, a, strata)
        mh = mh_fn(ones)
        yb, ab = y[det_mask], a[det_mask]
        neg_b = rt.log_odds_ratio(*rt.weighted_table(yb, ab, np.ones(len(yb))))
        rec.update(baseline=baseline, matched=matched, stratified=mh, neg_b_oracle=neg_b,
                   pair_conditional=pair_conditional(a, pos_idx, neg_idx),
                   detmatch_crude=rt.log_odds_ratio(*rt.weighted_table(d_y, a[d_index], np.ones(len(d_y))))
                   if len(dpos_idx) >= 10 else np.nan,
                   detmatch_pair=pair_conditional(a, dpos_idx, dneg_idx))
        if theor is not None:
            rec.update(artifact1_estimators(y, a, groups, theor))
        if want_ci and len(pos_idx) >= 10:
            b_ci, _ = rt.bootstrap_interval(lambda w: rt.log_odds_ratio(*table(w)), groups,
                                            replicates=n_boot, seed=rt.SEED)
            m_ci, _ = rt.bootstrap_interval(lambda w: rt.log_odds_ratio(*m_table(w)), m_groups,
                                            replicates=n_boot, seed=rt.SEED)
            s_ci, _ = rt.bootstrap_interval(mh_fn, groups, replicates=n_boot, seed=rt.SEED)
            r_point, r_ci, _, _ = rt.random_control_effect(y, a, groups, pos_idx, len(neg_idx),
                                                           replicates=n_boot)
            primary, _ = rt.classify(baseline, matched, m_ci, r_ci, blocked)
            secondary, _ = rt.classify(baseline, mh, s_ci, r_ci, blocked)
            disagreement = primary != secondary
            final = rt.DOWNGRADE.get(primary, primary) if disagreement else primary
            rec.update(random_control=r_point,
                       baseline_lo=b_ci[0], baseline_hi=b_ci[1], matched_lo=m_ci[0], matched_hi=m_ci[1],
                       stratified_lo=s_ci[0], stratified_hi=s_ci[1], random_lo=r_ci[0], random_hi=r_ci[1],
                       verdict_primary=primary, verdict_secondary=secondary,
                       verdict_final=final,
                       verdict_final_guarded=direction_guard(final, baseline, b_ci),
                       verdict_primary_guarded=direction_guard(primary, baseline, b_ci))
            if run_tool:
                rec.update(_tool_retest(y, a, cov, groups, tool_proteins, tool_positions, n_boot))
        else:
            rec["random_control"] = random_control_point(y, a, pos_idx, len(neg_idx))
        out[name] = rec
    return out


def _tool_retest(y, a, cov, groups, proteins, positions, n_boot):
    """Released Cys-Audit v0.2.2 claim retest, --claim-covariates VIS10, background proteome."""
    extra = {nm: cov[:, j] for j, nm in enumerate(VIS10_NAMES)}
    extra["claim_attr"] = a.astype(float)
    ds = Dataset(protein=list(proteins), position=np.asarray(positions), label=np.asarray(y),
                 detected=None, n_cys_peptide=None, abundance=None, peptide_mass=None,
                 cluster=list(groups), extra=extra, fasta=None, notes=[])
    cfg = {"claim_feature": "claim_attr", "claim_direction": "positive", "background": "proteome",
           "claim_covariates": VIS10_NAMES, "seed": tool_constants.DEFAULT_SEED, "reps": n_boot}
    res = tool_claim_retest(ds, cfg)
    d = res.get("details", {})
    ctrl = d.get("controlled") or {}
    ci = ctrl.get("ci") or [np.nan, np.nan]
    return {"tool_status": res.get("status"), "tool_verdict": res.get("verdict"),
            "tool_controlled": ctrl.get("estimate", np.nan),
            "tool_controlled_lo": ci[0], "tool_controlled_hi": ci[1],
            "tool_random_kind": (d.get("random_control") or {}).get("kind")}


SCEN_CODES = {"S0": 0, "S1": 1, "S1grid": 2, "S0rate": 3, "S1rate": 4, "S0perm": 5, "S1perm": 6,
              "S0boot5000": 7, "S0steep": 8, "S1steep": 9, "S1steepgrid": 10, "calib": 11,
              "S0art1": 12, "S1art1": 13, "S0art1steep": 14, "S1art1steep": 15,
              # revision after verification (post hoc)
              "S1permgrid": 16, "S0steepvis": 17, "S1steepvis": 18, "S0steepnc": 19, "S1steepnc": 20,
              "calibalt": 21, "S0art1steep4": 22, "S1art1steep4": 23,
              # revision round 2 (post hoc): proxy-chemistry model (chemistry on a continuous K/R-count feature)
              "P1grid": 24, "P1steepgrid": 25, "P1permgrid": 26, "P1kr5grid": 27, "P1akrv10grid": 28,
              "P1ci": 29, "P1permci": 30,
              # revision round 3 (post hoc): sequence-only discriminator (count vs offsets association)
              "DM1real": 31, "DM1steep": 32, "DM1perm": 33, "DM2real": 34, "DM2steep": 35, "DM2perm": 36,
              "DM2kr5": 37, "DM2akrv10": 38, "DMIXgrid": 39, "DMIX": 40, "D0r16": 41, "DM1r16": 42, "DM2r16": 43,
              # revision round 4 (post hoc): detectability-only selection on tryptic-peptide length
              "L0grid": 44, "L0cal": 45, "L0cal16": 46, "L0ci": 47}


def run_replicate(task):
    """task: dict(scenario, attribute (planted attribute or None), effect, rate, rep, want_ci,
    run_tool, detection_mode, analyse, n_boot). Returns flat records, one per analysed attribute."""
    U = universe()
    scen = task["scenario"]
    attr_code = ATTR_CODES[task["attribute"]] if task["attribute"] else 99
    seed_scen = SCEN_CODES[task.get("seed_scenario", scen)]
    seed_parts = [seed_scen, attr_code, int(round(task["effect"] * 100)), int(round(task["rate"] * 1000)), task["rep"]]
    if task.get("attribute2"):                         # revision round 3: mixture tasks only
        seed_parts += [ATTR_CODES[task["attribute2"]], int(round(float(task["effect2"]) * 100))]
    if task.get("detection_mode") == "lensel7":        # revision round 4: length-selection tasks only
        seed_parts += [LENSEL_SEED_TAG, int(round(float(task["gamma"]) * 1000))]
    ss = seed_for(*seed_parts)
    rng = np.random.default_rng(ss)
    n_all = len(U["detected"])
    latent = latent_attribute(rng, n_all)
    attrs_all = dict(U["attr"])
    attrs_all[LATENT] = latent
    attrs_all.update(U.get("extra_attr", {}))
    attrs_all.update(U.get("count_attr", {}))          # revision round 2 (no random draw)
    detected = U["detected"]
    det_logit_all = U["det_logit"]
    mode = task.get("detection_mode", "real")
    if mode == "permuted_within_protein":
        detected = detected.copy()
        codes = U["prot_codes"]
        order = np.argsort(codes, kind="stable")
        bounds = np.flatnonzero(np.diff(codes[order])) + 1
        for block in np.split(order, bounds):
            detected[block] = detected[block][rng.permutation(len(block))]
        det_logit_all = rt.propensity_score(U["vis10"], detected)[1]
    planted = task["attribute"]
    gamma = float(task.get("gamma", 0.0)) if mode in STEEP_MODES else 0.0
    steep_z = U[STEEP_MODES[mode]] if mode in STEEP_MODES else None
    if planted and planted not in attrs_all:           # revision round 2: proxy-chemistry feature
        planted_vec = U["chem_feat"][planted]
    else:
        planted_vec = attrs_all[planted] if planted else None
    extra_all = None
    planted2 = task.get("attribute2") or ""
    if planted2:                                       # revision round 3: second planted term (mixture)
        vec2 = attrs_all[planted2] if planted2 in attrs_all else U["chem_feat"][planted2]
        extra_all = float(task["effect2"]) * np.log(2.0) * vec2
    y_all = draw_positives(detected, scen, planted_vec,
                           task["effect"], task["rate"], rng, steep_z=steep_z, gamma=gamma, extra_all=extra_all)
    det = detected == 1
    idx = cohort_index(U["prot_codes"], y_all)
    y = y_all[idx]
    how = task.get("analyse")
    if task.get("attrs"):                              # revision round 2: explicit list
        names = list(task["attrs"])
    elif how == "planted_only" and planted:
        names = [planted]
    elif how == "extra_only":
        names = list(EXTRA_ATTRIBUTES)
    elif how == "with_extra":
        names = ALL_ATTRIBUTES + EXTRA_ATTRIBUTES
    elif how == "count_only":                          # revision round 2
        names = list(COUNT_ATTRIBUTES)
    elif how == "with_count":                          # revision round 2
        names = ALL_ATTRIBUTES + EXTRA_ATTRIBUTES + COUNT_ATTRIBUTES
    else:
        names = ALL_ATTRIBUTES
    attrs = {nm: attrs_all[nm][idx] for nm in names}
    res = analyse(y, attrs, U["vis10"][idx], U["protein"][idx], detected[idx], det_logit_all[idx],
                  task["want_ci"], n_boot=task.get("n_boot", N_BOOT), run_tool=task.get("run_tool", False),
                  tool_proteins=U["protein"][idx], tool_positions=U["position"][idx],
                  theor=U["theor"][idx] if task.get("artifact1") else None)
    rows = []
    for nm, rec in res.items():
        a_all = attrs_all[nm]
        realised = rt.log_odds_ratio(*rt.weighted_table(y_all[det], a_all[det], np.ones(int(det.sum()))))
        truth = float(task["effect"]) if (planted == nm and not scen.startswith("S0")) else 0.0
        row = {"scenario": scen, "detection_mode": mode, "planted_attribute": planted or "",
               "attribute": nm, "effect_log2_or": float(task["effect"]), "truth_log2_or": truth,
               "rate": float(task["rate"]), "rep": int(task["rep"]), "gamma": gamma,
               "realised_within_detected": realised, "with_ci": bool(task["want_ci"]),
               "n_boot": int(task.get("n_boot", N_BOOT)) if task["want_ci"] else 0}
        if planted2:                                   # revision round 3: mixture tasks only
            row.update(planted_attribute2=planted2, effect2_log2_or=float(task["effect2"]))
        row.update(rec)
        rows.append(row)
    return rows
