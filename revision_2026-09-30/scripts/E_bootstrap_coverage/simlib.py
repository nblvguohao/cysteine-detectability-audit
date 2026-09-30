# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- POST HOC revision analysis (2026-09-30). Not registered, not pre-specified.

Core library for the null / non-null simulation of the protein-clustered percentile bootstrap that
decides every claim verdict (Haldane log2 OR, 5,000 replicates, 2.5th-97.5th percentiles).

Fidelity:
  * percentile interval: cys_audit.stats.boot_log2_or (released tool, imported unchanged);
  * the tool's replicate generator (cys_audit.stats.multiplicities) is also used for BCa and for the
    Mantel-Haenszel (secondary caliber) bootstrap;
  * log2 OR, MH statistic, quantile strata, random control, verdict rule: compiled VERBATIM from the
    re-test pipeline's source files by AST extraction (no re-implementation; see `pipeline()`),
    and checked against the vectorised versions used for speed (validate.py).

Nothing here writes outside the caller's own output folder.
"""
from __future__ import annotations

import ast
import hashlib
import math
import os
import sys

sys.dont_write_bytecode = True

import numpy as np
from scipy import optimize, special
from scipy import stats as sps

TOOL_SRC = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/cys-audit/src"
REPO_SCRIPTS = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/repo/scripts"
PIPE_MAIN = os.path.join(REPO_SCRIPTS, "run_phase2_claims_under_detectability_control.py")
PIPE_2B = os.path.join(REPO_SCRIPTS, "run_phase2b_claims_backfill.py")
TOOL_STATS = os.path.join(TOOL_SRC, "cys_audit", "stats.py")
if TOOL_SRC not in sys.path:
    sys.path.insert(0, TOOL_SRC)
from cys_audit import stats as tool_stats  # noqa: E402

LEVEL = 0.95
Z975 = float(sps.norm.ppf(0.975))
LN2 = math.log(2.0)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------------------------
# verbatim pipeline functions
# ---------------------------------------------------------------------------------------------
MAIN_NAMES = ["SEED", "REPLICATES", "N_STRATA", "RANDOM_DRAWS", "SURVIVE_RETENTION", "HALDANE",
              "log_odds_ratio", "weighted_table", "make_weighted_table", "make_mh_log_odds_ratio",
              "cluster_weight_matrix", "bootstrap_interval", "quantile_strata",
              "random_control_effect", "classify", "DOWNGRADE"]
B2_NAMES = ["NULL_VERDICT_BUCKET", "crosses_zero", "classify_null"]


def _extract(path, names):
    src = open(path, encoding="utf-8").read()
    tree = ast.parse(src)
    keep = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            keep.append(node)
        elif isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in names
                                                  for t in node.targets):
            keep.append(node)
    ns = {"np": np, "__name__": "pipeline_verbatim"}
    exec(compile(ast.Module(body=keep, type_ignores=[]), path, "exec"), ns)
    missing = [n for n in names if n not in ns]
    if missing:
        raise RuntimeError(f"not found in {path}: {missing}")
    return ns


_PIPE = None


def pipeline():
    global _PIPE
    if _PIPE is None:
        ns = _extract(PIPE_MAIN, MAIN_NAMES)
        ns.update(_extract(PIPE_2B, B2_NAMES))
        _PIPE = ns
    return _PIPE


# ---------------------------------------------------------------------------------------------
# generative model
# ---------------------------------------------------------------------------------------------
GH_X, GH_W = np.polynomial.hermite_e.hermegauss(80)
GH_W = GH_W / GH_W.sum()

NB_SIZE, NB_MEAN_EXTRA = 2.0, 7.5          # cluster size = 1 + NB(size 2, mean 7.5): mean 8.5
LN_SIGMA = 0.94                            # anchored on the PXD063463 NEG_A cohort (log-scale SD)


def _disc_lognormal_mean(mu, sigma=LN_SIGMA, kmax=20000):
    k = np.arange(1, kmax + 1, dtype=float)
    hi = sps.norm.cdf((np.log(k + 0.5) - mu) / sigma)
    lo = np.where(k == 1, 0.0, sps.norm.cdf((np.log(np.maximum(k - 0.5, 1e-12)) - mu) / sigma))
    return float((k * (hi - lo)).sum())


LN_MU = optimize.brentq(lambda mu: _disc_lognormal_mean(mu) - 8.5, 0.0, 4.0)


def cluster_sizes(rng, K, kind, mean_size=8.5):
    """kind 'nb': 1 + NB(size 2); 'lognormal': max(1, round(LN(mu, 0.94))) (heavy right tail).
    mean_size rescales either family (default 8.5, the task's 8-9)."""
    if kind == "nb":
        extra = mean_size - 1.0
        return 1 + rng.negative_binomial(NB_SIZE, NB_SIZE / (NB_SIZE + extra), size=K)
    if kind == "lognormal":
        mu = LN_MU if mean_size == 8.5 else optimize.brentq(
            lambda m: _disc_lognormal_mean(m) - mean_size, -3.0, 6.0)
        return np.maximum(1, np.rint(rng.lognormal(mu, LN_SIGMA, size=K))).astype(np.int64)
    if kind == "singleton":
        return np.ones(K, dtype=np.int64)
    raise ValueError(kind)


def mean_expit(a, s):
    return float(np.sum(GH_W * special.expit(a + s * GH_X)))


def solve_intercept(target, s):
    return optimize.brentq(lambda a: mean_expit(a, s) - target, -30, 30)


def population_table(alpha0, beta0, su, sv, rho, gamma):
    """Expected per-site 2x2 probabilities (attr&pos, attr&neg, noattr&pos, noattr&neg)."""
    z1 = GH_X[:, None]
    z2 = GH_X[None, :]
    w = GH_W[:, None] * GH_W[None, :]
    u = su * z1 + 0.0 * z2
    v = sv * (rho * z1 + math.sqrt(max(0.0, 1 - rho * rho)) * z2)
    q = special.expit(beta0 + v)
    p1 = special.expit(alpha0 + u + gamma)
    p0 = special.expit(alpha0 + u)
    e11 = float((w * q * p1).sum())
    e10 = float((w * q * (1 - p1)).sum())
    e01 = float((w * (1 - q) * p0).sum())
    e00 = float((w * (1 - q) * (1 - p0)).sum())
    return e11, e10, e01, e00


def marginal_log2or(alpha0, beta0, su, sv, rho, gamma):
    e11, e10, e01, e00 = population_table(alpha0, beta0, su, sv, rho, gamma)
    return math.log2(e11 * e00 / (e10 * e01)), e11 + e01


def calibrate(pi_y, pi_x, su, sv, rho, target_log2or=None, gamma=0.0):
    """Returns alpha0, beta0, gamma and the population (pooled) log2 OR, holding the overall
    positive rate at pi_y and the attribute prevalence at pi_x. If target_log2or is given, gamma is
    solved so that the POOLED (marginal) log2 OR equals it (rho must then be 0 for a clean target)."""
    beta0 = solve_intercept(pi_x, sv)

    def alpha_for(g):
        def f(a):
            return marginal_log2or(a, beta0, su, sv, rho, g)[1] - pi_y
        return optimize.brentq(f, -30, 30)

    if target_log2or is not None:
        def h(g):
            return marginal_log2or(alpha_for(g), beta0, su, sv, rho, g)[0] - target_log2or
        gamma = optimize.brentq(h, -10, 10)
    alpha0 = alpha_for(gamma)
    theta, pos = marginal_log2or(alpha0, beta0, su, sv, rho, gamma)
    return {"alpha0": alpha0, "beta0": beta0, "gamma": gamma, "gamma_log2": gamma / LN2,
            "theta_pooled_log2": theta, "positive_rate": pos}


def generate(rng, K, size_kind, par, su, sv, rho, mean_size=8.5, strat_shift=0.55):
    m = cluster_sizes(rng, K, size_kind, mean_size)
    N = int(m.sum())
    z1 = rng.standard_normal(K)
    z2 = rng.standard_normal(K)
    u = su * z1
    v = sv * (rho * z1 + math.sqrt(max(0.0, 1 - rho * rho)) * z2)
    cid = np.repeat(np.arange(K), m)
    x = rng.random(N) < special.expit(par["beta0"] + v[cid])
    y = rng.random(N) < special.expit(par["alpha0"] + u[cid] + par["gamma"] * x)
    s = strat_shift * y + rng.standard_normal(N)  # propensity-like covariate, AUC ~0.65, x indep. given y
    return {"cid": cid, "x": x, "y": y, "s": s, "K": K}


# ---------------------------------------------------------------------------------------------
# per-cluster tables
# ---------------------------------------------------------------------------------------------
def tool_counts(cid, y, x, K):
    """cys_audit order: pos&flag, pos&~flag, bg&flag, bg&~flag."""
    code = np.where(y, np.where(x, 0, 1), np.where(x, 2, 3))
    out = np.zeros((K, 4), dtype=np.int64)
    np.add.at(out, (cid, code), 1)
    return out


def reindex(cid):
    uniq, inv = np.unique(cid, return_inverse=True)
    return inv, len(uniq)


# ---------------------------------------------------------------------------------------------
# interval methods on a K x 4 table (tool order)
# ---------------------------------------------------------------------------------------------
def log2_or_rows(t):
    return tool_stats.log2_or(t[..., 0], t[..., 1], t[..., 2], t[..., 3])


def percentile_and_bca(counts, reps, seed):
    """Percentile (the paper's, via the released tool) and BCa from the same replicates."""
    point, (lo, hi), _, vals = tool_stats.boot_log2_or(counts, reps, seed, LEVEL)
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    B = vals.size
    prop = (np.sum(vals < point) + 0.5 * np.sum(vals == point)) / B
    prop = min(max(prop, 0.5 / B), 1 - 0.5 / B)
    z0 = sps.norm.ppf(prop)
    tot = counts.sum(axis=0)
    jack = log2_or_rows(tot[None, :] - counts)
    d = jack.mean() - jack
    den = 6.0 * (np.sum(d * d) ** 1.5)
    a = float(np.sum(d ** 3) / den) if den > 0 else 0.0
    out = {"point": float(point), "pct_lo": lo, "pct_hi": hi, "boot_sd": float(vals.std(ddof=1)),
           "bca_z0": float(z0), "bca_a": a}
    zl, zh = sps.norm.ppf(0.025), sps.norm.ppf(0.975)
    try:
        a1 = sps.norm.cdf(z0 + (z0 + zl) / (1 - a * (z0 + zl)))
        a2 = sps.norm.cdf(z0 + (z0 + zh) / (1 - a * (z0 + zh)))
        out["bca_lo"] = float(np.quantile(vals, a1))
        out["bca_hi"] = float(np.quantile(vals, a2))
    except Exception:  # pragma: no cover
        out["bca_lo"] = out["bca_hi"] = float("nan")
    return out


def wald_sandwich(counts, haldane=True):
    """Cluster-robust (GEE independence, logit link) Wald intervals for the attribute log OR,
    centred on the Haldane log2 OR (haldane=True) with the sandwich evaluated at the
    Haldane-adjusted fitted proportions. Closed forms for a single binary covariate:
      influence of cluster i on the slope:  r1_i/(N1 w1) - r0_i/(N0 w0)
      Mancl-DeRouen (bias-reduced) version: r1_i/((N1-n1_i) w1) - r0_i/((N0-n0_i) w0)
    CR0: normal quantile; CR1: x G/(G-1) with t(G-1); MD: t(G-1)."""
    c = np.asarray(counts, float)
    A_, B_, C_, D_ = c[:, 0], c[:, 1], c[:, 2], c[:, 3]   # pos&attr, pos&~attr, neg&attr, neg&~attr
    n1, n0 = A_ + C_, B_ + D_
    N1, N0 = n1.sum(), n0.sum()
    h = 0.5 if haldane else 0.0
    p1 = (A_.sum() + h) / (N1 + 2 * h)
    p0 = (B_.sum() + h) / (N0 + 2 * h)
    beta = (special.logit(p1) - special.logit(p0))
    w1, w0 = p1 * (1 - p1), p0 * (1 - p0)
    r1 = A_ - p1 * n1
    r0 = B_ - p0 * n0
    G = c.shape[0]
    inf0 = r1 / (N1 * w1) - r0 / (N0 * w0)
    var0 = float(np.sum(inf0 ** 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        inf_md = r1 / ((N1 - n1) * w1) - r0 / ((N0 - n0) * w0)
    var_md = float(np.sum(inf_md ** 2)) if np.all(np.isfinite(inf_md)) else float("nan")
    tq = float(sps.t.ppf(0.975, max(G - 1, 1)))
    b2 = beta / LN2
    se0 = math.sqrt(var0) / LN2
    se1 = math.sqrt(var0 * G / max(G - 1, 1)) / LN2
    semd = math.sqrt(var_md) / LN2 if np.isfinite(var_md) else float("nan")
    return {"wald_point": b2, "se_cr0": se0, "se_cr1": se1, "se_md": semd,
            "cr0_lo": b2 - Z975 * se0, "cr0_hi": b2 + Z975 * se0,
            "cr1t_lo": b2 - tq * se1, "cr1t_hi": b2 + tq * se1,
            "mdt_lo": b2 - tq * semd, "mdt_hi": b2 + tq * semd}


def analyse_table(counts, reps, seed):
    out = percentile_and_bca(counts, reps, seed)
    out.update(wald_sandwich(counts))
    out["K"] = int(counts.shape[0])
    out["n"] = int(counts.sum())
    out["n_pos"] = int(counts[:, 0].sum() + counts[:, 1].sum())
    out["n_attr"] = int(counts[:, 0].sum() + counts[:, 2].sum())
    out["n11"] = int(counts[:, 0].sum())
    return out


# ---------------------------------------------------------------------------------------------
# Mantel-Haenszel over propensity-like quintile strata (secondary caliber), tool replicates
# ---------------------------------------------------------------------------------------------
def mh_stat(tab, n_strata=5):
    """Vectorised copy of the statistic in the pipeline's make_mh_log_odds_ratio (checked in
    validate.py). tab: R x (n_strata*4), pipeline cell order n11, n10, n01, n00 per stratum."""
    t = tab.reshape(tab.shape[0], n_strata, 4)
    n11, n10, n01, n00 = t[..., 0], t[..., 1], t[..., 2], t[..., 3]
    total = n11 + n10 + n01 + n00
    usable = (total > 0) & (n11 + n01 > 0) & (n10 + n00 > 0) & (n11 + n10 > 0) & (n01 + n00 > 0)
    safe = np.where(usable, total, 1.0)
    num = np.sum(np.where(usable, n11 * n00 / safe, 0.0), axis=1)
    den = np.sum(np.where(usable, n10 * n01 / safe, 0.0), axis=1)
    ok = usable.any(axis=1) & (num > 0) & (den > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(ok, np.log2(num / den), np.nan)


def mh_percentile(cid, K, y, x, s, reps, seed):
    P = pipeline()
    strata = P["quantile_strata"](s)
    cell_code = np.where(x & y, 0, np.where(x & ~y, 1, np.where(~x & y, 2, 3)))
    code = strata * 4 + cell_code
    cell = np.zeros((K, 20), dtype=np.int64)
    np.add.at(cell, (cid, code), 1)
    point = float(mh_stat(cell.sum(axis=0)[None, :].astype(float))[0])
    cellf = cell.astype(float)
    vals = [mh_stat(w.astype(float) @ cellf) for w in tool_stats.multiplicities(K, reps, seed)]
    vals = np.concatenate(vals)
    lo, hi = tool_stats.percentile_interval(vals, LEVEL)
    return {"mh_point": point, "mh_lo": lo, "mh_hi": hi, "mh_finite": int(np.isfinite(vals).sum())}


# ---------------------------------------------------------------------------------------------
# designs
# ---------------------------------------------------------------------------------------------
def subsample_1to1(rng, y):
    """All positives + an equal number of random negatives (the size of the matched set and of
    one random-control draw; matching on a covariate unrelated to the attribute is equivalent)."""
    pos = np.flatnonzero(y)
    neg = np.flatnonzero(~y)
    n_neg = min(len(pos), len(neg))
    picked = rng.choice(neg, size=n_neg, replace=False)
    return np.concatenate([pos, picked])


def random_control_vectorised(cid, y, x, reps, seed, draws=20):
    """Same algorithm as the pipeline's random_control_effect (20 fixed-seed draws of as many
    random negatives as positives; each bootstrap replicate picks one draw at random and resamples
    that draw's clusters), vectorised per draw. Point = mean of the 20 draw estimates."""
    pos = np.flatnonzero(y)
    neg = np.flatnonzero(~y)
    n_neg = min(len(pos), len(neg))
    tabs, points = [], []
    for i in range(draws):
        r = np.random.default_rng(seed + i)
        picked = r.choice(neg, size=n_neg, replace=False)
        idx = np.concatenate([pos, picked])
        inv, Kd = reindex(cid[idx])
        t = tool_counts(inv, y[idx], x[idx], Kd)
        tabs.append(t)
        points.append(float(log2_or_rows(t.sum(axis=0))))
    rng = np.random.default_rng(seed)
    pick = rng.integers(0, draws, size=reps)
    vals = []
    for d in range(draws):
        nd = int((pick == d).sum())
        if nd == 0:
            continue
        Kd = tabs[d].shape[0]
        # multinomial(Kd, 1/Kd) multiplicities == counts of Kd uniform draws with replacement
        dr = rng.integers(0, Kd, size=(nd, Kd))
        w = np.bincount((np.arange(nd)[:, None] * Kd + dr).ravel(), minlength=nd * Kd).reshape(nd, Kd)
        s = w @ tabs[d]
        vals.append(log2_or_rows(s))
    vals = np.concatenate(vals)
    return {"rc_point": float(np.mean(points)), "rc_lo": float(np.quantile(vals, 0.025)),
            "rc_hi": float(np.quantile(vals, 0.975))}


# ---------------------------------------------------------------------------------------------
# summaries
# ---------------------------------------------------------------------------------------------
METHODS = {"pct": ("pct_lo", "pct_hi"), "bca": ("bca_lo", "bca_hi"), "cr0": ("cr0_lo", "cr0_hi"),
           "cr1t": ("cr1t_lo", "cr1t_hi"), "mdt": ("mdt_lo", "mdt_hi")}


def ratio_to_target(point, lo, hi, target=0.0):
    """Multiplier of the half-width (measured from the point estimate toward `target`) at which
    that bound reaches `target`: > 1 means the interval excludes it. When `target` is the truth,
    the 95th percentile of this ratio over simulated datasets is the calibration multiplier w0
    (1 for an exactly calibrated interval; > 1 means the interval must be widened)."""
    point, lo, hi = (np.asarray(v, float) for v in (point, lo, hi))
    toward_low = point >= target
    hw = np.where(toward_low, point - lo, hi - point)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.abs(point - target) / hw
    return np.where(hw > 0, r, np.inf)


def ratio_to_zero(point, lo, hi):
    return ratio_to_target(point, lo, hi, 0.0)


# ---------------------------------------------------------------------------------------------
# NEG_A-type cohorts: every protein carries >= 1 positive (all cysteines of site-carrying proteins)
# ---------------------------------------------------------------------------------------------
def generate_selected(rng, K, size_kind, par, su, sv, rho, mean_size_pre, strat_shift=0.55):
    """Draw proteins from the same model and keep the first K that carry at least one positive.
    Parameters (alpha0, pre-selection mean size) are chosen by calibrate_selected() so that the
    SELECTED cohort has the target positive rate and mean cluster size."""
    keep_c, keep_x, keep_y, keep_s = [], [], [], []
    got = 0
    while got < K:
        d = generate(rng, 2 * K, size_kind, par, su, sv, rho, mean_size=mean_size_pre, strat_shift=strat_shift)
        npos = np.bincount(d["cid"], weights=d["y"].astype(float), minlength=2 * K)
        ok = np.flatnonzero(npos > 0)[: K - got]
        sel = np.isin(d["cid"], ok)
        remap = np.full(2 * K, -1)
        remap[ok] = np.arange(got, got + len(ok))
        keep_c.append(remap[d["cid"][sel]])
        keep_x.append(d["x"][sel])
        keep_y.append(d["y"][sel])
        keep_s.append(d["s"][sel])
        got += len(ok)
    cid = np.concatenate(keep_c)
    order = np.argsort(cid, kind="stable")
    return {"cid": cid[order], "x": np.concatenate(keep_x)[order], "y": np.concatenate(keep_y)[order],
            "s": np.concatenate(keep_s)[order], "K": K}


def calibrate_selected(pos_rate, attr_rate, mean_size, su, sv, size_kind="nb", n_prot=20000, seed=12345,
                       iters=8):
    """Fixed-seed simulation calibration of (alpha0, pre-selection mean size) for generate_selected."""
    beta0 = solve_intercept(attr_rate, sv)
    m_pre = mean_size
    alpha0 = solve_intercept(pos_rate, su)
    for _ in range(iters):
        def post(a):
            rng = np.random.default_rng(seed)
            d = generate_selected(rng, n_prot, size_kind, {"alpha0": a, "beta0": beta0, "gamma": 0.0},
                                  su, sv, 0.0, m_pre)
            return float(d["y"].mean()), len(d["y"]) / n_prot
        lo_a, hi_a = alpha0 - 3.0, alpha0 + 1.0
        for _ in range(25):
            mid = 0.5 * (lo_a + hi_a)
            if post(mid)[0] < pos_rate:
                lo_a = mid
            else:
                hi_a = mid
        alpha0 = 0.5 * (lo_a + hi_a)
        pr, ms = post(alpha0)
        m_pre = max(1.0 + 1e-6, m_pre * mean_size / ms) if abs(ms - mean_size) > 0.01 else m_pre
        if abs(ms - mean_size) <= 0.01:
            break
    pr, ms = post(alpha0)
    return {"alpha0": alpha0, "beta0": beta0, "gamma": 0.0, "mean_size_pre": m_pre,
            "post_positive_rate": pr, "post_mean_size": ms}
