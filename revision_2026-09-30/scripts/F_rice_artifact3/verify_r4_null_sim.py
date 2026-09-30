"""VERIFIER round 4 (adversarial; POST HOC), item F_rice_artifact3, part (b): independent re-simulation of the
calibrated nulls that the proposed Results text quotes (registered set A under N4 and N2; site-level set S
under N4), with the verifier's own null fits (Newton, converged; spline basis without intercept) and own
statistic code (verify_r4_lib), and own seeds (SeedSequence(990000 + combo, spawn_key=(chunk,))).
Also recomputes the 14 observed statistics on the 6,990 groups with >= 1 detectable Cys.

Usage: python verify_r4_null_sim.py [reps_A_N4 reps_A_N2 reps_S_N4]   (defaults 5000 2000 3000); 4 workers.
Outputs: verify_r4/null_sim_r4.json, verify_r4/null_sim_r4_draws.csv.gz
"""
import gzip
import io
import json
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import rankdata

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from verify_r4_lib import deciles, dummies, fit_binary, greedy_match, spline_basis  # noqa

W = pathlib.Path("C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30")
RES = W / "results" / "F_rice_artifact3"
OUT = RES / "verify_r4"
STATS14 = ["r_dec_len", "r_dec_cys", "r_dt_len", "r_dt_cys", "r_m_len", "r_m_cys", "r_m_undet", "b_undet_M3u",
           "b_len_M3u", "r_mc_det", "r_mc_len", "dAIC_D_minus_T", "b_det_B", "d_det_minus_undet_L"]


def rb(a, b):
    if len(a) < 2 or len(b) < 2:
        return np.nan
    rk = rankdata(np.concatenate([a, b]))
    U = rk[:len(a)].sum() - len(a) * (len(a) + 1) / 2.0
    return 2.0 * U / (len(a) * len(b)) - 1.0


def newton_logit(y, X, it=30):
    b = np.zeros(X.shape[1]); ll_old = -np.inf
    for _ in range(it):
        eta = X @ b; mu = expit(eta)
        ll = float(np.sum(y * eta - np.logaddexp(0, eta)))
        if abs(ll - ll_old) < 1e-9:
            break
        ll_old = ll
        w = mu * (1 - mu)
        b = b + np.linalg.solve(X.T @ (X * w[:, None]), X.T @ (y - mu))
    eta = X @ b
    return b, float(np.sum(y * eta - np.logaddexp(0, eta)))


class Des:
    def __init__(self, H):
        self.x = H.log10_ibaq.to_numpy(float); self.L = H.length.to_numpy(float)
        self.nc = H.n_cys.to_numpy(float); self.nd = H.n_det.to_numpy(float); self.nu = self.nc - self.nd
        self.dec = deciles(self.x)
        tert = np.where(self.nd <= 3, 0, np.where(self.nd <= 6, 1, 2))
        self.dt = self.dec * 10 + tert
        self.opp = np.minimum(self.nd, 20).astype(int); self.cys = np.minimum(self.nc, 25).astype(int)
        self.dec_idx = [np.flatnonzero(self.dec == c) for c in np.unique(self.dec)]
        self.dt_idx = [np.flatnonzero(self.dt == c) for c in np.unique(self.dt)]
        one = np.ones(len(self.x)); l2len = np.log2(self.L); l2d = np.log2(1 + self.nd); l2c = np.log2(1 + self.nc)
        self.X_m3u = np.column_stack([one, self.x, l2d, l2len, self.nu])
        self.X_D = np.column_stack([one, self.x, l2len, l2d]); self.X_T = np.column_stack([one, self.x, l2len, l2c])
        self.X_B = np.column_stack([one, self.x, l2len, l2d, l2c])
        self.X_L = np.column_stack([one, self.x, l2len, self.nd, self.nu])

    def strat(self, y, v, idxs):
        rs = []
        for idx in idxs:
            yy = y[idx]
            a = v[idx[yy == 1]]; b = v[idx[yy == 0]]
            if len(a) >= 10 and len(b) >= 10:
                rs.append(rb(a, b))
        return float(np.mean(rs)) if rs else np.nan

    def stats(self, y):
        y = np.asarray(y, int); yf = y.astype(float)
        o = {"n_pos": int(y.sum())}
        o["r_dec_len"] = self.strat(y, self.L, self.dec_idx); o["r_dec_cys"] = self.strat(y, self.nc, self.dec_idx)
        o["r_dt_len"] = self.strat(y, self.L, self.dt_idx); o["r_dt_cys"] = self.strat(y, self.nc, self.dt_idx)
        P, N = greedy_match(y, self.x, self.opp)
        o["r_m_len"] = rb(self.L[P], self.L[N]); o["r_m_cys"] = rb(self.nc[P], self.nc[N]); o["r_m_undet"] = rb(self.nu[P], self.nu[N])
        P2, N2 = greedy_match(y, self.x, self.cys)
        o["r_mc_det"] = rb(self.nd[P2], self.nd[N2]); o["r_mc_len"] = rb(self.L[P2], self.L[N2])
        b, _ = newton_logit(yf, self.X_m3u); o["b_len_M3u"] = b[3]; o["b_undet_M3u"] = b[4]
        _, lD = newton_logit(yf, self.X_D); _, lT = newton_logit(yf, self.X_T)
        o["dAIC_D_minus_T"] = 2 * (lT - lD)
        bB, _ = newton_logit(yf, self.X_B); o["b_det_B"] = bB[3]
        bL, _ = newton_logit(yf, self.X_L); o["d_det_minus_undet_L"] = bL[3] - bL[4]
        return o


def worker(args):
    H, p, reps, ent, key = args
    D = Des(H)
    rng = np.random.default_rng(np.random.SeedSequence(ent, spawn_key=key))
    return [D.stats((rng.random(len(p)) < p).astype(int)) for _ in range(reps)]


def null_probs(H, lab, null):
    y = H[lab].to_numpy(float)
    B = spline_basis(H.log10_ibaq.to_numpy(float))
    if null in ("N2", "N4"):
        v = np.minimum(H.n_det if null == "N2" else H.n_cys, 20).to_numpy(int)
        Dm, lev = dummies(v)
        X = np.hstack([B, Dm])
        f = fit_binary(y, X, "logit")
    else:
        v = np.log(H.n_det if null == "N1" else H.n_cys).to_numpy(float)
        X = np.column_stack([np.ones(len(H)), H.log10_ibaq, v])
        f = fit_binary(y, X, "cloglog")
    return f["fitted"], f


def main():
    reps = [int(a) for a in sys.argv[1:4]] if len(sys.argv) >= 4 else [5000, 2000, 3000]
    G = pd.read_csv(OUT / "groups_r4.csv")
    H = G[G.n_det >= 1].reset_index(drop=True)
    assert len(H) == 6990
    cols = ["log10_ibaq", "length", "n_cys", "n_det"]
    D = Des(H)
    observed = {lab: D.stats(H[lab].to_numpy(int)) for lab in ("A", "S")}
    item_obs = json.loads((RES / "s12_null_pooled_summary.json").read_text(encoding="utf-8"))["observed"]
    obs_diff = {lab: {s: float(observed[lab][s] - item_obs[lab][s]) for s in STATS14} for lab in ("A", "S")}
    combos = [("A", "N4", reps[0]), ("A", "N2", reps[1]), ("S", "N4", reps[2])]
    CH = 250
    tasks, meta = [], {}
    for ci, (lab, null, R) in enumerate(combos):
        p, f = null_probs(H, lab, null)
        # compare with the item's (non-converged IRLS) null fit: expected positives
        meta[f"{lab}|{null}"] = {"reps": R, "expected_pos": float(p.sum()), "observed_pos": int(H[lab].sum()),
                                 "converged": bool(f["converged"]), "llf": float(f["llf"]), "k": int(f["k"])}
        for c in range(R // CH):
            tasks.append(((lab, null), (H[cols], p, CH, 990000 + ci, (c,))))
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as ex:
        results = list(ex.map(worker, [t[1] for t in tasks]))
    rows = []
    for (key, _), res in zip(tasks, results):
        for r in res:
            rows.append({"label": key[0], "null": key[1], **r})
    DR = pd.DataFrame(rows)
    buf = io.StringIO(); DR.to_csv(buf, index=False, lineterminator="\n")
    with open(OUT / "null_sim_r4_draws.csv.gz", "wb") as fh:
        with gzip.GzipFile(filename="null_sim_r4_draws.csv", mode="wb", fileobj=fh, mtime=0) as gz:
            gz.write(buf.getvalue().encode("utf-8"))
    summ = {}
    for (lab, null, R) in combos:
        sub = DR[(DR.label == lab) & (DR.null == null)]
        rec = {}
        for s in STATS14:
            v = sub[s].to_numpy(float); v = v[np.isfinite(v)]
            o = observed[lab][s]
            pge = (1 + np.sum(v >= o)) / (1 + len(v)); ple = (1 + np.sum(v <= o)) / (1 + len(v))
            lo, hi = np.percentile(v, [2.5, 97.5])
            rec[s] = {"obs": float(o), "lo": float(lo), "hi": float(hi), "p2": float(min(1, 2 * min(pge, ple))),
                      "inside": bool(lo <= o <= hi)}
        summ[f"{lab}|{null}"] = {"per_stat": rec, "n_outside": sum(not rec[s]["inside"] for s in STATS14),
                                 "outside": [s for s in STATS14 if not rec[s]["inside"]],
                                 "min_p_inside": min(rec[s]["p2"] for s in STATS14 if rec[s]["inside"])}
    out = {"label": "VERIFIER round 4 (POST HOC): independent null re-simulation",
           "elapsed_s": round(time.time() - t0, 1), "meta": meta, "observed": observed,
           "observed_minus_item_observed": obs_diff, "summary": summ}
    (OUT / "null_sim_r4.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print("max |observed - item observed|:", {lab: max(abs(v) for v in d.values()) for lab, d in obs_diff.items()})
    print(json.dumps(meta, indent=1))
    for k, v in summ.items():
        print(k, "outside", v["n_outside"], v["outside"], "min P inside", round(v["min_p_inside"], 4))
        print("   ", {s: (round(v["per_stat"][s]["obs"], 4), round(v["per_stat"][s]["lo"], 4), round(v["per_stat"][s]["hi"], 4),
                     round(v["per_stat"][s]["p2"], 4)) for s in STATS14})
    print("elapsed", out["elapsed_s"])


if __name__ == "__main__":
    main()
