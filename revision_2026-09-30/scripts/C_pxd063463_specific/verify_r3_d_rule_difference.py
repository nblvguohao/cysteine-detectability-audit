"""verify_r3_d_rule_difference.py -- POST HOC adversarial verification (round 3), item C_pxd063463_specific.

The proposed text reads 'the own-rule distal signature persists in both strata whereas the trypsin rule shows no distal
association' as support for 'the signature follows the protease'. That is a comparison of two intervals. Here the
difference itself is estimated: own-rule minus trypsin-rule distal Haldane log2 OR for the same positives and the same
background, with one protein-clustered bootstrap resample shared by both rules (own SFC64 RNG), 95% and 97.5%
intervals. Designs: stored sets unstratified (proteome), split by tryptic-HydN identification with stratum-matched
background (proteome), and the pooled strict sets within their own identification stratum.
Output: verify_r3/vr3_d_rule_difference.csv
"""
import numpy as np
import pandas as pd

import verify_r3_common as V

seed = V.SEED0 + 40000
P = {a: V.load(a, "HydP") for a in ("AspN", "CT", "GluC")}
N = {a: V.load(a, "HydN") for a in V.ARMS}
keyN = {a: set(zip(N[a].protein, N[a].position)) for a in V.ARMS}
camN = {a: set(zip(N[a].protein[N[a].label == 1], N[a].position[N[a].label == 1])) for a in V.ARMS}
key_pool = set().union(*keyN.values())
cam_pool = set().union(*camN.values())


def diff_boot(prot, pos, f1, f2, seed, reps=V.REPS):
    labs, codes = np.unique(np.asarray(prot), return_inverse=True)
    K = len(labs)
    c1 = np.zeros((K, 4), np.int64)
    c2 = np.zeros((K, 4), np.int64)
    for j, (s1, s2) in enumerate(((pos & f1, pos & f2), (pos & ~f1, pos & ~f2), (~pos & f1, ~pos & f2),
                                  (~pos & ~f1, ~pos & ~f2))):
        c1[:, j] = np.bincount(codes[s1], minlength=K)
        c2[:, j] = np.bincount(codes[s2], minlength=K)
    t1, t2 = c1.sum(0), c2.sum(0)
    e1, e2 = float(V.hlog2or(*t1)), float(V.hlog2or(*t2))
    vals = []
    for w in V._mult(K, reps, seed):
        s1, s2 = w @ c1, w @ c2
        vals.append(V.hlog2or(s1[:, 0], s1[:, 1], s1[:, 2], s1[:, 3]) - V.hlog2or(s2[:, 0], s2[:, 1], s2[:, 2], s2[:, 3]))
    vals = np.concatenate(vals)
    q = lambda a: float(np.quantile(vals, a))
    return {"own": e1, "trypsin": e2, "difference": e1 - e2, "ci95_low": q(0.025), "ci95_high": q(0.975),
            "ci975_low": q(0.0125), "ci975_high": q(0.9875), "share_reps_le_0": float((vals <= 0).mean()),
            "n_pos": int(pos.sum()), "n_bg": int((~pos).sum())}


rows = []
for arm in ("AspN", "CT", "GluC"):
    own = V.OWN[arm]
    e = V.flags(V.expand(P[arm]), [own, "trypsin"])
    k = list(zip(e.protein, e.position))
    inT = np.array([x in keyN["Trypsin"] for x in k])
    inP = np.array([x in key_pool for x in k])
    cmP = np.array([x in cam_pool for x in k])
    cam = e.label.values == 1
    fo = e[f"f_{own}_distal_6_12"].values
    ft = e["f_trypsin_distal_6_12"].values
    strict_p = cam & inP & ~cmP
    designs = [("stored_unstratified", np.ones(len(e), bool), cam),
               ("in_trypsin_hydn_matched", inT, cam), ("not_in_trypsin_hydn_matched", ~inT, cam),
               ("strict_pooled_in_own_stratum", inP & ~(cam & ~strict_p), strict_p)]
    for name, m, ps in designs:
        pos = m & ps
        bg = m & ~cam
        keep = pos | bg
        r = diff_boot(e.protein.values[keep], pos[keep], fo[keep], ft[keep], seed)
        r.update({"arm": arm, "design": name, "seed": seed})
        seed += 1
        rows.append(r)
    print("done", arm, flush=True)
out = pd.DataFrame(rows)
out.to_csv(f"{V.OUT}/vr3_d_rule_difference.csv", index=False)
pd.set_option("display.width", 250, "display.max_columns", 30)
print(out[["arm", "design", "n_pos", "n_bg", "own", "trypsin", "difference", "ci95_low", "ci95_high", "ci975_low",
           "ci975_high", "share_reps_le_0"]].round(3).to_string(index=False))
V.dump_prov("verify_r3_d_rule_difference", {"reps": V.REPS, "seeds": f"{V.SEED0 + 40000}..{seed - 1}"})
