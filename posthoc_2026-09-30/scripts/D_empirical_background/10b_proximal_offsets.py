"""Step 10b (POST HOC revision analysis, item D_empirical_background), written 2026-09-30 after verification round 2
and after the step-10 results had been seen. Descriptive only; not used for any estimate or reading rule.

Question: the combined restriction (bc) reaches the observed class (e) in the distal K/R band but moves the proximal
(1-3) band away from it (+0.84). Where within the proximal band do the unenriched arm's observed cysteines differ from
the enriched arm's? Per offset -3..-1 and +1..+3, the share of cysteines with a trypsin-competent K/R (Cys-Audit
competent_mask, rule 'trypsin': K/R not before P) at that offset, for sites and for the backgrounds (b), (bc), (e) and
(e) restricted to theoretically detectable cysteines. A K/R at -1 makes the cysteine the N-terminal residue of the
tryptic peptide (after cleavage following that K/R).

Output: r2_proximal_offsets.csv
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

from common import BOOT_REPS, IN_FASTA_GZ, MASTER_SEED, RESULTS, read_fasta_gz
from cys_audit import proteases, stats  # noqa: E402

OFFSETS = [-3, -2, -1, 1, 2, 3]


def main():
    T = pd.read_csv(f"{RESULTS}/cys_scores_trypsin_arm.csv", dtype={"protein": str})
    seqs = read_fasta_gz(IN_FASTA_GZ)
    masks = {}
    rows = np.zeros((len(T), len(OFFSETS)), dtype=bool)
    for k, (p, pos1) in enumerate(zip(T.protein.to_numpy(), T.position.to_numpy())):
        s = seqs[p]
        if p not in masks:
            masks[p] = proteases.competent_mask(s, "trypsin")
        m = masks[p]
        assert s[pos1 - 1] == "C"
        for j, o in enumerate(OFFSETS):
            q = pos1 - 1 + o
            rows[k, j] = 0 <= q < len(s) and bool(m[q])
    lab, det = T.label.to_numpy(), T.detected.to_numpy() == 1
    theo, emp = T.theo_trypsin.to_numpy() == 1, T.emp_obs.to_numpy() == 1
    pool, pos = lab == 0, lab == 1
    # the band flag must equal "any offset" (self-check against the stored kr_prox flag)
    assert (rows.any(axis=1) == T.kr_prox.to_numpy().astype(bool)).all()
    groups = {"sites": pos, "(b) theoretically detectable": pool & theo,
              "(bc) theoretically detectable and observed in the unenriched arm": pool & theo & emp,
              "theoretically detectable, not observed in the unenriched arm": pool & theo & ~emp,
              "(e) observed in the ABE arm": pool & det, "(e) restricted to theoretically detectable": pool & det & theo,
              "sites observed in the unenriched arm": pos & emp, "sites not observed in the unenriched arm": pos & ~emp}
    no_m1 = rows[:, [j for j, o in enumerate(OFFSETS) if o != -1]].any(axis=1)
    out = []
    for g, m in groups.items():
        r = {"group": g, "n": int(m.sum()), "share_any_1_3": float(rows[m].any(axis=1).mean()),
             "share_any_1_3_excluding_offset_-1": float(no_m1[m].mean())}
        for j, o in enumerate(OFFSETS):
            r[f"share_offset_{o:+d}"] = float(rows[m, j].mean())
        out.append(r)
    # point log2 odds ratios (Haldane), sites against each background: whole band, band without offset -1, offset -1
    m1 = rows[:, OFFSETS.index(-1)]
    for g in ("(b) theoretically detectable", "(bc) theoretically detectable and observed in the unenriched arm",
              "(e) observed in the ABE arm", "(e) restricted to theoretically detectable"):
        bgm = groups[g]
        r = {"group": f"log2 OR sites vs {g} (point, Haldane)", "n": int(bgm.sum())}
        for nm, fl in (("share_any_1_3", rows.any(axis=1)), ("share_any_1_3_excluding_offset_-1", no_m1),
                       ("share_offset_-1", m1)):
            a, b = fl[pos].sum(), (~fl[pos]).sum()
            c, d = fl[bgm].sum(), (~fl[bgm]).sum()
            r[nm] = float(np.log2(((a + .5) / (b + .5)) / ((c + .5) / (d + .5))))
        out.append(r)
    # paired protein bootstrap (primary multiplicities of steps 3/9: all proteins, seed 20260930 + 7, 5,000 replicates)
    codes, labels = stats.cluster_index(list(T.protein.to_numpy()))
    K = len(labels)
    feats = {"band_1_3": rows.any(axis=1), "band_1_3_excluding_offset_-1": no_m1, "offset_-1": m1}
    bgs = {"bc": pool & theo & emp, "e": pool & det, "b": pool & theo}
    blocks, keys = [], []
    for fn, fl in feats.items():
        for bn, bm in bgs.items():
            role = np.full(len(T), -1)
            role[bm] = 0
            role[pos] = 1
            blocks.append(stats.two_by_two_counts(codes, K, role == 1, fl, mask=role >= 0))
            keys.append((fn, bn))
    cnt = np.concatenate(blocks, axis=1)
    tot = cnt.sum(axis=0).reshape(len(keys), 4)
    pts = {k: float(stats.log2_or(*tot[i])) for i, k in enumerate(keys)}
    reps = []
    for Wm in stats.multiplicities(K, BOOT_REPS, MASTER_SEED + 7):
        s = (Wm @ cnt).reshape(Wm.shape[0], len(keys), 4)
        reps.append(stats.log2_or(s[..., 0], s[..., 1], s[..., 2], s[..., 3]))
    reps = np.concatenate(reps)
    for fn in feats:
        for bn in bgs:
            v = reps[:, keys.index((fn, bn))]
            lo, hi = stats.percentile_interval(v, 0.95)
            out.append({"group": f"paired: sites vs ({bn}), {fn}", "n": int(bgs[bn].sum()), "log2_or": pts[(fn, bn)],
                        "ci95_low": lo, "ci95_high": hi})
        for x, y in (("bc", "e"), ("b", "e")):
            v = reps[:, keys.index((fn, x))] - reps[:, keys.index((fn, y))]
            lo, hi = stats.percentile_interval(v, 0.95)
            out.append({"group": f"paired: ({x}) - ({y}), {fn}", "n": np.nan,
                        "log2_or": pts[(fn, x)] - pts[(fn, y)], "ci95_low": lo, "ci95_high": hi,
                        "boot_p": stats.bootstrap_p(v, 0.0)})
    O = pd.DataFrame(out)
    O.to_csv(f"{RESULTS}/r2_proximal_offsets.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 20):
        print(O.round(4).to_string())


if __name__ == "__main__":
    main()
