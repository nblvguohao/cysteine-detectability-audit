"""Step 12 (POST HOC; revision round 3, 2026-09-30): 4,000 further label sets for every null combination of
s10b, pooled with the 1,000 s10b sets (5,000 per combination).

Why this step exists. The round-3 verifier re-simulated the registered-set (label A) nulls with its own seeds
(1,000 sets per null) and counted 1 of the 14 reported contrasts outside the N4 range where s10b counts 0, and
10 of 14 outside the N1 range where s10b counts 11. Each difference is a single statistic near the edge of the
simulated range: under N4 the cysteine-count contrast within iBAQ decile x detectable-count tertile (observed
+0.186; s10b upper limit 0.187, two-sided P 0.060; the verifier's 4,000 further sets gave a pooled P of about
0.05), under N1 the length contrast within the same strata (s10b P 0.038). With 1,000 sets, a count of
'contrasts outside the 95% range' changes with the seed at such edges, and the proposed Results text quoted
such a count ("all 14"). This step draws 4,000 further label sets per combination from the SAME fitted null
probabilities (same fitting call as s10b; asserted through each fit's AIC and expected positives), with seeds
fixed here before the run, and pools them with the s10b sets.

Declared before the run (applied whatever the outcome):
  - the pooled 5,000-set summaries replace the 1,000-set summaries of s10b as the primary null results in the
    report and in every proposed text; the s10b (1,000-set) and extension-only (4,000-set) values are listed
    beside the pooled ones in s12_null_pooled_summary.csv;
  - 'outside' means outside the pooled 2.5-97.5 percentile range (np.percentile, linear interpolation), as in
    s10b; P values as in s10b: p_ge = (1 + #null >= observed) / (1 + reps), p_le likewise,
    p_two = min(1, 2 min(p_ge, p_le)); with 5,000 sets the smallest two-sided P is 2/5001 = 0.0004;
  - a contrast whose pooled two-sided P lies in [0.04, 0.06] is flagged 'edge' (about +-2 Monte Carlo
    standard errors around 0.05 with 5,000 sets); proposed texts give such a contrast's P instead of counting
    it silently either way.
Seeds: numpy SeedSequence(entropy = SEED + 840 + i, spawn_key = (chunk,)), i = index of the combination in
s10b's COMBOS list (0..11), chunks 0..15 of 250 sets; 4 worker processes. SEED + 840..851 are not used by any
other step of this item (s04 +1..12, s08 +801..804, s08b/s10b +810..821, s10 +831..834).
Outputs: s12_null_ext_draws.csv.gz, s12_null_pooled_summary.csv, s12_contrast_table_pooled.csv,
s12_null_pooled_summary.json; the report tables are printed to stdout (saved as s12_stdout.txt by the caller).
"""
import gzip
import io
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
from s10_lib import Design2, simulate_chunk2
from s10b_null_total import COMBOS, CONTRAST, NULLS, STATS

warnings.filterwarnings("ignore")
EXT_REPS, CHUNK_SIZE, WORKERS = 4000, 250, 4
N_CHUNKS = EXT_REPS // CHUNK_SIZE
EDGE = (0.04, 0.06)
LABELS = ("A", "S", "S_map")
SHORT = {n: n[:2] for n in NULLS}


def summarize(v, o):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    pge = float((1 + (v >= o).sum()) / (1 + len(v)))
    ple = float((1 + (v <= o).sum()) / (1 + len(v)))
    lo, hi = float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))
    return {"reps": int(len(v)), "mean": float(v.mean()), "lo": lo, "hi": hi,
            "p_two": float(min(1.0, 2 * min(pge, ple))), "inside": bool(lo <= o <= hi)}


def fmt_val(stat, x):
    if stat == "dAIC_D_minus_T":
        return f"{x:.1f}"
    if stat == "b_det_B":
        return f"{x:.2f}"
    return f"{x:.3f}"


def fmt_p(p):
    if p >= 0.095:
        return f"{p:.2f}"
    if p >= 0.0095:
        return f"{p:.3f}"
    return f"{p:.4f}"


def main():
    F = pd.read_csv(RES / "s08_group_features_windows.csv")
    H = F[F.n_det >= 1].reset_index(drop=True)
    H["log_det"] = np.log(H.n_det); H["log_cys"] = np.log(H.n_cys)
    H["opp20"] = np.minimum(H.n_det, 20); H["cys20"] = np.minimum(H.n_cys, 20)
    cols = ["log10_ibaq", "length", "n_cys", "n_det"]
    D = Design2(H[cols])
    old_meta = json.loads((RES / "s10b_null_summary.json").read_text(encoding="utf-8"))
    old_combo = {(c["label"], c["null"]): c for c in old_meta["combos"]}

    tasks, meta, observed = [], {}, {}
    for i, (lab, null, ent_old) in enumerate(COMBOS):
        link, rhs = NULLS[null]
        fam = sm.families.Binomial() if link == "logit" else sm.families.Binomial(link=sm.families.links.CLogLog())
        m = smf.glm(f"{lab} ~ {rhs}", data=H, family=fam).fit()   # the s10b call, unchanged
        p = m.fittedvalues.to_numpy(float)
        oc = old_combo[(lab, null)]
        assert abs(float(m.aic) - oc["aic"]) < 1e-6 and abs(float(p.sum()) - oc["expected_pos"]) < 1e-6, (lab, null)
        if lab not in observed:
            observed[lab] = D.stats_all(H[lab].to_numpy(int))
            for s in STATS:
                a, b = float(observed[lab][s]), float(old_meta["observed"][lab][s])
                assert (np.isnan(a) and np.isnan(b)) or abs(a - b) < 1e-9, (lab, s)
        ent = SEED + 840 + i
        meta[(lab, null)] = {"label": lab, "null": null, "combo_index": i, "entropy_s10b": ent_old,
                             "entropy_extension": ent, "aic": float(m.aic), "expected_pos": float(p.sum()),
                             "observed_pos": int(H[lab].sum()), "irls_converged": bool(m.converged),
                             "fit_identical_to_s10b": True}
        for c in range(N_CHUNKS):
            tasks.append(((lab, null), (H[cols], p, CHUNK_SIZE, ent, (c,))))
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        results = list(ex.map(simulate_chunk2, [t[1] for t in tasks]))
    rows = []
    for (key, _), res in zip(tasks, results):
        for r in res:
            rows.append({"label": key[0], "null": key[1], **r})
    EXT = pd.DataFrame(rows)
    EXT["draw"] = 1000 + EXT.groupby(["label", "null"]).cumcount()
    buf = io.StringIO()
    EXT.to_csv(buf, index=False, lineterminator="\n")
    with open(RES / "s12_null_ext_draws.csv.gz", "wb") as fh:          # mtime=0: byte-identical reruns
        with gzip.GzipFile(filename="s12_null_ext_draws.csv", mode="wb", fileobj=fh, mtime=0) as gz:
            gz.write(buf.getvalue().encode("utf-8"))

    OLD = pd.read_csv(RES / "s10b_null_draws.csv")
    assert OLD.groupby(["label", "null"]).size().eq(1000).all()
    summ = []
    for lab, null, _ in COMBOS:
        o_sub = OLD[(OLD.label == lab) & (OLD.null == null)]
        e_sub = EXT[(EXT.label == lab) & (EXT.null == null)]
        assert len(e_sub) == EXT_REPS
        for s in STATS:
            o = float(observed[lab][s])
            a = summarize(o_sub[s], o); b = summarize(e_sub[s], o)
            c = summarize(np.concatenate([o_sub[s].to_numpy(float), e_sub[s].to_numpy(float)]), o)
            summ.append({"label": lab, "null": null, "statistic": s, "observed": o,
                         "s10b_lo": a["lo"], "s10b_hi": a["hi"], "s10b_p_two": a["p_two"], "s10b_inside": a["inside"],
                         "ext_lo": b["lo"], "ext_hi": b["hi"], "ext_p_two": b["p_two"], "ext_inside": b["inside"],
                         "pooled_reps": c["reps"], "pooled_mean": c["mean"], "pooled_lo": c["lo"], "pooled_hi": c["hi"],
                         "pooled_p_two": c["p_two"], "pooled_inside": c["inside"],
                         "pooled_edge": bool(EDGE[0] <= c["p_two"] <= EDGE[1])})
    SUM = pd.DataFrame(summ)
    SUM.to_csv(RES / "s12_null_pooled_summary.csv", index=False, lineterminator="\n")

    ct, counts = [], {}
    for lab in LABELS:
        for s, desc in CONTRAST:
            rec = {"label": lab, "statistic": s, "description": desc, "observed": float(observed[lab][s])}
            for null in NULLS:
                q = SUM[(SUM.label == lab) & (SUM.null == null) & (SUM.statistic == s)].iloc[0]
                k = SHORT[null]
                rec[f"{k}_lo"] = float(q.pooled_lo); rec[f"{k}_hi"] = float(q.pooled_hi)
                rec[f"{k}_p_two"] = float(q.pooled_p_two); rec[f"{k}_inside"] = bool(q.pooled_inside)
                rec[f"{k}_edge"] = bool(q.pooled_edge)
                rec[f"{k}_s10b_inside"] = bool(q.s10b_inside)
            ct.append(rec)
        for null in NULLS:
            k = SHORT[null]
            sub = [r for r in ct if r["label"] == lab]
            counts[f"{lab}|{k}"] = {
                "reported_contrasts": len(sub),
                "outside_pooled": [r["statistic"] for r in sub if not r[f"{k}_inside"]],
                "edge_pooled": [(r["statistic"], round(r[f"{k}_p_two"], 4)) for r in sub if r[f"{k}_edge"]],
                "outside_s10b": [r["statistic"] for r in sub if not r[f"{k}_s10b_inside"]],
                "min_p_two_inside_pooled": min([r[f"{k}_p_two"] for r in sub if r[f"{k}_inside"]] or [float("nan")]),
                "min_p_two_inside_pooled_excluding_edge": min(
                    [r[f"{k}_p_two"] for r in sub if r[f"{k}_inside"] and not r[f"{k}_edge"]] or [float("nan")])}
            counts[f"{lab}|{k}"]["n_outside_pooled"] = len(counts[f"{lab}|{k}"]["outside_pooled"])
            counts[f"{lab}|{k}"]["n_outside_s10b"] = len(counts[f"{lab}|{k}"]["outside_s10b"])
            # the 1,000-set counts recomputed here must equal those s10b reported
            assert set(counts[f"{lab}|{k}"]["outside_s10b"]) == set(
                old_meta["contrasts_outside_null_range"][f"{lab}|{null}"]["outside_95_range"]), (lab, null)
    CT = pd.DataFrame(ct)
    CT.to_csv(RES / "s12_contrast_table_pooled.csv", index=False, lineterminator="\n")
    write_json(RES / "s12_null_pooled_summary.json", {
        "analysis_label": "POST HOC revision analysis 2026-09-30, round 3 (not registered)",
        "reason": "count of contrasts outside the 1,000-set null range unstable at the edge (round-3 verifier)",
        "declared_before_run": {"primary": "pooled 5,000-set summaries", "edge_band_two_sided_p": list(EDGE),
                                "outside": "outside pooled 2.5-97.5 percentile range"},
        "groups": int(len(H)), "reps_s10b": 1000, "reps_extension": EXT_REPS, "reps_pooled": 1000 + EXT_REPS,
        "workers": WORKERS, "combos": [meta[(l, n)] for l, n, _ in COMBOS], "counts": counts,
        "observed": observed})

    # ---- report tables (markdown), in the order of the report: N2, N4, N1, N3
    order = ["N2", "N4", "N1", "N3"]
    by_short = {SHORT[n]: n for n in NULLS}
    print("POST HOC round 3: pooled nulls, 5,000 label sets per combination")
    for lab in LABELS:
        print(f"\n### label {lab} (observed positives on the 6,990 groups: {int(H[lab].sum())})")
        print("| statistic | observed | " + " | ".join(order) + " |")
        print("|---|---|" + "---|" * len(order))
        for r in [r for r in ct if r["label"] == lab]:
            cells = []
            for k in order:
                txt = f"{fmt_val(r['statistic'], r[k + '_lo'])} to {fmt_val(r['statistic'], r[k + '_hi'])} ({fmt_p(r[k + '_p_two'])})"
                if not r[k + "_inside"]:
                    txt = f"**{txt}**"
                if r[k + "_edge"]:
                    txt += " edge"
                cells.append(txt)
            print(f"| {r['statistic']} | {fmt_val(r['statistic'], r['observed'])} | " + " | ".join(cells) + " |")
        print("| outside (pooled) | | " + " | ".join(
            f"{counts[f'{lab}|{k}']['n_outside_pooled']} of 14 (1,000 sets: {counts[f'{lab}|{k}']['n_outside_s10b']})"
            for k in order) + " |")
    print("\ncounts:")
    print(json.dumps(counts, indent=1))


if __name__ == "__main__":
    main()
