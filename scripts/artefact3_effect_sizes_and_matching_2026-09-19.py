"""Separate effect size from sample size in the Artefact 3 recomputation, and show what 1:1 matching does to z.

WHY THIS EXISTS
  The pre-registered run (recompute_artefact3_osa_abundance_2026-09-19.py) landed on branch B2:
  abundance does not account for the protein-level separation. That is the opposite of the claim the
  manuscript used to carry, so before it is reported the statistic itself has to be examined.

  A z conflates effect size with sample size. "Matching" in the original round is described as a
  procedure that made z collapse. Two different things can produce that collapse:
     (i) the effect really disappears within abundance strata, or
     (ii) 1:1 matching throws away most of the control arm, n falls, and z falls with it while the
          effect size is unchanged.
  Only (i) supports the manuscript's claim. This script measures both, so the two are distinguishable.

CRITERIA, WRITTEN BEFORE THE RUN
  Effect size: rank-biserial correlation r = 2*U/(n1*n2) - 1, which is a difference of probabilities
  and does not grow with n. Reported before matching, within each abundance decile, and as the
  decile-weighted mean.
  A 1:1 matched arm is built for comparison: each persulfidated group is paired with the unmodified
  group nearest in log10 iBAQ, without replacement, nearest-first. Its z and r are both reported.
  Rules from the main script are IMPORTED, not reimplemented (project notes 9.24.2), and that module's
  sha256 is recorded here.

  READING RULE, fixed before the run:
    R1 |r| after matching is close to |r| before (within 0.05) AND matched-arm z is much smaller
       -> the collapse in the original round is a sample-size artefact, not a vanished effect.
    R2 |r| after matching falls below 0.05
       -> the effect really does vanish within abundance strata; the original claim reproduces.
    R3 anything else -> report both and claim neither.

This script computes. It does not edit the manuscript.
"""
import csv, hashlib, importlib.util, json, pathlib, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
MAIN = ROOT / "scripts/recompute_artefact3_osa_abundance_2026-09-19.py"
OUT = ROOT / "results/artefact3_effect_sizes_2026-09-19.csv"
AUD = ROOT / "results/artefact3_effect_sizes_2026-09-19_audit.json"
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

spec = importlib.util.spec_from_file_location("a3", MAIN)
a3 = importlib.util.module_from_spec(spec)
_stdout = sys.stdout
class _Null:
    def write(self, *_): pass
    def flush(self): pass
sys.stdout = _Null()
try:
    spec.loader.exec_module(a3)          # it runs; we want its rows and helpers
except SystemExit:
    pass
finally:
    sys.stdout = _stdout

rows = a3.rows
have = [r for r in rows if r["length"] is not None]
mwu_z = a3.mwu_z

def rank_biserial(x, y):
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

def by_decile(key):
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

def matched_11(key):
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

res = []
for key, label in (("length", "protein_length"), ("cys_count", "cysteine_count")):
    a = [r[key] for r in have if r["persulfidated"]]
    b = [r[key] for r in have if not r["persulfidated"]]
    r_before, z_before = rank_biserial(a, b), mwu_z(a, b)
    r_after, used = by_decile(key)
    z_strat, us, dr = a3.stratified_z(have, key)
    A, B = matched_11(key)
    r_m, z_m = rank_biserial(A, B), mwu_z(A, B)
    res += [
     (f"{label}_r_before_matching", round(r_before, 4), f"n_pos={len(a)}, n_neg={len(b)}"),
     (f"{label}_z_before_matching", round(z_before, 4), ""),
     (f"{label}_r_decile_weighted_after", round(r_after, 4), f"deciles used {used}"),
     (f"{label}_z_stratified_after", round(z_strat, 4), f"Stouffer over {us} deciles, {dr} dropped"),
     (f"{label}_r_1to1_matched", round(r_m, 4), f"pairs={len(A)}"),
     (f"{label}_z_1to1_matched", round(z_m, 4), "same pairs"),
    ]

d = dict((k, v) for k, v, _ in res)
verdicts = {}
for label in ("protein_length", "cysteine_count"):
    rb, ra, zm, zb = (d[f"{label}_r_before_matching"], d[f"{label}_r_decile_weighted_after"],
                      d[f"{label}_z_1to1_matched"], d[f"{label}_z_before_matching"])
    if abs(ra) < 0.05: v = "R2_effect_vanishes_within_strata"
    elif abs(abs(ra) - abs(rb)) <= 0.05 and abs(zm) < abs(zb): v = "R1_collapse_would_be_a_sample_size_artefact"
    else: v = "R3_report_both_claim_neither"
    verdicts[label] = v

with open(OUT, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["quantity", "value", "note"]); w.writerows(res)
    for k, v in verdicts.items(): w.writerow([f"verdict_{k}", v, "reading rule fixed before the run"])

AUD.write_text(json.dumps({
 "script": "scripts/artefact3_effect_sizes_and_matching_2026-09-19.py", "script_sha256": sha(__file__),
 "imported_module": str(MAIN.relative_to(ROOT)), "imported_module_sha256": sha(MAIN),
 "interpreter": sys.version.split()[0], "verdicts": verdicts,
 "what_this_does_not_establish": [
   "It does not show what the original round actually did; it shows what two different readings of "
   "'matching' produce on this data.",
   "Rank-biserial r is a probability difference, not a biological effect size."],
}, ensure_ascii=False, indent=1, sort_keys=True))
for k, v, n in res: print(f"  {k:<44} {v:>10}   {n}")
print("\nverdicts:", verdicts)
