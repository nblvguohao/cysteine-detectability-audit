"""The decisive test: run Artefact 3 on the ORIGINAL round's own positive set.

WHAT WAS RECONSTRUCTED, AND HOW EXACTLY
  The original paragraph reported "479 versus 6,537 proteins". The positive arm is reproduced here
  EXACTLY: the 640 proteins whose persulfidation sites were inferred from single-cysteine peptides
  (inputs/PXD072089_mapped_rice.json, built in this project on 2026-09-12), intersected with the
  7,692 protein groups carrying a usable iBAQ in PXD072035, gives 479. Not approximately 479: 479.

  The control arm was NOT reproduced. Candidates measured, none equal to 6,537:
     quantified minus single-cysteine set          7,213
     quantified minus leading-razor set            6,495   <- closest, 42 away
     quantified minus all-proteins-column set      5,692
  The control arm's definition is therefore unknown and the test below is run against BOTH of the two
  defensible ones, not against a guess tuned to hit 6,537.

CRITERIA, WRITTEN BEFORE THE RUN
  Same statistics as the main round: median iBAQ ratio, and rank-biserial r for protein length and
  cysteine count, before matching and decile-weighted after. Rules imported, not reimplemented.

  READING RULE, fixed before the run:
    T1 with the original positive set, |r| after matching falls below 0.05 on both features
       -> the original claim reproduces, and the main round's disagreement is explained by its wider
          positive arm. The manuscript's Artefact 3 stands, scoped to that positive set.
    T2 |r| after matching stays at or above 0.05 on either feature, under BOTH control arms
       -> the claim does not reproduce even on its own positive set. The manuscript's Artefact 3
          sentence is not supported by this data.
    T3 the two control arms disagree with each other -> report both, conclude nothing.

This script computes. It does not edit the manuscript.
"""
import csv, hashlib, importlib.util, json, pathlib, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
MAIN = ROOT / "scripts/recompute_artefact3_osa_abundance_2026-09-19.py"
EFF = ROOT / "scripts/artefact3_effect_sizes_and_matching_2026-09-19.py"
D = ROOT / "external/pride_osa_abundance_20260919"
OUT = ROOT / "results/artefact3_reconstructed_original_arms_2026-09-19.csv"
AUD = ROOT / "results/artefact3_reconstructed_original_arms_2026-09-19_audit.json"
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
csv.field_size_limit(10**7)

class _N:
    def write(self, *_): pass
    def flush(self): pass
def load(p, n):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s)
    o = sys.stdout; sys.stdout = _N()
    try: s.loader.exec_module(m)
    except SystemExit: pass
    finally: sys.stdout = o
    return m
a3 = load(MAIN, "a3"); eff = load(EFF, "a3eff")
rb = eff.rank_biserial

single = {p["acc"] for p in json.load(open(ROOT / "inputs/PXD072089_mapped_rice.json"))
          ["osa_persulfidation"]["proteins"]}
with open(D / "SS-all-peptides.tsv") as f:
    peps = list(csv.DictReader(f, delimiter="\t"))
defA = {(p.get("Leading razor protein") or "").strip() for p in peps} - {""}
with open(D / "20240326_061934_xyj_proteome_Report.tsv") as f:
    prot = list(csv.DictReader(f, delimiter="\t"))
accmap = {}
for p in prot:
    a = [x.strip() for x in (p.get("PG.ProteinAccessions") or "").split(";") if x.strip()]
    if a: accmap[a[0]] = a

base = a3.rows
for r in base:
    r["_accs"] = accmap.get(r["group_id"], [r["group_id"]])
    r["_is_single"] = any(a in single for a in r["_accs"])
    r["_is_defA"] = any(a in defA for a in r["_accs"])

POS = [r for r in base if r["_is_single"]]
CTRLS = {"control_quant_minus_single_cys_set": [r for r in base if not r["_is_single"]],
         "control_quant_minus_leading_razor_set": [r for r in base if not r["_is_defA"]]}

def stats(pos, neg):
    out = {"n_pos": len(pos), "n_neg": len(neg)}
    out["median_iBAQ_ratio"] = float(np.median([r["ibaq_median"] for r in pos]) /
                                     np.median([r["ibaq_median"] for r in neg]))
    hp = [r for r in pos if r["length"] is not None]
    hn = [r for r in neg if r["length"] is not None]
    allr = hp + hn
    vals = np.array([r["ibaq_median"] for r in allr], float)
    edges = np.percentile(vals, np.linspace(0, 100, 11)); edges[0] -= 1e-9; edges[-1] += 1e-9
    for key, lab in (("length", "length"), ("cys_count", "cys")):
        out[f"r_before_{lab}"] = rb([r[key] for r in hp], [r[key] for r in hn])
        rs, used = [], 0
        for k in range(10):
            a = [r[key] for r in hp if edges[k] < r["ibaq_median"] <= edges[k + 1]]
            b = [r[key] for r in hn if edges[k] < r["ibaq_median"] <= edges[k + 1]]
            if len(a) < 10 or len(b) < 10: continue
            rs.append(rb(a, b)); used += 1
        out[f"r_after_{lab}"] = float(np.mean(rs)) if rs else float("nan")
        out[f"deciles_{lab}"] = used
    return out

res = {name: stats(POS, neg) for name, neg in CTRLS.items()}
vanish = all(abs(v[f"r_after_{k}"]) < 0.05 for v in res.values() for k in ("length", "cys"))
persist = all(abs(v[f"r_after_{k}"]) >= 0.05 for v in res.values() for k in ("length", "cys"))
verdict = ("T1_original_claim_reproduces" if vanish else
           "T2_claim_does_not_reproduce_on_its_own_positive_set" if persist else
           "T3_control_arms_disagree")

hdr = ["control_arm", "n_pos", "n_neg", "median_iBAQ_ratio", "r_before_length", "r_after_length",
       "r_before_cys", "r_after_cys", "deciles_used"]
rows = [hdr]
for name, v in res.items():
    rows.append([name, v["n_pos"], v["n_neg"], round(v["median_iBAQ_ratio"], 4),
                 round(v["r_before_length"], 4), round(v["r_after_length"], 4),
                 round(v["r_before_cys"], 4), round(v["r_after_cys"], 4), v["deciles_length"]])
rows.append(["verdict", verdict] + [""] * 7)
with open(OUT, "w", newline="") as f: csv.writer(f).writerows(rows)

AUD.write_text(json.dumps({
 "script": "scripts/artefact3_reconstructed_original_arms_2026-09-19.py", "script_sha256": sha(__file__),
 "imported": {str(MAIN.relative_to(ROOT)): sha(MAIN), str(EFF.relative_to(ROOT)): sha(EFF)},
 "positive_arm_reconstruction": {
   "rule": "single-cysteine-peptide protein set (640) intersect quantified groups (7,692)",
   "n": len(POS), "matches_original_479": len(POS) == 479},
 "control_arm_reconstruction": {
   "reproduced": False, "original_reported": 6537,
   "candidates": {"quant_minus_single": 7213, "quant_minus_leading_razor": 6495,
                  "quant_minus_all_proteins_column": 5692}},
 "results": {k: {kk: (None if isinstance(vv, float) and np.isnan(vv) else vv)
                 for kk, vv in v.items()} for k, v in res.items()},
 "verdict": verdict,
 "what_this_does_not_establish": [
   "Reproducing 479 exactly identifies the positive arm's RULE; it does not prove the original round "
   "used that rule, only that this rule lands on the same number.",
   "The control arm was not reproduced, so a like-for-like replication of the original statistic is "
   "still not possible.",
   "It says nothing about the annotation-term and keyword parts of the same paragraph."],
}, ensure_ascii=False, indent=1, sort_keys=True))
for r in rows: print("  " + "  ".join(str(x).rjust(10) for x in r))
print("\nverdict:", verdict)
