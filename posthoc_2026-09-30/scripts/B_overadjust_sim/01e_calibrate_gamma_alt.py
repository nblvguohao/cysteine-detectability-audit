# -*- coding: utf-8 -*-
"""B_overadjust_sim step 1e (POST HOC, revision after adversarial verification): selection steepness of
the two alternative selection models added in the revision.

  steep_vis     positives selected on the VIS10 detection-model logit (VIS10 represents the selection)
  steep_nocomp  positives selected on the peptide detection model without the basic/acidic fractions

Rule (fixed before the runs below were summarised): for each model take the gamma in {1, 2, 3, 4, 6}
whose mean propensity AUC (claim label on VIS10, null tables) is closest to that of the original steep
model (mean over its 300 null tables in sim_summary.csv, 0.722), so that the three selection models are
compared at the same selection strength as the matching sees it. 40 null tables per gamma, point
estimates, seed scenario 'calibalt'. Writes results/B_overadjust_sim/gamma_calibration_alt.csv and
gamma_choice_alt.json.
Amendment (post hoc, after the first run): for steep_vis the first grid did not bracket the target
(AUC 0.784 already at gamma 1), so gammas 0.25, 0.5 and 0.75 were added for that model; the rule
(closest AUC) is unchanged. Seeds depend only on the table number, so the first run's values are
reproduced exactly.
"""
import json
import multiprocessing as mp
import os
import sys

sys.dont_write_bytecode = True
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402
import pandas as pd  # noqa: E402

GAMMAS = {"steep_vis": [0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 6.0],   # 0.25-0.75 added (amendment)
          "steep_nocomp": [1.0, 2.0, 3.0, 4.0, 6.0]}
MODES = ["steep_vis", "steep_nocomp"]
N_TABLES = 40


def _work(task):
    import sim_engine as se
    return se.run_replicate(task)


def main():
    s = pd.read_csv(os.path.join(cb.RESULTS, "sim_summary.csv"))
    target = float(s[(s.scenario == "S0steep")]["mean_prop_auc_before"].iloc[0])
    tasks = []
    for mode in MODES:
        for g in GAMMAS[mode]:
            for r in range(N_TABLES):
                # seed depends on rep only: the same random numbers are reused across gamma and mode
                tasks.append(dict(scenario="S0" + mode, seed_scenario="calibalt", detection_mode=mode, gamma=g,
                                  attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=False))
    with mp.get_context("spawn").Pool(4) as pool:
        rows = [row for res in pool.imap_unordered(_work, tasks, chunksize=10) for row in res]
    df = pd.DataFrame(rows)
    out = (df.groupby(["detection_mode", "gamma", "attribute"])[["prop_auc_before", "n_obs", "n_pos", "pos_rate",
                                                                 "baseline", "matched"]]
           .mean().reset_index())
    out.to_csv(os.path.join(cb.RESULTS, "gamma_calibration_alt.csv"), index=False)
    auc = df.groupby(["detection_mode", "gamma"])["prop_auc_before"].mean().reset_index()
    choice = {"target_prop_auc_original_steep": target, "rule": "gamma closest to the target AUC", "chosen": {}}
    for mode in MODES:
        a = auc[auc.detection_mode == mode].copy()
        a["dist"] = (a["prop_auc_before"] - target).abs()
        best = a.sort_values("dist").iloc[0]
        choice["chosen"][mode] = {"gamma": float(best["gamma"]), "prop_auc": float(best["prop_auc_before"])}
        choice[mode + "_auc_by_gamma"] = dict(zip(a["gamma"].astype(float), a["prop_auc_before"].astype(float)))
    with open(os.path.join(cb.RESULTS, "gamma_choice_alt.json"), "w", encoding="utf-8") as fh:
        json.dump(choice, fh, indent=1)
    print(json.dumps(choice, indent=1))
    print(out[out.attribute.isin(["a1_SFE006_KR", "a3_SNO006_DE3", "a3b_SFE002_E"])].round(3).to_string())


if __name__ == "__main__":
    main()
