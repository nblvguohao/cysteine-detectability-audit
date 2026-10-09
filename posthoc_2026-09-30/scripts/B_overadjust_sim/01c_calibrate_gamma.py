# -*- coding: utf-8 -*-
"""B_overadjust_sim step 1c (POST HOC): choose the selection steepness gamma of the 'steep' sensitivity
model. 40 null tables per gamma (point estimates, seed scenario 'calib'); the propensity AUC of the claim
label on VIS10 is compared with the 0.72-0.76 of the published site-level re-tests. gamma = 2.0 was taken
(AUC about 0.72); larger values were not explored because the AUC rises slowly and the model is a
stress test, not a fit. Writes results/B_overadjust_sim/gamma_calibration.csv."""
import os
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402
import pandas as pd  # noqa: E402
import sim_engine as se  # noqa: E402

rows = []
for g in [0.0, 0.5, 1.0, 1.5, 2.0]:
    for r in range(40):
        rows += se.run_replicate(dict(scenario="S0steep", seed_scenario="calib", detection_mode="steep", gamma=g,
                                      attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=False))
df = pd.DataFrame(rows)
out = df.groupby(["gamma", "attribute"])[["prop_auc_before", "n_obs", "n_pos", "pos_rate", "baseline", "matched"]].mean().reset_index()
out.to_csv(os.path.join(cb.RESULTS, "gamma_calibration.csv"), index=False)
print(out.groupby("gamma")[["prop_auc_before", "pos_rate"]].mean())
