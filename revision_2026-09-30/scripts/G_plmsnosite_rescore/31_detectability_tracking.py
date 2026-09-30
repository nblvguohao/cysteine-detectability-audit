"""G_plmsnosite_rescore, step 3b (system python): does the pLMSNOSite score track peptide
detectability within a label class?

POST HOC revision analysis (2026-09-30); not registered. Additional to tasks 4(a)-(d).
Within negatives and within positives separately (so the label cannot drive the association), the
AUROC with which each pLMSNOSite score separates cysteines on a theoretically detectable tryptic
peptide from the rest is computed, for two indicators from the DIG25 block:
  pep_detectable_both                 fully cleaved peptide, length 7-30 and mass 700-3,500 Da
  pep_detectable_any_missed_cleavage  the same window reached with up to two missed cleavages
AUROC 0.5 = no tracking. Protein-clustered bootstrap (5,000 resamples, seed 20260930).
Also reported: base rates of the indicators among positives and negatives.
Output: results/G_plmsnosite_rescore/detectability_tracking.csv
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

SEED, REPS = 20260930, 5000


def main():
    feats = pd.read_csv(C.OUT / "detectability_test_features.csv")
    plm = pd.read_csv(C.OUT / "plmsnosite_test_scores.csv")
    det = pd.read_csv(C.OUT / "detectability_test_scores.csv")
    assert (feats.UniProt.values == plm.UniProt.values).all() and (feats.Position.values == plm.Position.values).all()
    y = plm.Target.to_numpy()
    scores = {"pLMSNOSite": plm.plmsnosite_prob.to_numpy(), "ProtT5_arm": plm.prott5_base_prob.to_numpy(),
              "Embedding_arm": plm.embedding_base_prob.to_numpy(), "DIG25": det.dig25_score.to_numpy(),
              "VIS10": det.vis10_score.to_numpy()}
    rows = []
    for ind in ("pep_detectable_both", "pep_detectable_any_missed_cleavage"):
        z = feats[ind].to_numpy().astype(int)
        rows.append({"analysis": "base rate", "indicator": ind, "class": "positives", "score": "",
                     "n": int((y == 1).sum()), "n_indicator_1": int(z[y == 1].sum()), "value": float(z[y == 1].mean())})
        rows.append({"analysis": "base rate", "indicator": ind, "class": "negatives", "score": "",
                     "n": int((y == 0).sum()), "n_indicator_1": int(z[y == 0].sum()), "value": float(z[y == 0].mean())})
        for cls in (0, 1):
            sel = y == cls
            sub_codes, labels = C.protein_codes(plm.UniProt[sel].tolist())
            K = len(labels)
            for name, s in scores.items():
                rs = C.RankedScore(s[sel], z[sel])
                vals = [rs.auroc(m[sub_codes]) for m in C.multiplicities(K, REPS, SEED)]
                lo, hi = C.pctl(vals)
                rows.append({"analysis": "AUROC of score for indicator, within class", "indicator": ind,
                             "class": "positives" if cls else "negatives", "score": name, "n": int(sel.sum()),
                             "n_indicator_1": int(z[sel].sum()), "value": rs.auroc(), "ci_low": lo, "ci_high": hi,
                             "n_proteins": K,
                             "mean_score_indicator_1": float(s[sel][z[sel] == 1].mean()),
                             "mean_score_indicator_0": float(s[sel][z[sel] == 0].mean())})
    out = pd.DataFrame(rows)
    out.to_csv(C.OUT / "detectability_tracking.csv", index=False, float_format="%.10g")
    pd.set_option("display.width", 250)
    print(out.round(4).to_string())


if __name__ == "__main__":
    main()
