"""G_plmsnosite_rescore, step 1b (system python): diagnose why the VIS10 and label-permuted
intervals and the grouped-CV fold AUROCs of Supplemental Note 5 are not reproduced exactly,
although all point estimates and the DIG25 interval are.

POST HOC revision diagnostic (2026-09-30). A small, pre-listed set of plausible implementation
variants is tried once each and every outcome is reported (see reproduction_diagnostics.json);
none is used to change a reported analysis.
Hypotheses
  H1  one bootstrap generator (seed 20260924) was shared by the models, run in sequence
      (orders tried: DIG25->VIS10->NC, DIG25->NC->VIS10, VIS10->DIG25->NC)
  H2  grouped CV used shuffling: GroupKFold(shuffle=True, random_state=20260924) or
      StratifiedGroupKFold(shuffle=True, random_state=20260924) or StratifiedGroupKFold(no shuffle)
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
import importlib.util
import itertools
import json
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

SEED = 20260924
REPS = 5000
REPO_PROBE = r"/path/to/local/_cys_repo_work/repo/scripts/run_cross_protease_detectability_probe.py"


def main():
    sc = pd.read_csv(C.OUT / "detectability_test_scores.csv")
    y = sc["Target"].to_numpy()
    codes, labels = C.protein_codes(sc["UniProt"].tolist())
    K = len(labels)
    # label-permuted DIG25 scores (same construction as step 1)
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("xprobe", REPO_PROBE)
    probe = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(probe)
    train = C.read_blob_csv("data/train/sequence_train.csv")
    test = C.read_blob_csv("data/test/sequence_test.csv")
    rule = probe.PROTEASES["Trypsin"]

    def mat(df):
        out = []
        for seq, pos in zip(df["sequences"], df["Position"]):
            out.append(probe.features(seq, int(pos) - 1, probe.boundaries_for(seq, rule), rule))
        return np.asarray(out, float)

    x_tr, x_te = mat(train), mat(test)
    y_tr = train["Target"].to_numpy()
    perm = np.random.default_rng(SEED).permutation(y_tr)
    s_nc = HistGradientBoostingClassifier(**probe.HGB_KWARGS).fit(x_tr, perm).predict_proba(x_te)[:, 1]
    scores = {"VIS10": sc["vis10_score"].to_numpy(), "DIG25": sc["dig25_score"].to_numpy(), "NC": s_nc}
    target = {"VIS10": C.MANUSCRIPT["VIS10_ci"], "DIG25": C.MANUSCRIPT["DIG25_ci"], "NC": C.MANUSCRIPT["NC_ci"]}
    ranked = {k: C.RankedScore(v, y) for k, v in scores.items()}

    out = {"label": "POST HOC reproduction diagnostic; not registered", "H1_shared_generator": {}, "H2_cv_variants": {}}
    for order in (("DIG25", "VIS10", "NC"), ("DIG25", "NC", "VIS10"), ("VIS10", "DIG25", "NC")):
        rng = np.random.default_rng(SEED)
        got = {}
        for name in order:
            vals = []
            for _ in range(REPS):
                w = np.bincount(rng.integers(0, K, K), minlength=K)[codes]
                vals.append(ranked[name].auroc(w))
            got[name] = [round(v, 4) for v in C.pctl(vals)]
        out["H1_shared_generator"]["->".join(order)] = {
            n: {"ci": got[n], "manuscript": list(target[n]), "match": tuple(got[n]) == tuple(target[n])} for n in order}
    # H2: fold schemes
    groups = train["UniProt"].to_numpy()
    schemes = {
        "GroupKFold(shuffle=True, rs=20260924)": GroupKFold(n_splits=5, shuffle=True, random_state=SEED),
        "StratifiedGroupKFold(shuffle=True, rs=20260924)": StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED),
        "StratifiedGroupKFold(no shuffle)": StratifiedGroupKFold(n_splits=5),
    }
    for name, cv in schemes.items():
        folds = []
        for a, b in cv.split(x_tr, y_tr, groups=groups):
            m = HistGradientBoostingClassifier(**probe.HGB_KWARGS).fit(x_tr[a], y_tr[a])
            folds.append(round(float(roc_auc_score(y_tr[b], m.predict_proba(x_tr[b])[:, 1])), 4))
        out["H2_cv_variants"][name] = {"folds": folds, "mean": round(float(np.mean(folds)), 4),
                                       "manuscript": list(C.MANUSCRIPT["PC_folds"]),
                                       "match": tuple(folds) == tuple(C.MANUSCRIPT["PC_folds"])}
    (C.OUT / "reproduction_diagnostics.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
