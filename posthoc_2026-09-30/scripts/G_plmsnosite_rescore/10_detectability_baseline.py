"""G_plmsnosite_rescore, step 1 (system python): reproduce the VIS10 / DIG25 detectability-only
baselines of Supplemental Note 5 on the pLMSNOSite benchmark and write per-site test scores.

POST HOC revision analysis (2026-09-30). The original analysis script and its registration JSON
(detectability_baseline_plmsnosite_preregistration_2026-09-24.json) are not present in the
repository, so the baseline is re-implemented from Supplemental Note 5: the 25 tryptic-digest
features and the fixed HistGradientBoosting configuration are taken, by import, from
repo/scripts/run_cross_protease_detectability_probe.py ("the fixed configuration of the
cross-protease detectability probe", Note 5). Nothing is tuned to hit the manuscript numbers;
where a number does not reproduce, both values are reported.

Outputs (results/G_plmsnosite_rescore/):
  detectability_test_scores.csv        per-site VIS10 / DIG25 test scores (full precision)
  detectability_test_features.csv      the 25 features of every test site
  detectability_reproduction.csv       reproduced vs manuscript values
  detectability_reproduction.json      details (bootstrap RNG variants, controls, checks)
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")

import importlib.util
import json
import sys
import time

import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

REPO_PROBE = r"/path/to/local/_cys_repo_work/repo/scripts/run_cross_protease_detectability_probe.py"
SEED_REPRO = 20260924   # bootstrap seed stated in Supplemental Note 5
REPS = 5000


def load_probe_module():
    sys.dont_write_bytecode = True          # never write __pycache__ into the read-only repo
    spec = importlib.util.spec_from_file_location("xprobe", REPO_PROBE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_matrix(df, probe):
    rule = probe.PROTEASES["Trypsin"]
    rows, bad = [], 0
    cache = {}
    for uid, pos, seq in zip(df["UniProt"], df["Position"], df["sequences"]):
        idx = int(pos) - 1
        if seq[idx] != "C":
            bad += 1
        key = (uid, seq)
        if key not in cache:
            cache[key] = probe.boundaries_for(seq, rule)
        rows.append(probe.features(seq, idx, cache[key], rule))
    return np.asarray(rows, dtype=np.float64), bad


def boot_ci(y, scores, codes, n_clusters, seed, mode, reps=REPS):
    rs = C.RankedScore(scores, y)
    au, ap = [], []
    for m in C.multiplicities(n_clusters, reps, seed, mode):
        w = m[codes]
        au.append(rs.auroc(w))
        ap.append(rs.auprc(w))
    return C.pctl(au), C.pctl(ap)


def main():
    t0 = time.time()
    C.OUT.mkdir(parents=True, exist_ok=True)
    manifest = {}
    head = C.git_head()
    assert head == C.EXPECTED_COMMIT, head
    train = C.read_blob_csv("data/train/sequence_train.csv", manifest)
    test = C.read_blob_csv("data/test/sequence_test.csv", manifest)
    manifest["repo/scripts/run_cross_protease_detectability_probe.py"] = C.sha256_file(REPO_PROBE)
    probe = load_probe_module()

    names = list(probe.FEATURE_NAMES)
    vis_idx = [i for i, n in enumerate(names) if n in probe.VISIBILITY_ONLY]
    x_tr, bad_tr = build_matrix(train, probe)
    x_te, bad_te = build_matrix(test, probe)
    y_tr = train["Target"].to_numpy().astype(int)
    y_te = test["Target"].to_numpy().astype(int)

    checks = {
        "git_head": head,
        "n_train": int(len(train)), "n_train_pos": int(y_tr.sum()), "n_train_proteins": int(train.UniProt.nunique()),
        "n_test": int(len(test)), "n_test_pos": int(y_te.sum()), "n_test_proteins": int(test.UniProt.nunique()),
        "proteins_in_both": int(len(set(train.UniProt) & set(test.UniProt))),
        "non_cysteine_positions_train": int(bad_tr), "non_cysteine_positions_test": int(bad_te),
        "hgb_kwargs": probe.HGB_KWARGS,
        "feature_names_DIG25": names,
        "feature_names_VIS10": [names[i] for i in vis_idx],
        "sklearn": sklearn.__version__, "numpy": np.__version__, "pandas": pd.__version__,
        "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
    }
    print(json.dumps({k: v for k, v in checks.items() if not k.startswith("feature_names")}, indent=1))

    def fit_predict(cols, labels=y_tr, random_state=None):
        kw = dict(probe.HGB_KWARGS)
        if random_state is not None:
            kw["random_state"] = random_state
        model = HistGradientBoostingClassifier(**kw)
        model.fit(x_tr[:, cols], labels)
        return model.predict_proba(x_te[:, cols])[:, 1]

    all_idx = list(range(len(names)))
    s_vis = fit_predict(vis_idx)
    s_dig = fit_predict(all_idx)
    # random_state is inert here (no early stopping, < 2e5 rows): verify.
    s_dig_rs = fit_predict(all_idx, random_state=SEED_REPRO)
    s_vis_rs = fit_predict(vis_idx, random_state=SEED_REPRO)
    checks["random_state_inert_DIG25_max_abs_diff"] = float(np.max(np.abs(s_dig - s_dig_rs)))
    checks["random_state_inert_VIS10_max_abs_diff"] = float(np.max(np.abs(s_vis - s_vis_rs)))

    codes, labels = C.protein_codes(test["UniProt"].tolist())
    K = len(labels)
    res = {}
    for name, s in (("VIS10", s_vis), ("DIG25", s_dig)):
        rs = C.RankedScore(s, y_te)
        res[name] = {
            "auroc_sklearn": float(roc_auc_score(y_te, s)), "auprc_sklearn": float(average_precision_score(y_te, s)),
            "auroc_weighted_impl": rs.auroc(), "auprc_weighted_impl": rs.auprc(),
            "n_distinct_scores": int(len(np.unique(s))),
        }
        for mode in ("per_replicate", "block", "chunk250"):
            (alo, ahi), (plo, phi) = boot_ci(y_te, s, codes, K, SEED_REPRO, mode)
            res[name][f"ci_seed{SEED_REPRO}_{mode}"] = {"auroc": [alo, ahi], "auprc": [plo, phi]}
        print(name, json.dumps(res[name], indent=1))

    # Negative control (label-permuted DIG25). The permutation RNG of the original run is unknown;
    # one natural choice is used and reported without trying alternatives.
    perm = np.random.default_rng(SEED_REPRO).permutation(y_tr)
    s_nc = fit_predict(all_idx, labels=perm)
    (nlo, nhi), _ = boot_ci(y_te, s_nc, codes, K, SEED_REPRO, "per_replicate")
    res["DIG25_label_permuted"] = {"permutation": f"numpy default_rng({SEED_REPRO}).permutation(y_train)",
                                   "auroc": float(roc_auc_score(y_te, s_nc)), "ci": [nlo, nhi]}

    # Positive control: five-fold CV on training, grouped by protein (GroupKFold, no shuffle).
    folds = []
    for tr_i, va_i in GroupKFold(n_splits=5).split(x_tr, y_tr, groups=train["UniProt"].to_numpy()):
        mdl = HistGradientBoostingClassifier(**probe.HGB_KWARGS).fit(x_tr[tr_i], y_tr[tr_i])
        folds.append(float(roc_auc_score(y_tr[va_i], mdl.predict_proba(x_tr[va_i])[:, 1])))
    res["DIG25_groupkfold5_train"] = {"fold_auroc": folds, "mean": float(np.mean(folds))}
    print("NC", res["DIG25_label_permuted"], "\nPC", res["DIG25_groupkfold5_train"])

    # ---- outputs
    scores = pd.DataFrame({
        "row": np.arange(len(test)), "UniProt": test["UniProt"], "Position": test["Position"],
        "Target": y_te, "vis10_score": s_vis, "dig25_score": s_dig,
    })
    scores.to_csv(C.OUT / "detectability_test_scores.csv", index=False, float_format="%.17g")
    feats = pd.DataFrame(x_te, columns=names)
    feats.insert(0, "Position", test["Position"].to_numpy())
    feats.insert(0, "UniProt", test["UniProt"].to_numpy())
    feats.to_csv(C.OUT / "detectability_test_features.csv", index=False, float_format="%.17g")

    M = C.MANUSCRIPT
    rows = []
    for name in ("VIS10", "DIG25"):
        r = res[name]
        ci = r[f"ci_seed{SEED_REPRO}_per_replicate"]["auroc"]
        rows += [
            {"quantity": f"{name} test AUROC", "reproduced": r["auroc_sklearn"], "manuscript": M[f"{name}_auroc"],
             "match_4dp": round(r["auroc_sklearn"], 4) == M[f"{name}_auroc"]},
            {"quantity": f"{name} test AUPRC", "reproduced": r["auprc_sklearn"], "manuscript": M[f"{name}_auprc"],
             "match_4dp": round(r["auprc_sklearn"], 4) == M[f"{name}_auprc"]},
            {"quantity": f"{name} AUROC CI low (seed {SEED_REPRO}, per-replicate draws)", "reproduced": ci[0],
             "manuscript": M[f"{name}_ci"][0], "match_4dp": round(ci[0], 4) == M[f"{name}_ci"][0]},
            {"quantity": f"{name} AUROC CI high (seed {SEED_REPRO}, per-replicate draws)", "reproduced": ci[1],
             "manuscript": M[f"{name}_ci"][1], "match_4dp": round(ci[1], 4) == M[f"{name}_ci"][1]},
        ]
    nc = res["DIG25_label_permuted"]
    rows.append({"quantity": "DIG25 label-permuted test AUROC", "reproduced": nc["auroc"], "manuscript": M["NC_auroc"],
                 "match_4dp": round(nc["auroc"], 4) == M["NC_auroc"]})
    for k, (a, b) in enumerate(zip(res["DIG25_groupkfold5_train"]["fold_auroc"], M["PC_folds"])):
        rows.append({"quantity": f"DIG25 grouped 5-fold CV AUROC fold {k + 1}", "reproduced": a, "manuscript": b,
                     "match_4dp": round(a, 4) == b})
    pd.DataFrame(rows).to_csv(C.OUT / "detectability_reproduction.csv", index=False, float_format="%.17g")
    (C.OUT / "detectability_reproduction.json").write_text(json.dumps(
        {"label": "POST HOC revision analysis 2026-09-30 (G_plmsnosite_rescore); not registered",
         "checks": checks, "results": res, "inputs_sha256": manifest,
         "runtime_s": round(time.time() - t0, 1)}, indent=1), encoding="utf-8")
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main()
