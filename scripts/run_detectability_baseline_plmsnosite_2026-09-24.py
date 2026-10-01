"""Detectability-only baseline on the pLMSNOSite S-nitrosylation benchmark.

Every rule is fixed in protocols/detectability_baseline_plmsnosite_preregistration_2026-09-24.json, registered
(SHA-256 in results/detectability_baseline_plmsnosite_registered_2026-09-24.sha256) before this script existed. G0
refuses if the protocol or the data files changed. Features and model configuration are imported unchanged from
scripts/run_cross_protease_detectability_probe.py (FEATURE_NAMES, features(), boundaries_for(), PROTEASES,
HGB_KWARGS); nothing is tuned.

GATES
  G0 protocol and data SHA-256 as registered
  G1 every row's Position is a cysteine; 25 features per row; no missing value
  NC label-permuted model: test AUROC interval covers 0.5 (else nothing is used)
Interpreter: /Users/lyuguohao/opt/persulfidation_v2/venv/bin/python (Python 3.9.6, scikit-learn 1.6.1).
"""
import csv
import hashlib
import importlib.util
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
PROTOCOL = ROOT / "protocols/detectability_baseline_plmsnosite_preregistration_2026-09-24.json"
REG = ROOT / "results/detectability_baseline_plmsnosite_registered_2026-09-24.sha256"
DATA = ROOT / "external/benchmarks_2026-09-24/pLMSNOSite"
PROBE = ROOT / "scripts/run_cross_protease_detectability_probe.py"
OUT = ROOT / "results/detectability_baseline_plmsnosite_2026-09-24.csv"
AUD = ROOT / "results/detectability_baseline_plmsnosite_2026-09-24_audit.json"
SEED, BOOT = 20260924, 5000
PUBLISHED = {"pLMSNOSite": 0.754, "PreSNO": 0.756, "DeepNitro": 0.731}


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def load(name):
    return list(csv.DictReader(open(DATA / name, encoding="utf-8")))


def main():
    if AUD.exists():
        sys.exit("REFUSE: audit exists")
    proto = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if sha(PROTOCOL) != REG.read_text().split()[0]:
        sys.exit("REFUSE: G0 protocol changed after registration")
    for fn, h in proto["data"]["sha256"].items():
        if sha(DATA / fn) != h:
            sys.exit("REFUSE: G0 data file %s changed" % fn)
    spec = importlib.util.spec_from_file_location("probe", PROBE)
    P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import average_precision_score, roc_auc_score
    from sklearn.model_selection import GroupKFold

    rule = P.PROTEASES["Trypsin"]
    names = list(P.FEATURE_NAMES)
    vis = [i for i, n in enumerate(names) if n in P.VISIBILITY_ONLY]

    def featurise(rows):
        X, y, g, cache = [], [], [], {}
        for r in rows:
            seq, pos = r["sequences"], int(r["Position"])
            if seq[pos - 1] != "C":
                sys.exit("REFUSE: G1 non-cysteine position %s:%d" % (r["UniProt"], pos))
            if seq not in cache:
                cache[seq] = P.boundaries_for(seq, rule)
            X.append(P.features(seq, pos - 1, cache[seq], rule))
            y.append(int(r["Target"]))
            g.append(r["UniProt"])
        X = np.asarray(X, dtype=float)
        if X.shape[1] != 25 or not np.isfinite(X).all():
            sys.exit("REFUSE: G1 feature matrix malformed")
        return X, np.asarray(y), np.asarray(g)

    Xtr, ytr, gtr = featurise(load("sequence_train.csv"))
    Xte, yte, gte = featurise(load("sequence_test.csv"))
    shared = sorted(set(gtr) & set(gte))

    def fit(X, y):
        return HistGradientBoostingClassifier(**P.HGB_KWARGS).fit(X, y)

    rng = np.random.default_rng(SEED)
    uniq = np.unique(gte)
    idx_of = {p: np.flatnonzero(gte == p) for p in uniq}

    def boot_ci(scores):
        vals = []
        for _ in range(BOOT):
            pick = rng.choice(uniq, size=len(uniq), replace=True)
            ii = np.concatenate([idx_of[p] for p in pick])
            if len(np.unique(yte[ii])) < 2:
                continue
            vals.append(roc_auc_score(yte[ii], scores[ii]))
        lo, hi = np.percentile(vals, [2.5, 97.5])
        return float(lo), float(hi), len(vals)

    results = {}
    for label, cols in (("DIG25", list(range(25))), ("VIS10", vis)):
        m = fit(Xtr[:, cols], ytr)
        s = m.predict_proba(Xte[:, cols])[:, 1]
        auc = float(roc_auc_score(yte, s))
        lo, hi, nb = boot_ci(s)
        results[label] = dict(n_features=len(cols), test_auroc=round(auc, 4), ci_low=round(lo, 4), ci_high=round(hi, 4),
                              boot_replicates_used=nb, test_auprc=round(float(average_precision_score(yte, s)), 4))
    # NC: label-permuted training
    yperm = np.random.default_rng(SEED).permutation(ytr)
    snc = fit(Xtr, yperm).predict_proba(Xte)[:, 1]
    nc_auc = float(roc_auc_score(yte, snc))
    nc_lo, nc_hi, _ = boot_ci(snc)
    nc_ok = nc_lo <= 0.5 <= nc_hi
    # PC: 5-fold CV on the training set, grouped by protein
    cv = []
    for tr, va in GroupKFold(n_splits=5).split(Xtr, ytr, gtr):
        cv.append(roc_auc_score(ytr[va], fit(Xtr[tr], ytr[tr]).predict_proba(Xtr[va])[:, 1]))
    d = results["DIG25"]
    R = {k: round((d["test_auroc"] - 0.5) / (v - 0.5), 4) for k, v in PUBLISHED.items()}
    r_main = R["pLMSNOSite"]
    if not nc_ok:
        branch = "NC_FAILED_no_result_used"
    elif r_main >= 0.75 and d["ci_low"] > 0.5:
        branch = "B1"
    elif 0.25 <= r_main < 0.75 and d["ci_low"] > 0.5:
        branch = "B2"
    else:
        branch = "B3"
    rows = [dict(feature_set=k, **v) for k, v in results.items()]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    AUD.write_text(json.dumps(dict(
        script=pathlib.Path(__file__).name, script_sha256=sha(__file__), protocol=PROTOCOL.name,
        protocol_sha256=sha(PROTOCOL), probe_sha256=sha(PROBE), n_train=int(len(ytr)), n_test=int(len(yte)),
        train_pos=int(ytr.sum()), test_pos=int(yte.sum()), proteins_train=int(len(set(gtr))),
        proteins_test=int(len(uniq)), proteins_shared_train_test=len(shared), results=results,
        recovery_ratio=R, branch=branch,
        NC_label_permuted=dict(test_auroc=round(nc_auc, 4), ci=[round(nc_lo, 4), round(nc_hi, 4)], covers_half=nc_ok),
        PC_train_cv_auroc=dict(folds=[round(float(x), 4) for x in cv], mean=round(float(np.mean(cv)), 4)),
        published_reference=PUBLISHED, all_pass=nc_ok), indent=2), encoding="utf-8")
    print(json.dumps(dict(results=results, recovery_ratio=R, branch=branch,
                          NC=dict(auc=round(nc_auc, 4), ci=[round(nc_lo, 4), round(nc_hi, 4)]),
                          PC_cv_mean=round(float(np.mean(cv)), 4), proteins_shared=len(shared)), indent=1))


if __name__ == "__main__":
    main()
