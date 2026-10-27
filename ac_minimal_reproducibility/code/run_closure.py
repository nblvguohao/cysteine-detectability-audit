"""Low-concurrency, resumable, independently versioned training execution."""
import argparse
import datetime
import hashlib
import json
import os
import platform
import sys
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from threadpoolctl import threadpool_limits

from feature_definitions import DEFINITIONS, FEATURE_NAMES, FEATURE_SETS, matrix
from prepare_inputs import COHORTS, ROOT, fasta, sha
from public_probe_snapshot import HGB_KWARGS


def write_json(path, obj):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    tmp.replace(path)


def event(kind, **values):
    record = dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(), kind=kind, **values)
    print(json.dumps(record, ensure_ascii=False), flush=True)
    with (ROOT / "logs/events.jsonl").open("a") as h:
        h.write(json.dumps(record, ensure_ascii=False) + "\n")


class WeightedAUC:
    """Pre-sort ties once; unit weights change for each paired bootstrap draw."""
    def __init__(self, y, score, units, mask):
        self.units = units[mask]
        self.y = y[mask].astype(bool)
        _, self.bins = np.unique(score[mask], return_inverse=True)
        self.nbins = int(self.bins.max()) + 1

    def __call__(self, unit_weights):
        weights = unit_weights[self.units]
        pos = np.bincount(self.bins[self.y], weights=weights[self.y], minlength=self.nbins)
        neg = np.bincount(self.bins[~self.y], weights=weights[~self.y], minlength=self.nbins)
        denominator = pos.sum() * neg.sum()
        if not denominator:
            return np.nan
        return float(np.dot(pos, np.cumsum(neg) - 0.5 * neg) / denominator)


def evaluate(sites, scores, reps):
    y = sites.is_pos.to_numpy()
    _, unit_indices = np.unique(sites.unit, return_inverse=True)
    n_units = unit_indices.max() + 1
    evals = {b: (y == 1) | sites["in_" + b].to_numpy().astype(bool) for b in ("A", "B")}
    calculators = {a+b: WeightedAUC(y, scores[a], unit_indices, evals[b])
                   for a in ("A", "B") for b in ("A", "B")}
    point = {k: v(np.ones(n_units)) for k, v in calculators.items()}
    for k, value in point.items():
        mask = evals[k[1]]
        assert abs(value - roc_auc_score(y[mask], scores[k[0]][mask])) < 1e-12
    point["delta_fix_A"] = point["AB"] - point["AA"]
    point["delta_fix_B"] = point["BB"] - point["BA"]
    point["training_effect_on_B"] = point["BB"] - point["AB"]
    rng = np.random.default_rng(20261013)
    draws = {key: [] for key in point}
    for _ in range(reps):
        counts = np.bincount(rng.integers(0, n_units, n_units), minlength=n_units)
        values = {key: calc(counts) for key, calc in calculators.items()}
        values["delta_fix_A"] = values["AB"] - values["AA"]
        values["delta_fix_B"] = values["BB"] - values["BA"]
        values["training_effect_on_B"] = values["BB"] - values["AB"]
        for key, value in values.items():
            draws[key].append(value)
    intervals = {}
    for key, values in draws.items():
        a = np.asarray(values)
        finite = a[np.isfinite(a)]
        assert len(finite) >= .95 * reps
        intervals[key] = dict(estimate=point[key], lower=float(np.quantile(finite, .025)),
                              upper=float(np.quantile(finite, .975)), valid_draws=len(finite),
                              invalid_draws=int(reps-len(finite)))
    return dict(metrics=intervals, bootstrap_reps=reps, seed=20261013,
                uncertainty="paired unit percentile bootstrap conditional on new fitted OOF scores")


def environment():
    return dict(python=sys.version, platform=platform.platform(), numpy=np.__version__,
                scipy=scipy.__version__, sklearn=sklearn.__version__, pandas=pd.__version__,
                joblib=joblib.__version__, threads={k:os.environ.get(k) for k in
                ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "CUDA_VISIBLE_DEVICES"]})


def fingerprint():
    paths = [ROOT / "PROTOCOL.md", ROOT / "run_low_priority.sh", ROOT / "requirements.txt",
             ROOT / "inputs/source_manifest.json", *sorted((ROOT / "code").glob("*.py")),
             *sorted((ROOT / "inputs/prepared").glob("*/*.tsv")),
             *sorted((ROOT / "inputs/prepared").glob("*/*.fasta"))]
    manifest = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    identity = dict(files=manifest, environment=environment(), model=HGB_KWARGS,
                    feature_sets=FEATURE_SETS, definitions=DEFINITIONS)
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    return digest, identity


def run_task(cohort, definition, feature_set, identity, reps):
    started = time.monotonic()
    folder = ROOT / "results" / definition / cohort / feature_set
    folder.mkdir(parents=True, exist_ok=True)
    complete = folder / "complete.json"
    task_metadata = dict(global_identity=identity, cohort=cohort, definition=definition,
                         feature_set=feature_set, bootstrap_reps=reps)
    task_identity = hashlib.sha256(json.dumps(task_metadata, sort_keys=True).encode()).hexdigest()
    required_outputs = {"features.tsv.gz", "feature_schema.json", "folds.tsv.gz", "new_oof.tsv.gz", "evaluation.json"}
    required_outputs |= {"fold_%d_%s.%s" % (k, b, ext)
                         for k in range(5) for b in ("A", "B") for ext in ("json", "npz", "joblib")}
    if complete.exists():
        record = json.loads(complete.read_text())
        assert record["identity"] == task_identity, "Checkpoint identity changed"
        assert record["task_metadata"] == task_metadata
        assert set(record["outputs"]) == required_outputs, "Incomplete output inventory"
        assert all(sha(folder / name) == digest for name, digest in record["outputs"].items())
        saved_evaluation = json.loads((folder / "evaluation.json").read_text())
        assert saved_evaluation["identity"] == task_identity
        assert all(saved_evaluation[k] == v for k, v in task_metadata.items() if k != "global_identity")
        event("task_reused", cohort=cohort, definition=definition, feature_set=feature_set)
        return
    inp = ROOT / "inputs/prepared" / cohort
    sites = pd.read_csv(inp / "sites.tsv", sep="\t")
    assert list(sites.columns) == ["acc", "site", "unit", "is_pos", "in_A", "in_B"]
    assert not sites.duplicated(["acc", "site"]).any()
    seqs = fasta(inp / "reference.fasta")
    X_all = matrix(sites[["acc", "site"]], seqs, definition)
    keep = [FEATURE_NAMES.index(name) for name in FEATURE_SETS[feature_set]]
    X = X_all[:, keep]
    assert np.isfinite(X).all() and X.shape == (len(sites), len(keep))
    pd.DataFrame(X, columns=FEATURE_SETS[feature_set]).to_csv(folder / "features.tsv.gz", sep="\t", index=False)
    nonstandard = {a: sorted(set(seq)-set("ACDEFGHIKLMNPQRSTVWY")) for a, seq in seqs.items()
                   if set(seq)-set("ACDEFGHIKLMNPQRSTVWY")}
    write_json(folder / "feature_schema.json", dict(definition=definition, features=FEATURE_SETS[feature_set],
               nonstandard_residues_by_accession=nonstandard, features_sha256=sha(folder / "features.tsv.gz")))
    y, groups = sites.is_pos.to_numpy(), sites.unit.to_numpy()
    splitter = GroupKFold(n_splits=5)
    splits = list(splitter.split(X, y, groups))
    folds = np.full(len(y), -1)
    for k, (_, test) in enumerate(splits):
        assert np.all(folds[test] == -1)
        folds[test] = k
    assert np.all(folds >= 0)
    membership = sites.copy()
    membership["new_fold"] = folds
    membership.to_csv(folder / "folds.tsv.gz", sep="\t", index=False)
    scores = {b: np.full(len(y), np.nan) for b in ("A", "B")}
    fold_audits = []
    event("task_start", cohort=cohort, definition=definition, feature_set=feature_set, rows=len(sites))
    for k, (train, test) in enumerate(splits):
        assert not set(groups[train]) & set(groups[test])
        for b in ("A", "B"):
            selected = train[(y[train] == 1) | sites["in_" + b].to_numpy()[train].astype(bool)]
            assert len(np.unique(y[selected])) == 2
            prefix = folder / ("fold_%d_%s" % (k, b))
            audit_path = prefix.with_suffix(".json")
            prediction_path = prefix.with_suffix(".npz")
            model_path = prefix.with_suffix(".joblib")
            if audit_path.exists():
                audit = json.loads(audit_path.read_text())
                assert audit["identity"] == task_identity
                assert audit["fold"] == k and audit["background"] == b
                assert audit["train_rows"] == len(selected) and audit["test_rows"] == len(test)
                assert audit["train_positive"] == int(y[selected].sum())
                assert audit["test_positive"] == int(y[test].sum())
                assert audit["train_test_unit_overlap"] == 0
                assert sha(prediction_path) == audit["prediction_sha256"]
                assert sha(model_path) == audit["model_sha256"]
                saved = np.load(prediction_path)
                assert np.array_equal(saved["test_rows"], test)
                assert np.array_equal(saved["train_rows"], selected)
                prediction = saved["scores"]
            else:
                t0 = time.monotonic()
                model = HistGradientBoostingClassifier(**HGB_KWARGS)
                with threadpool_limits(limits=2, user_api="openmp"):
                    model.fit(X[selected], y[selected])
                    prediction = model.predict_proba(X[test])[:, 1]
                joblib.dump(model, model_path, compress=3)
                np.savez_compressed(prediction_path, test_rows=test, train_rows=selected, scores=prediction)
                assert np.isfinite(prediction).all() and ((0 <= prediction) & (prediction <= 1)).all()
                audit = dict(identity=task_identity, fold=k, background=b,
                             train_rows=len(selected), train_positive=int(y[selected].sum()),
                             train_units=int(len(set(groups[selected]))), test_rows=len(test),
                             test_positive=int(y[test].sum()), test_units=int(len(set(groups[test]))),
                             train_test_unit_overlap=0, seconds=time.monotonic()-t0,
                             prediction_sha256=sha(prediction_path), model_sha256=sha(model_path))
                write_json(audit_path, audit)
                event("fold_complete", cohort=cohort, definition=definition, feature_set=feature_set,
                      fold=k, background=b, seconds=audit["seconds"])
            assert prediction.shape == (len(test),) and np.isfinite(prediction).all()
            assert ((0 <= prediction) & (prediction <= 1)).all()
            scores[b][test] = prediction
            fold_audits.append(audit)
    assert all(np.isfinite(s).all() for s in scores.values())
    membership["new_score_A"] = scores["A"]
    membership["new_score_B"] = scores["B"]
    membership.to_csv(folder / "new_oof.tsv.gz", sep="\t", index=False)
    evaluation = evaluate(sites, scores, reps)
    evaluation.update(cohort=cohort, definition=definition, feature_set=feature_set,
                      fold_audits=fold_audits, identity=task_identity, seconds=time.monotonic()-started)
    write_json(folder / "evaluation.json", evaluation)
    outputs = {name:sha(folder / name) for name in sorted(required_outputs)}
    write_json(complete, dict(identity=task_identity, task_metadata=task_metadata, outputs=outputs,
               completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    event("task_complete", cohort=cohort, definition=definition, feature_set=feature_set,
          seconds=evaluation["seconds"], delta_fix_A=evaluation["metrics"]["delta_fix_A"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-only", action="store_true")
    parser.add_argument("--bootstrap-reps", type=int, default=2000)
    args = parser.parse_args()
    identity, record = fingerprint()
    existing = ROOT / "logs/execution_identity.json"
    if existing.exists():
        assert json.loads(existing.read_text())["identity"] == identity, "Run identity changed; use a new version directory"
    write_json(ROOT / "logs/execution_identity.json", dict(identity=identity, **record))
    cohorts = COHORTS[:1] if args.pilot_only else COHORTS
    feature_sets = ["VIS10"] if args.pilot_only else ["VIS10", "DIG25", "CPL15"]
    event("run_start", identity=identity, pilot_only=args.pilot_only, planned_fits=len(cohorts)*len(feature_sets)*2*2*5)
    for feature_set in feature_sets:
        for definition in DEFINITIONS:
            for cohort in cohorts:
                run_task(cohort, definition, feature_set, identity, args.bootstrap_reps)
    write_json(ROOT / "logs/run_complete.json", dict(identity=identity, pilot_only=args.pilot_only,
               tasks=len(cohorts)*len(feature_sets)*2, bootstrap_reps=args.bootstrap_reps))
    event("run_complete", identity=identity)


if __name__ == "__main__":
    try:
        main()
    except BaseException as error:
        event("run_failed", exception=repr(error))
        raise
