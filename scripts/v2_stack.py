"""Nested, homology-isolated model stack for Cys annotation ranking (v2).

Split discipline, identical for every member and for the blender:

* outer fold ``k`` is the frozen evaluation fold and is never read during
  preprocessing, hyperparameter search, epoch selection, calibration or
  blend-weight estimation;
* inner validation fold ``v = k % 5 + 1`` is one of the remaining frozen
  component folds, matching the v1 corrected-MLP rule;
* stage 1 fits on ``subtrain = folds not in {k, v}`` and selects everything on
  ``v``;
* stage 2 refits with the stage-1 selection on ``outer_train = folds != k`` and
  predicts fold ``k``.

Calibration and blend weights are estimated from stage-1 out-of-sample
predictions on ``v``, whose labels belong to the outer-training partition.
"""
from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression

from common import ROOT, load_frozen, validation_objective

FEATURES = ROOT / "features"
PLM = ROOT / "plm"


# --------------------------------------------------------------------------- #
# data


def load_v2_matrix():
    archive = np.load(FEATURES / "v2_features.npz", allow_pickle=False)
    X = archive["X"]
    names = archive["names"].astype(str)
    X_frozen, y, folds, proteins, components, positions, meta = load_frozen()
    assert np.array_equal(archive["proteins"].astype(str), proteins)
    assert np.array_equal(archive["positions"].astype(int), positions)
    assert np.array_equal(archive["folds"].astype(int), folds)
    assert np.array_equal(archive["y"].astype(int), y)
    assert np.array_equal(archive["components"].astype(str), components)
    assert np.allclose(X[:, :23], X_frozen.astype(np.float32), equal_nan=True)
    return X, names, y, folds, proteins, components, positions, meta


def load_embeddings(proteins, positions):
    archive = np.load(PLM / "protbert_cys_embeddings.npz", allow_pickle=False)
    assert np.array_equal(archive["proteins"].astype(str), proteins)
    assert np.array_equal(archive["positions"].astype(int), positions)
    raw = archive["embeddings"].astype(np.float32)
    centred = np.empty_like(raw)
    order = np.argsort(proteins, kind="stable")
    starts = np.flatnonzero(np.r_[True, proteins[order][1:] != proteins[order][:-1]])
    for group in np.split(order, starts[1:]):
        block = raw[group]
        centred[group] = block - block.mean(axis=0, keepdims=True)
    return raw, centred


def fold_partitions(folds, outer_fold, proteins, components):
    validation_fold = outer_fold % 5 + 1
    subtrain = np.flatnonzero((folds != outer_fold) & (folds != validation_fold))
    validation = np.flatnonzero(folds == validation_fold)
    outer_train = np.flatnonzero(folds != outer_fold)
    test = np.flatnonzero(folds == outer_fold)
    for left, right in (
        (subtrain, validation),
        (subtrain, test),
        (validation, test),
        (outer_train, test),
    ):
        assert not set(components[left]).intersection(components[right])
        assert not set(proteins[left]).intersection(proteins[right])
    return subtrain, validation, outer_train, test, validation_fold


# --------------------------------------------------------------------------- #
# helpers


def group_sorted(index, proteins):
    """Sort an index array by protein and return LightGBM group sizes."""
    order = np.argsort(proteins[index], kind="stable")
    ordered = index[order]
    _, counts = np.unique(proteins[ordered], return_counts=True)
    boundaries = np.flatnonzero(
        np.r_[True, proteins[ordered][1:] != proteins[ordered][:-1]]
    )
    sizes = np.diff(np.r_[boundaries, len(ordered)])
    assert sizes.sum() == len(ordered)
    return ordered, sizes


def objective_on(y, scores, proteins, components, index):
    return validation_objective(
        y[index], scores, proteins[index], components[index]
    )


def platt(scores, labels):
    """Platt scaling on inner-fold out-of-sample scores.

    The raw score is first standardised with the inner-fold mean and standard
    deviation.  Those two constants are carried to the test fold unchanged, so
    no test-fold statistic enters the calibration.
    """
    finite = np.isfinite(scores)
    if finite.sum() < 2:
        raise ValueError("Calibration needs at least two finite scores")
    centre = float(np.mean(scores[finite]))
    spread = float(np.std(scores[finite]))
    if not np.isfinite(spread) or spread < 1e-12:
        spread = 1.0
    standardised = ((scores[finite] - centre) / spread).reshape(-1, 1)
    model = LogisticRegression(C=1.0, solver="lbfgs", max_iter=2000)
    # Apple Accelerate raises spurious floating-point flags inside matmul during
    # the L-BFGS iterations; finiteness of the fitted coefficients is asserted
    # explicitly below instead of relying on those flags.
    with np.errstate(all="ignore"):
        model.fit(standardised, labels[finite])
    if not np.isfinite(model.coef_).all() or not np.isfinite(model.intercept_).all():
        raise ValueError("Calibration did not converge to finite coefficients")
    return {
        "model": model,
        "centre": centre,
        "spread": spread,
        "coefficient": float(model.coef_[0][0]),
        "intercept": float(model.intercept_[0]),
    }


def apply_platt(calibrator, scores):
    standardised = ((scores - calibrator["centre"]) / calibrator["spread"]).reshape(-1, 1)
    probabilities = calibrator["model"].predict_proba(standardised)[:, 1]
    if not np.isfinite(probabilities).all():
        raise ValueError("Calibrated probabilities contain non-finite values")
    return probabilities


def simplex_grid(n_members, step=0.1):
    """All non-negative weight vectors on a simplex grid."""
    levels = int(round(1.0 / step))

    def recurse(remaining, slots):
        if slots == 1:
            yield (remaining,)
            return
        for take in range(remaining + 1):
            for rest in recurse(remaining - take, slots - 1):
                yield (take,) + rest

    for combination in recurse(levels, n_members):
        yield np.asarray(combination, dtype=float) / levels
