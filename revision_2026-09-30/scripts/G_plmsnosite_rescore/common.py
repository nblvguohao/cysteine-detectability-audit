"""Shared helpers for item G_plmsnosite_rescore.

POST HOC revision analysis (2026-09-30) in response to pre-submission review criticism.
Nothing here is registered or pre-specified.

Contents
- paths and the git-blob reader (the Windows working tree of the pLMSNOSite clone has CRLF line
  endings, so every pLMSNOSite file is read from the git object store with `git show HEAD:<path>`,
  which returns the bytes whose SHA-256 values are registered in Supplemental Note 5);
- the protein-clustered bootstrap in the house style of the repository
  (sorted protein labels -> codes 0..K-1; numpy default_rng(seed); one `rng.integers(0, K, K)`
  draw per replicate; 95% percentile interval by numpy.quantile);
- tie-corrected AUROC and step-wise average precision (identical to scikit-learn's
  roc_auc_score / average_precision_score) computed from per-row integer weights, so that a
  bootstrap replicate (a vector of protein multiplicities) is evaluated without copying rows.
"""
from __future__ import annotations

import hashlib
import io
import pathlib
import subprocess

import numpy as np

W = pathlib.Path(r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30")
EXT = W / "external" / "pLMSNOSite"
OUT = W / "results" / "G_plmsnosite_rescore"
SCRIPTS = W / "scripts" / "G_plmsnosite_rescore"
EXPECTED_COMMIT = "e9158af06418f3e50c9bbbbeca77ea692e6507ce"

# Registered in Supplemental Note 5 (git-blob bytes).
REGISTERED_SHA256 = {
    "data/train/sequence_train.csv": "a3b4a1fe3858ed3e4a0aa62d47e52c8e1c118cceb06c4edcf75c0b82ff96bd0f",
    "data/test/sequence_test.csv": "d4498e6f2afa28b89d150dacf44b1de6bbfd791e772f1fefe120aa4aabcc5551",
}

# Values quoted in the manuscript / Supplemental Note 5 (for comparison only).
MANUSCRIPT = {
    "VIS10_auroc": 0.7095, "VIS10_ci": (0.6817, 0.7367), "VIS10_auprc": 0.1677,
    "DIG25_auroc": 0.7765, "DIG25_ci": (0.7231, 0.8274), "DIG25_auprc": 0.2195,
    "NC_auroc": 0.5124, "NC_ci": (0.471, 0.5546),
    "PC_folds": (0.7669, 0.8344, 0.7412, 0.7371, 0.7467), "PC_mean": 0.7653,
    "published_auroc_pLMSNOSite": 0.754, "published_auroc_PreSNO": 0.756,
    "published_auroc_DeepNitro": 0.731,
    "published_sn_pLMSNOSite": 0.735, "published_sp_pLMSNOSite": 0.773,
    "recovery_VIS10": {"pLMSNOSite": 0.8248, "PreSNO": 0.8184, "DeepNitro": 0.9069},
    "recovery_DIG25": {"pLMSNOSite": 1.0886, "PreSNO": 1.0801, "DeepNitro": 1.197},
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_head() -> str:
    out = subprocess.run(["git", "-C", str(EXT), "rev-parse", "HEAD"], capture_output=True, check=True)
    return out.stdout.decode().strip()


def git_blob(path: str) -> bytes:
    """Raw blob bytes of `path` at HEAD of the pLMSNOSite clone (no line-ending conversion)."""
    out = subprocess.run(["git", "-C", str(EXT), "show", f"HEAD:{path}"], capture_output=True, check=True)
    return out.stdout


def read_blob_csv(path: str, manifest: dict | None = None):
    import pandas as pd
    data = git_blob(path)
    digest = sha256_bytes(data)
    if path in REGISTERED_SHA256 and digest != REGISTERED_SHA256[path]:
        raise SystemExit(f"{path}: blob sha256 {digest} differs from the registered value")
    if manifest is not None:
        manifest[f"pLMSNOSite@{EXPECTED_COMMIT[:7]}:{path}"] = digest
    return pd.read_csv(io.BytesIO(data))


# ---------------------------------------------------------------- clustered bootstrap machinery
def protein_codes(groups):
    labels = sorted(set(groups))
    pos = {c: i for i, c in enumerate(labels)}
    return np.fromiter((pos[g] for g in groups), dtype=np.int64, count=len(groups)), labels


def multiplicities(n_clusters: int, reps: int, seed: int, mode: str = "per_replicate"):
    """Yield one length-K vector of protein multiplicities per replicate.

    per_replicate: rng.integers(0, K, K) once per replicate (run_qtrp_intervals_and_cross_label.py).
    block:         rng.integers(0, K, size=(reps, K)) drawn at once (run_phase2_synthetic_benchmark).
    chunk250:      cys_audit.stats.multiplicities (blocks of 250 replicates).
    """
    rng = np.random.default_rng(seed)
    if mode == "per_replicate":
        for _ in range(reps):
            yield np.bincount(rng.integers(0, n_clusters, n_clusters), minlength=n_clusters)
    elif mode == "block":
        draws = rng.integers(0, n_clusters, size=(reps, n_clusters))
        for r in range(reps):
            yield np.bincount(draws[r], minlength=n_clusters)
    elif mode == "chunk250":
        done = 0
        while done < reps:
            m = min(250, reps - done)
            draws = rng.integers(0, n_clusters, size=(m, n_clusters))
            for r in range(m):
                yield np.bincount(draws[r], minlength=n_clusters)
            done += m
    else:
        raise ValueError(mode)


class RankedScore:
    """Pre-sorted score so that weighted AUROC / AP are O(n) per replicate."""

    def __init__(self, scores, labels):
        scores = np.asarray(scores, dtype=float)
        self.labels = np.asarray(labels).astype(bool)
        uniq, self.grp = np.unique(scores, return_inverse=True)  # ascending distinct scores
        self.n_groups = len(uniq)

    def auroc(self, w=None):
        w = np.ones(len(self.labels)) if w is None else np.asarray(w, dtype=float)
        pg = np.bincount(self.grp, weights=w * self.labels, minlength=self.n_groups)
        ng = np.bincount(self.grp, weights=w * ~self.labels, minlength=self.n_groups)
        wp, wn = pg.sum(), ng.sum()
        if wp == 0 or wn == 0:
            return float("nan")
        below = np.concatenate(([0.0], np.cumsum(ng)[:-1]))
        return float((pg * (below + 0.5 * ng)).sum() / (wp * wn))

    def auprc(self, w=None):
        """Average precision as in sklearn.metrics.average_precision_score (step-wise)."""
        w = np.ones(len(self.labels)) if w is None else np.asarray(w, dtype=float)
        pg = np.bincount(self.grp, weights=w * self.labels, minlength=self.n_groups)[::-1]
        ng = np.bincount(self.grp, weights=w * ~self.labels, minlength=self.n_groups)[::-1]
        tp, fp = np.cumsum(pg), np.cumsum(ng)
        if tp[-1] == 0:
            return float("nan")
        with np.errstate(invalid="ignore", divide="ignore"):
            precision = np.where(tp + fp > 0, tp / (tp + fp), 0.0)
        recall = tp / tp[-1]
        d_recall = np.diff(np.concatenate(([0.0], recall)))
        return float((d_recall * precision).sum())


def pctl(values, level=0.95):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    a = (1.0 - level) / 2.0
    return float(np.quantile(v, a)), float(np.quantile(v, 1.0 - a))


def boot_p_le(values, null=0.0):
    """Share of replicates at or below `null` (one-sided, descriptive only)."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    return float(np.mean(v <= null))


def logit(p, eps=1e-6):
    p = np.clip(np.asarray(p, dtype=float), eps, 1.0 - eps)
    return np.log(p / (1.0 - p))
