"""G_plmsnosite_rescore, step 2b (system python; no TensorFlow): staging-free, TensorFlow-free
re-computation of the released pLMSNOSite ensemble on the released test set.

POST HOC revision analysis (2026-09-30), added after adversarial verification round 1; not registered,
not pre-specified.

Why this step exists
- The verbatim TensorFlow run (20_run_plmsnosite.py) had to stage the released files in an ASCII-only
  temporary directory outside the item folders, because TensorFlow's Windows file layer cannot open
  non-ASCII paths. This step reads the same git blobs into memory (h5py on io.BytesIO) and writes
  nothing outside results/G_plmsnosite_rescore.
- The comparison combines TensorFlow 2.15.1 inference (venv) with statistics computed in the system
  Python. This step recomputes the probabilities in the system Python without TensorFlow, so the
  cross-environment agreement can be stated.

Procedure
- Model configs and weights are read from models/{ProtT5,Embedding,pLMSNOSite}.h5 (git blobs at
  e9158af) with h5py; every layer type and hyper-parameter that the forward pass relies on is asserted.
- Forward pass in numpy float64: Normalization (x - mean) / max(sqrt(variance), 1e-7) (Keras 2.x
  definition, epsilon 1e-7); Dense; ReLU; Dropout = identity at inference; Embedding lookup; the
  serialized Lambda = expand_dims(x, 3) (decoded in 01_inspect_lambda.py; replaced by an equivalent
  Reshape in the TensorFlow cross-check of 20_run_plmsnosite.py); Conv2D 'valid', stride 1;
  MaxPooling2D 'valid'; Flatten in C order (channels_last); sigmoid.
- Inputs are parsed here: data/test/protT5_test.csv columns 3:-1 (as evaluate_model.py) and
  data/test/embedding_test.fasta encoded with the alphabet of evaluate_model.py
  ('ARNDCQEGHILKMFPSTWYVUX-').
- Stacked features = penultimate-layer outputs of the two base models (Dense 4 and Dense 16), in the
  order ProtT5 then Embedding, as in evaluate_model.py.

Outputs (results/G_plmsnosite_rescore/)
- numpy_forward_pass_scores.csv   per-site probabilities (stack and both arms) from numpy
- numpy_forward_pass_check.json    max |numpy - TensorFlow| differences, threshold calls, AUROC/AUPRC
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")

import io
import json
import platform
import sys

import h5py
import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import average_precision_score, confusion_matrix, roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

ALPHABET = "ARNDCQEGHILKMFPSTWYVUX-"  # evaluate_model.py, get_input_for_embedding
EPS = 1e-7                             # keras.backend.epsilon()


def _s(x):
    return x.decode() if isinstance(x, bytes) else x


def load_h5_blob(rel: str, manifest: dict):
    """Config (layer list), weights {layer: {name: array}} and keras_version of an .h5 git blob, in memory."""
    blob = C.git_blob(rel)
    manifest[f"pLMSNOSite@{C.EXPECTED_COMMIT[:7]}:{rel}"] = C.sha256_bytes(blob)
    with h5py.File(io.BytesIO(blob), "r") as f:
        cfg = json.loads(_s(f.attrs["model_config"]))
        kv = _s(f.attrs.get("keras_version"))
        mw = f["model_weights"]
        weights = {}
        for lname in (_s(n) for n in mw.attrs["layer_names"]):
            g = mw[lname]
            wd = {}
            for wn in (_s(n) for n in g.attrs["weight_names"]):
                key = wn.split("/")[-1].split(":")[0]
                wd[key] = np.asarray(g[wn][()], dtype=np.float64)
            weights[lname] = wd
    assert cfg["class_name"] == "Sequential", cfg["class_name"]
    return cfg["config"]["layers"], weights, kv


def activation(h, spec):
    if isinstance(spec, dict):
        assert spec["class_name"] == "ReLU", spec
        c = spec["config"]
        assert c.get("max_value") is None and float(c.get("negative_slope", 0)) == 0.0 \
            and float(c.get("threshold", 0)) == 0.0, c
        return np.maximum(h, 0.0)
    if spec == "relu":
        return np.maximum(h, 0.0)
    if spec == "sigmoid":
        return 1.0 / (1.0 + np.exp(-h))
    if spec in ("linear", None):
        return h
    raise ValueError(spec)


def conv2d_valid(h, kernel, bias, chunk=400):
    """h (N, H, W, Cin), kernel (kh, kw, Cin, Cout); stride 1, 'valid', channels_last."""
    kh, kw, _, _ = kernel.shape
    out = []
    for s in range(0, len(h), chunk):
        win = np.lib.stride_tricks.sliding_window_view(h[s:s + chunk], (kh, kw), axis=(1, 2))
        out.append(np.einsum("nijcuv,uvcd->nijd", win, kernel, optimize=True) + bias)
    return np.concatenate(out, axis=0)


def maxpool2d_valid(h, pool, strides):
    ph, pw = pool
    sh, sw = strides
    n, hh, ww, cc = h.shape
    oh, ow = (hh - ph) // sh + 1, (ww - pw) // sw + 1
    out = np.full((n, oh, ow, cc), -np.inf)
    for a in range(ph):
        for b in range(pw):
            out = np.maximum(out, h[:, a:a + sh * (oh - 1) + 1:sh, b:b + sw * (ow - 1) + 1:sw, :])
    return out


def forward(layers, weights, x, audit: list, model_name: str):
    """Outputs of every non-input layer (list), in order; asserts the hyper-parameters used."""
    acts, h = [], x
    for layer in layers:
        cls, c = layer["class_name"], layer["config"]
        name = c["name"]
        if cls == "InputLayer":
            continue
        if cls == "Normalization":
            w = weights[name]
            assert c.get("axis", -1) in (-1, [-1]) and not c.get("invert", False), c
            h = (h - w["mean"].reshape(1, -1)) / np.maximum(np.sqrt(w["variance"].reshape(1, -1)), EPS)
        elif cls == "Dense":
            w = weights[name]
            assert c.get("use_bias", True)
            h = activation(h @ w["kernel"] + w["bias"], c["activation"])
        elif cls == "Dropout":
            pass  # identity at inference
        elif cls == "Embedding":
            assert int(c["input_dim"]) == len(ALPHABET), c["input_dim"]
            h = weights[name]["embeddings"][h]
        elif cls == "Lambda":
            h = h[..., None]  # K.expand_dims(x, 3) on a (N, 37, 4) tensor
        elif cls == "Conv2D":
            assert c["padding"] == "valid" and list(c["strides"]) == [1, 1] and list(c["dilation_rate"]) == [1, 1]
            assert c.get("data_format", "channels_last") == "channels_last" and c.get("groups", 1) == 1
            assert c.get("use_bias", True)
            w = weights[name]
            h = activation(conv2d_valid(h, w["kernel"], w["bias"]), c["activation"])
        elif cls == "MaxPooling2D":
            assert c["padding"] == "valid" and c.get("data_format", "channels_last") == "channels_last"
            pool = tuple(c["pool_size"])
            strides = tuple(c["strides"]) if c.get("strides") is not None else pool
            h = maxpool2d_valid(h, pool, strides)
        elif cls == "Flatten":
            assert c.get("data_format", "channels_last") == "channels_last"
            h = h.reshape(len(h), -1)
        else:
            raise ValueError(f"unsupported layer {cls}")
        audit.append({"model": model_name, "layer": name, "class": cls, "output_shape": list(h.shape[1:])})
        acts.append(h)
    return acts


def parse_fasta(text):
    ids, seqs = [], []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            ids.append(line[1:])
            seqs.append("")
        else:
            seqs[-1] += line
    return ids, seqs


def main():
    manifest = {}
    seq = C.read_blob_csv("data/test/sequence_test.csv", manifest)
    pt5 = C.read_blob_csv("data/test/protT5_test.csv", manifest)
    fasta_blob = C.git_blob("data/test/embedding_test.fasta")
    manifest[f"pLMSNOSite@{C.EXPECTED_COMMIT[:7]}:data/test/embedding_test.fasta"] = C.sha256_bytes(fasta_blob)
    ids, wins = parse_fasta(fasta_blob.decode("ascii"))
    assert ids == [f"{u}|{p}" for u, p in zip(seq.UniProt, seq.Position)]
    assert (pt5.UniProt.values == seq.UniProt.values).all() and (pt5.Position.values == seq.Position.values).all()
    assert all(len(w) == 37 and w[18] == "C" for w in wins)
    assert not ({ch for w in wins for ch in w} - set(ALPHABET))
    cmap = {ch: i for i, ch in enumerate(ALPHABET)}
    x_emb = np.array([[cmap[ch] for ch in w] for w in wins], dtype=np.int64)
    x_pt5 = pt5.iloc[:, 3:-1].to_numpy(dtype=np.float32).astype(np.float64)  # as evaluate_model.py
    assert x_pt5.shape[1] == 1024

    audit = []
    L1, W1, kv1 = load_h5_blob("models/ProtT5.h5", manifest)
    L2, W2, kv2 = load_h5_blob("models/Embedding.h5", manifest)
    L3, W3, kv3 = load_h5_blob("models/pLMSNOSite.h5", manifest)
    a1 = forward(L1, W1, x_pt5, audit, "ProtT5")
    a2 = forward(L2, W2, x_emb, audit, "Embedding")
    stacked = np.concatenate([a1[-2], a2[-2]], axis=1)
    a3 = forward(L3, W3, stacked, audit, "pLMSNOSite")
    prob = a3[-1].ravel()
    p_pt5, p_emb = a1[-1].ravel(), a2[-1].ravel()

    tf_scores = pd.read_csv(C.OUT / "plmsnosite_test_scores.csv")
    tf_stack = pd.read_csv(C.OUT / "plmsnosite_stacked_features_test.csv")
    assert (tf_scores.UniProt.values == seq.UniProt.values).all() and (tf_scores.Position.values == seq.Position.values).all()
    y = seq.Target.to_numpy().astype(int)
    tfp = tf_scores.plmsnosite_prob.to_numpy()
    calls_np, calls_tf = prob > 0.5, tfp > 0.5
    cm = confusion_matrix(y, calls_np)
    sn, sp = cm[1, 1] / cm[1].sum(), cm[0, 0] / cm[0].sum()
    check = {
        "label": "POST HOC revision analysis 2026-09-30 (G_plmsnosite_rescore), added after verification round 1; "
                 "not registered",
        "what": "TensorFlow-free numpy forward pass of the released models, read in memory from the git blobs "
                "(no staging copy on disk)",
        "keras_version_in_h5": {"ProtT5": kv1, "Embedding": kv2, "pLMSNOSite": kv3},
        "layers": audit,
        "max_abs_diff_vs_tensorflow": {
            "stacked_probability": float(np.max(np.abs(prob - tfp))),
            "ProtT5_arm_probability": float(np.max(np.abs(p_pt5 - tf_scores.prott5_base_prob.to_numpy()))),
            "Embedding_arm_probability": float(np.max(np.abs(p_emb - tf_scores.embedding_base_prob.to_numpy()))),
            "stacked_features_vs_saved_float32_%.9g": float(np.max(np.abs(stacked - tf_stack.iloc[:, 2:].to_numpy()))),
        },
        "threshold_0.5_calls_identical_to_tensorflow": bool(np.array_equal(calls_np, calls_tf)),
        "min_abs_distance_of_probability_from_0.5": float(np.min(np.abs(prob - 0.5))),
        "confusion_matrix_numpy": cm.tolist(),
        "sensitivity": float(sn), "specificity": float(sp), "balanced_accuracy_(=AUROC_of_calls)": float((sn + sp) / 2),
        "auroc_numpy": float(roc_auc_score(y, prob)), "auroc_tensorflow": float(roc_auc_score(y, tfp)),
        "auprc_numpy": float(average_precision_score(y, prob)), "auprc_tensorflow": float(average_precision_score(y, tfp)),
        "auroc_ProtT5_arm_numpy": float(roc_auc_score(y, p_pt5)),
        "auroc_Embedding_arm_numpy": float(roc_auc_score(y, p_emb)),
        "inputs_sha256": manifest,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                        "h5py": h5py.__version__, "scikit-learn": sklearn.__version__, "tensorflow": "not imported"},
    }
    out = pd.DataFrame({"row": np.arange(len(y)), "UniProt": seq.UniProt, "Position": seq.Position, "Target": y,
                        "plmsnosite_prob_numpy": prob, "prott5_arm_prob_numpy": p_pt5, "embedding_arm_prob_numpy": p_emb})
    out.to_csv(C.OUT / "numpy_forward_pass_scores.csv", index=False, float_format="%.17g")
    (C.OUT / "numpy_forward_pass_check.json").write_text(json.dumps(check, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in check.items() if k not in ("layers", "inputs_sha256")}, indent=1))


if __name__ == "__main__":
    main()
