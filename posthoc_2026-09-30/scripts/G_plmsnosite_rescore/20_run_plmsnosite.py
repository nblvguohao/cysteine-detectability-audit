"""G_plmsnosite_rescore, step 2 (TensorFlow venv: W/.venv_tf/Scripts/python.exe):
run the released pLMSNOSite models on the released independent test set exactly as the authors'
evaluate_model.py does, and keep the per-site probabilities that the script thresholds away.

POST HOC revision analysis (2026-09-30); not registered.

Procedure
1. A staging copy of the authors' repository layout is written to an ASCII-only temporary
   directory (see STAGE below) from the git object store at commit e9158af
   (evaluate_model.py, models/*.h5, data/test/*) so that the data bytes are the LF blobs whose
   hashes Supplemental Note 5 registers (the Windows checkout converts line endings).
2. evaluate_model.py is executed verbatim with runpy (cwd = staging). Its imports of tqdm,
   matplotlib, seaborn and mpl_toolkits are unused by its computation and are not installed in the
   venv; inert stub modules stand in for them. Its console output is saved.
3. From the script's own namespace the stacked test features (X_stacked_test) and the stacking model
   (model_final) are taken, and model_final.predict(X_stacked_test) gives the probabilities the
   script thresholds at 0.5. Agreement of (probability > 0.5) with the script's y_pred is asserted.
4. Independent cross-check without the marshal-serialized Lambda layer: the Embedding model is
   rebuilt from its stored config with the Lambda (K.expand_dims(x, 3)) replaced by an equivalent
   Reshape((37, 4, 1)); weights are loaded by layer name; the test inputs are parsed by this script,
   not by evaluate_model.py; the whole stack is recomputed and compared with step 3.
5. Determinism: step 3 prediction repeated; max |difference| reported.
"""
from __future__ import annotations

import os

os.environ.setdefault("TF_DETERMINISTIC_OPS", "1")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import contextlib
import copy
import importlib.util
import io
import json
import platform
import runpy
import shutil
import sys
import time
import types

import numpy as np
import pandas as pd
import tensorflow as tf

tf.config.threading.set_intra_op_parallelism_threads(4)
tf.config.threading.set_inter_op_parallelism_threads(2)
tf.random.set_seed(20260930)
np.random.seed(20260930)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

# TensorFlow's C++ file layer on Windows cannot open paths (or a working directory) containing
# non-ASCII characters, and every path under W contains "小论文" (no 8.3 alias exists). The staging
# copy is therefore written to an ASCII-only temporary directory (default: a fixed temporary folder;
# override with G_PLMS_STAGE) and deleted after the run. All outputs go to results/G_plmsnosite_rescore.
import pathlib  # noqa: E402

STAGE = pathlib.Path(os.environ.get(
    "G_PLMS_STAGE",
    r"/path/to/tmp/"
    r"g_plmsnosite_stage"))
assert str(STAGE).isascii(), "staging path must be ASCII-only for TensorFlow file I/O"
FILES =["evaluate_model.py", "predict.py", "models/pLMSNOSite.h5", "models/ProtT5.h5", "models/Embedding.h5",
         "data/test/protT5_test.csv", "data/test/sequence_test.csv", "data/test/embedding_test.fasta"]
ALPHABET_EVAL = "ARNDCQEGHILKMFPSTWYVUX-"   # evaluate_model.py (23 symbols; Embedding input_dim = 23)
ALPHABET_PRED = "ARNDCQEGHILKMFPSTWYVX"     # predict.py (21 symbols) -- recorded, not used


def stage():
    if STAGE.exists():
        shutil.rmtree(STAGE)
    hashes = {}
    for rel in FILES:
        data = C.git_blob(rel)
        dest = STAGE / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
        hashes[f"pLMSNOSite@{C.EXPECTED_COMMIT[:7]}:{rel}"] = C.sha256_bytes(data)
    assert hashes[f"pLMSNOSite@{C.EXPECTED_COMMIT[:7]}:data/test/sequence_test.csv"] == \
        C.REGISTERED_SHA256["data/test/sequence_test.csv"]
    return hashes


def install_stubs():
    """Inert stand-ins for modules evaluate_model.py imports but never uses in its computation."""
    stubbed = []

    def stub(name, **attrs):
        mod = types.ModuleType(name)
        mod.__dict__.update(attrs)
        sys.modules[name] = mod
        stubbed.append(name)
        return mod

    if importlib.util.find_spec("tqdm") is None:
        stub("tqdm")
    if importlib.util.find_spec("matplotlib") is None:
        mpl = stub("matplotlib")
        mpl.pyplot = stub("matplotlib.pyplot")
        mpl.colors = stub("matplotlib.colors", ListedColormap=object)
    if importlib.util.find_spec("seaborn") is None:
        stub("seaborn")
    if importlib.util.find_spec("mpl_toolkits") is None:
        mt = stub("mpl_toolkits")
        mt.mplot3d = stub("mpl_toolkits.mplot3d", Axes3D=object)
    return stubbed


class Tee(io.TextIOBase):
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)
        return len(s)

    def flush(self):
        for st in self.streams:
            st.flush()


def run_verbatim():
    buf = io.StringIO()
    cwd = os.getcwd()
    os.chdir(STAGE)
    try:
        with contextlib.redirect_stdout(Tee(sys.__stdout__, buf)):
            g = runpy.run_path(str(STAGE / "evaluate_model.py"), run_name="__main__")
    finally:
        os.chdir(cwd)
    return g, buf.getvalue()


def penultimate(model, x):
    layer = model.layers[len(model.layers) - 2]
    sub = tf.keras.Model(inputs=model.input, outputs=layer.output)
    return sub.predict(x, verbose=0), layer.name


def rebuild_without_lambda(path):
    import h5py
    with h5py.File(path, "r") as f:
        cfg = json.loads(f.attrs["model_config"])
    cfg2 = copy.deepcopy(cfg)
    replaced = []
    for i, layer in enumerate(cfg2["config"]["layers"]):
        if layer["class_name"] == "Lambda":
            replaced.append(layer["config"]["name"])
            cfg2["config"]["layers"][i] = {"class_name": "Reshape",
                                           "config": {"name": layer["config"]["name"], "trainable": True,
                                                      "dtype": "float32", "target_shape": [37, 4, 1]}}
    model = tf.keras.models.model_from_json(json.dumps(cfg2))
    model.load_weights(str(path), by_name=True)
    return model, replaced


def parse_fasta(text):
    ids, seqs, cur = [], [], None
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
    t0 = time.time()
    C.OUT.mkdir(parents=True, exist_ok=True)
    hashes = stage()
    stubbed = install_stubs()

    # ---------------- step 2-3: verbatim run
    g, log = run_verbatim()
    (C.OUT / "evaluate_model_verbatim_stdout.txt").write_text(log, encoding="utf-8")
    model_final = g["model_final"]
    X_stacked = np.asarray(g["X_stacked_test"], dtype=np.float32)
    y = np.asarray(g["y_stacked_test"]).astype(int)
    prob = model_final.predict(X_stacked, verbose=0).ravel().astype(np.float64)
    prob_again = model_final.predict(X_stacked, verbose=0).ravel().astype(np.float64)
    ypred_script = np.asarray(g["y_pred"]).ravel().astype(bool)
    assert np.array_equal(prob > 0.5, ypred_script), "thresholded probabilities differ from the script's y_pred"
    members = g["members"]
    base_prott5 = members[0].predict(g["X_test_pt5"], verbose=0).ravel().astype(np.float64)
    base_emb = members[1].predict(g["X_test_embedding"], verbose=0).ravel().astype(np.float64)
    layer_names = {"ProtT5": members[0].layers[len(members[0].layers) - 2].name,
                   "Embedding": members[1].layers[len(members[1].layers) - 2].name}

    # ---------------- step 4: independent cross-check without the Lambda bytecode
    seq_df = pd.read_csv(STAGE / "data/test/sequence_test.csv")
    pt5 = pd.read_csv(STAGE / "data/test/protT5_test.csv")
    ids, wins = parse_fasta((STAGE / "data/test/embedding_test.fasta").read_text())
    align = {
        "fasta_ids_match_sequence_csv": bool(ids == [f"{u}|{p}" for u, p in zip(seq_df.UniProt, seq_df.Position)]),
        "prott5_rows_match_sequence_csv": bool((pt5.UniProt.values == seq_df.UniProt.values).all()
                                               and (pt5.Position.values == seq_df.Position.values).all()),
        "prott5_target_equals_sequence_target": bool((pt5.Target.values == seq_df.Target.values).all()),
        "window_lengths": sorted({len(w) for w in wins}),
        "centre_is_C": int(sum(w[18] == "C" for w in wins)),
        "windows_with_dash_padding": int(sum("-" in w for w in wins)),
        "windows_with_X": int(sum("X" in w for w in wins)),
        "windows_with_U": int(sum("U" in w for w in wins)),
        "symbols_outside_eval_alphabet": sorted({c for w in wins for c in w} - set(ALPHABET_EVAL)),
    }
    cmap = {c: i for i, c in enumerate(ALPHABET_EVAL)}
    x_emb = np.array([[cmap[c] for c in w] for w in wins], dtype=np.int64)
    x_pt5 = pt5.iloc[:, 3:-1].to_numpy(dtype=np.float32)
    m_pt5 = tf.keras.models.load_model(str(STAGE / "models/ProtT5.h5"), compile=False)
    m_emb, replaced = rebuild_without_lambda(STAGE / "models/Embedding.h5")
    m_fin = tf.keras.models.load_model(str(STAGE / "models/pLMSNOSite.h5"), compile=False)
    f1, n1 = penultimate(m_pt5, x_pt5)
    f2, n2 = penultimate(m_emb, x_emb)
    stacked2 = np.concatenate([f1, f2], axis=1).astype(np.float32)
    prob2 = m_fin.predict(stacked2, verbose=0).ravel().astype(np.float64)
    cross = {
        "lambda_layers_replaced_by_reshape": replaced,
        "penultimate_layers": {"ProtT5": n1, "Embedding": n2},
        "max_abs_diff_stacked_features": float(np.max(np.abs(stacked2 - X_stacked))),
        "max_abs_diff_probability": float(np.max(np.abs(prob2 - prob))),
        "max_abs_diff_repeat_prediction": float(np.max(np.abs(prob_again - prob))),
        "identical_threshold_calls": bool(np.array_equal(prob2 > 0.5, prob > 0.5)),
    }

    # ---------------- outputs
    out = pd.DataFrame({
        "row": np.arange(len(y)), "UniProt": seq_df.UniProt, "Position": seq_df.Position, "Target": y,
        "plmsnosite_prob": prob, "plmsnosite_pred_0.5": (prob > 0.5).astype(int),
        "prott5_base_prob": base_prott5, "embedding_base_prob": base_emb,
        "plmsnosite_prob_rebuilt_check": prob2,
    })
    out.to_csv(C.OUT / "plmsnosite_test_scores.csv", index=False, float_format="%.17g")
    feat = pd.DataFrame(X_stacked, columns=[f"prott5_h{i}" for i in range(4)] + [f"emb_h{i}" for i in range(16)])
    feat.insert(0, "Position", seq_df.Position.values)
    feat.insert(0, "UniProt", seq_df.UniProt.values)
    feat.to_csv(C.OUT / "plmsnosite_stacked_features_test.csv", index=False, float_format="%.9g")

    import h5py
    import keras
    import sklearn
    info = {
        "label": "POST HOC revision analysis 2026-09-30 (G_plmsnosite_rescore); not registered",
        "git_head": C.git_head(),
        "inputs_sha256_git_blobs": hashes,
        "stub_modules": stubbed,
        "evaluate_model_penultimate_layers": layer_names,
        "alignment_checks": align,
        "cross_check": cross,
        "n_test": int(len(y)), "n_pos": int(y.sum()),
        "environment": {"python": platform.python_version(), "tensorflow": tf.__version__, "keras": keras.__version__,
                        "numpy": np.__version__, "pandas": pd.__version__, "h5py": h5py.__version__,
                        "scikit-learn": sklearn.__version__, "platform": platform.platform(),
                        "TF_DETERMINISTIC_OPS": os.environ.get("TF_DETERMINISTIC_OPS"),
                        "intra_op_threads": 4, "inter_op_threads": 2},
        "saved_model_keras_version": "2.8.0 (h5 attribute keras_version)",
        "predict_py_note": ("predict.py pads windows with 'X' and encodes with a 21-symbol alphabet (X->20), "
                            "whereas evaluate_model.py and the 23-symbol Embedding layer use '-' padding "
                            "(-> 22, X -> 21, U -> 20); the benchmark is scored with evaluate_model.py."),
        "runtime_s": round(time.time() - t0, 1),
    }
    (C.OUT / "plmsnosite_run.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    print(json.dumps({k: info[k] for k in ("alignment_checks", "cross_check", "evaluate_model_penultimate_layers")},
                     indent=1))
    try:
        shutil.rmtree(STAGE)
    except OSError as exc:  # keep going; the staging copy is reproducible from the git blobs
        print("could not remove staging directory:", exc)


if __name__ == "__main__":
    main()
