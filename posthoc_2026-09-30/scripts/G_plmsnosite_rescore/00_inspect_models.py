"""POST HOC revision analysis (G_plmsnosite_rescore): inspect the released pLMSNOSite .h5 model configs.
Read-only; prints layer summaries. Run with the TF venv python."""
import json, pathlib, h5py
D = pathlib.Path(r"/path/to/local/_cys_repo_work/public/revision_2026-09-30/external/pLMSNOSite/models")
for name in ["ProtT5", "Embedding", "pLMSNOSite"]:
    with h5py.File(D / f"{name}.h5", "r") as f:
        attrs = dict(f.attrs)
        cfg = attrs.get("model_config")
        if isinstance(cfg, bytes):
            cfg = cfg.decode()
        cfg = json.loads(cfg)
        print("=" * 30, name, "keras_version", attrs.get("keras_version"), "backend", attrs.get("backend"))
        print("class", cfg["class_name"])
        for layer in cfg["config"]["layers"]:
            c = layer["config"]
            extra = {k: c.get(k) for k in ("units", "activation", "input_dim", "output_dim", "filters", "kernel_size", "rate", "batch_input_shape", "function", "target_shape") if k in c}
            if "function" in extra and extra["function"] is not None:
                extra["function"] = str(extra["function"])[:80]
            print(" ", layer["class_name"], c.get("name"), extra)
        print("training_config" in attrs)
