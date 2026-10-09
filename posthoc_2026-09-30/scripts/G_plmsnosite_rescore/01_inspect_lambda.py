"""POST HOC revision analysis (G_plmsnosite_rescore): decode the marshal-serialized Lambda function
stored in models/Embedding.h5 and disassemble it under the running Python (3.10). Read-only."""
import json, pathlib, h5py, codecs, marshal, dis, sys
D = pathlib.Path(r"/path/to/local/_cys_repo_work/public/revision_2026-09-30/external/pLMSNOSite/models")
with h5py.File(D / "Embedding.h5", "r") as f:
    cfg = json.loads(f.attrs["model_config"])
lam = [l for l in cfg["config"]["layers"] if l["class_name"] == "Lambda"][0]["config"]
print("Lambda config keys:", {k: (v if k != "function" else "...") for k, v in lam.items()})
fn = lam["function"]
print("function field type:", type(fn), len(fn))
code_b64 = fn[0] if isinstance(fn, list) else fn
raw = codecs.decode(code_b64.encode("ascii"), "base64")
print("raw marshal bytes:", len(raw), raw[:16])
code = marshal.loads(raw)
print("co_names", code.co_names, "co_consts", code.co_consts, "co_varnames", code.co_varnames, "argcount", code.co_argcount)
print("python", sys.version)
dis.dis(code)
