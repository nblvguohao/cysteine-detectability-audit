"""Step 0 (POST HOC revision analysis, item D_empirical_background): reproduce the stored Fig. 2 trypsin-arm
cleavage-geometry values with the released Cys-Audit CLI (v0.2.2; the cleavage, statistics, protease and I/O
modules are byte-identical to the archived v0.1.0 that produced the stored audits, apart from line endings).

The stored audits used the deposit-era FASTA 'Mouse_2025_UP000000589_10090.fasta' (sha256 033ac9c4...), which is
not available locally; this run uses UniProt 2026_03 (external/UP000000589_10090.fasta.gz), decompressed into
results/D_empirical_background/_work/ for the CLI and deleted afterwards. Differences from the stored values are
expected from proteins removed or changed between releases and are reported, not tuned away.
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
from common import (CYS_AUDIT_SRC, IN_FASTA_GZ, IN_HYDP, RESULTS, STORED_AUDIT_OBSERVED, STORED_AUDIT_PROTEOME,
                    TOOL_SEED, dump_json, sha256_file)


def main():
    work = f"{RESULTS}/_work"
    os.makedirs(work, exist_ok=True)
    fasta = f"{work}/UP000000589_10090_2026_03.fasta"
    with gzip.open(IN_FASTA_GZ, "rb") as src, open(fasta, "wb") as dst:
        shutil.copyfileobj(src, dst)
    fasta_sha = sha256_file(fasta)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONPATH=CYS_AUDIT_SRC, PYTHONHASHSEED="0")
    out = {"fasta_2026_03_decompressed_sha256": fasta_sha, "runs": {}}
    for bg, stored_path in (("proteome", STORED_AUDIT_PROTEOME), ("observed", STORED_AUDIT_OBSERVED)):
        od = f"{RESULTS}/cli_repro_trypsin_{bg}"
        cmd = [sys.executable, "-m", "cys_audit", "audit", "--input", IN_HYDP, "--fasta", fasta, "--expand-background",
               "--protease", "trypsin", "--background", bg, "--modification", "S-palmitoylation",
               "--identity-readout", "chemistry_inferred", "--seed", str(TOOL_SEED), "--output", od,
               "--dataset", f"pxd063463_Trypsin_HydP_trypsin_{bg}_fasta2026_03"]
        r = subprocess.run(cmd, env=env, capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit(r.stderr)
        new = json.load(open(f"{od}/audit.json", encoding="utf-8"))
        old = json.load(open(stored_path, encoding="utf-8"))
        rec = {"cli_stdout": r.stdout.splitlines(), "input_new": new["input"], "input_stored": old["input"]}
        for band in ("proximal_1_3", "distal_6_12"):
            nb = new["tests"]["cleavage_geometry"]["details"]["bands"][band]
            ob = old["tests"]["cleavage_geometry"]["details"]["bands"][band]
            rec[band] = {"stored": {k: ob[k] for k in ("estimate", "ci", "share_flagged_positive", "share_flagged_background")},
                         "reproduced_fasta_2026_03": {k: nb[k] for k in ("estimate", "ci", "share_flagged_positive",
                                                                          "share_flagged_background")},
                         "ci_level": new["tests"]["cleavage_geometry"]["ci_level"]}
        rec["n_positive"] = {"stored": old["tests"]["cleavage_geometry"]["n_positive"],
                             "reproduced": new["tests"]["cleavage_geometry"]["n_positive"]}
        rec["n_background"] = {"stored": old["tests"]["cleavage_geometry"]["n_background"],
                               "reproduced": new["tests"]["cleavage_geometry"]["n_background"]}
        out["runs"][bg] = rec
        print(bg, json.dumps({k: rec[k] for k in ("proximal_1_3", "distal_6_12", "n_positive", "n_background")}, indent=1))
    dump_json(out, f"{RESULTS}/cli_reproduction.json")
    os.remove(fasta)
    os.rmdir(work)


if __name__ == "__main__":
    main()
