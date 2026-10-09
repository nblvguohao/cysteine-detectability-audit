"""Command-line interface.

  cys-audit audit --input sites.tsv --protease trypsin --modification persulfidation \
      --background observed --output audit_report [--fasta proteome.fasta] [--expand-background] \
      [--precursor-ppm 4.5 --search-mods Sulfide,Oxidation] [--identity-readout direct_mass] \
      [--claim-feature COLUMN --claim-direction positive] [--seed 20260922] [--dataset NAME]
  cys-audit schema          print the input schema
  cys-audit simulate ...    write a synthetic input with a planted artefact (see synthetic.py)
"""
from __future__ import annotations

import argparse
import os
import sys

from . import constants as C
from .audit import run_audit, sha256_file
from .io import SCHEMA, InputError, read_sites
from .proteases import RULES
from .report import write_all


def _audit(a):
    try:
        ds = read_sites(a.input, fasta=a.fasta, expand_background=a.expand_background)
    except InputError as exc:
        print(f"input error: {exc}", file=sys.stderr)
        return 2
    cfg = {"dataset": a.dataset or os.path.splitext(os.path.basename(a.input))[0], "protease": a.protease.lower(),
           "modification": a.modification, "background": a.background, "seed": a.seed, "reps": a.reps,
           "precursor_ppm": a.precursor_ppm, "search_mods": [m for m in (a.search_mods or "").split(",") if m],
           "identity_readout": a.identity_readout, "claim_feature": a.claim_feature,
           "claim_direction": a.claim_direction,
           "claim_covariates": [c for c in (a.claim_covariates or "").split(",") if c],
           "mass_range": tuple(float(x) for x in a.mass_range.split(",")) if a.mass_range else None}
    inputs = {"sites_file": os.path.basename(a.input), "sites_sha256": sha256_file(a.input)}
    if a.fasta:
        inputs.update(fasta_file=os.path.basename(a.fasta), fasta_sha256=sha256_file(a.fasta))
    rec = run_audit(ds, cfg, inputs)
    command = ["cys-audit", "audit"] + [x if not os.path.isabs(x) else os.path.basename(x) for x in a._argv]
    write_all(rec, a.output, [a.input, a.fasta], command)
    for name, t in rec["tests"].items():
        print(f"{name:28s} {t['status']:12s} {t.get('estimate')}")
    print(f"{'overall_claim_status':28s} {rec['overall_claim_status']}")
    return 0


def _schema(_a):
    for c, (req, typ, desc) in SCHEMA.items():
        print(f"{c:14s} {req:9s} {typ:6s} {desc}")
    return 0


def _simulate(a):
    from .synthetic import simulate
    simulate(a.artefact, a.strength, a.seed, a.out_dir, n_proteins=a.n_proteins)
    print(f"wrote {a.out_dir}/sites.tsv and {a.out_dir}/proteome.fasta")
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    p = argparse.ArgumentParser(prog="cys-audit", description="Audit a cysteine-modification site table for "
                                "measurement-driven artefacts.")
    sub = p.add_subparsers(dest="cmd", required=True)
    au = sub.add_parser("audit", help="run the five artefact tests")
    au.add_argument("--input", required=True)
    au.add_argument("--output", required=True)
    au.add_argument("--protease", required=True, choices=sorted(RULES), type=str.lower)
    au.add_argument("--modification", default=None)
    au.add_argument("--background", default="observed", choices=["observed", "proteome"])
    au.add_argument("--fasta", default=None)
    au.add_argument("--expand-background", action="store_true")
    au.add_argument("--precursor-ppm", type=float, default=None)
    au.add_argument("--search-mods", default=None, help="comma-separated variable/fixed modification names")
    au.add_argument("--identity-readout", default="direct_mass", choices=["direct_mass", "chemistry_inferred"])
    au.add_argument("--mass-range", default=None, help="lo,hi peptide mass range (Da) if masses are not in the input")
    au.add_argument("--claim-feature", default=None)
    au.add_argument("--claim-covariates", default=None,
                    help="comma-separated numeric columns; the claim control stratifies on propensity quintiles of these")
    au.add_argument("--claim-direction", default="positive", choices=["positive", "negative", "null"])
    au.add_argument("--seed", type=int, default=C.DEFAULT_SEED)
    au.add_argument("--reps", type=int, default=C.BOOTSTRAP_REPS)
    au.add_argument("--dataset", default=None)
    sc = sub.add_parser("schema", help="print the input schema")
    si = sub.add_parser("simulate", help="write a synthetic input with a planted artefact")
    si.add_argument("--artefact", required=True, choices=["cleavage", "abundance", "multicys", "overlap", "none"])
    si.add_argument("--strength", type=float, required=True)
    si.add_argument("--seed", type=int, default=C.DEFAULT_SEED)
    si.add_argument("--n-proteins", type=int, default=600)
    si.add_argument("--out-dir", required=True)
    a = p.parse_args(argv)
    a._argv = argv[1:]
    return {"audit": _audit, "schema": _schema, "simulate": _simulate}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
