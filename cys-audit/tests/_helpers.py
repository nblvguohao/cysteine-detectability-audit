import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))

from cys_audit import synthetic  # noqa: E402
from cys_audit.io import read_sites  # noqa: E402

REPS = 400  # unit tests use fewer bootstrap replicates than the 5000 default; decisions here are not borderline


def cfg(**kw):
    base = {"protease": "trypsin", "background": "proteome", "seed": 7, "reps": REPS, "modification": "persulfidation",
            "precursor_ppm": 4.5, "search_mods": ["Sulfide"], "identity_readout": "direct_mass"}
    base.update(kw)
    return base


def synth_dataset(artefact, strength, seed=11, n_proteins=500, expand=False):
    d = tempfile.mkdtemp(prefix="cysaudit_test_")
    synthetic.simulate(artefact, strength, seed, d, n_proteins=n_proteins)
    return read_sites(os.path.join(d, "sites.tsv"), fasta=os.path.join(d, "proteome.fasta"), expand_background=expand), d


def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
