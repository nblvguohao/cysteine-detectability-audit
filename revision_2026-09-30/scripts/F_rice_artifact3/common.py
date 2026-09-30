"""Shared paths and loaders for item F_rice_artifact3 (POST HOC revision analysis, 2026-09-30).

Everything in this folder is a post hoc analysis written in response to pre-submission review
criticism of Artifact 3 and of the rice cohort size. Nothing here is registered or pre-specified.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import pathlib
import re
import sys

W = pathlib.Path("C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30")
REPO = pathlib.Path("C:/Users/admin/Desktop/小论文/_cys_repo_work/repo")
CYS_AUDIT_SRC = pathlib.Path("C:/Users/admin/Desktop/小论文/_cys_repo_work/public/cys-audit/src")
ITEM = "F_rice_artifact3"
SCRIPTS = W / "scripts" / ITEM
RES = W / "results" / ITEM
FETCHED = RES / "fetched"

PEP = W / "external" / "SS-all-peptides.tsv"
PRO = W / "external" / "20240326_061934_xyj_proteome_Report.tsv"
FASTA_REF = W / "external" / "UP000059680_39947.fasta.gz"
STORED_A3 = W / "inputs" / "repo_results" / "artefact3_osa_abundance_recomputed_2026-09-19.csv"
STORED_A3_SENS = W / "inputs" / "repo_results" / "artefact3_positive_set_sensitivity_2026-09-19.csv"
# read-only repo products used for the reconciliation (the self-audit cohort's stored site table)
REPO_SITE_SCORES = REPO / "external" / "v2_retrospective_site_scores.csv"
REPO_EFFECT_SIZES = REPO / "results" / "artefact3_effect_sizes_2026-09-19.csv"
REPO_RECON_ARMS = REPO / "results" / "artefact3_reconstructed_original_arms_2026-09-19.csv"
REPO_A3_SCRIPT = REPO / "scripts" / "recompute_artefact3_osa_abundance_2026-09-19.py"

SEED = 20260930          # new seed for every post hoc resampling in this item
SEED_ORIGINAL = 20260919  # seed of the stored Artifact 3 run, used only to reproduce it

csv.field_size_limit(10 ** 8)


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    # the stored script opened these with the platform default; they are ASCII in practice
    with open(path, encoding="utf-8", errors="strict", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta_any(path):
    """Return {accession: sequence} from a UniProt FASTA (gz or plain). Accession = 2nd '|' field."""
    op = gzip.open if str(path).endswith(".gz") else open
    out, acc, buf = {}, None, []
    with op(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc is not None:
                    out[acc] = "".join(buf)
                m = re.match(r">\w+\|([^|]+)\|", line)
                acc = m.group(1) if m else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc is not None:
        out[acc] = "".join(buf)
    return out


def split_accs(s):
    return [a.strip() for a in (s or "").split(";") if a.strip()]


def is_decoy_or_contaminant(acc: str) -> bool:
    return acc.startswith("REV__") or acc.startswith("CON__")


def write_json(path, obj):
    pathlib.Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=1, sort_keys=True,
                                             default=str), encoding="utf-8")


def load_all_sequences():
    """Reference proteome plus the sequences fetched from UniProt REST in fetch_sequences.py.

    Returns (seqs, source) where source[acc] in {'UP000059680_2026_03', 'uniprot_rest_active',
    'unisave_last_version'}. A fetched entry is stored under the REQUESTED accession.
    """
    seqs = read_fasta_any(FASTA_REF)
    source = {a: "UP000059680_2026_03" for a in seqs}
    fmap = FETCHED / "fetched_sequences.tsv"
    if fmap.exists():
        with open(fmap, encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["sequence"] and r["requested_accession"] not in seqs:
                    seqs[r["requested_accession"]] = r["sequence"]
                    source[r["requested_accession"]] = r["source"]
    return seqs, source
