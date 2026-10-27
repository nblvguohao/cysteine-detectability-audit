"""Rebuild members from retained census/FASTA, and groups from fixed components.

No saved fold or prediction is written into training inputs. This is not a new
workbook ingestion or de novo homology clustering.
"""
import csv
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
COHORTS = ["qtrp_S1_ph5", "qtrp_S2_ph5", "qpers_sid_tierB",
           "cysboost2019_human_sno_hela", "cysboost2019_human_sno_shsy5y",
           "abiotech2025_ath_sno", "fps2020_ath_sulfenyl"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fasta(path):
    result, name, chunks = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                result[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        result[name] = "".join(chunks)
    return result


def prepare():
    inp = ROOT / "inputs"
    for name, record in json.loads((inp / "source_manifest.json").read_text()).items():
        assert sha(inp / name) == record["sha256"], name
    human = fasta(inp / "human_additional.fasta")
    for acc, sequence in fasta(inp / "human_qtrp_qpers.fasta").items():
        assert acc not in human or human[acc] == sequence
        human[acc] = sequence
    ath = fasta(inp / "arabidopsis_tair.fasta")
    with (inp / "homology_components.csv").open(encoding="utf-8-sig") as h:
        components = {r["accession"]: r["component"] for r in csv.DictReader(h)}
    census = pd.read_csv(inp / "census.csv")
    legacy = pd.read_csv(inp / "legacy_comparison_ONLY.tsv.gz", sep="\t",
                         usecols=["cohort", "acc", "site", "unit", "is_pos", "in_A", "in_B"])
    target = inp / "prepared"
    target.mkdir(exist_ok=True)
    audits = []
    for cohort in COHORTS:
        c = census[census.dataset_id.eq(cohort)]
        assert len(c) and c.proteome_key.nunique() == 1
        species = c.proteome_key.iloc[0]
        reference = human if species == "hsa" else ath
        c = c[c.role.isin(["positive", "observed_unmodified"])]
        # Census role spelling is checked explicitly; never silently drop B.
        all_roles = set(census.loc[census.dataset_id.eq(cohort), "role"])
        assert all_roles <= {"positive", "observed_unmodified"}, (cohort, all_roles)
        pos = set(zip(c.loc[c.role.eq("positive"), "accession"], c.loc[c.role.eq("positive"), "site"]))
        neg_b = set(zip(c.loc[c.role.eq("observed_unmodified"), "accession"], c.loc[c.role.eq("observed_unmodified"), "site"])) - pos
        def seq_for(acc):
            return reference.get(acc, reference.get(acc.split("-")[0]))
        def valid(key):
            sequence = seq_for(key[0])
            return sequence is not None and 1 <= int(key[1]) <= len(sequence) and sequence[int(key[1])-1] == "C"
        pos = {k for k in pos if valid(k)}
        neg_b = {k for k in neg_b if valid(k)}
        neg_a = {(acc, i+1) for acc in sorted({a for a, _ in pos})
                 for i, residue in enumerate(seq_for(acc)) if residue == "C"} - pos
        rows = []
        for acc, site in sorted(pos | neg_a | neg_b):
            if species == "hsa":
                component = components.get(acc, components.get(acc.split("-")[0]))
                unit = "hsa_" + component if component is not None else "acc_" + acc
            else:
                unit = "acc_" + acc
            key = acc, site
            rows.append(dict(acc=acc, site=int(site), unit=unit, is_pos=int(key in pos),
                             in_A=int(key in neg_a), in_B=int(key in neg_b)))
        sites = pd.DataFrame(rows)
        old = legacy[legacy.cohort.eq(cohort)].sort_values(["acc", "site"]).reset_index(drop=True)
        assert sites[old[sites.columns].columns].equals(old[sites.columns]), cohort
        folder = target / cohort
        folder.mkdir(exist_ok=True)
        sites.to_csv(folder / "sites.tsv", sep="\t", index=False)
        with (folder / "reference.fasta").open("w") as h:
            for acc in sorted(sites.acc.unique()):
                h.write(">" + acc + "\n" + seq_for(acc) + "\n")
        audit = dict(cohort=cohort, rows=len(sites), positives=len(pos), neg_A=len(neg_a),
                     neg_B=len(neg_b), proteins=sites.acc.nunique(), units=sites.unit.nunique(),
                     species=species, membership_and_units_equal_legacy=True,
                     accession_fallback_proteins=int(sites.loc[sites.unit.str.startswith('acc_'), 'acc'].nunique()),
                     sites_sha256=sha(folder / "sites.tsv"), fasta_sha256=sha(folder / "reference.fasta"))
        audits.append(audit)
        print(json.dumps(audit), flush=True)
    (target / "preparation_audit.json").write_text(json.dumps(audits, indent=2) + "\n")


if __name__ == "__main__":
    prepare()
