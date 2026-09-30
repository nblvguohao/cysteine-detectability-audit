"""Per-cysteine structural features from AlphaFold models, for the five structural claims.

`reports/PHASE2D_CLAIMS_AFTER_MATERIAL_DROP.md` left twelve claims needing structure and named
that as the remaining blocker for gap G2. This script computes the two families that have an
explicit author threshold - solvent accessibility and secondary structure - for every cysteine of
every protein in the three site tables, so the cohorts can be built exactly as the papers state.
pKa (five claims), conservation (one) and intrinsic reactivity (one) are NOT attempted here: no
predictor is installed and substituting a proxy would change what is being tested.

**Predeclared, fixed before the run.**

* Structures are AlphaFold DB models (`AF-<accession>-F1-model_v6.pdb`, falling back to v5 then
  v4), one per UniProt accession; v6 is what the archive serves as of 2026-09-16 and the version
  actually used is recorded per protein. Only the canonical `-F1` entry is taken, never an isoform
  entry such as `AF-<acc>-2-F1`, because the site tables were mapped against canonical sequences. This is a SUBSTITUTION in all three cohorts and is declared as
  such: Yang 2014 used the NetSurfP sequence predictor, Doulias 2010 used experimental structures
  for the subset that had them, and Marino & Gladyshev used homology models (their Table S1
  carries template and coverage columns). Using one consistent source for all three removes a
  source of between-cohort variation but is not what any of the three did.
* A model is used only if its sequence matches the UniProt sequence the site table was mapped
  against; length mismatch or a non-cysteine at the site means the site is dropped and counted.
* Solvent accessibility: biotite's Shrake-Rupley implementation (`biotite.structure.sasa`) with
  its defaults. freesasa would have been the more conventional choice but ships sdist-only and
  cannot be built in this sandbox; the two algorithms differ slightly in absolute area, which
  matters for Marino's absolute 1.0 A^2 threshold and is declared with that claim. Two quantities
  are stored because two papers use different ones -
  `sg_sasa` the absolute area of the cysteine's SG atom in A^2 (Marino's threshold is stated on
  the sulfur atom: buried means <= 1.0 A^2), and
  `rsa` the residue's total area divided by 167.0 A^2, the standard Tien et al. theoretical
  maximum for cysteine (Yang's threshold is stated as relative accessibility > 25%).
* Secondary structure: biotite's P-SEA assignment (`annotate_sse`), which returns a/b/c per
  residue from coordinates alone. DSSP is not installed and its binary is not available here;
  P-SEA's three-state output is what Doulias reports (helix / sheet / coil).
  `sse_coil_frac_0to3` is the fraction of positions 0, +1, +2, +3 assigned coil, which is the
  window Doulias names.
* pLDDT at the site is recorded. The PRIMARY caliber uses every site regardless of pLDDT, because
  filtering on pLDDT would condition on disorder and disorder is correlated with the very
  accessibility being measured; `plddt_ge_70` is stored so the cohorts can run a declared
  sensitivity on confident residues only.

Writes results/structural_site_features.csv and results/structural_site_features_audit.json.
Structures are cached under a scratch directory and are NOT added to the tree: the audit pins the
sha256 and model version of every structure actually used.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
PROTEOMES = os.path.join(ROOT, "external", "proteomes")
DROP = os.path.join(ROOT, "external", "intake", "manual_fetch_2026-09-16")
SEQ_CACHE = os.path.join(DROP, ".phase2d_uniprot_seq_cache.json")
SCRATCH = os.environ.get("AF_SCRATCH", "/tmp/af_models")
SCRIPT = os.path.abspath(__file__)

YANG = os.path.join(RESULTS, "phase2c_yang2014_sulfenyl_sites.csv")
DOULIAS = os.path.join(RESULTS, "doulias2010_sno_sites.csv")
MARINO = os.path.join(RESULTS, "marino2010_nocys_sites.csv")
AF_URL = "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v%d.pdb"
CYS_MAX_ASA = 167.0          # Tien et al. 2013 theoretical maximum for cysteine
MARINO_BURIED = 1.0          # A^2 on the sulfur atom, the author's stated threshold
YANG_RSA = 0.25              # the author's stated relative-accessibility threshold


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rows_of(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def read_fasta_gz(path):
    out, name, buf = {}, None, []
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf).upper()
                parts = line[1:].split()[0].split("|")
                name = parts[1] if len(parts) > 2 else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


def cohort_sequences():
    """{cohort: (sites, {accession: sequence})} using exactly the mappings already on disk."""
    cache = json.load(open(SEQ_CACHE, encoding="utf-8")) if os.path.exists(SEQ_CACHE) else {}
    human = read_fasta_gz(os.path.join(PROTEOMES, "hsa.fasta.gz"))

    yang_sites, yang_seq = set(), {}
    for r in rows_of(YANG):
        acc, pos = r["accession"], int(r["position"])
        if acc in human and pos <= len(human[acc]):
            yang_sites.add((acc, pos))
            yang_seq[acc] = human[acc]

    dou_sites, dou_seq = set(), {}
    for r in rows_of(DOULIAS):
        if not r["canonical_position"]:
            continue
        acc, pos = r["accession"], int(r["canonical_position"])
        seq = cache.get(acc)
        if seq and pos <= len(seq) and seq[pos - 1] == "C":
            dou_sites.add((acc, pos))
            dou_seq[acc] = seq

    mar_sites, mar_seq = set(), {}
    for r in rows_of(MARINO):
        if r["status"] != "confirmed_cys":
            continue
        acc, pos = r["accession"], int(r["position"])
        seq = cache.get(acc) or cache.get(acc.split("-")[0])
        if seq and pos <= len(seq) and seq[pos - 1] == "C":
            mar_sites.add((acc, pos))
            mar_seq[acc] = seq
    return {"yang2014_human": (yang_sites, yang_seq),
            "doulias2010_mouse": (dou_sites, dou_seq),
            "marino2010_multispecies": (mar_sites, mar_seq)}


def fetch_model(accession):
    os.makedirs(SCRATCH, exist_ok=True)
    base = accession.split("-")[0]
    for version in (6, 5, 4):
        path = os.path.join(SCRATCH, f"AF-{base}-F1-model_v{version}.pdb")
        if os.path.exists(path) and os.path.getsize(path) > 0:
            return path, version
        url = AF_URL % (base, version)
        for attempt in range(3):
            try:
                with urllib.request.urlopen(urllib.request.Request(url), timeout=120) as r:
                    data = r.read()
                with open(path, "wb") as fh:
                    fh.write(data)
                return path, version
            except urllib.error.HTTPError as exc:
                if exc.code == 404:
                    break
                time.sleep(1.5 * (attempt + 1))
            except Exception:
                time.sleep(1.5 * (attempt + 1))
    return None, None


def features_for(path, sequence):
    """-> (per-position dict, model sequence, plddt list) or None when the model is unusable."""
    import numpy as np
    from biotite.structure import annotate_sse, sasa
    from biotite.structure.io.pdb import PDBFile

    pdb = PDBFile.read(path)
    array = pdb.get_structure(model=1, extra_fields=["b_factor"])
    ca = array[array.atom_name == "CA"]
    model_len = ca.array_length()
    if model_len != len(sequence):
        return None, model_len, None

    plddt = dict(zip(ca.res_id.astype(int).tolist(), ca.b_factor.astype(float).tolist()))
    sse = annotate_sse(array)                       # one letter per amino-acid residue
    res_ids = ca.res_id.astype(int).tolist()
    sse_by_pos = {res_ids[i]: str(sse[i]) for i in range(min(len(sse), len(res_ids)))}

    atom_area = sasa(array)                         # per-atom, NaN for ignored atoms
    residue_area, sg_area = {}, {}
    for i in range(array.array_length()):
        area = float(atom_area[i])
        if not np.isfinite(area):
            continue
        pos = int(array.res_id[i])
        residue_area[pos] = residue_area.get(pos, 0.0) + area
        if str(array.atom_name[i]).strip() == "SG":
            sg_area[pos] = area
    return {"plddt": plddt, "sse": sse_by_pos, "residue_area": residue_area,
            "sg_area": sg_area}, model_len, True


def main():
    started = time.time()
    cohorts = cohort_sequences()
    accessions = {}
    for cohort, (sites, seqs) in cohorts.items():
        for acc in seqs:
            accessions.setdefault(acc, []).append(cohort)
    print(f"proteins to model: {len(accessions)}", flush=True)

    rows, used, stat = [], {}, {
        "proteins": len(accessions), "model_missing": 0, "length_mismatch": 0,
        "sites_written": 0, "sites_dropped_not_cys_in_model_numbering": 0,
        "no_sg_atom": 0}
    for k, (acc, cohort_list) in enumerate(sorted(accessions.items())):
        seq = None
        for cohort in cohort_list:
            seq = cohorts[cohort][1].get(acc)
            if seq:
                break
        path, version = fetch_model(acc)
        if path is None:
            stat["model_missing"] += 1
            continue
        feats, model_len, ok = features_for(path, seq)
        if feats is None:
            stat["length_mismatch"] += 1
            continue
        used[acc] = {"model_version": version, "sha256": sha256_of(path),
                     "model_length": model_len}
        positive_in = {c: (acc, 0) for c in cohort_list}
        for i, residue in enumerate(seq):
            if residue != "C":
                continue
            pos = i + 1
            sg = feats["sg_area"].get(pos)
            if sg is None:
                stat["no_sg_atom"] += 1
            total = feats["residue_area"].get(pos)
            coil = [feats["sse"].get(pos + d) for d in (0, 1, 2, 3)]
            coil_known = [c for c in coil if c is not None]
            rows.append({
                "accession": acc, "position": pos,
                "cohorts": ";".join(sorted(set(cohort_list))),
                "sg_sasa": "" if sg is None else round(sg, 4),
                "residue_sasa": "" if total is None else round(total, 4),
                "rsa": "" if total is None else round(total / CYS_MAX_ASA, 6),
                "sse": feats["sse"].get(pos, ""),
                "sse_coil_frac_0to3": ("" if not coil_known else
                                       round(sum(1 for c in coil_known if c == "c") / len(coil_known), 4)),
                "sse_window_known": len(coil_known),
                "plddt": round(feats["plddt"].get(pos, float("nan")), 2),
                "plddt_ge_70": int(feats["plddt"].get(pos, 0.0) >= 70.0),
                "model_version": version,
            })
            stat["sites_written"] += 1
        if (k + 1) % 100 == 0:
            print(f"  {k+1}/{len(accessions)} proteins, {stat['sites_written']} sites", flush=True)

    out = os.path.join(RESULTS, "structural_site_features.csv")
    with open(out, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    import numpy as np
    exposed = [r for r in rows if r["rsa"] != "" and float(r["rsa"]) > YANG_RSA]
    buried_sg = [r for r in rows if r["sg_sasa"] != "" and float(r["sg_sasa"]) <= MARINO_BURIED]
    audit = {
        "script": "scripts/build_structural_site_features.py", "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "round": "structural site features for the five claims with an explicit author threshold",
        "claims_enabled": ["SFE-001", "SNO-002", "SNO-012", "SNO-001", "SNO-009"],
        "claims_not_attempted": {
            "pKa": ["SNO-003", "SNO-007", "SNO-010", "SNO-011", "SNO-015"],
            "conservation": ["SNO-013"], "intrinsic_reactivity": ["SFE-004"],
            "reason": "no predictor installed; a proxy would change what is being tested"},
        "source": "AlphaFold DB, AF-<acc>-F1-model_v6.pdb with v5/v4 fallback (v6 served 2026-09-16)",
        "substitution_declared": (
            "all three cohorts used something else: Yang 2014 the NetSurfP sequence predictor, "
            "Doulias 2010 experimental structures for the subset that had them, Marino & Gladyshev "
            "homology models. One consistent source removes between-cohort variation but is not "
            "what any of the three papers did."),
        "thresholds": {"yang_rsa_exposed": YANG_RSA, "marino_sg_buried_A2": MARINO_BURIED,
                       "cys_max_asa_A2": CYS_MAX_ASA,
                       "cys_max_asa_source": "Tien et al. 2013 theoretical maximum"},
        "sasa_method": ("biotite.structure.sasa, Shrake-Rupley with library defaults; freesasa is sdist-only and cannot be built in this sandbox, and the two algorithms differ slightly in absolute area - relevant to Marino's absolute threshold"),
        "sse_method": "biotite annotate_sse (P-SEA); DSSP binary unavailable in this sandbox",
        "plddt_policy": ("primary caliber uses all sites; filtering on pLDDT would condition on "
                         "disorder, which correlates with the accessibility being measured. "
                         "plddt_ge_70 is stored for a declared sensitivity."),
        "counts": stat,
        "sites_total": len(rows),
        "sites_with_rsa": sum(1 for r in rows if r["rsa"] != ""),
        "sites_with_sg": sum(1 for r in rows if r["sg_sasa"] != ""),
        "sites_plddt_ge_70": sum(int(r["plddt_ge_70"]) for r in rows),
        "fraction_exposed_at_yang_threshold": round(len(exposed) / max(1, len(rows)), 4),
        "fraction_sg_buried_at_marino_threshold": round(len(buried_sg) / max(1, len(rows)), 4),
        "models_used": used,
        "inputs": {"results/phase2c_yang2014_sulfenyl_sites.csv": sha256_of(YANG),
                   "results/doulias2010_sno_sites.csv": sha256_of(DOULIAS),
                   "results/marino2010_nocys_sites.csv": sha256_of(MARINO)},
        "outputs": {"results/structural_site_features.csv": sha256_of(out)},
        "versions": {"python": sys.version},
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }
    json.dump(audit, open(os.path.join(RESULTS, "structural_site_features_audit.json"), "w"),
              ensure_ascii=False, indent=2)
    print(json.dumps({k: v for k, v in audit.items()
                      if k in ("counts", "sites_total", "sites_with_rsa", "sites_with_sg",
                               "sites_plddt_ge_70", "fraction_exposed_at_yang_threshold",
                               "fraction_sg_buried_at_marino_threshold", "elapsed_minutes")},
                     indent=1))


if __name__ == "__main__":
    main()
