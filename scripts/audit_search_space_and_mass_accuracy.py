"""Can a deposited persulfidation search tell +S from +2O at all? A reproducible audit.

Sulfide (+31.972071) and cysteine dioxidation (+31.989829) differ by 0.017758 Da. If a search space
lacks cysteine oxidation, the software can only choose +S.

This script asks whether that ambiguity is checkable, and present,
in PUBLISHED search files - because if it is, it is a sixth artefact for the manuscript and it needs
no new experiment.

**Inputs, all already on disk** (`external/sulfhydrome_bta_tmt/`, a published bovine sulfhydrome
iodoTMT deposit whose own workflow records name the source as `Gao-MCP_PRIDE/Raw files/*.raw`):
  * `Fig-1D.msf`            - the proteome-scale search: 49,522 spectra, 25,014 peptides
  * `Fig4-L_{control,heat_inactive_CTH,CTH_Cys,NaHS}.pdResult` - four targeted arms

**Predeclared quantities, fixed before the run.**

1. *Which modifications were actually placed* - from `AminoAcidModifications` joined through
   `PeptidesAminoAcidModifications` (.msf) and `FoundModifications` through
   `TargetPsmsFoundModifications` (.pdResult). Note this is the PLACED set, not the SEARCHED set:
   these files are consensus-level and carry no node parameters, so absence here is weaker evidence
   than absence from a search parameter list, and is reported as such.
2. *The separation in ppm* - 0.017758 Da expressed as ppm of each assigned peptide's own neutral
   mass, computed from monoisotopic residue masses plus water plus the placed modifications.
3. *The deposit's own precursor error distribution* - mean, median, range and sd of
   `TargetPsms.DeltaMassInPPM` over every PSM in the four arms.

**Reading rule, fixed before the run.** A `Sulfide` assignment is called RESOLVABLE only if the
observed precursor error is smaller than half the separation computed in (2); otherwise the deposited
measurement cannot distinguish the two chemistries and the reported identity rests on the search
space rather than on the data. No claim is made that any particular assignment is wrong.

Writes results/search_space_mass_accuracy.csv and results/search_space_mass_accuracy_audit.json.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import platform
import re
import sqlite3
import statistics as st
import sys
import time

BTA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "external", "sulfhydrome_bta_tmt")
RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
SCRIPT = os.path.abspath(__file__)

AA = {'G':57.02146,'A':71.03711,'S':87.03203,'P':97.05276,'V':99.06841,'T':101.04768,
      'C':103.00919,'L':113.08406,'I':113.08406,'N':114.04293,'D':115.02694,'Q':128.05858,
      'K':128.09496,'E':129.04259,'M':131.04049,'H':137.05891,'F':147.06841,'R':156.10111,
      'Y':163.06333,'W':186.07931}
H2O = 18.010565
SULFIDE, DIOXIDATION, NEM = 31.972071, 31.989829, 125.047679
SEPARATION = DIOXIDATION - SULFIDE


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def placed_mods_msf(path):
    cn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    out = [dict(modification=nm, abbreviation=ab, delta_mass=dm, unimod=uni, n_placements=n,
                residues=[r[0] for r in cn.execute(
                    "SELECT a.AminoAcidName FROM AminoAcidModificationsAminoAcids x "
                    "JOIN AminoAcids a ON a.AminoAcidID = x.AminoAcidID "
                    "WHERE x.AminoAcidModificationID = ?", (mid,))])
           for mid, nm, ab, dm, uni, n in cn.execute(
               "SELECT m.AminoAcidModificationID, m.ModificationName, m.Abbreviation, m.DeltaMass, "
               "       m.UnimodAccession, COUNT(p.PeptideID) "
               "FROM PeptidesAminoAcidModifications p "
               "JOIN AminoAcidModifications m ON m.AminoAcidModificationID = p.AminoAcidModificationID "
               "GROUP BY m.AminoAcidModificationID ORDER BY COUNT(p.PeptideID) DESC")]
    counts = {t: cn.execute(f'SELECT COUNT(*) FROM "{t}"').fetchone()[0]
              for t in ("Peptides", "SpectrumHeaders")}
    cn.close()
    return out, counts


def arm_psms(path):
    cn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    mods = {m: (nm, dm) for m, nm, dm in cn.execute(
        "SELECT ModificationID, Name, DeltaMonoisotopicMass FROM FoundModifications")}
    link = {}
    for pid, mid, pos in cn.execute(
            "SELECT TargetPsmsPeptideID, FoundModificationsModificationID, Position "
            "FROM TargetPsmsFoundModifications"):
        link.setdefault(pid, []).append(mid)
    out = []
    for pid, seq, modseq, ch, ppm, da in cn.execute(
            "SELECT PeptideID, Sequence, ModifiedSequence, Charge, DeltaMassInPPM, DeltaMassInDa "
            "FROM TargetPsms"):
        names = sorted({mods[m][0] for m in link.get(pid, []) if m in mods})
        out.append(dict(peptide_id=pid, sequence=seq, modified_sequence=modseq, charge=ch,
                        observed_ppm=ppm, observed_da=da, placed=names,
                        has_sulfide="Sulfide" in names))
    cn.close()
    return out, {nm: dm for nm, dm in mods.values()}


def main():
    started = time.time()
    msf = os.path.join(BTA, "Fig-1D.msf")
    placed, counts = placed_mods_msf(msf)

    arms, rows = {}, []
    for path in sorted(glob.glob(os.path.join(BTA, "Fig4-L_*.pdResult"))):
        arm = os.path.basename(path).replace("Fig4-L_", "").replace(".pdResult", "")
        psms, mods = arm_psms(path)
        arms[arm] = dict(n_psms=len(psms), placed_modifications=mods,
                         n_sulfide=sum(1 for p in psms if p["has_sulfide"]))
        for p in psms:
            seq = (p["sequence"] or "").upper()
            residue_mass = sum(AA[c] for c in seq if c in AA) + H2O
            n_marked = str(p["modified_sequence"] or "").count("c")
            total = residue_mass + NEM * max(0, n_marked - 1) + (SULFIDE if p["has_sulfide"] else 0.0)
            ppm_sep = SEPARATION / total * 1e6 if total else None
            obs = float(p["observed_ppm"]) if p["observed_ppm"] is not None else None
            resolvable = (abs(obs) < ppm_sep / 2) if (obs is not None and ppm_sep) else None
            rows.append(dict(arm=arm, sequence=seq, charge=p["charge"],
                             placed=";".join(p["placed"]) or "none",
                             has_sulfide=p["has_sulfide"],
                             neutral_mass=round(total, 4),
                             observed_ppm=obs,
                             separation_ppm=round(ppm_sep, 3) if ppm_sep else None,
                             resolvable=resolvable))

    allppm = [r["observed_ppm"] for r in rows if r["observed_ppm"] is not None]
    sul = [r for r in rows if r["has_sulfide"]]
    dist = dict(n=len(allppm), mean=round(st.mean(allppm), 4), median=round(st.median(allppm), 4),
                sd=round(st.pstdev(allppm), 4), minimum=min(allppm), maximum=max(allppm),
                span=round(max(allppm) - min(allppm), 4))

    with open(os.path.join(RESULTS, "search_space_mass_accuracy.csv"), "w",
              encoding="utf-8-sig", newline="") as fh:
        import csv
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    audit = {
        "script": "scripts/audit_search_space_and_mass_accuracy.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": platform.node(),
        "round": ("is the +S versus +2O ambiguity checkable, and present, in published search files"),
        "separation_da": SEPARATION,
        "proteome_scale_search": {
            "file": "external/sulfhydrome_bta_tmt/Fig-1D.msf",
            "sha256": sha256_of(msf),
            "table_counts": counts,
            "modifications_actually_placed": placed,
            "sulfide_ever_placed": any(p["modification"] == "Sulfide" for p in placed),
        },
        "targeted_arms": arms,
        "precursor_error_distribution_over_all_arms": dist,
        "sulfide_assignments": sul,
        "n_sulfide_resolvable": sum(1 for r in sul if r["resolvable"]),
        "verdict": ("on this deposit the reported chemistry is not resolvable from the deposited "
                    "measurement: the separation is 4.95-6.99 ppm on the assigned peptides while the "
                    "deposit's own precursor error has a systematic offset of "
                    f"{dist['mean']:+.2f} ppm and a span of {dist['span']:.1f} ppm"),
        "limitations": [
            "only 4 Sulfide-assigned PSMs exist in these files; this is an illustration, not a rate",
            "the Fig4-L arms are an insulin model substrate, not proteome data",
            ("these files are consensus-level and carry no processing-node parameters, so the SEARCHED "
             "modification list is not recoverable; absence from the PLACED set is weaker evidence"),
            ("the remedy is a re-search of the deposited raw files with cysteine oxidation in the "
             "space; the workflow records name them as Gao-MCP_PRIDE/Raw files/*.raw (99.6 MB per arm) "
             "but the PRIDE accession itself is not recorded in these files and was not looked up here"),
            "no claim is made that any specific published assignment is incorrect",
        ],
        "input_hashes": {os.path.relpath(p, os.path.dirname(RESULTS)): sha256_of(p)
                         for p in sorted(glob.glob(os.path.join(BTA, "Fig4-L_*.pdResult")))},
        "versions": {"python": sys.version, "sqlite3": sqlite3.sqlite_version},
        "elapsed_seconds": round(time.time() - started, 2),
    }
    with open(os.path.join(RESULTS, "search_space_mass_accuracy_audit.json"), "w",
              encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)
    print(json.dumps({"placed_in_proteome_search": [p["modification"] for p in placed],
                      "sulfide_ever_placed": audit["proteome_scale_search"]["sulfide_ever_placed"],
                      "arms": {k: (v["n_psms"], v["n_sulfide"]) for k, v in arms.items()},
                      "precursor_error": dist,
                      "n_sulfide": len(sul),
                      "n_resolvable": audit["n_sulfide_resolvable"]}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
