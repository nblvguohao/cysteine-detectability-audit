"""H_artifact5_docs step 3 - POST HOC: what precursor accuracy can decide, and the phospho/sulfo check.

POST HOC revision analysis in response to a pre-submission review; not registered, not pre-specified.

(a) Two different criteria have been used under the name 'resolvable' / 'resolving mass':
    * window criterion (Cys-Audit search_space check): the alternative falls outside the precursor window
      when its separation in ppm exceeds the tolerance, i.e. for M < M* = delta / (tol x 1e-6);
    * accuracy criterion (PXD015307 re-search): with a realised error SD s about the file's own
      baseline, the two hypotheses are separated by at least 2 SD on each side of the midpoint when
      delta/M x 1e6 >= 4 s, i.e. for M <= delta x 1e6 / (4 s).
    Both are computed here from stored values, per PXD015307 file (re-search baselines) and for the
    deposit's own targeted-arm PSMs, and applied to the theoretical cysteine peptides of the mouse
    reference proteome (UniProt 2026_03 canonical; tryptic, 7-30 residues, up to 2 missed cleavages),
    with one +32 modification and either no other cysteine modification or iodoTMT6plex on every
    other cysteine.
(b) The phosphorylation/sulfation extension is re-verified: delta from the released Cys-Audit constants,
    M* on the registered grid, and the PXD071110 share recomputed from the stored tool input with the
    released search_space check; the other three deposits are read from the stored table.
No random numbers are drawn.

Revision after verification (round 2): the window criterion is labelled as assuming an unbiased (zero-offset)
measurement, and offset-aware window limits are added for every file's own baseline median; the 2-SD accuracy
criterion is kept (M <= delta x 1e6 / (4 SD)) and is post hoc.
"""
from __future__ import annotations

import os
import sys
import types

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5_lib as L  # noqa: E402

sys.path.insert(0, L.CYS_AUDIT_SRC)
from cys_audit import constants as CA  # noqa: E402
from cys_audit.checks import searchspace as SS  # noqa: E402

AUDIT = os.path.join(L.W, "inputs", "repo_results", "pxd015307_research_audit.json")
FASTA = os.path.join(L.W, "external", "UP000000589_10090.fasta.gz")
PH_INPUT = os.path.join(L.REPO, "results", "phospho_sulfo_tool_input_2026-09-23.tsv")
PH_ROUND1 = os.path.join(L.REPO, "results", "phospho_sulfo_search_space_2026-09-23.csv")
PH_ROUND1_AUDIT = os.path.join(L.REPO, "results", "phospho_sulfo_search_space_2026-09-23_audit.json")
PH_EIGHT = os.path.join(L.REPO, "results", "phospho_sulfo_eight_deposits_2026-09-23.csv")
PH_EIGHT_AUDIT = os.path.join(L.REPO, "results", "phospho_sulfo_eight_deposits_2026-09-23_audit.json")
PH_SOURCE = os.path.join(L.W, "inputs", "mcp_package", "Source_Data_text_phospho_sulfo_eight_deposits.csv")
DEP_BASE = os.path.join(L.OUT, "deposit_targeted_arm_baselines.csv")   # written by step 2


def main():
    for f in ("constants.py", os.path.join("checks", "searchspace.py")):
        L.register(os.path.join(L.CYS_AUDIT_SRC, "cys_audit", f), "released Cys-Audit source", text=True)
    audit = L.read_json(AUDIT, "stored re-search audit")
    delta_cys = 2 * CA.MASS_O - CA.MASS_S
    delta_ps = abs(SS.DELTAS["phosphorylation"] - SS.DELTAS["sulfation"])

    # ---------- (a) theoretical mouse cysteine peptides ----------
    fasta = L.read_fasta_gz(FASTA, "UniProt mouse reference proteome 2026_03 (canonical)")
    seen = set()
    m1, m2 = [], []
    for _, s in fasta:
        for pep, _mc in L.tryptic_peptides(s):
            if "C" not in pep or pep in seen or any(c not in L.AA for c in pep):
                continue
            seen.add(pep)
            m = L.peptide_mass(pep)
            nc = pep.count("C")
            m1.append(m + L.SULFIDE)
            m2.append(m + L.SULFIDE + L.IODOTMT * (nc - 1))
    m1, m2 = np.array(m1), np.array(m2)

    rows = []

    def add(context, kind, param, mmax):
        rows.append({"context": context, "criterion": kind, "parameter": param, "mass_limit_da": mmax,
                     "fraction_theoretical_cys_peptides_at_or_below_plus32_only": float(np.mean(m1 <= mmax)),
                     "fraction_theoretical_cys_peptides_at_or_below_plus32_and_iodoTMT_on_other_cys": float(np.mean(m2 <= mmax))})

    for tol in (20.0, 10.0, 4.5, 2.0):
        add("any file, unbiased measurement (zero offset)", "window: alternative outside the precursor window (M < M*)",
            f"tolerance {tol} ppm", delta_cys / (tol * 1e-6))
    dep = pd.read_csv(DEP_BASE)
    # revision after verification: offset-aware window limits. With a systematic offset b (ppm) the error of
    # the true interpretation sits at b; the alternative sits at b - s (true sulfide, dioxidation candidate)
    # or b + s (true dioxidation, sulfide candidate), s = delta/M x 1e6. It is inside a +/-t window when
    # s <= t + b (true sulfide) or s <= t - b (true dioxidation); below the corresponding mass it is excluded.
    tol = 10.0
    offsets = [(f"PXD015307 re-search {arm}", a["baseline_median_ppm"]) for arm, a in audit["per_arm"].items()]
    offsets += [(f"PXD015307 deposit's own search, {b.arm} ({b.subset})", b.median_ppm_non_sulfide)
                for _, b in dep.iterrows()]
    for ctx, b in offsets:
        for truth, lim in (("true sulfide: dioxidation alternative inside the window", tol + b),
                           ("true dioxidation: sulfide alternative inside the window", tol - b)):
            add(f"{ctx}, offset {b:+.3f} ppm", f"window, offset-aware ({truth} above the limit)",
                f"tolerance {tol} ppm, offset {b:+.3f} ppm", delta_cys * 1e6 / lim if lim > 0 else np.inf)
    for z in (2, 3):
        crit = f"accuracy: separation >= {2 * z} SD ({z} SD each side of the midpoint)"
        for arm, a in audit["per_arm"].items():
            sd = a["baseline_sd_ppm"]
            add(f"PXD015307 re-search {arm} (n baseline PSMs {a['n_baseline_psms_no_plus32']})",
                crit, f"SD {sd} ppm", delta_cys * 1e6 / (2 * z * sd))
        for _, b in dep.iterrows():
            if not np.isfinite(b.sd_ppm_non_sulfide):
                continue
            add(f"PXD015307 deposit's own search, {b.arm} ({b.subset}; n = {b.n_non_sulfide_psms})",
                crit, f"SD {b.sd_ppm_non_sulfide:.4f} ppm", delta_cys * 1e6 / (2 * z * b.sd_ppm_non_sulfide))
    acc = pd.DataFrame(rows)

    # ---------- (b) phospho / sulfo ----------
    ph = pd.read_csv(L.register(PH_INPUT, "stored Cys-Audit input for PXD071110 phosphopeptides"), sep="\t")
    r1 = pd.read_csv(L.register(PH_ROUND1, "stored PXD071110 phospho/sulfo result"))
    r1a = L.read_json(PH_ROUND1_AUDIT, "stored PXD071110 phospho/sulfo audit")
    r8 = pd.read_csv(L.register(PH_EIGHT, "stored eight-deposit phospho/sulfo result"))
    r8a = L.read_json(PH_EIGHT_AUDIT, "stored eight-deposit phospho/sulfo audit")
    src = pd.read_csv(L.register(PH_SOURCE, "manuscript Source Data (phospho/sulfo, MCP package copy)"))
    masses = ph.loc[ph.label == 1, "peptide_mass"].to_numpy(dtype=float)
    masses = masses[np.isfinite(masses)]
    prow = []
    for tol in (1.0, 2.0, 4.5, 10.0, 20.0):
        mstar = delta_ps / (tol * 1e-6)
        ds = types.SimpleNamespace(peptide_mass=ph.peptide_mass.to_numpy(dtype=float),
                                   label=ph.label.to_numpy())
        res = SS.run(ds, {"modification": "phosphorylation", "precursor_ppm": tol,
                          "search_mods": ["phosphorylation"], "identity_readout": "direct_mass"})
        stored = r1.loc[r1.tolerance_ppm == tol, "share_at_or_above_resolving_mass_primary"]
        prow.append({"deposit": "PXD071110", "tolerance_ppm": tol, "resolving_mass_da": mstar,
                     "n_phosphopeptides": int(masses.size),
                     "share_at_or_above_recomputed": float(np.mean(masses >= mstar)),
                     "share_from_released_tool": res["details"]["alternatives"][0]["share_positive_peptides_at_or_above_resolving_mass"],
                     "tool_status": res["status"],
                     "share_stored": float(stored.iloc[0]) if len(stored) else None,
                     "note": "tolerance 1 ppm is outside the registered grid (sensitivity only)" if tol == 1.0 else ""})
    for _, r in r8.iterrows():
        if r.accession == "PXD071110":
            continue
        prow.append({"deposit": r.accession, "tolerance_ppm": r.main_search_ppm, "resolving_mass_da":
                     (delta_ps / (r.main_search_ppm * 1e-6)) if pd.notna(r.main_search_ppm) else None,
                     "n_phosphopeptides": r.n_primary, "share_at_or_above_recomputed": None,
                     "share_from_released_tool": None, "tool_status": r.verdict,
                     "share_stored": r.share_at_main_tolerance,
                     "note": "stored value only; peptide masses not in the local copy"})
    phdf = pd.DataFrame(prow)
    mods_recorded = {acc: d.get("modifications_searched") for acc, d in r8a["per_deposit"].items()}
    # PXD071110's list is recorded in the first-round audit, not in the eight-deposit audit
    mods_recorded["PXD071110"] = r1a["search_configuration_read_from_the_deposit"]["modifications_searched"]

    summary = {
        "label": L.POSTHOC_LABEL,
        "delta_cys_da": delta_cys, "delta_phospho_sulfo_da": delta_ps,
        "n_theoretical_mouse_cys_peptides": int(m1.size),
        "median_mass_plus32_only": float(np.median(m1)),
        "median_mass_plus32_iodoTMT_other_cys": float(np.median(m2)),
        "pxd071110_share_4p5_recomputed": float(phdf[(phdf.deposit == "PXD071110") & (phdf.tolerance_ppm == 4.5)].share_at_or_above_recomputed.iloc[0]),
        "pxd071110_share_4p5_stored": float(r1a["result"]["share_at_or_above_resolving_mass"]),
        "source_data_matches_stored_eight": bool((src.share_at_main_tolerance.fillna(-1).round(6).values
                                                  == r8.share_at_main_tolerance.fillna(-1).round(6).values).all()),
        "modification_lists_recorded_in_eight_deposit_audit": mods_recorded,
        "sulfation_absent_where_recorded": all(("Sulfation" not in " ".join(v or [])) for v in mods_recorded.values()),
        "n_deposits_with_recorded_modification_list": sum(1 for v in mods_recorded.values() if v),
        "heaviest_peptide_da": {acc: d.get("heaviest_peptide_da") for acc, d in r8a["per_deposit"].items()},
    }
    L.write_csv(acc, "precursor_accuracy_decidability.csv")
    L.write_csv(phdf, "phospho_sulfo_reverification.csv")
    L.write_json(summary, "mass_accuracy_phospho_summary.json")
    return summary, acc, phdf


if __name__ == "__main__":
    import json
    s, acc, phdf = main()
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.max_colwidth", 70)
    print(acc.to_string())
    print(phdf.to_string())
    print(json.dumps(s, indent=1, default=str))
