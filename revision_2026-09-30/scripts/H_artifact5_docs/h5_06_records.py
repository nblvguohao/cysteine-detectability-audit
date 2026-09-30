"""H_artifact5_docs step 6 - POST HOC: records consulted for the revised text (2026-09-30, round 2).

POST HOC revision analysis in response to a pre-submission review; not registered, not pre-specified.
Added in the revision after verification. Nothing is computed; the step extracts, from stored internal
records, the facts that the revised text and recommendation cite, so that each is traceable and hashed:

(a) where the phosphorylation/sulfation configuration was read from (the verifier found that the round-1
    proposal named parameters.txt as a source, which the stored records contradict);
(b) the unpublished collaborator re-search of recombinant PYL3 (the only instance in the records in which
    adding dioxidation to the search changed +32 assignments), its limits and its consent/use status;
(c) what the registration fixed before any search and which runs were discarded.
"""
from __future__ import annotations

import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5_lib as L  # noqa: E402

PH8_SCRIPT = os.path.join(L.REPO, "scripts", "analyse_phospho_sulfo_eight_deposits_2026-09-23.py")
PH1_SCRIPT = os.path.join(L.REPO, "scripts", "analyse_phospho_sulfo_search_space_2026-09-23.py")
PH1_AUDIT = os.path.join(L.REPO, "results", "phospho_sulfo_search_space_2026-09-23_audit.json")
PYL3_AUDIT = os.path.join(L.REPO, "results", "ms_research_persulfide_vs_dioxidation_audit.json")
PYL3_TABLE = os.path.join(L.REPO, "results", "ms_research_persulfide_vs_dioxidation.csv")
UNPUB = os.path.join(L.REPO, "results", "unpublished_data_provenance_2026-09-19.csv")
VALUE = os.path.join(L.REPO, "results", "private_data_value_assessment_2026-09-19.csv")
CONSENT = os.path.join(L.REPO, "results", "consent_grant_register_2026-09-17.csv")
REG = os.path.join(L.W, "inputs", "repo_reports", "pxd015307_research_preregistration_2026-09-17.json")
RPT = os.path.join(L.W, "inputs", "repo_reports", "PXD015307_RESEARCH_SEARCH_SPACE.md")


def main():
    ph8 = L.read_text(PH8_SCRIPT, "repository script, eight-deposit phospho/sulfo round (docstring: config source)")
    ph1 = L.read_text(PH1_SCRIPT, "repository script, PXD071110 phospho/sulfo round (docstring: config source)")
    ph1a = L.read_json(PH1_AUDIT, "stored PXD071110 phospho/sulfo audit")
    pyl = L.read_json(PYL3_AUDIT, "stored audit of the collaborator PYL3 re-search (transcription)")
    pyt = pd.read_csv(L.register(PYL3_TABLE, "stored PYL3 site table"), encoding="utf-8-sig")
    unp = L.read_text(UNPUB, "unpublished-data provenance register")
    val = L.read_text(VALUE, "private-data value assessment")
    con = L.read_text(CONSENT, "consent register")
    reg = L.read_json(REG, "re-search registration JSON (with amendments)")
    rpt = L.read_text(RPT, "re-search round report (Chinese)")

    phospho = {
        "PXD071110_sources": ph1a["search_configuration_read_from_the_deposit"]["source"],
        "PXD071110_note": ph1a["search_configuration_read_from_the_deposit"]["note"],
        "eight_deposit_round_config_source_mqpar_only": "mqpar.xml only" in ph8,
        "eight_deposit_round_parameters_txt_lacked_mods_and_tolerance":
            "confirmed to lack the modification list and tolerance" in ph8,
        "PXD071110_script_says_parameters_txt_states_neither":
            "parameters.txt states neither" in ph1 or "states neither the variable modifications" in ph1,
    }

    pyl_rows = pyt[pyt.protein == "PYL3"]
    pyl3 = {
        "engine": pyl["original_run"]["engine"],
        "raw_file": pyl["original_run"]["raw_file"],
        "instrument": pyl["instrument"]["model"], "fragmentation": pyl["instrument"]["fragmentation"],
        "research_parameters": {k: pyl["search_parameters_verified_on_a100"][k] for k in
                                ("variable_mods", "peptide_mass_tolerance_ppm", "fragment_bin_tol",
                                 "fragment_bin_offset", "theoretical_fragment_ions", "isotope_error")},
        "vendor_parameters": {k: pyl["vendor_parameters_verified_from_pdf"][k] for k in
                              ("software", "peptide_mass_tolerance", "fragment_mass_tolerance", "dynamic_mods")},
        "pyl3_site_psms_recount": {"persulfide_31_97": int(pyl_rows.recount_persulfide_31_97.sum()),
                                   "dioxidation_31_99": int(pyl_rows.recount_dioxidation_31_99.sum()),
                                   "alkyl_57_02": int(pyl_rows.recount_alkyl_57_02.sum()),
                                   "free": int(pyl_rows.recount_free.sum())},
        "gst_site_psms_recount": {"persulfide_31_97": int(pyt[pyt.protein == "GST tag"].recount_persulfide_31_97.sum()),
                                  "dioxidation_31_99": int(pyt[pyt.protein == "GST tag"].recount_dioxidation_31_99.sum())},
        "mass_accuracy": pyl["mass_accuracy"],
        "same_spectrum_comparisons_xcorr_ppm": next(c["effect"] for c in pyl["corrections"] if c["id"] == "C2"),
        "limits": pyl["limits"],
        "provenance_register_summary": next((ln for ln in unp.splitlines() if ln.startswith("U-03")), ""),
        "value_assessment": next((ln for ln in val.splitlines() if "PYL3" in ln), ""),
        "consent_register": next((ln for ln in con.splitlines() if ln.startswith("INV-19")), ""),
    }

    registration = {
        "registered_utc_label": reg["registered_utc"],
        "registered_before_any_search_was_run": reg["registered_before_any_search_was_run"],
        "fixed_before_any_search": ["question", "primary_statistic_and_reading_rule", "predeclared_outcomes",
                                    "search_parameters_fixed_before_the_run"],
        "twenty_psm_branch": reg["predeclared_outcomes"]["4_too_few_psms"],
        "amendments": [{"utc_label": a["utc"], "what": a["what"]} for a in reg["amendments"]],
        "static_cam_in_registration": any("add_C" in str(v) for v in reg["search_parameters_fixed_before_the_run"].values()),
        "run5_in_round_report": bool(re.search(r"\|\s*5（靶向臂）", rpt)),
        "run5_in_registration": any("run 5" in a["why"] for a in reg["amendments"]),
        "wide_window_diagnostic_in_registration": any("wide-window" in a["why"] for a in reg["amendments"]),
    }
    out = {"label": L.POSTHOC_LABEL, "phospho_config_sources": phospho, "pyl3_collaborator_research": pyl3,
           "registration": registration}
    L.write_json(out, "records_consulted_summary.json")
    return out


if __name__ == "__main__":
    import json
    print(json.dumps(main(), ensure_ascii=False, indent=1)[:6000])
