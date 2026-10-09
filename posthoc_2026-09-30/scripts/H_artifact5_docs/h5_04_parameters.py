"""H_artifact5_docs step 4 - POST HOC extraction and cross-check of every PXD015307 re-search parameter.

POST HOC revision analysis in response to a pre-submission review; not registered, not pre-specified.

Every parameter an MCP reviewer expects for a database search is listed with (i) the value to report,
(ii) the record(s) it comes from, (iii) whether independent records agree, and (iv) whether it could be
checked against a stored machine-readable output or against the search output itself (inference from
the stored Comet output lines, step 1). Parameters that only a comet.params file would settle are
flagged: the two parameter files named in the round report sit in external/sulfhydrome_bta_tmt/research/
on the analysis host and are not in the local repository copy.

Revision after verification (round 2): contaminants are marked '[to confirm]' (no record states them); the
static carbamidomethyl term is recorded as a configuration error of run 1 rather than a plan amendment; a
row states what was registered before any search; the fifth discarded run (insulin files without insulin,
round report only) is added.
"""
from __future__ import annotations

import json
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5_lib as L  # noqa: E402

REG = os.path.join(L.W, "inputs", "repo_reports", "pxd015307_research_preregistration_2026-09-17.json")
RPT = os.path.join(L.W, "inputs", "repo_reports", "PXD015307_RESEARCH_SEARCH_SPACE.md")
AUD = os.path.join(L.W, "inputs", "repo_results", "pxd015307_research_audit.json")
SUMCSV = os.path.join(L.W, "inputs", "repo_results", "pxd015307_research_summary.csv")
PSMCSV = os.path.join(L.W, "inputs", "repo_results", "pxd015307_research_psms.csv")
SCR = os.path.join(L.REPO, "scripts", "run_pxd015307_research_mass_accuracy.py")
PHS = os.path.join(L.REPO, "scripts", "posthoc_pxd015307_score_ties.py")
TIES = os.path.join(L.W, "inputs", "repo_results", "pxd015307_posthoc_score_ties.csv")
MS = os.path.join(L.MCP, "01_manuscript.tex")
TIES_SUMMARY = os.path.join(L.OUT, "ties_summary.json")
DEP_SUMMARY = os.path.join(L.OUT, "deposit_plus32_summary.json")


def main():
    reg = L.read_json(REG, "re-search registration JSON (with amendments)")
    rpt = L.read_text(RPT, "re-search round report (Chinese)")
    aud = L.read_json(AUD, "stored re-search audit JSON")
    sm = pd.read_csv(L.register(SUMCSV, "stored re-search per-file summary"))
    psm_txt = L.read_text(PSMCSV, "stored +32 PSM table (header only)")
    scr = L.read_text(SCR, "repository reading script (docstring)")
    phs = L.read_text(PHS, "repository post hoc tie script (docstring)")
    ms = L.read_text(MS, "manuscript (read-only)")
    ts = json.load(open(TIES_SUMMARY, encoding="utf-8"))
    ds = json.load(open(DEP_SUMMARY, encoding="utf-8"))
    ties = pd.read_csv(TIES, dtype={"e_value_best": str})

    amend = " ".join(a["what"] + " " + a["why"] for a in reg["amendments"])
    ep = ms[ms.find(r"\subsection*{Search space and mass accuracy}"):ms.find(r"\subsection*{Software}")]
    pa = aud["per_arm"]
    vlines = {k: v["comet_version_line"] for k, v in pa.items()}
    comet_ok = all(v.startswith("CometVersion 2026.01 rev. 1 (e4f767c)") for v in vlines.values())
    db_fig1d = "db/mmu_UP000000589_sp_iso.fasta" in vlines["Fig-1D"]
    db_targ = all("db/mmu_plus_bovine_insulin.fasta" in v for k, v in vlines.items() if k != "Fig-1D")
    run_times = {k: re.search(r"(\d\d/\d\d/\d{4}, [\d:]+ [AP]M)", v).group(1) for k, v in vlines.items()}
    big = ties[(ties.ppm_min < -10) | (ties.ppm_max > 10)]

    def in_ep(*words):
        return "yes" if all(w in ep for w in words) else "no"

    rows = []

    def row(section, parameter, value, applies, source, corroboration, status, in_manuscript, note=""):
        rows.append({"section": section, "parameter": parameter, "value_to_report": value,
                     "applies_to": applies, "primary_source": source, "corroborating_source": corroboration,
                     "verification_status": status, "in_current_manuscript_EP": in_manuscript, "note": note})

    files = reg["materials"]["files_and_expected_bytes"]
    bytes_in_report = all(f"{b:,}" in rpt for b in files.values())
    row("data", "deposit", "PRIDE PXD015307 ('Quantitative profiling of protein sulfhydrome by the BTA-TMT'); original study Gao et al. 2020 [52]",
        "all", "registration materials", "audit 'accession'", "verified_consistent", in_ep("PXD015307") if "PXD015307" in ep else "no (Results only)")
    row("data", "raw files and sizes", "; ".join(f"{k} {v:,} bytes" for k, v in files.items()),
        "all", "registration materials.files_and_expected_bytes", "round report table",
        "verified_consistent" if bytes_in_report else "conflicting_records", "no",
        "report states the sizes match the PRIDE file list")
    row("data", "raw file sha256", "NOT RECORDED", "all", "registration: 'sha256 recorded in the round's audit JSON'",
        "audit JSON has no raw-file hash", "missing", "no", "compute from the PRIDE files and add to Supplemental Data 12")
    row("data", "download source", reg["materials"]["source"], "all", "registration", "-", "prose_only", "no")
    row("data", "MS2 analyser per file",
        "Fig-1D: 31,450 FTMS + 3 ITMS MS2 scans (high-resolution MS2); Fig4-L_NaHS: 1,080 FTMS + 5,004 ITMS (ion-trap MS2); other three Fig4-L files ion trap per report (~1,100 FTMS / ~5,000 ITMS, not counted in a stored output)",
        "per file", "registration amendment (2026-09-17T02:00)", "round report", "partially_recorded", "partly ('mass analyzer actually used')",
        "survey scans FTMS in all files")
    row("data", "instrument model", "NOT RECORDED", "all", "-", "-", "missing", "no", "read from raw-file headers or PRIDE metadata")
    row("data", "raw-to-mzML conversion", "ThermoRawFileParser (bioconda, micromamba env ms2026); version and peak-picking options NOT RECORDED",
        "all", "registration search_parameters + amendment 1", "-", "missing_version", "no")
    row("engine", "search engine and version", "Comet 2026.01 rev. 1 (e4f767c)", "all five files",
        "audit comet_version_line (all five)", "report; manuscript EP", "verified_stored_output" if comet_ok else "conflicting_records",
        in_ep("2026.01"))
    row("engine", "run date/time (host local time) and host", "; ".join(f"{k} {v}" for k, v in run_times.items()) + "; host amax, 24 threads, CPU only",
        "per file", "audit comet_version_line", "registration amendment 1", "verified_stored_output", "no",
        "audit generated_utc 2026-09-16T18:20:20Z implies host local time = UTC+8")
    row("engine", "comet.params files", "two files, external/sulfhydrome_bta_tmt/research/ (analysis host)", "all",
        "round report section 8", "-", "not_in_local_copy", "no",
        "fragment_bin_offset for Fig-1D, theoretical_fragment_ions for Fig-1D, ion series, charge range, peptide length/mass range, num_enzyme_termini, spectrum preprocessing cannot be verified without them")
    row("database", "target database, Fig-1D", "mmu_UP000000589_sp_iso.fasta: UniProtKB/Swiss-Prot Mus musculus reference proteome UP000000589, canonical + isoforms, 25,750 entries",
        "Fig-1D", "audit comet_version_line (file name)", "report; registration amendment (2026-09-17T01:50)",
        "file_name_verified_entries_prose_only" if db_fig1d else "conflicting_records", "no",
        "UniProt release, download date and FASTA sha256 NOT RECORDED")
    row("database", "target database, Fig4-L files", "mmu_plus_bovine_insulin.fasta: the mouse database plus bovine insulin P01317 (25,751 entries)",
        "four Fig4-L files", "audit comet_version_line (file name)", "report",
        "file_name_verified_entries_prose_only" if db_targ else "conflicting_records", "no",
        "sha256 NOT RECORDED. CONFLICT: the deposit's own search of these files assigned HUMAN insulin chains (" +
        ", ".join(ds["insulin_sequences_assigned_by_deposit"]) + "); human insulin (P01308) is absent from the re-search database")
    row("database", "species rationale", "deposit's own Fig-1D.msf: 19,303 protein annotations, all _MOUSE",
        "Fig-1D", "registration amendment (2026-09-17T01:50)", "report", "prose_only", "no")
    row("database", "database originally registered", "Bos taurus UP000009136 (canonical); replaced by mouse after run 2 (174 PSMs at 1% FDR from 49,522 spectra)",
        "all", "registration materials.database", "registration amendment; report section 5",
        "amended_after_failed_run", "no",
        "the reading script's docstring still describes the bovine database (UP000009136, 51,611 entries, 'only run 2 is read here') and the audit's search_ran_on says 'run 2': stale records" if ("Bos taurus UP000009136" in scr and "51,611" in scr) else "")
    row("database", "contaminant sequences", "[to confirm]: no record states whether contaminant sequences were appended", "all", "-",
        "-",
        "missing", "no", "confirm from the FASTA or comet.params on the analysis host before stating 'none' in the manuscript")
    ut = ts["unique_peptide_td_counts"]
    row("database", "decoys", "Comet internal decoys (reversed peptide, C-terminal residue retained), concatenated search; decoy prefix DECOY_",
        "all", "script docstring ('internal concatenated decoys'); audit decoy_prefix", "output inference: "
        f"{ut.get('decoy', 0)} of {sum(ut.values())} unique output peptides map to the mouse proteome only after reversing all but the C-terminal residue",
        "verified_output_inference", "partly ('target--decoy FDR' only)", "")
    ez = ts["enzyme_inference"]
    row("enzyme", "enzyme and specificity", "trypsin (cleavage C-terminal to K/R, not before P), fully tryptic",
        "all", "registration ('trypsin, up to 2 missed cleavages'); script docstring", "output inference: "
        f"{ez['fraction_with_a_fully_tryptic_occurrence']:.4f} of {ez['n_target_peptides_checked']} target peptides have a fully tryptic occurrence; "
        f"{ez['n_with_internal_KP_or_RP']} contain internal KP/RP and {ez['n_exceeding_2_missed_cleavages_if_trypsin_P']} would exceed 2 missed cleavages under Trypsin/P",
        "verified_output_inference", "no", "num_enzyme_termini not recorded; inferred = 2")
    row("enzyme", "missed cleavages", "up to 2", "all", "registration; script docstring",
        f"output inference: max {ez['max_missed_cleavages']} (distribution {ez['missed_cleavage_distribution']})",
        "verified_output_inference", "no")
    row("enzyme", "peptide length range", f"observed {ez['length_min_all_output_peptides']}-{ez['length_max_all_output_peptides']} residues (consistent with Comet's default 5-50); setting NOT RECORDED",
        "all", "output inference", "-", "inferred_only", "no")
    row("tolerance", "precursor tolerance", "+/-10 ppm, all files (FTMS survey scans)", "all",
        "registration search_parameters; amendment (2026-09-17T02:00)", "report; script assertion |ppm| <= 10.0001 on 1%-FDR PSMs",
        "verified_consistent", "no",
        f"output ppm range {ts['ppm_range_all_rows']}: {len(big)} output lines exceed 10 ppm on neutral mass by at most 0.02 ppm, consistent with the tolerance being evaluated on another mass reference (e.g. MH+); tolerance type not recorded")
    row("tolerance", "isotope error", "0 (no 13C offsets)", "all", "registration amendment (2026-09-17T01:15)", "report; script docstring",
        "verified_consistent", "no")
    row("tolerance", "fragment settings, Fig-1D", "fragment_bin_tol 0.02 (high-resolution); fragment_bin_offset and theoretical_fragment_ions NOT RECORDED",
        "Fig-1D", "registration search_parameters (0.02); amendment (2026-09-17T02:00)", "report", "partially_recorded", "partly")
    row("tolerance", "fragment settings, Fig4-L files", "fragment_bin_tol 1.0005, fragment_bin_offset 0.4, theoretical_fragment_ions 1 (ion trap)",
        "four Fig4-L files", "registration amendment (2026-09-17T02:00)", "report",
        "prose_only", "partly", "script docstring still says 'fragment bin 0.02' for all files (stale)" if "fragment bin 0.02" in scr else "")
    row("tolerance", "ion series / neutral losses", "NOT RECORDED (Comet default b and y)", "all", "-", "-", "missing", "no")
    row("tolerance", "precursor charge range", "NOT RECORDED", "all", "-", "-", "missing", "no")
    row("modifications", "static modifications", "none (the registration listed no static modification; run 1 was discarded because Comet's default static add_C_cysteine 57.021464 had been left on, a configuration error rather than a plan amendment)",
        "all", "registration (variable modifications only); amendment 6 (discarded runs); report", "script docstring; manuscript EP", "verified_consistent", in_ep("static"))
    row("modifications", "variable modifications",
        "cysteine: sulfide +31.972071, dioxidation +31.989829, iodoTMT6plex +329.226595, carbamidomethyl +57.021464; methionine: oxidation +15.9949",
        "all", "report; script docstring", "registration (oxidation +15.994915; also listed protein N-terminal acetyl, later dropped)",
        "verified_consistent", "partly ('both candidates')",
        "the +32 masses are verified in the output (script tolerance 1e-4 Da); oxidation mass quoted as 15.9949 vs 15.994915")
    row("modifications", "maximum variable modifications per peptide", "3", "all", "registration amendment (2026-09-17T01:15)",
        "report; script docstring", "verified_consistent", "no")
    row("modifications", "protein N-terminal acetylation", "not searched (registered, then dropped)", "all",
        "registration amendment (2026-09-17T01:15)", "report limitation 4", "amended", "no",
        "stated reason ('five variable-modification slots') should be checked: Comet parameter files provide more than five variable_mod entries")
    row("modifications", "N-ethylmaleimide (+125.047679, C)", "NOT searched", "four Fig4-L files", "absence from every record of the search space",
        f"deposit's own search placed NEM on {ds['n_targeted_psms_with_nem']} of {ds['n_targeted_psms']} targeted-arm PSMs",
        "conflict_with_deposit", "no", "all four deposit Sulfide PSMs lie on human insulin chains absent from the re-search database, and three of them also carry NEM; none could be re-matched")
    row("modifications", "iodoTMT6plex on D/E/H/K and oxidation on residues other than M", "not searched",
        "Fig-1D", "search records", "deposit's Fig-1D.msf placed iodoTMT6plex on C, D, E, H, K and oxidation on 11 residue types",
        "difference_from_deposit", "no")
    row("output", "candidates written per spectrum", "5 (num_output_lines)", "all", "registration amendment (2026-09-17T01:15)", "-", "prose_only", "no")
    row("fdr", "score, level and unit", "Comet E-value; rank-1 PSMs; PSM level, 1%, each file separately; no peptide- or protein-level FDR",
        "all", "script (fdr_threshold per arm)", "audit fdr 0.01", "verified_script", "partly ('1% target-decoy FDR')")
    row("fdr", "FDR estimator", "decoys/targets among rank-1 PSMs at or below an E-value; threshold = largest E-value with ratio <= 0.01; a PSM is a decoy when all its proteins carry DECOY_; no +1 correction",
        "all", "script", "audit decoy_prefix", "verified_script", "no",
        "with 14-26 accepted targets per Fig4-L file the estimator cannot resolve 1% (0 decoys accepted)")
    row("fdr", "rescoring / post-processing", "none (no Percolator); precursor error = (exp - calc)/calc x 1e6 on neutral masses",
        "all", "registration declared limits; script", "report limitation 4", "verified_consistent", "no")
    for _, r in sm.iterrows():
        a = pa[r.arm]
        row("results", f"per-file counts: {r.arm}",
            f"{a['n_rank1_psms']:,} rank-1 PSMs; E-value threshold {a['e_value_threshold_at_1pct_fdr']}; {r.n_psms_at_1pct_fdr} PSMs at 1% FDR; "
            f"{r.n_plus32_cys_psms} +32-class Cys PSMs; baseline median {r.baseline_median_ppm} ppm, SD {r.baseline_sd_ppm} ppm",
            r.arm, "audit per_arm", "summary CSV", "verified_stored_output", "no",
            f"Comet output sha256 {a['input_sha256']}")
    row("results", "+32 PSM table", "empty (header only)", "all", "results/pxd015307_research_psms.csv", "audit", "verified_stored_output" if psm_txt.strip() == "arm" else "check", "no")
    row("reading", "resolvable rule", "|ppm - file baseline median| < separation_ppm / 2", "all", "audit resolvable_rule; script",
        "registration says it is the rule of scripts/audit_search_space_and_mass_accuracy.py, which instead uses |ppm| < separation_ppm / 2",
        "conflicting_records", "yes", "the manuscript EP describes the baseline-corrected rule")
    row("provenance", "registration timestamps", f"registered_utc {reg['registered_utc']}; amendments 2026-09-17T01:15-02:00Z",
        "all", "registration JSON", f"audit generated_utc {aud['generated_utc']} (gmtime); Comet runs at 02:13-02:19 host local time",
        "conflicting_records", "yes ('fixed before the search')",
        "taken as UTC, the registration postdates the search and its reading by about 6.3 h; consistent only if the times are UTC+8 local times mislabelled 'Z'")
    row("provenance", "discarded runs", "five: runs 1-4 (Comet default static carbamidomethyl left on; the registered bovine database; high-resolution fragment settings on the ion-trap files; ion-trap settings on the high-resolution file) and run 5 (insulin files searched without insulin in the database; bovine insulin was then appended); plus a wide-window diagnostic search (precursor -200/+500 Da, no variable modifications) terminated before completion",
        "all", "registration amendment 6 (runs 1-4); round report section 5 (runs 1-5 and the diagnostic search)",
        "run 5 is in the round report only, not among the registered amendments", "not_in_local_copy", "no",
        "the audit's own amendment says the reading script's first execution stopped on an arm without baseline PSMs; every final file has baseline PSMs, so that execution read an output other than the final one")
    row("provenance", "what was registered before any search", "question; reading rule (baseline-corrected half-separation rule); primary and secondary readouts; the 20-PSM minimum and four outcome branches; the search plan (bovine database, 0.02-Da fragment bins for all files, protein N-terminal acetylation)",
        "all", f"registration JSON (registered_utc {reg['registered_utc']}, registered_before_any_search_was_run = {reg['registered_before_any_search_was_run']})",
        "round report section 1", "verified_consistent", "partly ('Outcome branches were fixed before the search')",
        "the settings, not the question, rule or branches, were amended after the failed runs: database species (amendment 4), per-file fragment settings (amendment 5); compute host, N-terminal acetylation, isotope error/maximum modifications and output lines (amendments 0-3) are dated 01:15, before the first recorded failure")
    df = pd.DataFrame(rows)
    L.write_csv(df, "search_parameters_verified.csv")
    counts = df.verification_status.value_counts().to_dict()
    L.write_json({"label": L.POSTHOC_LABEL, "n_parameters": int(len(df)), "status_counts": counts,
                  "comet_version_verified": comet_ok, "db_names_verified": [db_fig1d, db_targ],
                  "output_lines_beyond_10ppm": int(len(big)),
                  "stale_docstring_bovine": "Bos taurus UP000009136" in scr,
                  "stale_search_ran_on": "run 2" in aud["search_ran_on"],
                  "posthoc_script_mentions_human_vs_bovine": "HUMAN insulin" in phs},
                 "search_parameters_summary.json")
    return df


if __name__ == "__main__":
    d = main()
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 120)
    print(d[["parameter", "verification_status", "in_current_manuscript_EP"]].to_string())
