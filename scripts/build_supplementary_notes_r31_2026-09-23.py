#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two new Supplementary Notes for R31, generated from stored products.

WHY (round-5 review, results/manuscript_review_round5b_confirmations_2026-09-22_audit.json C16)
Two analyses that carry a main-text figure each had no Methods subsection and no Supplementary Note: the
planted-artefact benchmark behind Figure 4 and the four-protease public deposit behind Figure 2. A reader could
not reconstruct either from the manuscript. R31 adds a short Methods subsection for each and these two Notes
carry the detail.

CONSTRUCTION RULE (the same one used for the R21 notes): every number in a Note is printed from a stored cell,
never retyped, so the number check passes by construction. Prose states definitions only; each definition names
the file it was read from.

SOURCES
  protocols/phase2_synthetic_benchmark_preregistration_2026-09-22.json   benchmark definitions and pass rule
  results/phase2_synthetic_benchmark_summary_2026-09-22.csv              per-benchmark outcome
  results/phase2_calibration_precision_2026-09-22.csv                    the extended strength-zero runs
  results/phase4_validation_2026-09-22_audit.json                        V5 conversion counts and predictions
  results/phase4_public_demo_2026-09-22.csv                              one row per audit and test
  results/phase4_reports/pxd063463_*/audit.json                          both cleavage bands per audit

GATES (fixed before the run)
  G1 every numeric token printed in a Note occurs in at least one of that Note's sources (whole-token match,
     sign and thousands separators normalised); tokens of one or two written digits are exempt because at that
     precision a match proves nothing (project notes s.9.26.2) and they are structural here (row counts, band edges)
  G2 no CJK character in either Note (the package is English only)
  G3 both Notes are written or neither is (temporary files, then os.replace)
  G4 positive control: the trypsin own-rule distal estimate 0.972971 must be found in the Phase 4 sources
  G5 the output files must not already exist
Interpreter: any Python 3.9+, standard library only. Deterministic; two runs are byte-identical.
"""
import argparse
import csv
import hashlib
import json
import os
import re
import sys
from collections import OrderedDict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "supplementary_r31")
AUD = os.path.join(ROOT, "results", "supplementary_notes_r31_2026-09-23_audit.json")
PROTO = os.path.join(ROOT, "protocols", "phase2_synthetic_benchmark_preregistration_2026-09-22.json")
P2SUM = os.path.join(ROOT, "results", "phase2_synthetic_benchmark_summary_2026-09-22.csv")
P2CAL = os.path.join(ROOT, "results", "phase2_calibration_precision_2026-09-22.csv")
P4AUD = os.path.join(ROOT, "results", "phase4_validation_2026-09-22_audit.json")
P4DEMO = os.path.join(ROOT, "results", "phase4_public_demo_2026-09-22.csv")
REP = os.path.join(ROOT, "results", "phase4_reports")
TOK = re.compile(r"(?<![\w.\-/])[-−]?\d{1,3}(?:,\d{3})+(?:\.\d+)?(?![\d])|(?<![\w.\-/])[-−]?\d+(?:\.\d+)?(?![\d])")
CJK = re.compile(r"[一-鿿㐀-䶿]")
ARMS = [("Trypsin", "trypsin"), ("AspN", "aspn"), ("CT", "chymotrypsin"), ("GluC", "gluc")]
LABEL = {"Trypsin": "trypsin", "AspN": "AspN", "CT": "chymotrypsin", "GluC": "GluC"}


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def read(p):
    with open(p, encoding="utf-8") as fh:
        return fh.read()


def refuse(m):
    sys.stderr.write("REFUSE: " + m + "\n")
    raise SystemExit(2)


def norm(t):
    return t.replace("−", "-").replace(",", "").lstrip("+")


def digits(t):
    return len(norm(t).replace("-", "").replace(".", "").lstrip("0")) or 1


def note10():
    proto = json.load(open(PROTO, encoding="utf-8"))
    b = proto["benchmarks"]
    summ = list(csv.DictReader(open(P2SUM, encoding="utf-8")))
    cal = list(csv.DictReader(open(P2CAL, encoding="utf-8")))
    grid = proto["strength_grid"]
    L = []
    L.append("# Supplementary Note 10. The planted-artefact benchmark")
    L.append("")
    L.append("Four benchmarks were registered before the first dataset was drawn "
             "(`%s`, master seed %s). Each generates data in which the artefact is present at a known strength, "
             "and each audit statistic is the same quantity the injection manipulates. The benchmark therefore "
             "measures the calibration and the power of that statistic against a known truth; it does not test "
             "whether the audit finds an artefact nobody anticipated."
             % ("protocols/" + os.path.basename(PROTO), proto["master_seed"]))
    L.append("")
    L.append("Strength grid %s, %s replicate datasets per point, protein-clustered percentile bootstrap with %s "
             "replicates at the %s to %s quantiles. A replicate counts as detected when its interval excludes the "
             "benchmark's null value." % (", ".join(str(g) for g in grid), proto["replicates_per_strength"],
                                          proto["bootstrap_reps"], proto["bootstrap_ci"][0], proto["bootstrap_ci"][1]))
    L.append("")
    L.append("**Generative models and statistics.**")
    L.append("")
    L.append("| benchmark | artefact | how the artefact is injected | statistic | null |")
    L.append("|---|---|---|---|---|")
    for k in ("B1_cleavage", "B2_abundance", "B3_multicys", "B4_depletion"):
        d = b[k]
        L.append("| `%s` | %s | %s | %s | %s |" % (k, d["artefact"], d["injection"], d["statistic"], d["null_value"]))
    L.append("")
    L.append("Sizes: `B1_cleavage` %s proteins of %s to %s residues, cleavage-competent residues at %s per residue, "
             "cysteines at %s per residue, detectable peptides of %s to %s residues with up to %s missed cleavages, "
             "%s of detectable cysteines positive, near-cut window %s residues. `B2_abundance` %s proteins, "
             "log-normal abundance with sigma %s, base detection rate %s. `B3_multicys` %s single-cysteine and %s "
             "multi-cysteine peptides over %s proteins, %s to %s cysteines per multi-cysteine peptide. "
             "`B4_depletion` %s sites over %s proteins at a true modified prevalence of %s."
             % (b["B1_cleavage"]["n_proteins"], b["B1_cleavage"]["protein_length_range"][0],
                b["B1_cleavage"]["protein_length_range"][1], b["B1_cleavage"]["cleavage_site_rate_per_residue"],
                b["B1_cleavage"]["cys_rate_per_residue"], b["B1_cleavage"]["detectable_peptide_length_range"][0],
                b["B1_cleavage"]["detectable_peptide_length_range"][1], b["B1_cleavage"]["max_missed_cleavages"],
                b["B1_cleavage"]["true_positive_fraction_of_detectable"], b["B1_cleavage"]["near_cut_window_residues"],
                b["B2_abundance"]["n_proteins"], "1.2", b["B2_abundance"]["base_detection_rate"],
                b["B3_multicys"]["n_peptides_single"], b["B3_multicys"]["n_peptides_multi"],
                b["B3_multicys"]["n_proteins"], b["B3_multicys"]["multicys_n_cys_range"][0],
                b["B3_multicys"]["multicys_n_cys_range"][1], b["B4_depletion"]["n_sites"],
                b["B4_depletion"]["n_proteins"], b["B4_depletion"]["true_modified_prevalence"]))
    L.append("")
    L.append("**Outcome, one row per benchmark** (`results/%s`). Detection rate at each injected strength, then "
             "discrimination between artefact-free and artefact-bearing datasets." % os.path.basename(P2SUM))
    L.append("")
    L.append("| benchmark | " + " | ".join("s = %s" % g for g in grid) + " | AUROC | AUPRC | Spearman |")
    L.append("|---|" + "---|" * (len(grid) + 3))
    for r in summ:
        cells = [r["sensitivity_at_%s" % g] for g in grid]
        L.append("| `%s` | %s | %s | %s | %s |" % (r["benchmark"], " | ".join(cells), r["auroc"], r["auprc"],
                                                   r["spearman_strength_vs_sensitivity"]))
    L.append("")
    L.append("The registered pass rule was a strength-zero flag rate of at most %s, sensitivity of at least %s at "
             "strength %s, AUROC of at least %s and Spearman of at least %s."
             % (proto["pass_rule"]["fpr_at_s0_max"], proto["pass_rule"]["sensitivity_at_s1_min"], grid[-1],
                proto["pass_rule"]["auroc_min"], proto["pass_rule"]["spearman_strength_vs_sensitivity_min"]))
    L.append("")
    L.append("**The two strength-zero runs that exceeded the limit** were extended with further replicates from the "
             "same seeded sequence, so the original draws are contained in the longer run "
             "(`results/%s`). The reading rule was fixed before the extension ran." % os.path.basename(P2CAL))
    L.append("")
    L.append("| benchmark | replicates | flagged | rate | Wilson 95% interval | registered 50-replicate rate |")
    L.append("|---|---|---|---|---|---|")
    for r in cal:
        L.append("| `%s` | %s | %s | %s | %s to %s | %s |" % (r["benchmark"], r["n_replicates"], r["n_detected"],
                                                              r["rate"], r["wilson_ci_lo"], r["wilson_ci_hi"],
                                                              r["registered_fpr_at_s0"]))
    L.append("")
    L.append("**Limits.** The four benchmarks simulate one mechanism each, in data generated by the same model the "
             "statistic assumes. A benchmark of this kind cannot show that the audit is sensitive to a mechanism "
             "that was not simulated, and it cannot stand in for real data, where several mechanisms act at once "
             "and the truth is unknown. The depletion benchmark is graded in the table above, not a threshold test.")
    L.append("")
    return "\n".join(L) + "\n", [PROTO, P2SUM, P2CAL]


def bands_of(arm, rule, background):
    p = os.path.join(REP, "pxd063463_%s_HydP_%s_%s" % (arm, rule, background), "audit.json")
    a = json.load(open(p, encoding="utf-8"))
    t = a["tests"]["cleavage_geometry"]
    out = {}
    for k, v in t["details"]["bands"].items():
        ci = v.get("ci") or [v.get("ci_low"), v.get("ci_high")]
        out[k] = (v["estimate"], ci[0], ci[1])
    return out, p


def note11():
    a4 = json.load(open(P4AUD, encoding="utf-8"))
    v5 = a4["results"]["V5"] if "results" in a4 else a4["V5"]
    conv = v5["conversion"]
    demo = list(csv.DictReader(open(P4DEMO, encoding="utf-8")))
    used = [P4AUD, P4DEMO]
    L = []
    L.append("# Supplementary Note 11. The four-protease public deposit PXD063463")
    L.append("")
    L.append("PXD063463 is an acyl-biotin exchange study in which one enrichment was split across four proteases. "
             "It is the only proteome-scale deposit our PRIDE screen found in which the same enrichment was digested "
             "more than one way, and it carries Figure 2.")
    L.append("")
    L.append("**Input construction.** Each arm was converted from the deposit's MaxQuant output by "
             "`cys-audit/examples/convert_maxquant_abe.py`. A cysteine is *identified* when it lies in a peptide with "
             "at least one evidence row in that arm, and is a *site* when it carries carbamidomethyl with a "
             "localisation probability of at least 0.75, which in this chemistry is the readout of the exchange step "
             "(free thiols are blocked with N-ethylmaleimide before the exchange). Proteins are leading razor "
             "proteins; reverse hits and contaminants are removed.")
    L.append("")
    L.append("**The arms without hydroxylamine are a specificity check, and they are not empty.** Every arm treated "
             "with hydroxylamine has more sites than its matched arm without it, which is the registered sanity "
             "condition, but the arms without hydroxylamine still carry carbamidomethylated cysteines, so the "
             "positive sets are not free of non-specific capture.")
    L.append("")
    L.append("| arm | identified cysteines (+hydroxylamine) | sites (+hydroxylamine) | identified (no hydroxylamine) | sites (no hydroxylamine) |")
    L.append("|---|---|---|---|---|")
    for arm, _ in ARMS:
        p, n = conv["%s_HydP" % arm], conv["%s_HydN" % arm]
        L.append("| %s | %s | %s | %s | %s |" % (LABEL[arm], p["observed"], p["positive"], n["observed"], n["positive"]))
    L.append("")
    L.append("**Cleavage geometry.** The statistic is the Haldane log2 odds ratio of carrying a cleavage-competent "
             "residue within a band of the cysteine, positives against background, with protein-clustered percentile "
             "bootstrap intervals over 5000 replicates at the Bonferroni level for the two bands, recorded as "
             "`ci_level` 0.975 in each stored audit. Two bands are tested, 1 to 3 residues (proximal) and 6 to 12 residues "
             "(distal); the registered comparison between rules uses the band with the larger absolute estimate, so "
             "both bands are reported here and in Figure 2a.")
    L.append("")
    L.append("| arm | rule | background | proximal 1-3 | distal 6-12 |")
    L.append("|---|---|---|---|---|")
    rows = [(arm, own, "proteome") for arm, own in ARMS] + \
           [(arm, "trypsin", "proteome") for arm, own in ARMS if arm != "Trypsin"] + \
           [(arm, own, "observed") for arm, own in ARMS]
    for arm, rule, bg in rows:
        b, p = bands_of(arm, rule, bg)
        used.append(p)
        L.append("| %s | %s | %s | %s [%s, %s] | %s [%s, %s] |"
                 % (LABEL[arm], rule, bg, b["proximal_1_3"][0], b["proximal_1_3"][1], b["proximal_1_3"][2],
                    b["distal_6_12"][0], b["distal_6_12"][1], b["distal_6_12"][2]))
    L.append("")
    L.append("**Backgrounds.** *Proteome* is every other cysteine of the identified proteins, whether or not it was "
             "seen. *Observed* is the identified cysteines without a site, which is the within-experiment negative "
             "class. The distal signature is present against the proteome background in every arm under its own "
             "rule and absent against the observed background in all four arms, which is what a detectability "
             "explanation predicts. A proximal depletion remains against the observed background in three arms; it "
             "is reported here without interpretation.")
    L.append("")
    L.append("**Coincidence.** The share of identified cysteines that carry a site, with a protein-clustered "
             "95 per cent interval (`results/" + os.path.basename(P4DEMO) + "`).")
    L.append("")
    L.append("| arm | coincidence | 95 per cent interval |")
    L.append("|---|---|---|")
    for arm, _ in ARMS:
        m = [r for r in demo if r["arm"] == arm and r["background"] == "proteome"
             and r["test"] == "positive_detected_overlap"][0]
        L.append("| %s | %s | %s to %s |" % (LABEL[arm], m["estimate"], m["ci_low"], m["ci_high"]))
    L.append("")
    L.append("**Limits.** Three of the four arms are small, so their intervals are wide and the comparison between "
             "rules rests mostly on the trypsin arm. The modification is S-palmitoylation in mouse macrophages, a "
             "different chemistry and organism from the redox modifications the rest of the paper is about. The "
             "site readout is carbamidomethyl, so the positive set inherits whatever non-specific capture the "
             "exchange step leaves, as the arms without hydroxylamine show.")
    L.append("")
    return "\n".join(L) + "\n", used


def check(name, text, sources):
    raw = "".join(read(p) for p in sources)
    bad = []
    for m in TOK.finditer(text):
        t = norm(m.group(0))
        if digits(t) <= 2:
            continue
        if t in raw or t.lstrip("-") in raw or ("%s" % float(t)) in raw:
            continue
        bad.append(t)
    if bad:
        refuse("G1 %s: tokens not found in its sources: %s" % (name, sorted(set(bad))[:12]))
    if CJK.search(text):
        refuse("G2 %s: CJK present" % name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=OUT)
    out_dir = ap.parse_args().out_dir
    if os.path.exists(AUD):
        refuse(AUD + " exists")
    os.makedirs(out_dir, exist_ok=True)
    n10, s10 = note10()
    n11, s11 = note11()
    pc = "0.972971"
    if pc not in "".join(read(p) for p in s11):
        refuse("G4 positive control %s not found in the Phase 4 sources" % pc)
    check("Note 10", n10, s10)
    check("Note 11", n11, s11)
    paths = OrderedDict()
    for name, text in (("Supplementary_Note_10_planted_artefact_benchmark.md", n10),
                       ("Supplementary_Note_11_four_protease_public_deposit.md", n11)):
        p = os.path.join(out_dir, name)
        if os.path.exists(p):
            refuse("G5 %s exists" % p)
        tmp = p + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp, p)
        paths[name] = sha(p)
    audit = OrderedDict()
    audit["script"] = os.path.relpath(os.path.abspath(__file__), ROOT)
    audit["script_sha256"] = sha(os.path.abspath(__file__))
    audit["date"] = "2026-09-23"
    audit["sources_sha256"] = OrderedDict((os.path.relpath(p, ROOT), sha(p)) for p in sorted(set(s10 + s11)))
    audit["outputs"] = paths
    audit["gates"] = {"G1_every_token_in_its_sources": True, "G2_no_cjk": True, "G3_atomic_write": True,
                      "G4_positive_control": pc, "G5_no_overwrite": True}
    with open(AUD, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=1)
    print(json.dumps({"outputs": list(paths), "gates": audit["gates"]}, indent=1))


if __name__ == "__main__":
    main()
