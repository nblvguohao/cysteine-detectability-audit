"""VERIFIER round 4 (adversarial; POST HOC), item F_rice_artifact3: comparison table (implementer vs verifier)
and a provenance record for the verify_r4 outputs. Reads only verify_r4/*.json written by the verifier's own
scripts; writes verify_r4/comparison_r4.csv and verify_r4/provenance_verify_r4.json."""
import csv
import hashlib
import json
import pathlib
import platform
import sys

W = pathlib.Path("C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30")
OUT = W / "results" / "F_rice_artifact3" / "verify_r4"
SCR = W / "scripts" / "F_rice_artifact3"
b = json.loads((OUT / "build_r4.json").read_text(encoding="utf-8"))
m = json.loads((OUT / "models_r4.json").read_text(encoding="utf-8"))
nr = json.loads((OUT / "null_recount_r4.json").read_text(encoding="utf-8"))["summary"]
ns = json.loads((OUT / "null_sim_r4.json").read_text(encoding="utf-8"))
rc = json.loads((OUT / "rice_r4.json").read_text(encoding="utf-8"))
c = m["comparisons"]; ln = m["linear_symmetric"]
r = lambda x, k=3: round(float(x), k)
rows = [
    ("group table rebuilt from raw inputs: mismatches (labels, length, n_cys, n_det, 4 windows)", "0", str(sum(v for k, v in b["comparison_with_item"].items() if k.startswith("mismatch"))), True),
    ("labels A / S / S1 / S_map; A not S", "1197 / 750 / 692 / 725; 447", f"{b['labels']['A']} / {b['labels']['S']} / {b['labels']['S1']} / {b['labels']['S_map']}; {b['A_not_S']}", True),
    ("iBAQ ratio A / S / A-not-S; raw A / S", "10.36 / 7.21 / 17.8; 4.88 / 3.93",
     f"{r(m['ratios']['A_ibaq'][0],2)} / {r(m['ratios']['S_ibaq'][0],2)} / {r(m['ratios']['A_not_S_ibaq'][0],2)}; {r(m['ratios']['A_qty'][0],2)} / {r(m['ratios']['S_qty'][0],2)}", True),
    ("rank-biserial A, all 7,692: length / Cys (unmatched, deciles, 1:1)", "-0.026/+0.199/+0.150; +0.082/+0.286/+0.244",
     "/".join(f"{r(v,3):+}" for v in m["rank_biserial_A_all_groups"]["length"].values()) + "; " + "/".join(f"{r(v,3):+}" for v in m["rank_biserial_A_all_groups"]["n_cys"].values()), True),
    ("Spearman log iBAQ vs length / Cys", "-0.28 / -0.24", f"{r(m['spearman']['logibaq_length'])} / {r(m['spearman']['logibaq_ncys'])}", True),
    ("dAIC(D-T) five forms, A", "10.6, 19.0, 17.0, 19.8, 12.5", ", ".join(str(r(c[f'A|{f}']['dAIC_D_minus_T'], 2)) for f in ("logit_all", "logit_h", "cloglog_h", "splq_h", "splc_h")), True),
    ("dAIC(D-T) five forms, S", "19.5, 22.9, 23.4, 24.1, 22.3", ", ".join(str(r(c[f'S|{f}']['dAIC_D_minus_T'], 2)) for f in ("logit_all", "logit_h", "cloglog_h", "splq_h", "splc_h")), True),
    ("P(detectable | total), A, four forms", "0.26, 0.95, 0.82, 0.48", ", ".join(str(r(c[f'A|{f}']['p_det_given_total'], 3)) for f in ("logit_all", "logit_h", "cloglog_h", "splq_h")), True),
    ("P(detectable | total), S, four forms", "0.005, 0.10, 0.11, 0.32", ", ".join(str(r(c[f'S|{f}']['p_det_given_total'], 4)) for f in ("logit_all", "logit_h", "cloglog_h", "splq_h")), True),
    ("round 3 indicator test P(det | total, 1[n_det>=1]) A / S / S1 / S_map", "0.11 / 0.149 / 0.38 / 0.022",
     " / ".join(str(r(c[f'{l}|ind_all']['p_det_given_total'], 4)) for l in ("A", "S", "S1", "S_map")), True),
    ("round 3 indicator test P(total | det, ind), A (max over labels)", "<= 7e-4", f"{c['A|ind_all']['p_total_given_det']:.2g}", True),
    ("round 3 indicator test dAIC(DI-TI) A / S / S1 / S_map", "9.0 / 27.4 / 23.6 / 14.4",
     " / ".join(str(r(c[f'{l}|ind_all']['dAIC_DI_minus_TI'], 2)) for l in ("A", "S", "S1", "S_map")), True),
    ("indicator test without length, A (not in texts)", "(not quoted; texts say P >= 0.11 in every form)", str(r(c["A|ind_all"]["p_det_given_total_nolen"], 4)), False),
    ("per-Cys OR det / undet (P), A and S", "1.077/1.071 (0.80); 1.111/1.081 (0.22)",
     f"{r(ln['A|all']['OR_det'][0])}/{r(ln['A|all']['OR_undet'][0])} ({r(ln['A|all']['p_equal'],2)}); {r(ln['S|all']['OR_det'][0])}/{r(ln['S|all']['OR_undet'][0])} ({r(ln['S|all']['p_equal'],2)})", True),
    ("length OR per doubling at fixed total count, A; S; S1 linear", "1.07 [0.97,1.17]; 0.92 [0.82,1.03]; 1.21 [1.08,1.35]",
     f"{[r(x,2) for x in c['A|logit_all']['len_OR_per_doubling_at_fixed_l2c'][:3]]}; {[r(x,2) for x in c['S|logit_all']['len_OR_per_doubling_at_fixed_l2c'][:3]]}; {[r(x,2) for x in c['S1|logit_all']['len_OR_per_doubling_at_fixed_linear_ncys']]}", True),
    ("pooled nulls recounted from the item's draws (outside of 14): A N2/N4/N1/N3", "8 / 1 (edge P 0.049) / 11 / 3",
     f"{nr['A|N2']['n_outside']} / {nr['A|N4']['n_outside']} (edge {nr['A|N4']['edge']}) / {nr['A|N1']['n_outside']} / {nr['A|N3']['n_outside']}", True),
    ("pooled nulls: A N4 min P of the other 13; A N2 max P outside", ">= 0.46; < 0.01", f"{r(nr['A|N4']['min_p_inside'],3)}; {r(nr['A|N2']['max_p_outside'],4)}", True),
    ("pooled nulls: S N4 outside (length P); cysteine min P", "2 (0.022, 0.013); > 0.09",
     f"{nr['S|N4']['n_outside']} ({nr['S|N4']['len_stats']['r_m_len'][1]}, {nr['S|N4']['len_stats']['r_mc_len'][1]}); {r(nr['S|N4']['cys_stats_min_p'],3)}", True),
    ("pooled nulls: S_map N4/N2/N1/N3; N3 with 1,000 sets", "4/9/9/4; 5",
     f"{nr['S_map|N4']['n_outside']}/{nr['S_map|N2']['n_outside']}/{nr['S_map|N1']['n_outside']}/{nr['S_map|N3']['n_outside']}; {nr['S_map|N3']['n_outside_first1000']}", True),
    ("independent re-simulation A|N4 (10,000 sets, converged fit, own code/seeds): outside; edge-stat P; min P other 13",
     "1 of 14 at edge (P 0.049); >= 0.46",
     f"{ns['summary']['A|N4']['n_outside']} of 14; P {r(ns['summary']['A|N4']['per_stat']['r_dt_cys']['p2'],4)}; {r(min(v['p2'] for k, v in ns['summary']['A|N4']['per_stat'].items() if k != 'r_dt_cys'),3)}", True),
    ("independent re-simulation A|N2 (5,000 sets): outside; max P outside", "8 of 14; < 0.01",
     f"{ns['summary']['A|N2']['n_outside']} of 14; {r(max(v['p2'] for v in ns['summary']['A|N2']['per_stat'].values() if not v['inside']),4)}", True),
    ("independent re-simulation S|N4 (5,000 sets): outside (length); min P of the 9 cysteine statistics", "2; > 0.09",
     f"{ns['summary']['S|N4']['n_outside']}; {r(min(ns['summary']['S|N4']['per_stat'][s]['p2'] for s in ['r_dec_cys','r_dt_cys','r_m_cys','r_m_undet','b_undet_M3u','r_mc_det','dAIC_D_minus_T','b_det_B','d_det_minus_undet_L']),3)}", False),
    ("rice chain step 0 / 0b / 1 (sites/proteins)", "1796/1192; 1769/1170; 1753/1164",
     "; ".join(f"{rc['chain'][k]['sites']}/{rc['chain'][k]['proteins']}" for k in ("0_all_incl_decoys_contaminants", "0b_decoys_removed", "1_decoys_contaminants_removed")), True),
    ("rice chain steps 2-5 and equality with stored list", "1377; 1340/1040; 817/640/5598; identical",
     f"{rc['chain']['2_single_cys']['sites']}; {rc['chain']['4_ge2_cys']['sites']}/{rc['chain']['4_ge2_cys']['proteins']}; {rc['chain']['5_active_2026_03']['sites']}/{rc['chain']['5_active_2026_03']['proteins']}/{rc['chain']['5_active_2026_03']['cys_in_proteins']}; {rc['final_equals_stored']}", True),
    ("lost proteins / sites / Cys", "400 / 523 / 4112", f"{rc['lost']['proteins']} / {rc['lost']['sites']} / {rc['lost']['cys']}", True),
    ("'0/400 saved entry records carry a sequence' as evidence that UniProtKB returns no sequence for deleted entries",
     "0/400 (cited as support)", f"0/400 but query was '{rc['deleted_entry_sequence_evidence']['s01b_query_string']}' (uninformative); the s01 batch query with a sequence field returned {rc['deleted_entry_sequence_evidence']['lost_returned_by_s01_batch_query_with_sequence_field']}/400 (valid support)", False),
    ("self-audit exclusion -> kept", "594/80/59 -> 5004/737/581",
     f"{rc['self_audit_gate']['rows_in_table']-rc['self_audit_gate']['rows_kept']}/817-{rc['self_audit_gate']['positives']}/{rc['self_audit_gate']['proteins_in_table']-rc['self_audit_gate']['proteins_matched']} -> {rc['self_audit_gate']['rows_kept']}/{rc['self_audit_gate']['positives']}/{rc['self_audit_gate']['proteins_matched']}", True),
    ("kept vs lost: median length; median Cys; sites/protein; quantified", "392 vs 411; 7 vs 8; 1.28 vs 1.31; 75% vs 52%",
     f"{rc['kept_vs_lost']['median_length']}; {rc['kept_vs_lost']['median_cys']}; {[r(x,2) for x in rc['kept_vs_lost']['sites_per_protein']]}; {[r(x,2) for x in rc['kept_vs_lost']['frac_in_quantified_groups_7692']]}", True),
    ("peptides unique to one group, every C", "1692", str(rc["unique_groups_every_C"]), True),
    ("comparison arm 'absent from that table'", "6,495 groups absent from the peptide table",
     f"{b['negatives_with_an_accession_in_Proteins_column']} of {b['negatives']} have an accession in the table's Proteins column", False),
]
with open(OUT / "comparison_r4.csv", "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh); w.writerow(["name", "implementer_value", "verifier_value", "match"])
    for row in rows:
        w.writerow(row)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
prov = {"label": "VERIFIER round 4 (adversarial; POST HOC) for F_rice_artifact3, 2026-09-30",
        "scripts_sha256": {p.name: sha(p) for p in sorted(SCR.glob("verify_r4_*.py"))},
        "outputs_sha256": {p.name: sha(p) for p in sorted(OUT.glob("*")) if p.is_file() and p.name != "provenance_verify_r4.json"},
        "inputs": b["inputs_sha256"], "network": "none", "workers": "4 (verify_r4_null_sim.py), otherwise 1",
        "seeds": {"ratio bootstrap": 4040404, "null re-simulation": "SeedSequence(990000 + combo, spawn_key=(chunk,))"},
        "python": sys.version.split()[0], "platform": platform.platform(),
        "run_order": "verify_r4_build -> verify_r4_models -> verify_r4_null_recount -> verify_r4_null_sim (10000 5000 5000) -> verify_r4_rice -> verify_r4_compare"}
(OUT / "provenance_verify_r4.json").write_text(json.dumps(prov, indent=1), encoding="utf-8")
for row in rows:
    print(" | ".join(str(x) for x in row))
