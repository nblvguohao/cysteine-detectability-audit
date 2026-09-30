"""Step 7 (POST HOC): collect every number quoted in reports/F_rice_artifact3.md into one table,
read from the stored outputs of s02-s05, s08-s11 and (round 3) s12-s13 (nothing is recomputed here).
Output: results/F_rice_artifact3/s07_key_numbers.csv (full precision)
"""
import pandas as pd

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa

rows = []
def add(section, name, value, source, note=""):
    rows.append({"section": section, "quantity": name, "value": value, "source": source, "note": note})

rep = pd.read_csv(RES / "s02_reproduction.csv")
for _, r in rep.iterrows():
    add("1_reproduction", r["quantity"], f"stored={r['stored']}; SEQ_REF={r['reproduced_SEQ_REF']}; SEQ_FULL={r['reproduced_SEQ_FULL']}",
        "s02_reproduction.csv", "SEQ_REF = UP000059680 2026_03 only, original first-found rule; SEQ_FULL = all accessions resolvable")

s4 = json.loads((RES / "s04_opportunity_summary.json").read_text(encoding="utf-8"))
d = s4["descriptives"]
for k in ("A_positive_without_any_identified_Cys_peptide", "groups_with_zero_detectable_cys", "fraction_cys_detectable",
          "identified_sites_mapped", "identified_sites_on_detectable_cys",
          "observed_share_Cys_peptides_in_table_real",
          "insilico_share_Cys_peptides_fully_tryptic_7_30_identified_proteins_mean"):
    add("2_descriptives", k, d[k], "s04_opportunity_summary.json")
for k, v in d["labels_all_groups"].items():
    add("2_descriptives", f"n_positive_label_{k}", v, "s04_opportunity_summary.json")
for k, v in d["spearman"].items():
    add("2_descriptives", f"spearman_{k}", v, "s04_opportunity_summary.json")

rb = pd.read_csv(RES / "s04_rank_biserial_strata.csv")
for _, r in rb.iterrows():
    for c in ("r_before", "r_within_iBAQ_deciles", "r_within_rawintensity_deciles", "r_1to1_matched",
              "r_within_decile_x_opp_tertile", "r_within_decile_x_exact_n_det_cys"):
        add("3_rank_biserial", f"{r['label']}|{r['feature']}|{c}", r[c], "s04_rank_biserial_strata.csv")
bt = pd.read_csv(RES / "s04_rank_biserial_bootstrap.csv")
for _, r in bt.iterrows():
    add("3_rank_biserial_boot", f"{r['label']}|{r['feature']}|{r['strata']}", f"{r['r']} {r['ci95']}", "s04_rank_biserial_bootstrap.csv", "400 group-bootstrap resamples")
mb = json.loads((RES / "s04b_summary.json").read_text(encoding="utf-8"))
for r in mb["matched_on_abundance_and_opportunity"]:
    add("4_matched_abundance_opportunity", f"{r['label']}|r_length", f"{r['r_length']} {r['r_length_ci95_boot200']}", "s04b_summary.json")
    add("4_matched_abundance_opportunity", f"{r['label']}|r_cys", f"{r['r_cys']} {r['r_cys_ci95_boot200']}", "s04b_summary.json")
    add("4_matched_abundance_only", f"{r['label']}|r_length", r["abundance_only_1to1"]["r_length"], "s04b_summary.json")
    add("4_matched_abundance_only", f"{r['label']}|r_cys", r["abundance_only_1to1"]["r_cys"], "s04b_summary.json")
orr = pd.read_csv(RES / "s04_logistic_OR.csv")
for _, r in orr[orr.model.isin(["M1_iBAQ", "M2_iBAQ_len_cys", "M3_iBAQ_len_cys_opp", "Q2_raw_len_cys",
                                "Q3_raw_len_cys_opp"])].iterrows():
    add("5_logistic", f"{r['label']}|{r['model']}|{r['term']}", f"{r['OR']} [{r['ci_low']}, {r['ci_high']}] p={r['p']}", "s04_logistic_OR.csv")
cl = pd.read_csv(RES / "s04_cloglog_opportunity.csv")
for _, r in cl.iterrows():
    add("6_cloglog", f"{r['label']}|{r['model']}|{r['term']}", f"{r['coef']} [{r['ci_low']}, {r['ci_high']}] p={r['p']}", "s04_cloglog_opportunity.csv")
pc = pd.read_csv(RES / "s04_per_cysteine_logit.csv")
for _, r in pc.iterrows():
    add("7_per_cysteine", f"{r['model']}|{r['term']}", f"{r['OR']} [{r['ci_low']}, {r['ci_high']}] p={r['p']}", "s04_per_cysteine_logit.csv")
pcd = pd.read_csv(RES / "s04_per_cysteine_by_decile.csv")
for _, r in pcd.iterrows():
    add("7_per_cysteine_decile", f"iBAQ_decile_{int(r['iBAQ_decile'])}", f"{r['n_sites']}/{r['n_det_cys']}={r['p_site']}", "s04_per_cysteine_by_decile.csv")
for k, v in s4["pure_opportunity_simulation_label_S"].items():
    add("8_simulation", k, json.dumps(v), "s04_opportunity_summary.json", "label S_map observed vs 1,000 (or 200/100) simulated label sets")
lr = pd.read_csv(RES / "s04c_label_ratios.csv")
for _, r in lr.iterrows():
    add("9_ratios", f"{r['positive_set']}|{r['measure']}", f"{r['median_ratio']} [{r['ci_low']}, {r['ci_high']}] n={r['n_pos']}", "s04c_label_ratios.csv")
ch = pd.read_csv(RES / "s05_filter_chain.csv")
for _, r in ch.iterrows():
    add("10_rice_chain", r["step"], f"{r['n_sites']} sites / {r['n_proteins']} proteins / {r['n_cys_in_proteins']} Cys; equals stored: {r['equals_stored_sites']}", "s05_filter_chain.csv", r["note"])
gr = pd.read_csv(RES / "s05_rule_grid.csv")
near = gr.loc[(gr.n_sites - 1691).abs().sort_values().index[:4]]
for _, r in near.iterrows():
    add("11_rule_grid_nearest_1691", f"{r['protein_assignment']}|{r['peptides']}|{r['decoys_contaminants']}|{r['uniqueness']}|{r['evidence']}", r["n_sites"], "s05_rule_grid.csv", "exploratory grid of 72 rules")
# ---------------- revision after adversarial verification (s08, s08b, s09) ----------------
s8 = json.loads((RES / "s08_summary.json").read_text(encoding="utf-8"))
for k, v in s8["descriptives"].items():
    add("12_residual_descriptives", k, json.dumps(v) if isinstance(v, dict) else v, "s08_summary.json")
for k, v in s8["first_listed_member"].items():
    add("12_first_listed_member", k, json.dumps(v) if isinstance(v, dict) else v, "s08_summary.json")
lo = pd.read_csv(RES / "s08_logit_undet.csv")
for _, r in lo[lo.term.isin(["log2_len", "cys", "n_undet"])].iterrows():
    add("13_logit_undet", f"{r['label']}|{r['model']}|{r['term']}", f"{r['OR']} [{r['ci_low']}, {r['ci_high']}] p={r['p']} aic={r['aic']}",
        "s08_logit_undet.csv", "M3u = iBAQ + log2(1+n_det) + log2 length + n_undet; M5u = spline iBAQ + C(n_det) (1..20) + log2 length + n_undet")
for _, r in pd.read_csv(RES / "s08_logit_LR.csv").iterrows():
    add("13_logit_LR", f"{r['label']}|{r['small']}->{r['big']}", f"LR={r['LR']} df={r['df']} p={r['p']} dAIC={r['dAIC']}", "s08_logit_LR.csv")
for _, r in pd.read_csv(RES / "s08_matched_pairs.csv").iterrows():
    for k in ("r_len", "r_cys", "r_undet"):
        add("14_matched_iBAQ_exact_ndet", f"{r['label']}|{k}", f"{r[k]} {r[k + '_ci95']}", "s08_matched_pairs.csv", "400 group-bootstrap resamples, matching redone")
    add("14_matched_iBAQ_exact_ndet", f"{r['label']}|mean_undet_pos_neg", f"{r['mean_undet_pos']} vs {r['mean_undet_neg']}; Wilcoxon p={r['wilcoxon_p_undet']}", "s08_matched_pairs.csv")
for _, r in pd.read_csv(RES / "s08_windows_undet.csv").iterrows():
    add("15_windows_undet", f"{r['window']}|{r['label']}", f"M3u OR_undet={r['M3u_OR_undet']} {r['M3u_OR_undet_ci']}; M5u={r['M5u_OR_undet']} {r['M5u_OR_undet_ci']}; "
        f"M3u OR_len={r['M3u_OR_len']}; matched r_undet={r['matched_r_undet']} r_cys={r['matched_r_cys']} r_len={r['matched_r_len']}; frac Cys detectable={r['frac_cys_detectable']}",
        "s08_windows_undet.csv")
for _, r in pd.read_csv(RES / "s08_sensitivity_undet.csv").iterrows():
    add("16_sensitivity_undet", f"{r['label']}|{r['model']}|{r['term']}", f"{r['OR']} [{r.get('ci_low')}, {r.get('ci_high')}] n={r['n']} n_pos={r['n_pos']}", "s08_sensitivity_undet.csv")
for _, r in pd.read_csv(RES / "s08_cloglog_free.csv").iterrows():
    add("17_cloglog_free", f"{r['label']}|{r['model']}|{r['term']}", f"{r['coef']} [{r['ci_low']}, {r['ci_high']}] p={r['p']} aic={r['aic']} expected_pos={r['expected_pos']}", "s08_cloglog_free.csv")
for _, r in pd.read_csv(RES / "s08_cloglog_LR.csv").iterrows():
    add("17_cloglog_LR", f"{r['label']}|{r['small']}->{r['big']}", f"LR={r['LR']} df={r['df']} p={r['p']} dAIC={r['dAIC']}", "s08_cloglog_LR.csv")
for _, r in pd.read_csv(RES / "s08_per_cysteine_undet.csv").iterrows():
    add("18_per_cysteine_undet", f"{r['model']}|{r['term']}", f"{r['OR']} [{r['ci_low']}, {r['ci_high']}] p={r['p']}", "s08_per_cysteine_undet.csv", "cluster-robust by group")
for _, r in pd.read_csv(RES / "s08b_null_summary.csv").iterrows():
    add("19_calibrated_null", f"{r['label']}|{r['null']}|{r['statistic']}", f"observed={r['observed']} null_mean={r['null_mean']} "
        f"[{r['null_2.5']}, {r['null_97.5']}] p_ge={r['p_ge_observed']} p_le={r['p_le_observed']} reps={r['reps']}",
        "s08b_null_summary.csv", f"expected positives {r['expected_pos_null']} vs observed {r['observed_pos']}; groups with >=1 detectable Cys")
for _, r in pd.read_csv(RES / "s08c_undet_designs.csv").iterrows():
    add("19b_designs_all_features", f"{r['label']}|{r['feature']}",
        f"before={r['r_before']} deciles={r['r_within_iBAQ_deciles']} 1:1_iBAQ={r['r_1to1_iBAQ_only']} dec_x_tertile={r['r_within_decile_x_opp_tertile']}",
        "s08c_undet_designs.csv")
s9 = json.loads((RES / "s09_artifact2_trace.json").read_text(encoding="utf-8"))
add("20_artifact2_trace", "rules_hitting_any_target", json.dumps(s9["rules_with_any_hit"]), "s09_artifact2_trace.json", "targets 756/600/689/545")
add("20_artifact2_trace", "range_asym_sites_over_60_rules", json.dumps(s9["range_asym_sites"]), "s09_artifact2_trace.json")
for r in s9["stored_vs_complete_list"]:
    add("20_artifact2_trace", f"{r['list']}|{r['enzyme']}|mc{r['missed']}|{r['min_len']}-{r['max_len']}",
        f"asym {r['asym_sites']} sites/{r['asym_prots']} proteins; sym {r['sym_sites']}/{r['sym_prots']}", "s09_artifact2_stored_vs_complete.csv")
# ---------------- revision round 2: detectable versus total cysteine count (s10, s10b, s11) ----------------
s10 = json.loads((RES / "s10_summary.json").read_text(encoding="utf-8"))
for k in ("spearman_ncys_ndet_all", "spearman_ncys_ndet_ndet_ge1", "frac_cys_detectable", "frac_groups_no_undetectable"):
    add("21_count_descriptives", k, s10[k], "s10_summary.json")
for k, v in s10["top_decile_fraction_by_tertile"].items():
    add("21_count_descriptives", f"top_iBAQ_decile_fraction_positive_by_tertile|{k}", json.dumps(v), "s10_summary.json")
cmp_cols = ["dAIC_D_minus_T", "dAIC0_D_minus_T", "LR_det_given_total", "p_det_given_total", "LR_total_given_det", "p_total_given_det",
            "p_det_given_total_nolen", "p_total_given_det_nolen", "LR_len_given_T", "p_len_given_T", "LR_len_given_D", "p_len_given_D",
            "len_OR_per_doubling_T", "len_OR_per_doubling_D", "len_coef_T", "len_coef_D", "exponent_D0", "exponent_T0",
            "exponent_D", "exponent_T", "B_coef[l2d]", "B_coef[l2c]", "B_coef[log_det]", "B_coef[log_cys]", "frac_det_given_total_logit"]
for _, r in pd.read_csv(RES / "s10_comparisons.csv").iterrows():
    for c in cmp_cols:
        if c in r and not (isinstance(r[c], float) and pd.isna(r[c])):
            add("22_detectable_vs_total", f"{r['family']}|{r['label']}|{c}", r[c], "s10_comparisons.csv",
                "D = iBAQ + log2 length + detectable count; T = same with total count; B = both; *0 = without length")
for _, r in pd.read_csv(RES / "s10_linear_symmetric.csv").iterrows():
    for c in ("OR_per_detectable_cys", "OR_per_undetectable_cys", "p_equal", "OR_per_cys_total_only",
              "OR_len_per_doubling_total_only", "p_len_total_only", "AIC_split", "AIC_total"):
        add("23_linear_symmetric", f"{r['sample']}|{r['label']}|{c}", r[c], "s10_linear_symmetric.csv",
            "y ~ log10 iBAQ + log2 length + n_det + n_undet (split) or + n_cys (total)")
for _, r in pd.read_csv(RES / "s10_windows.csv").iterrows():
    add("24_windows_detectable_vs_total", f"{r['window']}|{r['label']}",
        f"dAIC(D-T)={r['dAIC_D_minus_T']}; P det|total={r['p_det_given_total']}; P total|det={r['p_total_given_det']}; "
        f"OR det={r['OR_per_detectable_cys']}; OR undet={r['OR_per_undetectable_cys']}; P equal={r['p_equal']}; "
        f"frac Cys detectable={r['frac_cys_detectable']}; Spearman={r['spearman_ncys_ndet']}", "s10_windows.csv")
for _, r in pd.read_csv(RES / "s10_mirror_matched.csv").iterrows():
    for k in ("r_det", "r_undet", "r_len"):
        add("25_mirror_matched_iBAQ_exact_ncys", f"{r['label']}|{k}", f"{r[k]} {r[k + '_ci95']}", "s10_mirror_matched.csv",
            "1:1 on log10 iBAQ within exact total Cys count (cap 25); 400 group-bootstrap resamples, matching redone")
    add("25_mirror_matched_iBAQ_exact_ncys", f"{r['label']}|mean_det_pos_neg",
        f"{r['mean_det_pos']} vs {r['mean_det_neg']}; Wilcoxon p={r['wilcoxon_p_det']}; pairs equal n_cys={r['frac_pairs_equal_ncys']}",
        "s10_mirror_matched.csv")
    add("25_mirror_matched_iBAQ_exact_ncys", f"{r['label']}|cells_dec_x_ncys",
        f"r_det={r['r_cells_dec_x_ncys_det']} r_len={r['r_cells_dec_x_ncys_len']} cells={r['n_cells']}", "s10_mirror_matched.csv")
for _, r in pd.read_csv(RES / "s10_per_cysteine.csv").iterrows():
    ors = {c: r[c] for c in r.index if c.startswith("OR[") and isinstance(r[c], str)}
    add("26_per_cysteine_split", r["model"], f"{json.dumps(ors)}; p_equal={r.get('p_equal')}; aic={r['aic']}; n={r['n']}",
        "s10_per_cysteine.csv", "unit = detectable cysteine; protein-clustered SEs")
for _, r in pd.read_csv(RES / "s10b_null_summary.csv").iterrows():
    add("27_nulls_detectable_and_total", f"{r['label']}|{r['null']}|{r['statistic']}",
        f"observed={r['observed']} null_mean={r['null_mean']} [{r['null_2.5']}, {r['null_97.5']}] p_ge={r['p_ge_observed']} "
        f"p_le={r['p_le_observed']} p_two={r['p_two_sided']} reps={r['reps']}", "s10b_null_summary.csv",
        f"expected positives {r['expected_pos_null']} vs observed {r['observed_pos']}; groups with >=1 detectable Cys")
s10b = json.loads((RES / "s10b_null_summary.json").read_text(encoding="utf-8"))
for k, v in s10b["contrasts_outside_null_range"].items():
    add("27_nulls_detectable_and_total", f"outside_95_range|{k}", f"{len(v['outside_95_range'])}/{v['reported_contrasts']}: "
        f"{', '.join(v['outside_95_range'])}", "s10b_null_summary.json")
add("27_nulls_detectable_and_total", "s08b_combinations_reproduced_exactly", s10b["s08b_combinations_reproduced_exactly"], "s10b_null_summary.json")
for _, r in pd.read_csv(RES / "s11_rice_chain.csv").iterrows():
    add("28_rice_chain_start_position", r["step"], f"{r['n_sites']} sites / {r['n_proteins']} proteins / {r['n_cys_in_proteins']} Cys; "
        f"equals stored: {r['equals_stored_sites']}", "s11_rice_chain.csv", r["note"])
s11 = json.loads((RES / "s11_rice_chain.json").read_text(encoding="utf-8"))
add("28_rice_chain_start_position", "lost_at_step5", json.dumps({k: v for k, v in s11["lost_at_step5"].items() if k != "entry_records"}),
    "s11_rice_chain.json", "entry records: " + json.dumps(s11["lost_at_step5"]["entry_records"]["deletedReason"]))
add("28_rice_chain_start_position", "self_audit_exclusion", json.dumps(s11["self_audit_exclusion_repository_audit"]), "s11_rice_chain.json",
    "repository results/self_audit_public_cohorts_2026-09-20_audit.json, gates/coverage_rice_PXD072089")
add("28_rice_chain_start_position", "difference_to_s05", s11["difference_between_conventions"], "s11_rice_chain.json")
# ---------------- revision round 3: pooled nulls (s12), indicator test (s13), decoy rows of the rice chain (s11) --------
if "round3_decoys_and_contaminants" in s11:
    add("28_rice_chain_start_position", "round3_decoys_and_contaminants", json.dumps(s11["round3_decoys_and_contaminants"]),
        "s11_rice_chain.json", "round 2's step 0 (1,769 / 1,170) skipped the decoy rows without a start position")
if "round3_lost_entry_records" in s11:
    add("28_rice_chain_start_position", "round3_lost_entry_records", json.dumps(s11["round3_lost_entry_records"]),
        "s11_rice_chain.json", "UniProtKB entry records (saved by s01b) of the 400 accessions lost at step 5")
for _, r in pd.read_csv(RES / "s12_null_pooled_summary.csv").iterrows():
    add("29_nulls_pooled_5000", f"{r['label']}|{r['null']}|{r['statistic']}",
        f"observed={r['observed']} pooled=[{r['pooled_lo']}, {r['pooled_hi']}] p_two={r['pooled_p_two']} "
        f"inside={r['pooled_inside']} edge={r['pooled_edge']}; s10b=[{r['s10b_lo']}, {r['s10b_hi']}] p_two={r['s10b_p_two']}; "
        f"extension=[{r['ext_lo']}, {r['ext_hi']}] p_two={r['ext_p_two']}", "s12_null_pooled_summary.csv",
        "5,000 label sets per combination = 1,000 (s10b) + 4,000 (s12, SEED + 840..851); groups with >=1 detectable Cys")
s12 = json.loads((RES / "s12_null_pooled_summary.json").read_text(encoding="utf-8"))
for k, v in s12["counts"].items():
    add("29_nulls_pooled_5000", f"outside_95_range|{k}",
        f"pooled {v['n_outside_pooled']}/{v['reported_contrasts']}: {', '.join(v['outside_pooled'])}; edge (P in [0.04, 0.06]): "
        f"{json.dumps(v['edge_pooled'])}; 1,000 sets: {v['n_outside_s10b']}; smallest P inside, edge excluded: "
        f"{v['min_p_two_inside_pooled_excluding_edge']}", "s12_null_pooled_summary.json")
for _, r in pd.read_csv(RES / "s13_indicator_test.csv").iterrows():
    add("30_indicator_test", r["label"],
        f"P(det | total, 1[n_det>=1])={r['p_det_given_total_ind']} (no length {r['p_det_given_total_ind_nolen']}); "
        f"P(total | det, 1[n_det>=1])={r['p_total_given_det_ind']}; dAIC(DI-TI)={r['dAIC_DI_minus_TI']}; "
        f"positives without detectable Cys {r['n_pos_without_detectable_cys']}/{r['n_groups_without_detectable_cys']}",
        "s13_indicator_test.csv", "logistic, all 7,692 groups, log10 iBAQ + log2 length")
K = pd.DataFrame(rows)
K.to_csv(RES / "s07_key_numbers.csv", index=False)
print(len(K), "rows")
