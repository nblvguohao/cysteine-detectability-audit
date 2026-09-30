"""Step 1 - absolute AUCs behind Figure 5a, verification of the plotted shares, derived
caliber/composition diagnostics, the NEG_B within-protein table (task 3) and a proposed
Source Data addition.

POST HOC revision analysis (2026-09-30), item I_fig5a_estimand. Not registered, not pre-specified.
From STORED OUTPUTS ONLY: no model is fitted and no site-level data are read. Every number
below is either copied from results/ptm_detectability_share.csv (4 decimals, as stored) or
an arithmetic function of those stored numbers; derived quantities without a stored
bootstrap carry no interval and are labelled as point estimates.

Run:  python -B s01_tables.py
"""
from __future__ import annotations

import csv
import io
import json
import math
import re

import numpy as np
import pandas as pd

import common as C

pd.set_option("display.width", 250)


def load_share() -> pd.DataFrame:
    df = pd.read_csv(C.SHARE_CSV, dtype={"dataset_id": str})
    df["neg"] = df["negative_construction"].map(C.NEG)
    df["in_fig5a"] = df["dataset_id"].isin(C.LABEL).astype(int)
    df["cohort_label"] = df["dataset_id"].map(C.LABEL).fillna("(not in Fig. 5a)")
    return df


def long_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        for cal, est, lo, hi in (("global", "auc_global", "auc_global_lo", "auc_global_hi"),
                                 ("within_protein", "auc_within_unit", "auc_within_unit_lo", "auc_within_unit_hi")):
            rows.append({
                "dataset_id": r.dataset_id, "cohort_label": r.cohort_label, "in_fig5a": r.in_fig5a,
                "negative_set": r.neg, "feature_set": r.feature_set, "n_features": r.n_features,
                "caliber": cal, "auc": r[est], "auc_lo": r[lo], "auc_hi": r[hi],
                "interval_excludes_0.5": (None if pd.isna(r[lo]) else bool(r[lo] > 0.5 or r[hi] < 0.5)),
                "n_positive": r.n_positive, "n_negative": r.n_negative, "n_units": r.n_units,
                "n_units_with_positive": r.n_units_with_positive,
                "n_units_both_classes": r.n_informative_units, "grouping_unit": r.grouping_unit,
                "usable": r.usable, "exclusion_reason": r.exclusion_reason,
                "source": "results/ptm_detectability_share.csv (stored; 4 decimals)",
                "interval_note": ("95% percentile cluster bootstrap over grouping units, 5,000 replicates, seed "
                                  "20260915, out-of-fold scores not refitted per replicate" if r.usable == 1 else ""),
            })
    return pd.DataFrame(rows)


def share_rounding_bound(v: float, d: float, delta: float = 5e-5) -> float:
    """Exact worst-case |share(v', d') - share(v, d)| over |v'-v|, |d'-d| <= delta, plus the
    rounding of the stored share itself (delta). share is monotone in v and in d on each side."""
    base = (v - 0.5) / (d - 0.5)
    worst = 0.0
    for dv in (-delta, delta):
        for dd in (-delta, delta):
            worst = max(worst, abs((v + dv - 0.5) / (d + dd - 0.5) - base))
    return worst + delta


def wide_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for did, lab in C.COHORTS:
        for neg in ("NEG_A", "NEG_B"):
            sub = df[(df.dataset_id == did) & (df.neg == neg)].set_index("feature_set")
            v, d, c = sub.loc["VIS10"], sub.loc["DIG25"], sub.loc["CPL15"]
            row = {"dataset_id": did, "cohort_label": lab, "negative_set": neg,
                   "usable": int(v.usable), "exclusion_reason": v.exclusion_reason if isinstance(v.exclusion_reason, str) else "",
                   "n_positive": int(v.n_positive), "n_negative": int(v.n_negative),
                   "n_units": int(v.n_units), "n_units_with_positive": int(v.n_units_with_positive),
                   "grouping_unit": v.grouping_unit}
            if v.usable != 1:
                rows.append(row)
                continue
            nb = int(v.n_informative_units)
            row.update({
                "n_units_both_classes": nb,
                "frac_units_both_classes": nb / v.n_units,
                "n_units_positive_only": int(v.n_units_with_positive - nb),
                "n_units_negative_only": int(v.n_units - v.n_units_with_positive),
                "frac_units_negative_only": (v.n_units - v.n_units_with_positive) / v.n_units,
            })
            for fs, s in (("VIS10", v), ("DIG25", d), ("CPL15", c)):
                row.update({f"{fs}_auc_global": s.auc_global, f"{fs}_auc_global_lo": s.auc_global_lo,
                            f"{fs}_auc_global_hi": s.auc_global_hi,
                            f"{fs}_auc_within": s.auc_within_unit, f"{fs}_auc_within_lo": s.auc_within_unit_lo,
                            f"{fs}_auc_within_hi": s.auc_within_unit_hi})
            # stored share (Fig. 5a) and its verification from the stored rounded AUCs
            share_rec = (v.auc_global - 0.5) / (d.auc_global - 0.5)
            bound = share_rounding_bound(v.auc_global, d.auc_global)
            row.update({
                "share_global_stored": v.share_visibility, "share_global_lo_stored": v.share_visibility_lo,
                "share_global_hi_stored": v.share_visibility_hi,
                "share_global_recomputed_from_rounded_aucs": share_rec,
                "share_abs_deviation": abs(share_rec - v.share_visibility),
                "share_rounding_bound": bound,
                "share_recomputation_within_rounding": bool(abs(share_rec - v.share_visibility) <= bound),
                "dig25_global_minus_0.5": d.auc_global - 0.5,
                "vis10_global_minus_0.5": v.auc_global - 0.5,
                "dig25_minus_vis10_global": d.auc_global - v.auc_global,
                "dig25_minus_vis10_within": d.auc_within_unit - v.auc_within_unit,
                # derived, point estimates only (no stored bootstrap for these)
                "complement_ratio_global_point": (c.auc_global - 0.5) / (d.auc_global - 0.5),
                "complement_ratio_within_point": (c.auc_within_unit - 0.5) / (d.auc_within_unit - 0.5),
                "share_within_point": ((v.auc_within_unit - 0.5) / (d.auc_within_unit - 0.5)
                                       if (d.auc_within_unit - 0.5) >= C.SHARE_DENOM_FLOOR else np.nan),
                "share_within_denominator": d.auc_within_unit - 0.5,
                "global_minus_within_VIS10": v.auc_global - v.auc_within_unit,
                "global_minus_within_DIG25": d.auc_global - d.auc_within_unit,
                "global_minus_within_CPL15": c.auc_global - c.auc_within_unit,
            })
            row["share_plus_complement_global_point"] = share_rec + row["complement_ratio_global_point"]
            # the share interval conditions on replicates whose DIG25 AUC >= 0.52 (floor 0.02);
            # if the stored 2.5th percentile of DIG25 is below 0.52, more than 2.5% of replicates were dropped
            row["share_interval_conditional_on_floor"] = bool(d.auc_global_lo - 0.5 < C.SHARE_DENOM_FLOOR)
            se = (d.auc_global_hi - d.auc_global_lo) / (2 * 1.959964)
            from scipy.stats import norm
            row["approx_frac_replicates_below_floor_normal"] = float(norm.cdf((0.5 + C.SHARE_DENOM_FLOOR - d.auc_global) / se))
            rows.append(row)
    return pd.DataFrame(rows)


def source_data_check(wide: pd.DataFrame) -> pd.DataFrame:
    sd = pd.read_csv(C.SOURCE_DATA_FIG5)
    sd = sd[(sd.figure == "Fig5") & (sd.panel == "a")]
    out = []
    for _, r in wide[wide.usable == 1].iterrows():
        key = f"{r.dataset_id}|{r.negative_set}"
        got = sd[sd.row == key].set_index("field")["value"].astype(float)
        for fld, stored in (("share_visibility", r.share_global_stored), ("lo", r.share_global_lo_stored),
                            ("hi", r.share_global_hi_stored)):
            val = got.get(fld, np.nan)
            out.append({"row": key, "field": fld, "source_data_value": val, "stored_share_csv_value": stored,
                        "match": bool(abs(val - stored) < 1e-9)})
    return pd.DataFrame(out)


def manuscript_number_check(wide: pd.DataFrame) -> dict:
    """Do the share values quoted in the Results paragraph match the stored table?"""
    tex = C.MANUSCRIPT.read_text(encoding="utf-8")
    para = next(line for line in tex.splitlines() if line.startswith("The same published site table can be scored"))
    nums = [float(x) for x in re.findall(r"\d+\.\d+", para)]
    missing = []
    for _, r in wide[wide.usable == 1].iterrows():
        for val in (r.share_global_stored, r.share_global_lo_stored, r.share_global_hi_stored):
            if not any(abs(abs(val) - n) < 5e-5 for n in nums):
                missing.append((r.dataset_id, r.negative_set, val))
    return {"paragraph_prefix": para[:90], "n_numbers_in_paragraph": len(nums),
            "stored_values_checked": int(3 * (wide.usable == 1).sum()), "stored_values_not_found": missing,
            "all_found": not missing}


def parse_snote7() -> pd.DataFrame:
    text = C.SNOTE7.read_text(encoding="utf-8")
    block = text[text.index("| label | comparability |"):]
    lines = [ln for ln in block.splitlines() if ln.startswith("|")]
    rows = [[c.strip() for c in ln.strip("|").split("|")] for ln in lines]
    hdr, body = rows[0], [r for r in rows[2:] if len(r) == len(rows[0])]
    df = pd.DataFrame(body, columns=hdr)
    inv = {v: k for k, v in C.LABEL.items()}
    df["dataset_id"] = df["label"].map(inv)
    return df


def sd7_crosscheck(wide: pd.DataFrame) -> pd.DataFrame:
    """Supplemental Data 7 lists the cohorts admitted to the share; check n and usability agree."""
    sd7 = pd.read_csv(C.SD7_CENSUS)
    sd7 = sd7[sd7.source_of_row == "computed_in_this_study"].copy()
    sd7["dataset_id"] = sd7.record_id.str.split(":").str[1]
    sd7["negative_set"] = sd7.record_id.str.split(":").str[2].map(C.NEG)
    out = []
    for _, r in wide.iterrows():
        m = sd7[(sd7.dataset_id == r.dataset_id) & (sd7.negative_set == r.negative_set)]
        row = {"dataset_id": r.dataset_id, "negative_set": r.negative_set, "in_sd7": len(m) == 1}
        if len(m) == 1:
            m = m.iloc[0]
            row.update({"sd7_method_name": m.method_name, "sd7_species": m.species,
                        "sd7_n_sites_reported": m.n_sites_reported, "share_csv_n_positive": r.n_positive,
                        "n_match": bool(int(m.n_sites_reported) == int(r.n_positive)),
                        "sd7_usable_for_share": int(m.usable_for_share), "share_csv_usable": int(r.usable),
                        "usable_match": bool(int(m.usable_for_share) == int(r.usable)),
                        "sd7_negative_semantics": m.negative_semantics})
        out.append(row)
    return pd.DataFrame(out)


def negA_vs_negB(wide: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for did, lab in C.COHORTS:
        a = wide[(wide.dataset_id == did) & (wide.negative_set == "NEG_A")].iloc[0]
        b = wide[(wide.dataset_id == did) & (wide.negative_set == "NEG_B")].iloc[0]
        if b.usable != 1:
            rows.append({"dataset_id": did, "cohort_label": lab, "neg_b_available": 0})
            continue
        row = {"dataset_id": did, "cohort_label": lab, "neg_b_available": 1}
        for col in ("VIS10_auc_global", "DIG25_auc_global", "CPL15_auc_global", "VIS10_auc_within",
                    "DIG25_auc_within", "CPL15_auc_within", "share_global_stored", "share_within_point"):
            row[f"{col}_NEG_A"] = a[col]
            row[f"{col}_NEG_B"] = b[col]
            row[f"{col}_B_minus_A"] = b[col] - a[col]
        row["neg_b_lower_share_global"] = bool(b.share_global_stored < a.share_global_stored)
        row["neg_b_lower_share_within_point"] = bool(b.share_within_point < a.share_within_point)
        rows.append(row)
    return pd.DataFrame(rows)


def task3_table(wide: pd.DataFrame, snote7: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for did, lab in C.COHORTS:
        b = wide[(wide.dataset_id == did) & (wide.negative_set == "NEG_B")].iloc[0]
        a = wide[(wide.dataset_id == did) & (wide.negative_set == "NEG_A")].iloc[0]
        s7 = snote7[snote7.dataset_id == did]
        base = {"dataset_id": did, "cohort_label": lab,
                "snote7_comparability": s7.comparability.iloc[0] if len(s7) else "",
                "snote7_observed_class_definition": s7.observed_class_definition.iloc[0] if len(s7) else ""}
        if b.usable != 1:
            rows.append({**base, "neg_b_available": 0, "note": b.exclusion_reason})
            continue
        rows.append({**base, "neg_b_available": 1,
                     "n_positive": b.n_positive, "n_negative_NEG_B": b.n_negative,
                     "n_units_both_classes": b.n_units_both_classes,
                     "DIG25_auc_within": b.DIG25_auc_within, "DIG25_auc_within_lo": b.DIG25_auc_within_lo,
                     "DIG25_auc_within_hi": b.DIG25_auc_within_hi,
                     "DIG25_within_interval_excludes_0.5": bool(b.DIG25_auc_within_lo > 0.5),
                     "DIG25_within_above_chance": b.DIG25_auc_within - 0.5,
                     "VIS10_auc_within": b.VIS10_auc_within, "VIS10_auc_within_lo": b.VIS10_auc_within_lo,
                     "VIS10_auc_within_hi": b.VIS10_auc_within_hi,
                     "CPL15_auc_within": b.CPL15_auc_within, "CPL15_auc_within_lo": b.CPL15_auc_within_lo,
                     "CPL15_auc_within_hi": b.CPL15_auc_within_hi,
                     "DIG25_auc_within_NEG_A": a.DIG25_auc_within,
                     "retained_fraction_of_NEG_A_within_above_chance_point": (b.DIG25_auc_within - 0.5) / (a.DIG25_auc_within - 0.5),
                     "DIG25_auc_global_NEG_B": b.DIG25_auc_global,
                     "note": "point ratio has no interval (per-unit bootstrap draws were not stored)"})
    return pd.DataFrame(rows)


def source_data_addition(wide: pd.DataFrame) -> pd.DataFrame:
    """Rows in the submitted Source Data schema. Values are written exactly as stored (4 decimals)."""
    src = "analysis output: ptm_detectability_share.csv"
    rows = []
    for _, r in wide.iterrows():
        key = f"{r.dataset_id}|{r.negative_set}"
        if r.usable != 1:
            continue
        for fld in ("n_positive", "n_negative", "n_units", "n_units_with_positive", "n_units_both_classes"):
            rows.append(("Fig5", "a", key, fld, str(int(r[fld])), src, "count; post hoc revision addition"))
        for fs in ("VIS10", "DIG25", "CPL15"):
            for cal, tag in (("global", "auc_global"), ("within", "auc_within_protein")):
                for suff, lab in (("", ""), ("_lo", "_lo"), ("_hi", "_hi")):
                    val = r[f"{fs}_auc_{cal}{suff}"]
                    rows.append(("Fig5", "a", key, f"{tag}_{fs}{lab}", f"{val:.4f}".rstrip("0").rstrip("."),
                                 src, "absolute out-of-fold AUC; post hoc revision addition"))
    return pd.DataFrame(rows, columns=["figure", "panel", "row", "field", "value", "source_table",
                                      "specification_label"])


def ranges(wide: pd.DataFrame, t3: pd.DataFrame) -> dict:
    u = wide[wide.usable == 1]
    A, B = u[u.negative_set == "NEG_A"], u[u.negative_set == "NEG_B"]
    notS2 = lambda x: x[x.dataset_id != "qtrp_S2_ph5"]
    five = lambda x: x[~x.dataset_id.isin(["qtrp_S2_ph5", "fps2020_ath_sulfenyl"])]
    rng = lambda s: [float(s.min()), float(s.max())]
    out = {
        "NEG_A": {
            "n_arms": int(len(A)),
            "VIS10_global_all": rng(A.VIS10_auc_global), "VIS10_global_excl_S2": rng(notS2(A).VIS10_auc_global),
            "VIS10_within_all": rng(A.VIS10_auc_within),
            "DIG25_global_all": rng(A.DIG25_auc_global), "DIG25_within_all": rng(A.DIG25_auc_within),
            "share_global_all": rng(A.share_global_stored), "share_within_point_all": rng(A.share_within_point),
            "complement_ratio_global_all": rng(A.complement_ratio_global_point),
            "share_plus_complement_global": rng(A.share_plus_complement_global_point),
            "frac_units_both_classes": rng(A.frac_units_both_classes),
            "global_minus_within_DIG25": rng(A.global_minus_within_DIG25),
        },
        "NEG_B": {
            "n_arms": int(len(B)),
            "VIS10_global_five": rng(five(B).VIS10_auc_global), "VIS10_within_five": rng(five(B).VIS10_auc_within),
            "VIS10_global_all": rng(B.VIS10_auc_global), "VIS10_within_all": rng(B.VIS10_auc_within),
            "DIG25_global_all": rng(B.DIG25_auc_global), "DIG25_within_all": rng(B.DIG25_auc_within),
            "DIG25_within_excl_S2": rng(notS2(B).DIG25_auc_within),
            "share_global_five": rng(five(B).share_global_stored),
            "share_within_point_excl_S2": rng(notS2(B).share_within_point),
            "complement_ratio_global_all": rng(B.complement_ratio_global_point),
            "share_plus_complement_global": rng(B.share_plus_complement_global_point),
            "frac_units_both_classes": rng(B.frac_units_both_classes),
            "frac_units_negative_only": rng(B.frac_units_negative_only),
            "global_minus_within_DIG25": rng(B.global_minus_within_DIG25),
        },
        "task3_DIG25_within_NEG_B_excl_S2": rng(t3[(t3.neg_b_available == 1) & (t3.dataset_id != "qtrp_S2_ph5")].DIG25_auc_within),
        "task3_all_intervals_above_0.5_excl_S2": bool(t3[(t3.neg_b_available == 1) & (t3.dataset_id != "qtrp_S2_ph5")]["DIG25_within_interval_excludes_0.5"].all()),
    }
    return out


def main() -> None:
    df = load_share()
    long = long_table(df)
    long.to_csv(C.OUT / "t1_absolute_auc_long.csv", index=False)
    wide = wide_table(df)
    wide.to_csv(C.OUT / "t1_absolute_auc_by_arm.csv", index=False)

    ver = wide[wide.usable == 1][["dataset_id", "cohort_label", "negative_set", "VIS10_auc_global", "DIG25_auc_global",
                                  "share_global_stored", "share_global_recomputed_from_rounded_aucs",
                                  "share_abs_deviation", "share_rounding_bound",
                                  "share_recomputation_within_rounding", "share_global_lo_stored",
                                  "share_global_hi_stored", "share_interval_conditional_on_floor",
                                  "approx_frac_replicates_below_floor_normal"]]
    sdc = source_data_check(wide)
    ver.to_csv(C.OUT / "t1_share_verification.csv", index=False)
    sdc.to_csv(C.OUT / "t1_source_data_match.csv", index=False)
    msc = manuscript_number_check(wide)

    s7 = parse_snote7()
    sd7 = sd7_crosscheck(wide)
    sd7.to_csv(C.OUT / "t1_sd7_crosscheck.csv", index=False)
    pair = negA_vs_negB(wide)
    pair.to_csv(C.OUT / "t2_negA_vs_negB_paired_points.csv", index=False)
    comp = wide[wide.usable == 1][["dataset_id", "cohort_label", "negative_set", "n_units", "n_units_both_classes",
                                   "frac_units_both_classes", "n_units_positive_only", "n_units_negative_only",
                                   "frac_units_negative_only", "global_minus_within_VIS10", "global_minus_within_DIG25",
                                   "global_minus_within_CPL15", "share_global_stored", "share_within_point",
                                   "share_within_denominator", "complement_ratio_global_point",
                                   "complement_ratio_within_point", "share_plus_complement_global_point"]]
    comp.to_csv(C.OUT / "t2_caliber_composition_and_overlap.csv", index=False)
    t3 = task3_table(wide, s7)
    t3.to_csv(C.OUT / "t3_negB_within_protein_digest.csv", index=False)

    add = source_data_addition(wide)
    add.to_csv(C.OUT / "Source_Data_Fig5a_absolute_auc_ADDITION_PROPOSED.csv", index=False, lineterminator="\n")
    # full proposed file = submitted file bytes + the addition rows, same column order and newline style
    raw = C.SOURCE_DATA_FIG5.read_bytes()
    nl = b"\r\n" if b"\r\n" in raw else b"\n"
    buf = io.StringIO()
    add.to_csv(buf, index=False, header=False, lineterminator=nl.decode())
    body = raw if raw.endswith(nl) else raw + nl
    (C.OUT / "Source_Data_Fig5_three_axes_PROPOSED.csv").write_bytes(body + buf.getvalue().encode("utf-8"))

    summary = {
        "label": "POST HOC revision analysis 2026-09-30, item I_fig5a_estimand; from stored outputs only",
        "n_arms_usable_fig5a": int((wide.usable == 1).sum()),
        "share_recomputation": {
            "n_checked": int(len(ver)), "n_within_rounding_bound": int(ver.share_recomputation_within_rounding.sum()),
            "max_abs_deviation": float(ver.share_abs_deviation.max()),
            "max_deviation_row": ver.loc[ver.share_abs_deviation.idxmax(), ["dataset_id", "negative_set"]].tolist()},
        "source_data_match": {"n_checked": int(len(sdc)), "n_match": int(sdc.match.sum())},
        "sd7_crosscheck": {"n_rows": int(len(sd7)), "all_in_sd7": bool(sd7.in_sd7.all()),
                           "all_n_match": bool(sd7.n_match.all()), "all_usable_match": bool(sd7.usable_match.all())},
        "manuscript_results_paragraph_numbers": msc,
        "ranges": ranges(wide, t3),
        "neg_b_lower_share_global": {r.dataset_id: r.neg_b_lower_share_global for _, r in pair.iterrows() if r.neg_b_available == 1},
        "neg_b_lower_share_within_point": {r.dataset_id: r.neg_b_lower_share_within_point for _, r in pair.iterrows() if r.neg_b_available == 1},
        "share_interval_conditional_arms": ver[ver.share_interval_conditional_on_floor][["dataset_id", "negative_set", "approx_frac_replicates_below_floor_normal"]].values.tolist(),
        "snote7_neg_b_definitions": s7[["label", "comparability", "observed_class_definition"]].values.tolist(),
    }
    (C.OUT / "summary_numbers.json").write_text(json.dumps(summary, indent=1, default=float, ensure_ascii=False),
                                                encoding="utf-8")
    C.record_inputs("s01_tables", [C.SHARE_CSV, C.SOURCE_DATA_FIG5, C.MANUSCRIPT, C.SNOTE7, C.SD7_CENSUS])

    show = wide[wide.usable == 1][["cohort_label", "negative_set", "n_positive", "n_negative", "n_units_both_classes",
                                   "VIS10_auc_global", "DIG25_auc_global", "VIS10_auc_within", "DIG25_auc_within",
                                   "share_global_stored", "share_within_point", "complement_ratio_global_point"]]
    print(show.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(json.dumps(summary["share_recomputation"], default=float))
    print(json.dumps(summary["source_data_match"]))
    print(json.dumps(msc, default=str)[:600])
    print(t3[["cohort_label", "DIG25_auc_within", "DIG25_auc_within_lo", "DIG25_auc_within_hi", "n_units_both_classes",
              "retained_fraction_of_NEG_A_within_above_chance_point"]].to_string(index=False))


if __name__ == "__main__":
    main()
