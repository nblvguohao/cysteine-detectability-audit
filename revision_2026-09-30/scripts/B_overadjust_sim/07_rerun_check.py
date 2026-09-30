# -*- coding: utf-8 -*-
"""B_overadjust_sim step 7 (POST HOC): determinism check. sim_engine.py gained the Artifact-1 estimators after
the main groups had run, and (revision after adversarial verification) two selection models, one exploratory
attribute and new scenario codes. This re-runs a sample of stored tables from every group with the FINAL engine
and compares every shared column with the stored rows (tolerance 1e-9; stored CSVs carry 16-17 significant
digits). It also checks that the art1extra null tables, re-generated from the S0art1 / S0art1steep seeds, are the
same tables as the stored ones (identical table-level columns for all 600 tables).
Revision round 2: the round-2 groups (countclaim, proxygrid, proxysens, proxyci) are added to the re-run sample,
and the countclaim null tables, re-generated from the S0 / S0steep / S0perm seeds, are compared with the stored
tables of those scenarios (table-level columns, all 1,000 tables).
Revision round 3: the round-3 groups (discrim, discrim2, proxysensext) are added to the re-run sample, and the six
detectability-only cells of the discrim group, re-generated from the S0 / S0steep / S0art1steep4 / S0steepvis /
S0steepnc / S0perm seeds, are compared with the stored tables of those scenarios (table-level columns, all 1,600).
Revision round 4: sim_engine.py gained the length-selection detection mode ('lensel7'); the whole sample above is
re-run with that engine (backward compatibility), and the round-4 groups (lenselgrid, lenselcal, lenselci, kr5ext,
kr5cal) are added to the sample, with the calibrated gamma / b read back from the stored rows.
Writes results/B_overadjust_sim/rerun_check.csv, rerun_check_art1extra_tables.csv,
rerun_check_countclaim_tables.csv and rerun_check_discrim_tables.csv."""
import os
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import sim_engine as se  # noqa: E402

# revision round 3: the attribute list of the discriminator groups (02_run_simulation.DISCRIM_ATTRS)
DISCRIM_ATTRS = ["a1_SFE006_KR", "x1_KRcount20_hi", "x3_AKRV10_hi", "x4_KR5_hi", "x2_KR5_ge3", "a1b_SNO021_K",
                 "a2b_PERS009_AKRV", "a3b_SFE002_E", "a4_LATENT30"]

checks = [
    ("reps_S0.csv.gz", [dict(scenario="S0", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=True, run_tool=True) for r in (0, 1)]),
    ("reps_S1.csv.gz", [dict(scenario="S1", attribute="a1_SFE006_KR", effect=1.0, rate=0.2, rep=r, want_ci=True,
                             run_tool=True, analyse="planted_only") for r in (0, 1)]),
    ("reps_S1grid.csv.gz", [dict(scenario="S1grid", attribute="a2b_PERS009_AKRV", effect=0.75, rate=0.2, rep=r,
                                 want_ci=False, analyse="planted_only") for r in (0, 1, 2)]),
    ("reps_steep.csv.gz", [dict(scenario="S1steep", attribute="a3_SNO006_DE3", effect=0.5, rate=0.2, rep=r, want_ci=True,
                                analyse="planted_only", detection_mode="steep", gamma=2.0) for r in (0,)]),
    # revision: groups that ran before the engine changes
    ("reps_perm.csv.gz", [dict(scenario="S1perm", attribute="a1_SFE006_KR", effect=1.0, rate=0.2, rep=r, want_ci=False,
                               analyse="planted_only", detection_mode="permuted_within_protein") for r in (0, 1)]),
    ("reps_artifact1.csv.gz", [dict(scenario="S0art1", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=False,
                                    artifact1=True) for r in (0,)]
     + [dict(scenario="S1art1steep", attribute="a1_SFE006_KR", effect=1.0, rate=0.2, rep=r, want_ci=False,
             analyse="planted_only", artifact1=True, detection_mode="steep", gamma=2.0) for r in (0,)]),
    # revision: new groups
    ("reps_permgrid.csv.gz", [dict(scenario="S1permgrid", attribute="a4_LATENT30", effect=3.5, rate=0.2, rep=r,
                                   want_ci=False, analyse="planted_only", detection_mode="permuted_within_protein")
                              for r in (0, 1)]),
    ("reps_steepalt.csv.gz", [dict(scenario="S1steepnc", attribute="a3b_SFE002_E", effect=1.0, rate=0.2, rep=r,
                                   want_ci=False, analyse="planted_only", detection_mode="steep_nocomp", gamma=3.0,
                                   artifact1=True) for r in (0,)]
     + [dict(scenario="S0steepvis", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=True,
             detection_mode="steep_vis", gamma=0.25, analyse="with_extra", artifact1=True) for r in (0,)]),
    ("reps_art1steep4.csv.gz", [dict(scenario="S1art1steep4", attribute="a1c_KRH_p5p8", effect=1.0, rate=0.2, rep=r,
                                     want_ci=False, analyse="planted_only", artifact1=True, detection_mode="steep",
                                     gamma=4.0) for r in (0,)]),
    # revision round 2
    ("reps_countclaim.csv.gz", [dict(scenario="S1permgrid", attribute="x1_KRcount20_hi", effect=1.0, rate=0.2, rep=r,
                                     want_ci=False, analyse="planted_only", detection_mode="permuted_within_protein")
                                for r in (0,)]
     + [dict(scenario="S1steep", attribute="x2_KR5_ge3", effect=0.5, rate=0.2, rep=r, want_ci=False,
             analyse="planted_only", detection_mode="steep", gamma=2.0) for r in (1,)]),
    ("reps_proxygrid.csv.gz", [dict(scenario="P1grid", attribute="F_KRcount20_z", effect=3.0, rate=0.2, rep=r,
                                    want_ci=False, analyse="with_count") for r in (0,)]
     + [dict(scenario="P1permgrid", attribute="F_KRcount20_z", effect=2.0, rate=0.2, rep=r, want_ci=False,
             analyse="with_count", detection_mode="permuted_within_protein") for r in (5,)]),
    ("reps_proxyart1.csv.gz", [dict(scenario="P1steepgridA1", seed_scenario="P1steepgrid", attribute="F_KRcount20_z",
                                    effect=2.0, rate=0.2, rep=r, want_ci=False, artifact1=True,
                                    detection_mode="steep", gamma=2.0,
                                    attrs=["a1_SFE006_KR", "a1c_KRH_p5p8", "x1_KRcount20_hi", "a3_SNO006_DE3",
                                           "a4_LATENT30"]) for r in (3,)]),
    ("reps_proxysens.csv.gz", [dict(scenario="P1akrv10grid", attribute="F_AKRV10_z", effect=2.0, rate=0.2, rep=r,
                                    want_ci=False, analyse="with_count") for r in (0,)]),
    # revision round 3
    ("reps_proxysensext.csv.gz", [dict(scenario="P1kr5grid", attribute="F_KR5_z", effect=5.0, rate=0.2, rep=r,
                                       want_ci=False, attrs=DISCRIM_ATTRS) for r in (0,)]),
    ("reps_discrim.csv.gz", [dict(scenario="DMIXgrid", attribute="F_KRcount20_z", effect=1.5,
                                  attribute2="a1_SFE006_KR", effect2=0.75, rate=0.2, rep=r, want_ci=False,
                                  attrs=DISCRIM_ATTRS) for r in (2,)]),
]


def _calibrated_tasks():
    """Cells whose planted strength is the full-precision calibrated value: read that value from the stored rows
    rather than re-deriving it, so the check compares the same cell (the seed uses round(100 * effect), so a
    rounded value would give the same table but a different recorded effect column)."""
    out = []
    for fname, scen, planted, rate in (("reps_discrim.csv.gz", "DM2real", "F_KRcount20_z", 0.2),
                                       ("reps_discrim.csv.gz", "DM1r16", "a1_SFE006_KR", 0.16),
                                       ("reps_discrim2.csv.gz", "DM2akrv10", "F_AKRV10_z", 0.2)):
        path = os.path.join(cb.RESULTS, fname)
        if not os.path.exists(path):
            continue
        st = pd.read_csv(path, low_memory=False)
        rows = st[st.scenario == scen]
        if not len(rows):
            continue
        e = float(rows["effect_log2_or"].iloc[0])
        out.append((fname, [dict(scenario=scen, attribute=planted, effect=e, rate=rate, rep=0, want_ci=False,
                                 attrs=DISCRIM_ATTRS)]))
    return out


def _proxyci_task():
    """One interval table of the proxyci group, with the effect it was run at (read from the stored rows)."""
    path = os.path.join(cb.RESULTS, "reps_proxyci.csv.gz")
    if not os.path.exists(path):
        return []
    st = pd.read_csv(path)
    e = float(st[st.scenario == "P1ci"]["effect_log2_or"].iloc[0])
    return [("reps_proxyci.csv.gz", [dict(scenario="P1ci", attribute="F_KRcount20_z", effect=e, rate=0.2, rep=0,
                                          want_ci=True, attrs=["a1_SFE006_KR", "a1b_SNO021_K", "a2b_PERS009_AKRV",
                                                               "a3b_SFE002_E", "a4_LATENT30"])])]


ROUND4_ATTRS = ["a1_SFE006_KR", "a1b_SNO021_K", "a2_SNO016_KRH", "a2b_PERS009_AKRV", "a3_SNO006_DE3", "a3b_SFE002_E",
                "a4_LATENT30", "x1_KRcount20_hi", "x2_KR5_ge3", "x3_AKRV10_hi", "x4_KR5_hi"]


def _round4_tasks():
    """Revision round 4: one or two tables from each round-4 group; calibrated gamma / b read from the stored rows."""
    out = []
    out.append(("reps_lenselgrid.csv.gz", [dict(scenario="L0grid", attribute=None, effect=0.0, rate=0.2, rep=r,
                                                want_ci=False, attrs=ROUND4_ATTRS, detection_mode="lensel7", gamma=3.4)
                                           for r in (0, 7)]))
    for fname, scen, rate, ci, attrs in (("reps_lenselcal.csv.gz", "L0cal", 0.2, False, ROUND4_ATTRS),
                                         ("reps_lenselcal.csv.gz", "L0cal16", 0.16, False, ROUND4_ATTRS),
                                         ("reps_lenselci.csv.gz", "L0ci", 0.2, True,
                                          ["a1_SFE006_KR", "a1b_SNO021_K", "a2b_PERS009_AKRV", "a3b_SFE002_E",
                                           "a4_LATENT30"])):
        path = os.path.join(cb.RESULTS, fname)
        if not os.path.exists(path):
            continue
        st = pd.read_csv(path, low_memory=False)
        g = float(st[st.scenario == scen]["gamma"].iloc[0])
        out.append((fname, [dict(scenario=scen, attribute=None, effect=0.0, rate=rate, rep=0, want_ci=ci, attrs=attrs,
                                 detection_mode="lensel7", gamma=g)]))
    out.append(("reps_kr5ext.csv.gz", [dict(scenario="P1kr5grid", attribute="F_KR5_z", effect=7.0, rate=0.2, rep=0,
                                            want_ci=False, attrs=DISCRIM_ATTRS)]))
    path = os.path.join(cb.RESULTS, "reps_kr5cal.csv.gz")
    if os.path.exists(path):
        st = pd.read_csv(path, low_memory=False)
        e = float(st[st.scenario == "DM2kr5"]["effect_log2_or"].iloc[0])
        out.append(("reps_kr5cal.csv.gz", [dict(scenario="DM2kr5", attribute="F_KR5_z", effect=e, rate=0.2, rep=r,
                                                want_ci=False, attrs=DISCRIM_ATTRS) for r in (0, 3)]))
    return out


checks = checks + _proxyci_task() + _calibrated_tasks() + _round4_tasks()


def compare(nr, q):
    maxdiff, mism = 0.0, []
    for c in nr.index:
        if c not in q.index:
            continue
        a, b = nr[c], q[c]
        if isinstance(a, (float, int, np.floating, np.integer)) and not isinstance(a, (bool, np.bool_)):
            if not (pd.isna(a) and pd.isna(b)):
                d = abs(float(a) - float(b))
                maxdiff = max(maxdiff, d)
                if d > 1e-9:
                    mism.append(c)
        else:
            empty_a = a is None or (isinstance(a, float) and np.isnan(a)) or a == ""
            empty_b = b is None or (isinstance(b, float) and np.isnan(b)) or b == ""
            if not (empty_a and empty_b) and str(a) != str(b):
                mism.append(c)
    return maxdiff, mism


rows = []
for fname, tasks in checks:
    path = os.path.join(cb.RESULTS, fname)
    if not os.path.exists(path):
        print("missing", fname)
        continue
    stored = pd.read_csv(path)
    for t in tasks:
        new = pd.DataFrame(se.run_replicate(t))
        for _, nr in new.iterrows():
            q = stored[(stored.scenario == nr.scenario) & (stored.attribute == nr.attribute) & (stored.rep == nr.rep)
                       & np.isclose(stored.effect_log2_or, nr.effect_log2_or)
                       & np.isclose(stored.gamma.fillna(0.0) if "gamma" in stored else 0.0, nr.gamma)]
            if "effect2_log2_or" in stored.columns:   # revision round 3: mixture rows
                q = q[np.isclose(q.effect2_log2_or.fillna(0.0), float(getattr(nr, "effect2_log2_or", 0.0) or 0.0))]
            assert len(q) == 1, (fname, nr.attribute, nr.rep, len(q))
            maxdiff, mism = compare(nr, q.iloc[0])
            rows.append({"file": fname, "scenario": nr.scenario, "attribute": nr.attribute, "rep": int(nr.rep),
                         "max_abs_numeric_diff": maxdiff, "mismatching_columns": ";".join(mism)})
out = pd.DataFrame(rows)
out.to_csv(os.path.join(cb.RESULTS, "rerun_check.csv"), index=False)
print(out.to_string())
print("all identical:", bool((out.mismatching_columns == "").all()))

# art1extra null tables == stored artifact1 null tables (table-level columns)
ex_path = os.path.join(cb.RESULTS, "reps_art1extra.csv.gz")
if os.path.exists(ex_path):
    ex = pd.read_csv(ex_path)
    st = pd.read_csv(os.path.join(cb.RESULTS, "reps_artifact1.csv.gz"))
    table_cols = ["n_obs", "n_pos", "n_proteins", "pos_rate", "n_pairs", "prop_auc_before", "prop_auc_after",
                  "matched_neg_detected_frac", "all_neg_detected_frac", "detmatch_n_pairs"]
    res = []
    for scen in ("S0art1", "S0art1steep"):
        a = ex[(ex.scenario == scen) & (ex.attribute == "a1c_KRH_p5p8")].set_index("rep")[table_cols]
        b = st[(st.scenario == scen) & (st.attribute == "a1_SFE006_KR")].set_index("rep")[table_cols]
        j = a.join(b, lsuffix="_new", rsuffix="_stored", how="inner")
        maxd = max(float((j[c + "_new"] - j[c + "_stored"]).abs().max()) for c in table_cols)
        res.append({"scenario": scen, "n_tables_compared": int(len(j)), "n_new": int(len(a)), "n_stored": int(len(b)),
                    "max_abs_diff_table_columns": maxd})
    r = pd.DataFrame(res)
    r.to_csv(os.path.join(cb.RESULTS, "rerun_check_art1extra_tables.csv"), index=False)
    print(r.to_string())

# revision round 2: countclaim null tables == stored S0 / S0steep / S0perm tables (table-level columns)
cc_path = os.path.join(cb.RESULTS, "reps_countclaim.csv.gz")
if os.path.exists(cc_path):
    cc = pd.read_csv(cc_path)
    table_cols = ["n_obs", "n_pos", "n_proteins", "pos_rate", "n_pairs", "prop_auc_before", "prop_auc_after",
                  "matched_neg_detected_frac", "all_neg_detected_frac", "detmatch_n_pairs"]
    res = []
    for scen, fname in (("S0", "reps_S0.csv.gz"), ("S0steep", "reps_steep.csv.gz"), ("S0perm", "reps_perm.csv.gz")):
        st = pd.read_csv(os.path.join(cb.RESULTS, fname))
        a = cc[(cc.scenario == scen) & (cc.attribute == "x1_KRcount20_hi")].set_index("rep")[table_cols]
        b = st[(st.scenario == scen) & (st.attribute == "a1_SFE006_KR")].set_index("rep")[table_cols]
        j = a.join(b, lsuffix="_new", rsuffix="_stored", how="inner")
        maxd = max(float((j[c + "_new"] - j[c + "_stored"]).abs().max()) for c in table_cols)
        res.append({"scenario": scen, "n_tables_compared": int(len(j)), "n_new": int(len(a)), "n_stored": int(len(b)),
                    "max_abs_diff_table_columns": maxd})
    # proxyart1 tables == proxygrid tables of the same seed scenario, effect and replicate
    pa_path, pg_path = os.path.join(cb.RESULTS, "reps_proxyart1.csv.gz"), os.path.join(cb.RESULTS, "reps_proxygrid.csv.gz")
    if os.path.exists(pa_path) and os.path.exists(pg_path):
        pa, pg = pd.read_csv(pa_path), pd.read_csv(pg_path)
        for scen, seed_scen in (("P1gridA1", "P1grid"), ("P1steepgridA1", "P1steepgrid")):
            a = pa[(pa.scenario == scen) & (pa.attribute == "a1_SFE006_KR")].set_index(["effect_log2_or", "rep"])
            b = pg[(pg.scenario == seed_scen) & (pg.attribute == "a1_SFE006_KR")].set_index(["effect_log2_or", "rep"])
            j = a[table_cols + ["baseline", "matched"]].join(b[table_cols + ["baseline", "matched"]], lsuffix="_new",
                                                              rsuffix="_stored", how="inner")
            maxd = max(float((j[c + "_new"] - j[c + "_stored"]).abs().max()) for c in table_cols + ["baseline", "matched"])
            res.append({"scenario": scen + " vs " + seed_scen, "n_tables_compared": int(len(j)), "n_new": int(len(a)),
                        "n_stored": int(len(b)), "max_abs_diff_table_columns": maxd})
    r = pd.DataFrame(res)
    r.to_csv(os.path.join(cb.RESULTS, "rerun_check_countclaim_tables.csv"), index=False)
    print(r.to_string())

# revision round 3: the detectability-only cells of the discrim group == the stored null tables of those models
dsc_path = os.path.join(cb.RESULTS, "reps_discrim.csv.gz")
if os.path.exists(dsc_path):
    dsc = pd.read_csv(dsc_path, low_memory=False)
    table_cols = ["n_obs", "n_pos", "n_proteins", "pos_rate", "n_pairs", "prop_auc_before", "prop_auc_after",
                  "matched_neg_detected_frac", "all_neg_detected_frac", "detmatch_n_pairs"]
    res = []
    for scen, stored_scen, fname in (("D0real", "S0", "reps_S0.csv.gz"),
                                     ("D0steep", "S0steep", "reps_steep.csv.gz"),
                                     ("D0steep4", "S0art1steep4", "reps_art1steep4.csv.gz"),
                                     ("D0steepvis", "S0steepvis", "reps_steepalt.csv.gz"),
                                     ("D0steepnc", "S0steepnc", "reps_steepalt.csv.gz"),
                                     ("D0perm", "S0perm", "reps_perm.csv.gz")):
        p = os.path.join(cb.RESULTS, fname)
        if not os.path.exists(p):
            continue
        st = pd.read_csv(p, low_memory=False)
        a = dsc[(dsc.scenario == scen) & (dsc.attribute == "a1_SFE006_KR")].set_index("rep")[table_cols]
        b = st[(st.scenario == stored_scen) & (st.attribute == "a1_SFE006_KR")].set_index("rep")[table_cols]
        j = a.join(b, lsuffix="_new", rsuffix="_stored", how="inner")
        maxd = max(float((j[c + "_new"] - j[c + "_stored"]).abs().max()) for c in table_cols)
        res.append({"scenario": "%s vs %s" % (scen, stored_scen), "n_tables_compared": int(len(j)),
                    "n_new": int(len(a)), "n_stored": int(len(b)), "max_abs_diff_table_columns": maxd})
    # the baselines of the SFE-006 attribute must also be identical, not only the table-level columns
    for scen, stored_scen, fname in (("D0real", "S0", "reps_S0.csv.gz"), ("D0perm", "S0perm", "reps_perm.csv.gz")):
        st = pd.read_csv(os.path.join(cb.RESULTS, fname), low_memory=False)
        a = dsc[(dsc.scenario == scen) & (dsc.attribute == "a1_SFE006_KR")].set_index("rep")[["baseline", "matched"]]
        b = st[(st.scenario == stored_scen) & (st.attribute == "a1_SFE006_KR")].set_index("rep")[["baseline", "matched"]]
        j = a.join(b, lsuffix="_new", rsuffix="_stored", how="inner")
        maxd = max(float((j[c + "_new"] - j[c + "_stored"]).abs().max()) for c in ("baseline", "matched"))
        res.append({"scenario": "%s vs %s (a1 baseline/matched)" % (scen, stored_scen),
                    "n_tables_compared": int(len(j)), "n_new": int(len(a)), "n_stored": int(len(b)),
                    "max_abs_diff_table_columns": maxd})
    r = pd.DataFrame(res)
    r.to_csv(os.path.join(cb.RESULTS, "rerun_check_discrim_tables.csv"), index=False)
    print(r.to_string())
