# -*- coding: utf-8 -*-
"""B_overadjust_sim step 2 (POST HOC): run the simulation groups with <= 4 worker processes.

Each group writes results/B_overadjust_sim/reps_<group>.csv.gz (one row per replicate x analysed
attribute, full precision). A group whose file exists is skipped, so the run can be resumed.

Groups (rate = P(positive | detected); every seed derives from MASTER_SEED 20260930 and the task
fields through numpy SeedSequence, see sim_engine.run_replicate):
  S0          null chemistry, real detection, 500 replicates, intervals (B = 1000) for all;
              released-tool cross-check on replicates 0-99
  S1          planted chemistry 0.5/1.0/1.5 log2 OR x 7 attributes x 200 replicates, intervals;
              tool cross-check on replicates 0-49
  S1grid      planted 0..3 log2 OR (10 values) x 7 attributes x 200 replicates, point estimates
  S0perm/S1perm   detection permuted within protein (no detection-attribute link), 200 replicates
  S0rate/S1rate   rate 0.10 and 0.40 instead of 0.20, 200 replicates
  S0steep/S1steep 'steep' selection of positives on peptide detectability (gamma 2.0), intervals
  S1steepgrid planted 0..3 under the steep model, point estimates
  S0boot5000  S0 replicates 0-49 re-analysed with B = 5000 (same seeds) to check B = 1000
  artifact1   the Artifact-1 control (within-protein z and protein-stratified MH log2 OR, all
              cysteines vs positives + theoretically detectable cysteines), real and steep models,
              point estimates (S0 300 replicates; S1 0.5/1.0/1.5 x 7 attributes x 200)

Groups added in the revision after adversarial verification (post hoc):
  permgrid    detection permuted within protein (no detection-attribute link: the whole baseline is
              planted chemistry and any loss under matching is over-adjustment), planted 0..3 log2 OR
              (latent also 3.5, 4.0) x 7 attributes x 200 replicates, point estimates; gives the
              conservative 'whole baseline genuine' band for the observed claims
  steepalt    two further selection models at the propensity AUC of the original steep model
              (01e_calibrate_gamma_alt.py): steep_vis (selection on the VIS10 detection logit,
              gamma 0.25) and steep_nocomp (selection on the peptide detection model without the
              basic/acidic fractions, gamma 3); S0 200 tables with intervals (B = 1000), all
              attributes plus a1c, Artifact-1 estimators on; S1 planted 1.0 x 7 attributes x 200,
              point estimates
  art1extra   exploratory Artifact-1 analogue a1c (K/R/H at +5..+8): the stored S0art1 and S0art1steep
              tables re-generated from their own seeds (identical tables) and analysed for a1c only;
              planted a1c 0.5/1.0/1.5 x 200 under the real and steep models
  art1steep4  the original steep model at gamma 4 (stronger selection), Artifact-1 estimators on:
              S0 200 tables (all attributes plus a1c), planted 1.0 into a1 and into a1c x 200

Groups added in revision round 2 (post hoc; the verifier's proxy-chemistry model and count-defined claim):
  countclaim  two HYPOTHETICAL claims defined by the local K/R count (x1: count within +/-20 above the median;
              x2: >= 3 K/R within +/-5). Null: the stored S0 (500), S0steep (300) and S0perm (200) tables are
              re-generated from their own seeds (identical tables; 07_rerun_check.py) and analysed for x1, x2
              only. Planted 0.5 / 1.0 / 1.5 log2 OR on x1 or x2 x 200 tables under real detection (scenario S1),
              the steep model (S1steep, gamma 2) and detection permuted within protein (S1permgrid). Point
              estimates; the latent reference cells of the same scenarios are the stored ones.
  proxygrid   PROXY CHEMISTRY: positives drawn with P = expit(a + b ln2 F_KRcount20_z), i.e. genuine chemistry
              that depends on the local K/R count (+/-20), b = log2 OR per SD of the count; b = 0, 0.5, ..., 4,
              5, 6 x 200 tables under real detection (P1grid), the steep model (P1steepgrid, gamma 2) and
              detection permuted within protein (P1permgrid). Every table is analysed for all published
              attributes, the latent attribute, a1c and the two count-defined attributes. Point estimates.
              (A first launch with the Artifact-1 estimators on for all ten attributes was stopped before any
              output was written: those estimators cost ~7 s per table; they run in proxyart1 instead.)
  proxyart1   Artifact-1 restriction under proxy chemistry: the proxygrid tables of P1grid and P1steepgrid at
              b = 0, 1, 2, 3, replicates 0-99, re-generated from their own seeds (seed_scenario; identical tables)
              under the scenario names P1gridA1 / P1steepgridA1 and analysed with the Artifact-1 estimators for
              a1, a1c, x1, a3 and the latent attribute.
  proxysens   sensitivity proxies under real detection: K/R count within +/-5 (P1kr5grid) and A/K/R/V count
              within +/-10 (P1akrv10grid), b = 1, 2, 3 x 200 tables; same analysis.
  proxyci     verdicts under proxy chemistry: 200 tables WITH intervals (B = 1000) at the grid b nearest to
              the b whose mean simulated SFE-006 baseline equals the observed 1.384 (rule applied to
              reps_proxygrid.csv.gz, see proxy_ci_effect()), under real detection (P1ci) and permuted detection
              (P1permci); analysed attributes a1, a1b, a2b, a3b, a4.

Groups added in revision round 3 (post hoc; the verifier's sequence-only discriminator, i.e. the association of the
label with a local count relative to its association with the SFE-006 offsets, in the same table). Every table is
analysed for DISCRIM_ATTRS (a1, x1, x3, x4, x2, a1b, a2b, a3b, a4); point estimates only.
  discrim     stage 1.
              (1) detectability only: the stored null tables of every detection model re-generated from their own
                  seeds (identical tables; 07_rerun_check.py): S0 (500), S0steep (300), S0art1steep4 (gamma 4, 200),
                  S0steepvis (200), S0steepnc (200), S0perm (200), under the scenario names D0real ... D0perm;
              (2) M1, chemistry on the SFE-006 offsets, and (3) M2, chemistry on the K/R count within +/-20, each at
                  the EXACT strength whose mean simulated SFE-006 baseline equals the observed 1.3840 (linear
                  interpolation along the stored calibration grid; exact_sfe006_effect()), real / steep / permuted
                  detection, 200 tables each (DM1*, DM2*);
              (4) cohort size: rate 0.16 (about 1,100 positives, as in the SFE-006 cohort) under real detection for
                  no chemistry, M1 and M2 at the rate-0.20 calibration (D0r16, DM1r16, DM2r16), 200 tables each;
              (5) mixture grid, real detection: chemistry on the +/-20 count (b = 0.5, 1.0, 1.5, 2.0 per SD) plus
                  chemistry on the offsets (c = 0.25, 0.5, 0.75, 1.0, 1.25 log2 OR), 100 tables per cell (DMIXgrid).
  proxysensext  the M2' sensitivity grids (P1kr5grid, P1akrv10grid; real detection) extended to b = 4, 5, 6
              (200 tables each), because neither reached the observed SFE-006 baseline at b <= 3 in round 2.
  discrim2    stage 2, strengths fixed by rule from the stage-1 outputs before these runs:
              (6) M2' at the exact b reproducing the observed SFE-006 baseline, where the extended grid reaches it
                  (DM2kr5, DM2akrv10), 200 tables each;
              (7) mixtures at the exact c (for each b) reproducing the observed SFE-006 baseline (DMIX), 200 each.

Groups added in revision round 4 (post hoc; the round-4 verifier's detectability-only model at the observed baseline,
re-implemented in this engine, and the K/R +/-5 grid extension). Every table is analysed for ROUND4_ATTRS (the seven
published/latent attributes plus the four count-defined ones) unless stated.
  lenselgrid  DETECTABILITY ONLY, selection on tryptic-peptide length among detected cysteines (detection mode
              'lensel7': P = expit(a + gamma * -ln max(pep_len, 7)); no chemistry), gamma = 3.0, 3.2, ..., 4.2,
              200 tables each (L0grid), point estimates.
  lenselcal   the same model at the EXACT gamma whose mean simulated SFE-006 baseline equals the observed 1.3840
              (linear interpolation along the lenselgrid output; exact_lensel_gamma(), rule fixed before this run):
              300 tables at rate 0.20 (L0cal) and 300 at rate 0.16 (L0cal16, about the SFE-006 cohort's 1,095
              positives), point estimates.
  lenselci    200 tables WITH intervals (B = 1000) at that gamma (L0ci), analysed for PROXY_CI_ATTRS: the verdicts
              the control returns for a detectability-only association of the observed size.
  kr5ext      chemistry on the K/R count within +/-5 (P1kr5grid) extended to b = 7, 8 per SD, 200 tables each
              (DISCRIM_ATTRS), because the round-3 grid (to b = 6) stopped just below the observed baseline.
  kr5cal      that chemistry at the EXACT b reproducing the observed SFE-006 baseline (exact_sfe006_effect over the
              proxysens, proxysensext and kr5ext grids; rule fixed before this run): DM2kr5, 200 tables.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import sys
import time

sys.dont_write_bytecode = True
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402

import pandas as pd  # noqa: E402

N_WORKERS = 4
EFFECTS_MAIN = [0.5, 1.0, 1.5]
EFFECTS_GRID = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]
GAMMA_STEEP = 2.0
# revision after verification: gammas chosen by 01e_calibrate_gamma_alt.py (gamma_choice_alt.json)
GAMMA_STEEP_VIS = 0.25
GAMMA_STEEP_NOCOMP = 3.0
GAMMA_STEEP_STRONG = 4.0
# revision round 2 (post hoc)
PROXY_GRID = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]
PROXY_SENS_GRID = [1.0, 2.0, 3.0]
PROXY_CI_ATTRS = ["a1_SFE006_KR", "a1b_SNO021_K", "a2b_PERS009_AKRV", "a3b_SFE002_E", "a4_LATENT30"]
OBSERVED_SFE006_BASELINE = 1.3840  # Source_Data_text_phase2c_sfe006_reproduction.json (checked in proxy_ci_effect)
# revision round 3 (post hoc)
DISCRIM_ATTRS = ["a1_SFE006_KR", "x1_KRcount20_hi", "x3_AKRV10_hi", "x4_KR5_hi", "x2_KR5_ge3", "a1b_SNO021_K",
                 "a2b_PERS009_AKRV", "a3b_SFE002_E", "a4_LATENT30"]
MIX_B = [0.5, 1.0, 1.5, 2.0]
MIX_C = [0.25, 0.5, 0.75, 1.0, 1.25]
PROXY_SENS_EXT = [4.0, 5.0, 6.0]
STEEP2 = dict(detection_mode="steep", gamma=2.0)
PERM = dict(detection_mode="permuted_within_protein")
# revision round 4 (post hoc)
LENSEL_GRID = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2]
KR5_EXT = [7.0, 8.0]
ROUND4_ATTRS = ["a1_SFE006_KR", "a1b_SNO021_K", "a2_SNO016_KRH", "a2b_PERS009_AKRV", "a3_SNO006_DE3", "a3b_SFE002_E",
                "a4_LATENT30", "x1_KRcount20_hi", "x2_KR5_ge3", "x3_AKRV10_hi", "x4_KR5_hi"]


def _observed_sfe006_baseline():
    import json
    j = json.load(open(cb.SFE006_JSON, encoding="utf-8"))
    assert abs(float(j["baseline_log2_or"]) - OBSERVED_SFE006_BASELINE) < 5e-5
    return float(j["baseline_log2_or"])


def _first_crossing(x, yb, target):
    """Linear interpolation at the first grid interval where the mean baseline crosses `target`; None if never."""
    for i in range(len(x) - 1):
        if (yb[i] - target) * (yb[i + 1] - target) <= 0 and yb[i + 1] != yb[i]:
            return float(x[i] + (target - yb[i]) * (x[i + 1] - x[i]) / (yb[i + 1] - yb[i]))
    return None


def exact_sfe006_effect(files, scenario, planted, allow_missing=False):
    """Revision round 3. EXACT planted strength (not the nearest grid value) at which the mean simulated baseline of
    the SFE-006 attribute equals the observed SFE-006 baseline, by linear interpolation along the stored grid of
    `scenario` (rows with planted_attribute == `planted`). Rule fixed before the round-3 runs."""
    target = _observed_sfe006_baseline()
    parts = []
    for f in files:
        p = os.path.join(cb.RESULTS, f)
        if os.path.exists(p):
            r = pd.read_csv(p, usecols=["scenario", "planted_attribute", "attribute", "effect_log2_or", "baseline"])
            parts.append(r[(r.scenario == scenario) & (r.attribute == "a1_SFE006_KR")
                           & (r.planted_attribute.fillna("") == planted)])
    r = pd.concat(parts)
    q = r.groupby("effect_log2_or")["baseline"].mean().sort_index()
    e = _first_crossing(q.index.to_numpy(), q.to_numpy(), target)
    if e is None and not allow_missing:
        raise RuntimeError("observed SFE-006 baseline not reached on the %s grid" % scenario)
    return e


def exact_mix_c(b):
    """Revision round 3. For chemistry on the +/-20 count at b (per SD), the strength c of added chemistry on the
    SFE-006 offsets at which the mean simulated SFE-006 baseline equals the observed one: linear interpolation along
    c = 0 (the stored P1grid cell at b) and the DMIXgrid cells (c = MIX_C). Rule fixed before the discrim2 runs."""
    import numpy as np
    target = _observed_sfe006_baseline()
    pg = pd.read_csv(os.path.join(cb.RESULTS, "reps_proxygrid.csv.gz"),
                     usecols=["scenario", "attribute", "effect_log2_or", "baseline"])
    c0 = float(pg[(pg.scenario == "P1grid") & (pg.attribute == "a1_SFE006_KR")
                  & np.isclose(pg.effect_log2_or, b)]["baseline"].mean())
    d = pd.read_csv(os.path.join(cb.RESULTS, "reps_discrim.csv.gz"),
                    usecols=["scenario", "attribute", "effect_log2_or", "effect2_log2_or", "baseline"])
    d = d[(d.scenario == "DMIXgrid") & (d.attribute == "a1_SFE006_KR") & np.isclose(d.effect_log2_or, b)]
    q = d.groupby("effect2_log2_or")["baseline"].mean().sort_index()
    x = np.concatenate([[0.0], q.index.to_numpy()])
    yb = np.concatenate([[c0], q.to_numpy()])
    c = _first_crossing(x, yb, target)
    if c is None:
        raise RuntimeError("observed SFE-006 baseline not reached on the mixture grid at b = %g" % b)
    return c


def exact_lensel_gamma():
    """Revision round 4. EXACT gamma of the length-selection model (detection mode 'lensel7', no chemistry) at which
    the mean simulated SFE-006 baseline equals the observed SFE-006 baseline, by linear interpolation along the
    stored lenselgrid output (scenario L0grid). Rule fixed before the lenselcal and lenselci runs."""
    target = _observed_sfe006_baseline()
    r = pd.read_csv(os.path.join(cb.RESULTS, "reps_lenselgrid.csv.gz"),
                    usecols=["scenario", "attribute", "gamma", "baseline"])
    q = (r[(r.scenario == "L0grid") & (r.attribute == "a1_SFE006_KR")].groupby("gamma")["baseline"].mean()
         .sort_index())
    g = _first_crossing(q.index.to_numpy(), q.to_numpy(), target)
    if g is None:
        raise RuntimeError("observed SFE-006 baseline not reached on the length-selection grid")
    return g


def proxy_ci_effect(scenario):
    """Grid b nearest to the b at which the mean simulated baseline of the SFE-006 attribute equals the observed
    SFE-006 baseline (linear interpolation along the proxy grid of `scenario`). Rule fixed before the interval
    runs; applied to the stored reps_proxygrid.csv.gz."""
    import json
    import numpy as np
    j = json.load(open(cb.SFE006_JSON, encoding="utf-8"))
    assert abs(float(j["baseline_log2_or"]) - OBSERVED_SFE006_BASELINE) < 5e-5
    r = pd.read_csv(os.path.join(cb.RESULTS, "reps_proxygrid.csv.gz"))
    q = (r[(r.scenario == scenario) & (r.attribute == "a1_SFE006_KR")].groupby("effect_log2_or")["baseline"]
         .mean().sort_index())
    x, yb = q.index.to_numpy(), q.to_numpy()
    target = float(j["baseline_log2_or"])
    for i in range(len(x) - 1):
        if (yb[i] - target) * (yb[i + 1] - target) <= 0 and yb[i + 1] != yb[i]:
            e = x[i] + (target - yb[i]) * (x[i + 1] - x[i]) / (yb[i + 1] - yb[i])
            return float(x[np.argmin(np.abs(x - e))]), float(e)
    raise RuntimeError("observed SFE-006 baseline not reached on the %s grid" % scenario)


def tasks_for(group):
    import sim_engine as se
    A = se.ALL_ATTRIBUTES
    T = []
    if group == "S0":
        for r in range(500):
            T.append(dict(scenario="S0", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=True,
                          run_tool=r < 100))
    elif group == "S1":
        for a in A:
            for e in EFFECTS_MAIN:
                for r in range(200):
                    T.append(dict(scenario="S1", attribute=a, effect=e, rate=0.2, rep=r, want_ci=True,
                                  run_tool=r < 50, analyse="planted_only"))
    elif group == "S1grid":
        for a in A:
            for e in EFFECTS_GRID:
                for r in range(200):
                    T.append(dict(scenario="S1grid", attribute=a, effect=e, rate=0.2, rep=r, want_ci=False,
                                  analyse="planted_only"))
    elif group == "perm":
        for r in range(200):
            T.append(dict(scenario="S0perm", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=False,
                          detection_mode="permuted_within_protein"))
        for a in A:
            for r in range(200):
                T.append(dict(scenario="S1perm", attribute=a, effect=1.0, rate=0.2, rep=r, want_ci=False,
                              analyse="planted_only", detection_mode="permuted_within_protein"))
    elif group == "rate":
        for rate in (0.10, 0.40):
            for r in range(200):
                T.append(dict(scenario="S0rate", attribute=None, effect=0.0, rate=rate, rep=r, want_ci=False))
            for a in A:
                for r in range(200):
                    T.append(dict(scenario="S1rate", attribute=a, effect=1.0, rate=rate, rep=r, want_ci=False,
                                  analyse="planted_only"))
    elif group == "steep":
        for r in range(300):
            T.append(dict(scenario="S0steep", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=True,
                          detection_mode="steep", gamma=GAMMA_STEEP))
        for a in A:
            for e in EFFECTS_MAIN:
                for r in range(100):
                    T.append(dict(scenario="S1steep", attribute=a, effect=e, rate=0.2, rep=r, want_ci=True,
                                  analyse="planted_only", detection_mode="steep", gamma=GAMMA_STEEP))
    elif group == "steepgrid":
        for a in A:
            for e in EFFECTS_GRID:
                for r in range(200):
                    T.append(dict(scenario="S1steepgrid", attribute=a, effect=e, rate=0.2, rep=r, want_ci=False,
                                  analyse="planted_only", detection_mode="steep", gamma=GAMMA_STEEP))
    elif group == "artifact1":
        # Artifact-1 control (background restricted to theoretically detectable cysteines) on the
        # same kind of tables; point estimates only
        for r in range(300):
            T.append(dict(scenario="S0art1", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=False,
                          artifact1=True))
            T.append(dict(scenario="S0art1steep", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=False,
                          artifact1=True, detection_mode="steep", gamma=GAMMA_STEEP))
        for a in A:
            for e in EFFECTS_MAIN:
                for r in range(200):
                    T.append(dict(scenario="S1art1", attribute=a, effect=e, rate=0.2, rep=r, want_ci=False,
                                  analyse="planted_only", artifact1=True))
                    T.append(dict(scenario="S1art1steep", attribute=a, effect=e, rate=0.2, rep=r,
                                  want_ci=False, analyse="planted_only", artifact1=True,
                                  detection_mode="steep", gamma=GAMMA_STEEP))
    elif group == "latentext":
        # extends the calibration grids of the latent attribute to 3.5 and 4.0 log2 OR so that the
        # structural claim SFE-001 (observed baseline 2.57) falls inside the simulated range
        for scen, mode in (("S1grid", "real"), ("S1steepgrid", "steep")):
            for e in (3.5, 4.0):
                for r in range(200):
                    t = dict(scenario=scen, attribute="a4_LATENT30", effect=e, rate=0.2, rep=r, want_ci=False,
                             analyse="planted_only")
                    if mode == "steep":
                        t.update(detection_mode="steep", gamma=GAMMA_STEEP)
                    T.append(t)
    elif group == "boot5000":
        for r in range(50):
            T.append(dict(scenario="S0boot5000", seed_scenario="S0", attribute=None, effect=0.0, rate=0.2,
                          rep=r, want_ci=True, n_boot=5000))
    # ---------------- revision after adversarial verification (post hoc) ----------------
    elif group == "permgrid":
        for a in A:
            effects = EFFECTS_GRID + ([3.5, 4.0] if a == "a4_LATENT30" else [])
            for e in effects:
                for r in range(200):
                    T.append(dict(scenario="S1permgrid", attribute=a, effect=e, rate=0.2, rep=r, want_ci=False,
                                  analyse="planted_only", detection_mode="permuted_within_protein"))
    elif group == "steepalt":
        for mode, s0, s1, gam in (("steep_vis", "S0steepvis", "S1steepvis", GAMMA_STEEP_VIS),
                                  ("steep_nocomp", "S0steepnc", "S1steepnc", GAMMA_STEEP_NOCOMP)):
            for r in range(200):
                T.append(dict(scenario=s0, attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=True,
                              detection_mode=mode, gamma=gam, analyse="with_extra", artifact1=True))
            for a in A:
                for r in range(200):
                    T.append(dict(scenario=s1, attribute=a, effect=1.0, rate=0.2, rep=r, want_ci=False,
                                  analyse="planted_only", detection_mode=mode, gamma=gam, artifact1=True))
    elif group == "art1extra":
        for r in range(300):
            T.append(dict(scenario="S0art1", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=False,
                          artifact1=True, analyse="extra_only"))
            T.append(dict(scenario="S0art1steep", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=False,
                          artifact1=True, analyse="extra_only", detection_mode="steep", gamma=GAMMA_STEEP))
        for e in EFFECTS_MAIN:
            for r in range(200):
                T.append(dict(scenario="S1art1", attribute="a1c_KRH_p5p8", effect=e, rate=0.2, rep=r,
                              want_ci=False, analyse="planted_only", artifact1=True))
                T.append(dict(scenario="S1art1steep", attribute="a1c_KRH_p5p8", effect=e, rate=0.2, rep=r,
                              want_ci=False, analyse="planted_only", artifact1=True,
                              detection_mode="steep", gamma=GAMMA_STEEP))
    elif group == "art1steep4":
        for r in range(200):
            T.append(dict(scenario="S0art1steep4", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=False,
                          artifact1=True, analyse="with_extra", detection_mode="steep", gamma=GAMMA_STEEP_STRONG))
        for a in ("a1_SFE006_KR", "a1c_KRH_p5p8"):
            for r in range(200):
                T.append(dict(scenario="S1art1steep4", attribute=a, effect=1.0, rate=0.2, rep=r, want_ci=False,
                              analyse="planted_only", artifact1=True, detection_mode="steep",
                              gamma=GAMMA_STEEP_STRONG))
    # ---------------- revision round 2 (post hoc) ----------------
    elif group == "countclaim":
        # null: stored tables re-generated from their own seeds, analysed for the count-defined attributes only
        for scen, n, extra in (("S0", 500, {}), ("S0steep", 300, dict(detection_mode="steep", gamma=GAMMA_STEEP)),
                               ("S0perm", 200, dict(detection_mode="permuted_within_protein"))):
            for r in range(n):
                t = dict(scenario=scen, attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=False,
                         analyse="count_only")
                t.update(extra)
                T.append(t)
        for scen, extra in (("S1", {}), ("S1steep", dict(detection_mode="steep", gamma=GAMMA_STEEP)),
                            ("S1permgrid", dict(detection_mode="permuted_within_protein"))):
            for a in ("x1_KRcount20_hi", "x2_KR5_ge3"):
                for e in EFFECTS_MAIN:
                    for r in range(200):
                        t = dict(scenario=scen, attribute=a, effect=e, rate=0.2, rep=r, want_ci=False,
                                 analyse="planted_only")
                        t.update(extra)
                        T.append(t)
    elif group == "proxygrid":
        for scen, extra in (("P1grid", {}), ("P1steepgrid", dict(detection_mode="steep", gamma=GAMMA_STEEP)),
                            ("P1permgrid", dict(detection_mode="permuted_within_protein"))):
            for e in PROXY_GRID:
                for r in range(200):
                    t = dict(scenario=scen, attribute="F_KRcount20_z", effect=e, rate=0.2, rep=r, want_ci=False,
                             analyse="with_count")
                    t.update(extra)
                    T.append(t)
    elif group == "proxyart1":
        for scen, seed_scen, extra in (("P1gridA1", "P1grid", {}),
                                       ("P1steepgridA1", "P1steepgrid", dict(detection_mode="steep", gamma=GAMMA_STEEP))):
            for e in (0.0, 1.0, 2.0, 3.0):
                for r in range(100):
                    t = dict(scenario=scen, seed_scenario=seed_scen, attribute="F_KRcount20_z", effect=e, rate=0.2,
                             rep=r, want_ci=False, artifact1=True,
                             attrs=["a1_SFE006_KR", "a1c_KRH_p5p8", "x1_KRcount20_hi", "a3_SNO006_DE3", "a4_LATENT30"])
                    t.update(extra)
                    T.append(t)
    elif group == "proxysens":
        for scen, feat in (("P1kr5grid", "F_KR5_z"), ("P1akrv10grid", "F_AKRV10_z")):
            for e in PROXY_SENS_GRID:
                for r in range(200):
                    T.append(dict(scenario=scen, attribute=feat, effect=e, rate=0.2, rep=r, want_ci=False,
                                  analyse="with_count"))
    elif group == "proxyci":
        for scen, grid_scen, extra in (("P1ci", "P1grid", {}),
                                       ("P1permci", "P1permgrid", dict(detection_mode="permuted_within_protein"))):
            e, e_exact = proxy_ci_effect(grid_scen)
            print("  proxyci %s: b reproducing the observed SFE-006 baseline %.3f -> nearest grid %.2f"
                  % (scen, e_exact, e), flush=True)
            for r in range(200):
                t = dict(scenario=scen, attribute="F_KRcount20_z", effect=e, rate=0.2, rep=r, want_ci=True,
                         attrs=list(PROXY_CI_ATTRS))
                t.update(extra)
                T.append(t)
    # ---------------- revision round 3 (post hoc): sequence-only discriminator ----------------
    elif group == "discrim":
        DA = list(DISCRIM_ATTRS)
        # (1) detectability only: stored null tables re-generated from their own seeds
        for scen, seed_scen, n, extra in (
                ("D0real", "S0", 500, {}),
                ("D0steep", "S0steep", 300, STEEP2),
                ("D0steep4", "S0art1steep4", 200, dict(detection_mode="steep", gamma=GAMMA_STEEP_STRONG)),
                ("D0steepvis", "S0steepvis", 200, dict(detection_mode="steep_vis", gamma=GAMMA_STEEP_VIS)),
                ("D0steepnc", "S0steepnc", 200, dict(detection_mode="steep_nocomp", gamma=GAMMA_STEEP_NOCOMP)),
                ("D0perm", "S0perm", 200, PERM)):
            for r in range(n):
                t = dict(scenario=scen, seed_scenario=seed_scen, attribute=None, effect=0.0, rate=0.2, rep=r,
                         want_ci=False, attrs=DA)
                t.update(extra)
                T.append(t)
        # (2) M1 and (3) M2 at the exact strength reproducing the observed SFE-006 baseline
        exact = {}
        for scen, fname, grid_scen, planted, extra in (
                ("DM1real", "reps_S1grid.csv.gz", "S1grid", "a1_SFE006_KR", {}),
                ("DM1steep", "reps_steepgrid.csv.gz", "S1steepgrid", "a1_SFE006_KR", STEEP2),
                ("DM1perm", "reps_permgrid.csv.gz", "S1permgrid", "a1_SFE006_KR", PERM),
                ("DM2real", "reps_proxygrid.csv.gz", "P1grid", "F_KRcount20_z", {}),
                ("DM2steep", "reps_proxygrid.csv.gz", "P1steepgrid", "F_KRcount20_z", STEEP2),
                ("DM2perm", "reps_proxygrid.csv.gz", "P1permgrid", "F_KRcount20_z", PERM)):
            e = exact_sfe006_effect([fname], grid_scen, planted)
            exact[scen] = e
            print("  discrim %s: strength reproducing the observed SFE-006 baseline on %s = %.4f"
                  % (scen, grid_scen, e), flush=True)
            for r in range(200):
                t = dict(scenario=scen, attribute=planted, effect=e, rate=0.2, rep=r, want_ci=False, attrs=DA)
                t.update(extra)
                T.append(t)
        # (4) cohort size: rate 0.16, real detection, at the rate-0.20 calibration
        for scen, planted, e in (("D0r16", None, 0.0), ("DM1r16", "a1_SFE006_KR", exact["DM1real"]),
                                 ("DM2r16", "F_KRcount20_z", exact["DM2real"])):
            for r in range(200):
                T.append(dict(scenario=scen, attribute=planted, effect=e, rate=0.16, rep=r, want_ci=False, attrs=DA))
        # (5) mixture grid, real detection
        for b in MIX_B:
            for c in MIX_C:
                for r in range(100):
                    T.append(dict(scenario="DMIXgrid", attribute="F_KRcount20_z", effect=b,
                                  attribute2="a1_SFE006_KR", effect2=c, rate=0.2, rep=r, want_ci=False, attrs=DA))
    elif group == "proxysensext":
        for scen, feat in (("P1kr5grid", "F_KR5_z"), ("P1akrv10grid", "F_AKRV10_z")):
            for e in PROXY_SENS_EXT:
                for r in range(200):
                    T.append(dict(scenario=scen, attribute=feat, effect=e, rate=0.2, rep=r, want_ci=False,
                                  attrs=list(DISCRIM_ATTRS)))
    elif group == "discrim2":
        DA = list(DISCRIM_ATTRS)
        # (6) M2' at the exact b reproducing the observed SFE-006 baseline, where the extended grid reaches it
        for scen, grid_scen, feat in (("DM2kr5", "P1kr5grid", "F_KR5_z"), ("DM2akrv10", "P1akrv10grid", "F_AKRV10_z")):
            e = exact_sfe006_effect(["reps_proxysens.csv.gz", "reps_proxysensext.csv.gz"], grid_scen, feat,
                                    allow_missing=True)
            if e is None:
                print("  discrim2 %s: observed SFE-006 baseline not reached for b <= 6; no cell" % scen, flush=True)
                continue
            print("  discrim2 %s: b reproducing the observed SFE-006 baseline = %.4f" % (scen, e), flush=True)
            for r in range(200):
                T.append(dict(scenario=scen, attribute=feat, effect=e, rate=0.2, rep=r, want_ci=False, attrs=DA))
        # (7) mixtures at the exact c reproducing the observed SFE-006 baseline
        for b in MIX_B:
            c = exact_mix_c(b)
            print("  discrim2 DMIX b=%.2f: c reproducing the observed SFE-006 baseline = %.4f" % (b, c), flush=True)
            for r in range(200):
                T.append(dict(scenario="DMIX", attribute="F_KRcount20_z", effect=b, attribute2="a1_SFE006_KR",
                              effect2=c, rate=0.2, rep=r, want_ci=False, attrs=DA))
    # ---------------- revision round 4 (post hoc): detectability-only model at the observed baseline ----------------
    elif group == "lenselgrid":
        for g in LENSEL_GRID:
            for r in range(200):
                T.append(dict(scenario="L0grid", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=False,
                              attrs=list(ROUND4_ATTRS), detection_mode="lensel7", gamma=g))
    elif group == "lenselcal":
        g = exact_lensel_gamma()
        print("  lenselcal: gamma reproducing the observed SFE-006 baseline on L0grid = %.4f" % g, flush=True)
        for scen, rate in (("L0cal", 0.2), ("L0cal16", 0.16)):
            for r in range(300):
                T.append(dict(scenario=scen, attribute=None, effect=0.0, rate=rate, rep=r, want_ci=False,
                              attrs=list(ROUND4_ATTRS), detection_mode="lensel7", gamma=g))
    elif group == "lenselci":
        g = exact_lensel_gamma()
        print("  lenselci: gamma reproducing the observed SFE-006 baseline on L0grid = %.4f" % g, flush=True)
        for r in range(200):
            T.append(dict(scenario="L0ci", attribute=None, effect=0.0, rate=0.2, rep=r, want_ci=True,
                          attrs=list(PROXY_CI_ATTRS), detection_mode="lensel7", gamma=g))
    elif group == "kr5ext":
        for e in KR5_EXT:
            for r in range(200):
                T.append(dict(scenario="P1kr5grid", attribute="F_KR5_z", effect=e, rate=0.2, rep=r, want_ci=False,
                              attrs=list(DISCRIM_ATTRS)))
    elif group == "kr5cal":
        e = exact_sfe006_effect(["reps_proxysens.csv.gz", "reps_proxysensext.csv.gz", "reps_kr5ext.csv.gz"],
                                "P1kr5grid", "F_KR5_z")
        print("  kr5cal DM2kr5: b reproducing the observed SFE-006 baseline = %.4f" % e, flush=True)
        for r in range(200):
            T.append(dict(scenario="DM2kr5", attribute="F_KR5_z", effect=e, rate=0.2, rep=r, want_ci=False,
                          attrs=list(DISCRIM_ATTRS)))
    else:
        raise ValueError(group)
    return T


def _work(task):
    import sim_engine as se
    return se.run_replicate(task)


def run_group(group, pool):
    out = os.path.join(cb.RESULTS, "reps_%s.csv.gz" % group)
    if os.path.exists(out):
        print("skip (exists)", out, flush=True)
        return
    tasks = tasks_for(group)
    t0 = time.time()
    rows = []
    heavy = any(t["want_ci"] for t in tasks)
    chunk = 1 if heavy else 20
    for i, res in enumerate(pool.imap_unordered(_work, tasks, chunksize=chunk), 1):
        rows.extend(res)
        if i % max(1, len(tasks) // 20) == 0:
            print("  %s %d/%d  %.0fs" % (group, i, len(tasks), time.time() - t0), flush=True)
    df = pd.DataFrame(rows)
    keys = ["scenario", "detection_mode", "planted_attribute", "attribute", "effect_log2_or", "rate", "rep"]
    if "effect2_log2_or" in df.columns:                # revision round 3: mixture rows (deterministic order)
        keys.insert(5, "effect2_log2_or")
    if df["gamma"].nunique() > 1:                      # revision round 4: gamma grid rows (deterministic order)
        keys.insert(keys.index("rate"), "gamma")
    df = df.sort_values(keys).reset_index(drop=True)
    df.to_csv(out + ".tmp", index=False, compression="gzip")
    os.replace(out + ".tmp", out)
    print("wrote %s rows=%d in %.0fs" % (out, len(df), time.time() - t0), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--groups", default="S0,S1,S1grid,perm,rate,steep,steepgrid,boot5000")
    args = ap.parse_args()
    ctx = mp.get_context("spawn")
    with ctx.Pool(N_WORKERS) as pool:
        for g in args.groups.split(","):
            print("[group]", g, flush=True)
            run_group(g, pool)


if __name__ == "__main__":
    main()
