### Type-I error, whole cohort (baseline interval): P(95% interval excludes 0) under no association (MC SE)

| method | K=44 null_iid | K=44 null_re | K=164 null_iid | K=164 null_re | K=776 null_iid | K=776 null_re | K=1475 null_iid | K=1475 null_re |
|---|---|---|---|---|---|---|---|---|
| percentile (paper) | 0.065 (0.002) | 0.068 (0.002) | 0.055 (0.002) | 0.056 (0.002) | 0.053 (0.004) | 0.051 (0.003) | 0.054 (0.004) | 0.046 (0.003) |
| BCa | 0.053 (0.002) | 0.061 (0.002) | 0.052 (0.002) | 0.056 (0.002) | 0.053 (0.004) | 0.052 (0.004) | 0.054 (0.004) | 0.046 (0.003) |
| sandwich CR0, z | 0.059 (0.002) | 0.073 (0.002) | 0.054 (0.002) | 0.056 (0.002) | 0.052 (0.004) | 0.051 (0.003) | 0.055 (0.004) | 0.046 (0.003) |
| sandwich CR1, t(G-1) | 0.050 (0.002) | 0.061 (0.002) | 0.052 (0.002) | 0.053 (0.002) | 0.051 (0.003) | 0.051 (0.003) | 0.054 (0.004) | 0.045 (0.003) |
| sandwich Mancl-DeRouen, t(G-1) | 0.043 (0.002) | 0.047 (0.002) | 0.050 (0.002) | 0.049 (0.002) | 0.051 (0.003) | 0.050 (0.003) | 0.054 (0.004) | 0.045 (0.003) |

### Type-I error, 1:1 subsample (matched-set proxy): P(95% interval excludes 0) under no association (MC SE)

| method | K=44 null_iid | K=44 null_re | K=164 null_iid | K=164 null_re | K=776 null_iid | K=776 null_re | K=1475 null_iid | K=1475 null_re |
|---|---|---|---|---|---|---|---|---|
| percentile (paper) | 0.063 (0.002) | 0.068 (0.002) | 0.055 (0.002) | 0.055 (0.002) | 0.051 (0.003) | 0.050 (0.003) | 0.059 (0.004) | 0.050 (0.003) |
| BCa | 0.046 (0.002) | 0.054 (0.002) | 0.051 (0.002) | 0.054 (0.002) | 0.051 (0.003) | 0.052 (0.004) | 0.059 (0.004) | 0.051 (0.003) |
| sandwich CR0, z | 0.052 (0.002) | 0.062 (0.002) | 0.052 (0.002) | 0.054 (0.002) | 0.051 (0.003) | 0.050 (0.003) | 0.058 (0.004) | 0.049 (0.003) |
| sandwich CR1, t(G-1) | 0.042 (0.002) | 0.051 (0.002) | 0.049 (0.002) | 0.050 (0.002) | 0.050 (0.003) | 0.049 (0.003) | 0.058 (0.004) | 0.048 (0.003) |
| sandwich Mancl-DeRouen, t(G-1) | 0.030 (0.001) | 0.035 (0.001) | 0.045 (0.002) | 0.045 (0.002) | 0.049 (0.003) | 0.048 (0.003) | 0.057 (0.004) | 0.047 (0.003) |

### Coverage of the pooled log2 OR 0.5 (alt05) and 1.0 (alt10), whole cohort

| method | alt05 K=44 | alt05 K=164 | alt05 K=776 | alt05 K=1475 | alt10 K=44 | alt10 K=164 | alt10 K=776 | alt10 K=1475 |
|---|---|---|---|---|---|---|---|---|
| percentile (paper) | 0.933 | 0.941 | 0.947 | 0.949 | 0.933 | 0.940 | 0.949 | 0.950 |
| BCa | 0.938 | 0.940 | 0.947 | 0.950 | 0.934 | 0.939 | 0.948 | 0.950 |
| sandwich CR0, z | 0.927 | 0.940 | 0.946 | 0.950 | 0.927 | 0.939 | 0.949 | 0.949 |
| sandwich CR1, t(G-1) | 0.937 | 0.942 | 0.947 | 0.950 | 0.936 | 0.942 | 0.949 | 0.949 |
| sandwich Mancl-DeRouen, t(G-1) | 0.951 | 0.948 | 0.949 | 0.950 | 0.950 | 0.947 | 0.951 | 0.950 |

### Coverage of the pooled log2 OR 0.5 (alt05) and 1.0 (alt10), 1:1 subsample

| method | alt05 K=44 | alt05 K=164 | alt05 K=776 | alt05 K=1475 | alt10 K=44 | alt10 K=164 | alt10 K=776 | alt10 K=1475 |
|---|---|---|---|---|---|---|---|---|
| percentile (paper) | 0.932 | 0.939 | 0.949 | 0.948 | 0.933 | 0.945 | 0.952 | 0.949 |
| BCa | 0.946 | 0.942 | 0.949 | 0.948 | 0.945 | 0.946 | 0.952 | 0.948 |
| sandwich CR0, z | 0.938 | 0.942 | 0.949 | 0.948 | 0.936 | 0.946 | 0.951 | 0.948 |
| sandwich CR1, t(G-1) | 0.950 | 0.944 | 0.950 | 0.948 | 0.948 | 0.949 | 0.952 | 0.948 |
| sandwich Mancl-DeRouen, t(G-1) | 0.965 | 0.950 | 0.951 | 0.949 | 0.965 | 0.955 | 0.953 | 0.949 |

### Percentile interval by design: pooled type-I error and calibration multiplier w0 (mean [range] over the 8 null_iid/null_re cells)

| design | K=44 | K=164 | K=776 | K=1475 |
|---|---|---|---|---|
| full | type-I 0.066; w0 1.073 [1.059-1.094] | type-I 0.056; w0 1.026 [1.005-1.064] | type-I 0.052; w0 1.008 [0.972-1.030] | type-I 0.050; w0 1.001 [0.971-1.073] |
| full_mh | type-I 0.066; w0 1.072 [1.049-1.100] | type-I 0.056; w0 1.026 [1.004-1.057] | type-I 0.052; w0 1.009 [0.969-1.034] | type-I 0.049; w0 0.992 [0.953-1.047] |
| sub | type-I 0.065; w0 1.078 [1.030-1.135] | type-I 0.055; w0 1.020 [1.009-1.040] | type-I 0.050; w0 0.998 [0.946-1.047] | type-I 0.054; w0 1.018 [0.989-1.057] |
| rc | type-I 0.006; w0 0.684 [0.657-0.711] | type-I 0.005; w0 0.685 [0.667-0.709] | type-I 0.004; w0 0.683 [0.655-0.715] | type-I 0.006; w0 0.687 [0.665-0.708] |

### Between-protein confounding (no within-protein association): P(excludes 0) / coverage of the pooled target, percentile

| scenario | design | K=44 | K=164 | K=776 | K=1475 | pooled target log2 OR |
|---|---|---|---|---|---|---|
| null_re | full | 0.068 / 0.932 | 0.056 / 0.944 | 0.051 / 0.949 | 0.046 / 0.954 | 0.000 |
| null_re | sub | 0.068 / 0.932 | 0.055 / 0.945 | 0.050 / 0.950 | 0.050 / 0.950 | 0.000 |
| null_conf15 | full | 0.078 / 0.931 | 0.091 / 0.945 | 0.233 / 0.945 | 0.398 / 0.944 | 0.160, 0.168 |
| null_conf15 | sub | 0.075 / 0.930 | 0.081 / 0.948 | 0.177 / 0.951 | 0.277 / 0.949 | 0.160, 0.168 |
| null_conf50 | full | 0.167 / 0.932 | 0.465 / 0.941 | 0.978 / 0.944 | 1.000 / 0.950 | 0.532, 0.557 |
| null_conf50 | sub | 0.137 / 0.934 | 0.333 / 0.940 | 0.899 / 0.945 | 0.989 / 0.946 | 0.532, 0.557 |

### Percentile, whole cohort, per cell (few-cluster detail)

| scenario | K | pi_x | size | datasets | type-I | miss below / above | w0 [95% CI] | mean width |
|---|---|---|---|---|---|---|---|---|
| null_iid | 44 | 0.15 | nb | 4000 | 0.068 (0.004) | 0.022 / 0.045 | 1.059 [1.035, 1.107] | 2.47 |
| null_iid | 44 | 0.15 | lognormal | 4000 | 0.065 (0.004) | 0.028 / 0.038 | 1.066 [1.032, 1.097] | 2.48 |
| null_iid | 44 | 0.35 | nb | 4000 | 0.060 (0.004) | 0.028 / 0.033 | 1.060 [1.022, 1.107] | 1.75 |
| null_iid | 44 | 0.35 | lognormal | 4000 | 0.066 (0.004) | 0.031 / 0.035 | 1.074 [1.037, 1.112] | 1.77 |
| null_iid | 164 | 0.15 | nb | 3000 | 0.062 (0.004) | 0.025 / 0.038 | 1.064 [1.015, 1.109] | 1.21 |
| null_iid | 164 | 0.15 | lognormal | 3000 | 0.053 (0.004) | 0.021 / 0.032 | 1.008 [0.983, 1.054] | 1.21 |
| null_iid | 164 | 0.35 | nb | 3000 | 0.051 (0.004) | 0.024 / 0.027 | 1.005 [0.973, 1.046] | 0.89 |
| null_iid | 164 | 0.35 | lognormal | 3000 | 0.055 (0.004) | 0.029 / 0.026 | 1.019 [0.989, 1.066] | 0.89 |
| null_re | 44 | 0.15 | nb | 4000 | 0.070 (0.004) | 0.025 / 0.044 | 1.087 [1.051, 1.132] | 2.66 |
| null_re | 44 | 0.15 | lognormal | 4000 | 0.068 (0.004) | 0.026 / 0.042 | 1.094 [1.044, 1.138] | 2.73 |
| null_re | 44 | 0.35 | nb | 4000 | 0.069 (0.004) | 0.033 / 0.036 | 1.072 [1.045, 1.097] | 1.93 |
| null_re | 44 | 0.35 | lognormal | 4000 | 0.065 (0.004) | 0.027 / 0.037 | 1.069 [1.036, 1.104] | 1.99 |
| null_re | 164 | 0.15 | nb | 3000 | 0.051 (0.004) | 0.019 / 0.032 | 1.007 [0.971, 1.047] | 1.29 |
| null_re | 164 | 0.15 | lognormal | 3000 | 0.054 (0.004) | 0.022 / 0.032 | 1.026 [0.987, 1.066] | 1.32 |
| null_re | 164 | 0.35 | nb | 3000 | 0.057 (0.004) | 0.023 / 0.033 | 1.024 [0.994, 1.052] | 0.98 |
| null_re | 164 | 0.35 | lognormal | 3000 | 0.061 (0.004) | 0.027 / 0.034 | 1.055 [1.016, 1.088] | 1.02 |

### Claim-matched null simulations (2,000 datasets each)

| claim | variant | design | n_datasets | type-I (MC SE) | w0 [95% CI] |
|---|---|---|---|---|---|
| PERS-003 | re | sub | 2000 | 0.049 (0.005) | 0.994 [0.958, 1.059] |
| PERS-003 | re | rc | 2000 | 0.017 (0.003) | 0.817 [0.782, 0.863] |
| PERS-006 | iid | sub | 2000 | 0.054 (0.005) | 1.046 [0.956, 1.090] |
| PERS-006 | iid | rc | 2000 | 0.006 (0.002) | 0.602 [0.576, 0.638] |
| SFE-006 | re | sub | 2000 | 0.053 (0.005) | 1.016 [0.956, 1.079] |
| SFE-006 | re | mh | 2000 | 0.054 (0.005) | 1.012 [0.978, 1.088] |
| SFE-008 | re | sub | 2000 | 0.059 (0.005) | 1.035 [0.986, 1.082] |
| SFE-008 | re | rc | 2000 | 0.003 (0.001) | 0.657 [0.630, 0.682] |
| SFE-011 | iid | sub | 2000 | 0.050 (0.005) | 0.999 [0.968, 1.046] |
| SFE-011 | iid | rc | 2000 | 0.005 (0.002) | 0.603 [0.581, 0.638] |
| SNO-001 | nega | sub | 2000 | 0.043 (0.005) | 0.974 [0.955, 1.016] |
| SNO-001 | re | sub | 2000 | 0.059 (0.005) | 1.029 [0.995, 1.087] |
| SNO-004 | nega | full | 2000 | 0.053 (0.005) | 1.010 [0.980, 1.053] |
| SNO-004 | nega | sub | 2000 | 0.056 (0.005) | 1.018 [0.981, 1.073] |
| SNO-004 | re | full | 2000 | 0.059 (0.005) | 1.031 [0.995, 1.082] |
| SNO-004 | re | sub | 2000 | 0.048 (0.005) | 0.986 [0.949, 1.059] |
| SNO-012 | nega | sub | 2000 | 0.058 (0.005) | 1.047 [1.000, 1.103] |
| SNO-012 | nega | mh | 2000 | 0.060 (0.005) | 1.045 [1.002, 1.092] |
| SNO-012 | re | sub | 2000 | 0.062 (0.005) | 1.060 [1.010, 1.112] |
| SNO-012 | re | mh | 2000 | 0.069 (0.006) | 1.056 [1.026, 1.108] |
| SNO-014 | nega | sub | 2000 | 0.065 (0.006) | 1.105 [1.027, 1.174] |
| SNO-014 | nega | rc | 2000 | 0.006 (0.002) | 0.674 [0.641, 0.702] |
| SNO-014 | re | sub | 2000 | 0.070 (0.006) | 1.124 [1.048, 1.184] |
| SNO-014 | re | rc | 2000 | 0.007 (0.002) | 0.715 [0.681, 0.762] |

### Decisive intervals with a bound near 0 (0.5 < w* < 1.6)

| claim_id | type | interval | w* | w0 (range) | flips at w0 | could flip in range | verdict_if_toggled | z_MC | w0_source |
|---|---|---|---|---|---|---|---|---|---|
| PERS-003 | matched | 0.382 [-0.015, 0.785] | 0.961 | 0.994 (0.958-1.059) | no | yes | attenuated | 2.0 | claim-matched sub re (range over re) |
| PERS-006 | matched | 1.558 [-0.139, 4.196] | 0.918 | 1.046 (0.956-1.090) | no | no | attenuated | 3.3 | claim-matched sub iid (range over iid) |
| SFE-006 | matched | 0.176 [-0.081, 0.438] | 0.684 | 1.016 (0.956-1.079) | no | no | attenuated | 16.3 | claim-matched sub re (range over re) |
| SFE-006 | stratified | 0.182 [-0.033, 0.402] | 0.846 | 1.012 (0.978-1.088) | no | no | undecidable | 7.9 | claim-matched mh re (range over re) |
| SFE-008 | random | 0.411 [-0.070, 0.870] | 0.855 | 0.657 (0.630-0.682) | yes | yes | vanishes | 7.7 | claim-matched rc re (range over re) |
| SFE-011 | matched | 1.168 [-0.116, 2.836] | 0.909 | 0.999 (0.968-1.046) | no | no | attenuated | 4.1 | claim-matched sub iid (range over iid) |
| SFE-011 | random | 0.768 [-0.663, 2.803] | 0.536 | 0.603 (0.581-0.638) | no | no | vanishes | 19.9 | claim-matched rc iid (range over iid) |
| SNO-001 | matched | 0.363 [-0.177, 0.915] | 0.672 | 0.974 (0.955-1.087) | no | no | attenuated | 16.8 | claim-matched sub nega (range over nega,re) |
| SNO-004 | baseline | 0.339 [-0.151, 0.914] | 0.692 | 1.010 (0.980-1.082) | no | no | undecidable | 14.7 | claim-matched full nega (range over nega,re) |
| SNO-004 | matched | 0.725 [0.062, 1.424] | 1.094 | 1.018 (0.949-1.073) | no | no | null_survives | 4.8 | claim-matched sub nega (range over nega,re) |
| SNO-012 | matched | 1.748 [0.528, 3.216] | 1.433 | 1.047 (1.000-1.112) | no | no | undecidable | 20.4 | claim-matched sub nega (range over nega,re) |
| SNO-012 | stratified | 0.580 [-0.316, 1.547] | 0.647 | 1.045 (1.002-1.108) | no | no | survives | 17.6 | claim-matched mh nega (range over nega,re) |
| SNO-014 | random | 0.994 [-0.665, 2.945] | 0.599 | 0.674 (0.641-0.762) | no | no | vanishes | 19.1 | claim-matched rc nega (range over nega,re) |

### Verdicts re-derived with calibrated intervals

| claim_id | unit | clusters | stored_verdict | calibrated_verdict_w0 | verdict_w0_low_end | verdict_w0_high_end |
|---|---|---|---|---|---|---|
| PERS-003 | site | 1089 | undecidable | undecidable | attenuated | undecidable |
| SFE-008 | site | 3662 | undecidable | vanishes | vanishes | vanishes |

### Verdicts under multiplicity-adjusted intervals (normal approximation)

| claim_id | unit | stored_verdict | stored_95 | bonferroni_F84 | bh_fcr_F84 | bonferroni_F112 | bh_fcr_F112 | bonferroni_F28_per_type | bh_fcr_F28_per_type |
|---|---|---|---|---|---|---|---|---|---|
| SNO-004 | site | null_broken_by_control | null_broken_by_control | null_survives | null_survives | null_survives | null_survives | null_survives | null_survives |
| SNO-006 | site | survives | survives | survives | survives | undecidable | survives | survives | survives |
| SNO-012 | site | attenuated | attenuated | undecidable | attenuated | undecidable | attenuated | undecidable | attenuated |

```
{
 "levels": {
  "F84": {
   "m": 84,
   "bonferroni_level": 0.9994047619047619,
   "bh_R": 46,
   "fcr_level": 0.9726190476190476
  },
  "F112": {
   "m": 112,
   "bonferroni_level": 0.9995535714285714,
   "bh_R": 60,
   "fcr_level": 0.9732142857142857
  },
  "F28_baseline": {
   "m": 28,
   "bonferroni_level": 0.9982142857142857,
   "bh_R": 16,
   "fcr_level": 0.9714285714285714
  },
  "F28_matched": {
   "m": 28,
   "bonferroni_level": 0.9982142857142857,
   "bh_R": 15,
   "fcr_level": 0.9732142857142857
  },
  "F28_random": {
   "m": 28,
   "bonferroni_level": 0.9982142857142857,
   "bh_R": 15,
   "fcr_level": 0.9732142857142857
  },
  "F28_stratified": {
   "m": 28,
   "bonferroni_level": 0.9982142857142857,
   "bh_R": 14,
   "fcr_level": 0.975
  }
 },
 "expected_chance_exclusions_global_null": {
  "nominal_F84": 4.2,
  "nominal_F112": 5.6000000000000005,
  "p_at_least_one_nominal_F84_independent": 0.9865481243423165,
  "simulated_F84": 3.184499999999999,
  "simulated_F112": 4.657338175932566,
  "p_at_least_one_simulated_F84_independent": 0.962031957669314
 },
 "observed": {
  "observed_exclusions_F84": 47,
  "observed_exclusions_F112": 62,
  "observed_exclusions_by_type": {
   "baseline": 16,
   "matched": 16,
   "random": 15,
   "stratified": 15
  }
 }
}
```

## Revision after verification (round 1), POST HOC

### R1-a. Random control: decomposition of the over-width (8 null cells per K; mean [range])

| K | mixing factor width(rc)/width(one draw) | SD(rc point)/SD(one-draw point) | SD(whole-cohort point)/SD(one-draw point) | w0_rc (nominal) | over-width vs own point (1/w0_rc) | w_single (single-draw yardstick) | over-width vs one calibrated draw | w_uncal | type-I rc |
|---|---|---|---|---|---|---|---|---|---|
| 44 | 1.188 [1.173-1.201] | 0.786 [0.765-0.811] | 0.774 | 0.684 | 1.46 | 0.908 [0.864-0.958] | 1.10 [1.04-1.16] | 0.842 | 0.0063 |
| 164 | 1.174 [1.159-1.186] | 0.798 [0.775-0.826] | 0.785 | 0.685 | 1.46 | 0.868 [0.851-0.880] | 1.15 [1.14-1.17] | 0.852 | 0.0047 |
| 776 | 1.168 [1.150-1.179] | 0.797 [0.757-0.832] | 0.785 | 0.683 | 1.46 | 0.855 [0.810-0.889] | 1.17 [1.12-1.24] | 0.856 | 0.0043 |
| 1475 | 1.169 [1.146-1.180] | 0.794 [0.759-0.822] | 0.782 | 0.687 | 1.46 | 0.871 [0.848-0.897] | 1.15 [1.11-1.18] | 0.856 | 0.0056 |

### R1-b. Random-control yardsticks at claim-matched cohort shapes (95% CI from 2,000 dataset resamples)

| claim | run | datasets | w0_rc (nominal) | w_single | w_uncal | mixing | point SD ratio |
|---|---|---|---|---|---|---|---|
| PERS-003 | first claim-matched run, re | 2000 | 0.817 [0.780, 0.860] | 0.907 [0.872, 0.964] | 0.912 | 1.096 | 0.879 |
| PERS-006 | first claim-matched run, iid | 2000 | 0.602 [0.576, 0.637] | 0.842 [0.769, 0.877] | 0.805 | 1.242 | 0.726 |
| SFE-008 | first claim-matched run, re | 2000 | 0.657 [0.630, 0.675] | 0.871 [0.831, 0.902] | 0.842 | 1.188 | 0.764 |
| SFE-011 | first claim-matched run, iid | 2000 | 0.603 [0.581, 0.631] | 0.817 [0.790, 0.851] | 0.818 | 1.223 | 0.773 |
| SNO-014 | first claim-matched run, nega | 2000 | 0.674 [0.641, 0.700] | 0.921 [0.855, 0.975] | 0.834 | 1.200 | 0.783 |
| SNO-014 | first claim-matched run, re | 2000 | 0.715 [0.681, 0.759] | 0.948 [0.883, 0.991] | 0.843 | 1.186 | 0.787 |
| SFE-008 | r1_sfe008_sim iid_nb077 (singletons 0.582) | 4000 | 0.649 [0.634, 0.672] | 0.839 [0.814, 0.864] | 0.842 | 1.188 | 0.761 |
| SFE-008 | r1_sfe008_sim re_nb077 (singletons 0.581) | 4000 | 0.654 [0.633, 0.676] | 0.841 [0.812, 0.873] | 0.843 | 1.186 | 0.764 |
| SFE-008 | r1_sfe008_sim re_nb2 (singletons 0.515) | 8000 | 0.649 [0.639, 0.664] | 0.868 [0.846, 0.883] | 0.843 | 1.187 | 0.753 |

### R1-c. Decisive random-control intervals under the single-draw yardstick (primary)

| claim | stored random control | w* | w_single (range) | bound at w_single (range) | bound at w_uncal | w0_rc (nominal) | flips at w_single | could flip in range | source |
|---|---|---|---|---|---|---|---|---|---|
| PERS-003 | 0.2018 [-0.2679, 0.6513] | 0.430 | 0.907 (0.872-0.964) | -0.224 (-0.251 to -0.208) | -0.227 | 0.817 | no | no | claim-matched re (range over re) |
| PERS-006 | 0.9425 [-0.9430, 3.7036] | 0.500 | 0.842 (0.769-0.877) | -0.646 (-0.711 to -0.508) | -0.576 | 0.602 | no | no | claim-matched iid (range over iid) |
| PERS-008 | -0.0069 [-1.2511, 1.3116] | 0.005 | 0.831 (0.771-0.933) | 1.089 (1.010 to 1.224) | 1.113 | 0.665 | no | no | grid null_iid K=776 pi=0.15 (range over null_iid,null_re x nb,lognormal) |
| SFE-006 | 1.1826 [0.8993, 1.4720] | 4.174 | 0.874 (0.804-0.957) | 0.935 (0.912 to 0.955) | 0.937 | 0.689 | no | no | grid null_re K=1475 pi=0.35 (range over null_iid,null_re x nb,lognormal) |
| SFE-007 | 0.0315 [-0.2218, 0.2862] | 0.124 | 0.849 (0.804-0.957) | -0.183 (-0.211 to -0.172) | -0.186 | 0.687 | no | no | grid null_re K=1475 pi=0.15 (range over null_iid,null_re x nb,lognormal) |
| SFE-008 | 0.4111 [-0.0698, 0.8696] | 0.855 | 0.841 (0.812-0.902) | 0.007 (-0.023 to 0.021) | 0.006 | 0.654 | yes | yes | claim-matched r1_sfe008 re_nb077 (range: CIs of re_nb2, re_nb077, iid_nb077 and the first run) |
| SFE-011 | 0.7676 [-0.6632, 2.8026] | 0.536 | 0.817 (0.790-0.851) | -0.401 (-0.450 to -0.363) | -0.402 | 0.603 | no | no | claim-matched iid (range over iid) |
| SNO-001 | 0.2377 [-0.4028, 0.9003] | 0.371 | 0.872 (0.831-0.925) | -0.321 (-0.355 to -0.295) | -0.314 | 0.705 | no | no | grid null_re K=164 pi=0.35 (range over null_iid,null_re x nb,lognormal) |
| SNO-002 | 0.0830 [-0.5452, 0.7083] | 0.132 | 0.872 (0.831-0.925) | -0.465 (-0.498 to -0.439) | -0.458 | 0.705 | no | no | grid null_re K=164 pi=0.35 (range over null_iid,null_re x nb,lognormal) |
| SNO-005 | 0.1542 [-0.8596, 1.3366] | 0.152 | 0.874 (0.823-0.919) | -0.732 (-0.778 to -0.680) | -0.713 | 0.690 | no | no | grid null_re K=164 pi=0.15 (range over null_iid,null_re x nb,lognormal) |
| SNO-009 | -0.0429 [-0.6295, 0.5400] | 0.074 | 0.872 (0.831-0.925) | 0.466 (0.441 to 0.496) | 0.459 | 0.705 | no | no | grid null_re K=164 pi=0.35 (range over null_iid,null_re x nb,lognormal) |
| SNO-014 | 0.9940 [-0.6646, 2.9449] | 0.599 | 0.921 (0.855-0.991) | -0.533 (-0.650 to -0.424) | -0.388 | 0.674 | no | no | claim-matched nega (range over nega,re) |

### R1-d. Verdicts that change under any calibration reading (primary = percentile intervals at w0, random control at the single-draw yardstick)

| claim_id | is_transfer | stored_verdict | reading | verdict_primary | verdict_primary_low_end | verdict_primary_high_end | verdict_rc_uncalibrated_single_draw_only | verdict_secondary_rc_nominal |
|---|---|---|---|---|---|---|---|---|
| PERS-003 | False | undecidable | boundary | undecidable | attenuated | undecidable | undecidable | undecidable |
| SFE-008 | True | undecidable | boundary | vanishes | vanishes | undecidable | vanishes | vanishes |

### R1-e. SNO-004: registered reading and post hoc sensitivity

stored verdict `null_broken_by_control`; registered bucket `reverses`; released tool: `reverses` (FAIL); stored matched null resolution 0.681; P (normal approx.) matched 0.0320, stratified 0.0409

| adjustment (post hoc, normal approximation) | interval | adjusted interval | excludes 0 | null resolution (half-width, log2) | largest OR inside |
|---|---|---|---|---|---|
| bonferroni_F84 | matched | [-0.436, 1.949] | no | 1.192 | 3.86 |
| bonferroni_F84 | stratified | [-0.367, 1.566] | no | 0.967 | 2.96 |
| bh_fcr_F84 | matched | [-0.021, 1.511] | no | 0.766 | 2.85 |
| bh_fcr_F84 | stratified | [-0.043, 1.200] | no | 0.621 | 2.30 |
| bonferroni_F112 | matched | [-0.462, 1.976] | no | 1.219 | 3.94 |
| bonferroni_F112 | stratified | [-0.388, 1.589] | no | 0.989 | 3.01 |
| bh_fcr_F112 | matched | [-0.024, 1.514] | no | 0.769 | 2.86 |
| bh_fcr_F112 | stratified | [-0.045, 1.202] | no | 0.624 | 2.30 |

### R1-f. Registered benchmark B1 (cleavage) at strength 0, regenerated with its own generator

validation: {"stored_replicates_compared": 50, "max_abs_diff_statistic": 9.71445146547012e-17, "max_abs_diff_ci_lo": 5.551115123125783e-17, "max_abs_diff_ci_hi": 1.1102230246251565e-16, "n_positions_equal": true, "n_positive_equal": true, "detected_first50_stored": 7, "detected_first50_recomputed": 7, "extension_stored_detected_of_300": 21, "extension_recomputed_detected_of_300": 21}

| set | datasets | flagged | rate [Wilson 95%] | below 0 / above 0 | w0 [95% CI] | clusters | positives | near-cut prevalence |
|---|---|---|---|---|---|---|---|---|
| registered_first50 | 50 | 7 | 0.1400 [0.0695, 0.2619] | 0.1000 / 0.0400 | 1.215 [1.028, 1.343] | 298.7 | 148.0 | 0.592 |
| registered_300 | 300 | 21 | 0.0700 [0.0462, 0.1046] | 0.0433 / 0.0267 | 1.055 [0.952, 1.258] | 298.4 | 147.6 | 0.591 |
| fresh_10000 | 10000 | 513 | 0.0513 [0.0471, 0.0558] | 0.0213 / 0.0300 | 1.004 [0.986, 1.020] | 298.3 | 147.9 | 0.591 |

### R1-g. SUPERSEDED IN ROUND 2 (parameters from the inaccurate 20 x 20 Gauss-Hermite fit; see R2-d). Anchored between-protein confounding: P(excludes 0) / coverage of the pooled target (2,000 datasets per cell)

| attribute (anchor) | label SD | attribute SD | rho | prevalence | pooled log2 OR (no within effect) | K=44 | K=164 | K=776 | K=1475 |
|---|---|---|---|---|---|---|---|---|---|
| A_distal | 0.95 | 0.48 | -0.012 | 0.787 | -0.0066 | 0.070 / 0.928 | 0.056 / 0.945 | 0.052 / 0.947 | 0.053 / 0.946 |
| A_acidic | 0.94 | 0.44 | -0.093 | 0.511 | -0.0474 | 0.062 / 0.941 | 0.053 / 0.950 | 0.071 / 0.950 | 0.082 / 0.958 |
| A_nbcys | 0.93 | 1.34 | -0.142 | 0.226 | -0.1802 | 0.079 / 0.936 | 0.102 / 0.946 | 0.257 / 0.943 | 0.435 / 0.948 |

### R1-h. SUPERSEDED IN ROUND 2 (see R2-a, R2-f). Pooled log2 OR with no within-protein association (Gauss-Hermite), anchors and rho +/- 1.96 SE

| attribute | rho_setting | rho | su | sv | prevalence | pooled log2 OR |
|---|---|---|---|---|---|---|
| A_distal | anchor | -0.012 | 0.946 | 0.478 | 0.787 | -0.0066 |
| A_distal | rho_minus_1.96se | -0.185 | 0.946 | 0.478 | 0.787 | -0.1043 |
| A_distal | rho_plus_1.96se | 0.162 | 0.946 | 0.478 | 0.787 | 0.0913 |
| A_distal | first_report_conf15 (SD 1, rho 0.15) | 0.15 | 1.0 | 1.0 | 0.787 | 0.1663 |
| A_acidic | anchor | -0.093 | 0.937 | 0.438 | 0.511 | -0.0474 |
| A_acidic | rho_minus_1.96se | -0.473 | 0.937 | 0.438 | 0.511 | -0.2409 |
| A_acidic | rho_plus_1.96se | 0.287 | 0.937 | 0.438 | 0.511 | 0.1460 |
| A_acidic | first_report_conf15 (SD 1, rho 0.15) | 0.15 | 1.0 | 1.0 | 0.511 | 0.1587 |
| A_nbcys | anchor | -0.142 | 0.928 | 1.339 | 0.226 | -0.1802 |
| A_nbcys | rho_minus_1.96se | -0.211 | 0.928 | 1.339 | 0.226 | -0.2686 |
| A_nbcys | rho_plus_1.96se | -0.073 | 0.928 | 1.339 | 0.226 | -0.0924 |
| A_nbcys | first_report_conf15 (SD 1, rho 0.15) | 0.15 | 1.0 | 1.0 | 0.226 | 0.1643 |

## Revision round 2, POST HOC

### R2-a. Anchoring model refitted with accurate integration (dense trapezoid grid; independent adaptive Gauss-Hermite fit)

| attribute | prevalence | label SD | attribute SD | rho [profile 95% CI] | rho Wald SE | LR test rho = 0: P | within log2 OR (model) | nll (dense) | rho, independent AGHQ-25 fit | pooled log2 OR, no within effect: MLE [CI ends] | round-1 GH20 rho (SE) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| neighbouring Cys | 0.226 | 0.928 | 1.337 | -0.1945 [-0.3251, -0.0588] | 0.0683 | 0.0051 | -0.847 | 8171.7119 | -0.1945 | -0.2472 [-0.4223, -0.0736] | -0.1420 (0.0352) |
| acidic | 0.511 | 0.937 | 0.440 | -0.1040 [-0.2817, 0.0749] | 0.0915 | 0.256 | 0.167 | 10188.2616 | -0.1040 | -0.0532 [-0.1419, 0.0388] | -0.0931 (0.1938) |
| distal cleavage | 0.787 | 0.946 | 0.476 | -0.0046 [-0.1935, 0.1842] | 0.0969 | 0.962 | 1.248 | 8493.0621 | -0.0046 | -0.0026 [-0.1071, 0.1041] | -0.0118 (0.0885) |

### R2-b. Negative log-likelihood under other integrators, minus the primary dense grid (same theta)

| attribute | theta | dense primary nll | dense_z9_h0.02 | dense_z7_h0.05 | aghq_Q25 | aghq_Q15 | aghq_Q40 | gh20_nonadaptive | gh60_nonadaptive | gh100_nonadaptive |
|---|---|---|---|---|---|---|---|---|---|---|
| neighbouring Cys | at_new_mle | 8171.711886 | +3.27e-11 | +2.70e-09 | +2.83e-07 | +1.36e-05 | +5.46e-12 | -1.01e+00 | -1.45e-01 | +5.58e-02 |
| neighbouring Cys | at_round1_gh20_optimum | 8172.011246 | +3.27e-11 | +2.75e-09 | +3.81e-07 | +1.93e-05 | -2.18e-11 | -2.38e+00 | +6.84e-01 | -2.20e-02 |
| acidic | at_new_mle | 10188.261595 | +3.27e-11 | +5.24e-09 | -5.31e-10 | -1.24e-06 | +0.00e+00 | +2.32e-01 | +1.06e-03 | +9.17e-06 |
| acidic | at_round1_gh20_optimum | 10188.272812 | +3.09e-11 | +5.27e-09 | -5.29e-10 | -1.24e-06 | -3.64e-12 | +2.16e-01 | +7.93e-04 | +6.04e-06 |
| distal cleavage | at_new_mle | 8493.062136 | +3.27e-11 | +4.33e-09 | -4.71e-10 | -1.06e-06 | +0.00e+00 | -6.21e-02 | +1.32e-04 | -1.30e-08 |
| distal cleavage | at_round1_gh20_optimum | 8493.066041 | +3.27e-11 | +4.34e-09 | -4.71e-10 | -1.06e-06 | +0.00e+00 | -7.05e-02 | +1.31e-04 | -2.27e-08 |

### R2-c. Profile likelihood of rho: 2 x (nll_profile - nll_min) on a fixed grid (95% cut-off 3.84)

| attribute | -0.40 | -0.36 | -0.32 | -0.28 | -0.24 | -0.20 | -0.16 | -0.12 | -0.08 | -0.04 | +0.00 | +0.04 | +0.08 | +0.12 | +0.16 | +0.20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| neighbouring Cys | 9.88 | 6.28 | 3.54 | 1.62 | 0.45 | 0.01 | 0.25 | 1.17 | 2.74 | 4.97 | 7.84 | 11.38 | 15.59 | 20.50 | 26.13 | 32.53 |
| acidic | 10.92 | 8.09 | 5.71 | 3.77 | 2.24 | 1.11 | 0.38 | 0.03 | 0.07 | 0.49 | 1.29 | 2.48 | 4.06 | 6.05 | 8.46 | 11.30 |
| distal cleavage | 17.50 | 14.00 | 10.93 | 8.27 | 6.00 | 4.11 | 2.59 | 1.42 | 0.61 | 0.13 | 0.00 | 0.21 | 0.76 | 1.66 | 2.91 | 4.52 |

### R2-d. Anchored between-protein confounding at the accurate estimates and across the profile-likelihood CI of rho: P(excludes 0) (MC SE) / coverage of the pooled target; whole cohort, log-normal cluster sizes mean 8.5

| attribute | setting | SDs (label, attribute) | rho | pooled target log2 OR | datasets per cell | K=44 | K=164 | K=776 | K=1475 |
|---|---|---|---|---|---|---|---|---|---|
| neighbouring Cys | mle | 0.928, 1.337 | -0.1945 | -0.2472 | 4000 | 0.092 (0.005) / 0.931 | 0.131 (0.005) / 0.944 | 0.412 (0.008) / 0.945 | 0.648 (0.008) / 0.951 |
| neighbouring Cys | profile_lo | 0.934, 1.360 | -0.3251 | -0.4223 | 2000 | 0.126 (0.007) / 0.927 | 0.292 (0.010) / 0.946 | 0.827 (0.008) / 0.952 | 0.976 (0.003) / 0.950 |
| neighbouring Cys | profile_hi | 0.926, 1.321 | -0.0588 | -0.0736 | 2000 | 0.068 (0.006) / 0.934 | 0.066 (0.006) / 0.940 | 0.096 (0.007) / 0.952 | 0.107 (0.007) / 0.953 |
| neighbouring Cys | round 1, GH20 parameters (superseded) | 0.928, 1.339 | -0.1420 | -0.1802 | 2000 | 0.079 (0.006) / 0.936 | 0.102 (0.007) / 0.946 | 0.257 (0.010) / 0.943 | 0.435 (0.011) / 0.948 |
| acidic | mle | 0.937, 0.440 | -0.1040 | -0.0532 | 4000 | 0.068 (0.004) / 0.934 | 0.061 (0.004) / 0.944 | 0.084 (0.004) / 0.949 | 0.106 (0.005) / 0.950 |
| acidic | profile_lo | 0.938, 0.432 | -0.2817 | -0.1419 | 2000 | 0.078 (0.006) / 0.933 | 0.094 (0.007) / 0.954 | 0.268 (0.010) / 0.941 | 0.475 (0.011) / 0.950 |
| acidic | profile_hi | 0.941, 0.444 | 0.0749 | 0.0388 | 2000 | 0.064 (0.005) / 0.934 | 0.054 (0.005) / 0.949 | 0.063 (0.005) / 0.953 | 0.082 (0.006) / 0.949 |
| acidic | round 1, GH20 parameters (superseded) | 0.937, 0.438 | -0.0931 | -0.0474 | 2000 | 0.062 (0.005) / 0.941 | 0.053 (0.005) / 0.950 | 0.071 (0.006) / 0.950 | 0.082 (0.006) / 0.958 |
| distal cleavage | mle | 0.946, 0.476 | -0.0046 | -0.0026 | 4000 | 0.066 (0.004) / 0.934 | 0.061 (0.004) / 0.940 | 0.050 (0.003) / 0.950 | 0.056 (0.004) / 0.945 |
| distal cleavage | profile_lo | 0.949, 0.467 | -0.1935 | -0.1071 | 2000 | 0.070 (0.006) / 0.934 | 0.079 (0.006) / 0.938 | 0.129 (0.007) / 0.945 | 0.210 (0.009) / 0.955 |
| distal cleavage | profile_hi | 0.947, 0.478 | 0.1842 | 0.1041 | 2000 | 0.082 (0.006) / 0.922 | 0.071 (0.006) / 0.950 | 0.128 (0.007) / 0.959 | 0.203 (0.009) / 0.945 |
| distal cleavage | round 1, GH20 parameters (superseded) | 0.946, 0.478 | -0.0118 | -0.0066 | 2000 | 0.070 (0.006) / 0.928 | 0.056 (0.005) / 0.945 | 0.052 (0.005) / 0.947 | 0.053 (0.005) / 0.946 |

### R2-e. Sensitivity: the anchoring cohort's own 766 cluster sizes (fixed design, mean 12.4, maximum 332)

| attribute | setting | rho | K | mean cysteines | pooled target log2 OR | datasets | P(excludes 0) (MC SE) | below / above 0 | coverage |
|---|---|---|---|---|---|---|---|---|---|
| neighbouring Cys | mle | -0.1945 | 766 | 9492 | -0.2472 | 2000 | 0.488 (0.011) | 0.488 / 0.000 | 0.929 |
| acidic | mle | -0.1040 | 766 | 9492 | -0.0532 | 2000 | 0.085 (0.006) | 0.081 / 0.004 | 0.952 |
| distal cleavage | mle | -0.0046 | 766 | 9492 | -0.0026 | 2000 | 0.051 (0.005) | 0.033 / 0.018 | 0.949 |
| neighbouring Cys | rho0_control | 0.0000 | 766 | 9492 | 0.0000 | 2000 | 0.054 (0.005) | 0.028 / 0.026 | 0.946 |

### R2-f. Pooled log2 OR with no within-protein association (simlib.calibrate, 80-node Gauss-Hermite on the population integrand)

| attribute | setting | rho | su | sv | prevalence | pooled log2 OR |
|---|---|---|---|---|---|---|
| A_nbcys | mle | -0.1945 | 0.9284 | 1.3374 | 0.2259 | -0.2472 |
| A_nbcys | profile_lo | -0.3251 | 0.9337 | 1.3604 | 0.2259 | -0.4223 |
| A_nbcys | profile_hi | -0.0588 | 0.9258 | 1.3206 | 0.2259 | -0.0736 |
| A_acidic | mle | -0.104 | 0.9373 | 0.4403 | 0.5107 | -0.0532 |
| A_acidic | profile_lo | -0.2817 | 0.9381 | 0.4324 | 0.5107 | -0.1419 |
| A_acidic | profile_hi | 0.0749 | 0.9405 | 0.4443 | 0.5107 | 0.0388 |
| A_distal | mle | -0.0046 | 0.946 | 0.4762 | 0.7873 | -0.0026 |
| A_distal | profile_lo | -0.1935 | 0.9493 | 0.4675 | 0.7873 | -0.1071 |
| A_distal | profile_hi | 0.1842 | 0.947 | 0.4779 | 0.7873 | 0.1041 |
| A_nbcys | grid scenario null_conf15 (SD 1, rho 0.15), for comparison | 0.15 | 1.0 | 1.0 | 0.2259 | 0.1643 |
| A_acidic | grid scenario null_conf15 (SD 1, rho 0.15), for comparison | 0.15 | 1.0 | 1.0 | 0.5107 | 0.1587 |
| A_distal | grid scenario null_conf15 (SD 1, rho 0.15), for comparison | 0.15 | 1.0 | 1.0 | 0.7873 | 0.1663 |
