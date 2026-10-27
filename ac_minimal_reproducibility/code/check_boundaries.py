"""Analytical boundary cases and independent weighted-AUC checks."""
import numpy as np
from sklearn.metrics import roc_auc_score
from feature_definitions import FEATURE_NAMES, feature_row
from public_probe_snapshot import boundaries_for, descriptors, PROTEASES
from run_closure import WeightedAUC, write_json
from prepare_inputs import ROOT


def check():
    assert boundaries_for("ACKPCCRAAC", PROTEASES["Trypsin"]) == [0, 7, 10]
    length, mass, gravy, count = descriptors("AC", 0, 2)
    assert length == 2 and count == 1
    assert abs(mass - 192.056865) < 1e-9 and abs(gravy - 2.15) < 1e-12
    def get(row, name):
        return row[FEATURE_NAMES.index(name)]
    row = feature_row("ACKCCKCCCC", 4, "strict_mc_v1")
    assert get(row, "pep_len") == 3 and get(row, "pep_offset_in_peptide") == 0
    strict = feature_row("GGGCGGG", 4, "strict_mc_v1")
    compat = feature_row("GGGCGGG", 4, "probe_compatible_v1")
    assert get(strict, "pep_detectable_length") == 1 and get(strict, "pep_detectable_mass") == 0
    assert get(strict, "pep_detectable_any_missed_cleavage") == 0
    assert get(compat, "pep_detectable_any_missed_cleavage") == 1
    strict = feature_row("CCCCCCC", 1, "strict_mc_v1")
    compat = feature_row("CCCCCCC", 1, "probe_compatible_v1")
    assert get(strict, "pep_detectable_both") == 1
    assert get(strict, "pep_mc1_detectable") == get(strict, "pep_mc2_detectable") == 0
    assert get(compat, "pep_mc1_detectable") == get(compat, "pep_mc2_detectable") == 1
    seq = "A"*49 + "K" + "CCKCCKCCK" + "A"*49 + "K"
    strict, compat = [feature_row(seq, 54, d) for d in ("strict_mc_v1", "probe_compatible_v1")]
    assert get(strict, "pep_mc1_detectable") == 0
    assert get(strict, "pep_mc2_detectable") == 1 and get(compat, "pep_mc2_detectable") == 0
    y = np.array([0, 1, 1, 0, 0, 1, 0, 1])
    s = np.array([.2, .2, .7, .4, .7, .8, .1, .4])
    units = np.array([0, 0, 1, 1, 2, 2, 3, 3])
    fast = WeightedAUC(y, s, units, np.ones(len(y), dtype=bool))
    rng = np.random.default_rng(111)
    for _ in range(100):
        weights = rng.integers(1, 5, 4)
        w = weights[units]
        slow = sum(w[i]*w[j]*(float(s[i]>s[j])+.5*float(s[i]==s[j]))
                   for i in np.flatnonzero(y) for j in np.flatnonzero(y == 0)) / (w[y == 1].sum()*w[y == 0].sum())
        assert abs(fast(weights)-slow) < 1e-12
        assert abs(fast(weights)-roc_auc_score(y, s, sample_weight=w)) < 1e-12
    write_json(ROOT / "logs/boundary_checks.json", dict(passed=True, auc_crosschecks=100,
               cases=["K/R-P", "mass-and-GRAVY", "cut-boundary-target", "low-mass-length-valid",
                      "protein-ends-and-exact-missed-count", "mixed-missed-cleavage", "weighted-AUC-ties"]))
    print("boundary checks passed", flush=True)


if __name__ == "__main__":
    check()
