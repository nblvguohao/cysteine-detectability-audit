"""Two explicitly versioned feature definitions; neither asserts historical identity.

probe_compatible_v1 preserves the recovered public probe exactly.
strict_mc_v1 enumerates all contiguous peptides containing the target and at
most two missed cuts, including a mixed left/right missed-cut peptide.
"""
import bisect
import numpy as np
import public_probe_snapshot as probe

DEFINITIONS = ("strict_mc_v1", "probe_compatible_v1")
FEATURE_NAMES = probe.FEATURE_NAMES
FEATURE_SETS = {
    "VIS10": [x for x in FEATURE_NAMES if x in probe.VISIBILITY_ONLY],
    "DIG25": list(FEATURE_NAMES),
    "CPL15": [x for x in FEATURE_NAMES if x not in probe.VISIBILITY_ONLY],
}


def feature_row(sequence, site, definition):
    assert definition in DEFINITIONS
    index = int(site) - 1
    assert 0 <= index < len(sequence) and sequence[index] == "C"
    rule = probe.PROTEASES["Trypsin"]
    bounds = probe.boundaries_for(sequence, rule)
    row = probe.features(sequence, index, bounds, rule)
    if definition == "probe_compatible_v1":
        return row
    k = bisect.bisect_right(bounds, index) - 1
    any_detectable = False
    for missed in (0, 1, 2):
        candidates = []
        for left_count in range(missed + 1):
            left, right = k - left_count, k + 1 + missed - left_count
            if left < 0 or right >= len(bounds):
                continue
            length, mass, _, _ = probe.descriptors(sequence, bounds[left], bounds[right])
            candidates.append(probe.DETECTABLE_LENGTH[0] <= length <= probe.DETECTABLE_LENGTH[1]
                              and probe.DETECTABLE_MASS[0] <= mass <= probe.DETECTABLE_MASS[1])
        detectable = any(candidates)
        if missed:
            row[FEATURE_NAMES.index("pep_mc%d_detectable" % missed)] = float(detectable)
        any_detectable = any_detectable or detectable
    row[FEATURE_NAMES.index("pep_detectable_any_missed_cleavage")] = float(any_detectable)
    return row


def matrix(sites, reference, definition):
    return np.asarray([feature_row(reference[r.acc], r.site, definition)
                       for r in sites.itertuples()], dtype=np.float64)
