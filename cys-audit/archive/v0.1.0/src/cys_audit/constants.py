"""Fixed decision constants. Every threshold that can change a status lives here and nowhere else.

These values were written into protocols/phase4_cys_audit_preregistration_2026-09-22.json (in the
analysis tree that produced this tool) BEFORE the tool was run on any real dataset, and the Phase 4
validation script refuses to run if the two disagree. Change them only with a new pre-registration.
"""

VERSION = "0.1.0"

# Status vocabulary for a single artefact test (brief s.12)
PASS, WARNING, FAIL, UNDECIDABLE = "PASS", "WARNING", "FAIL", "UNDECIDABLE"
STATUSES = (PASS, WARNING, FAIL, UNDECIDABLE)
# Severity order used only to pick the headline among sub-statistics and to summarise a dataset
SEVERITY = {FAIL: 3, WARNING: 2, UNDECIDABLE: 1, PASS: 0}

# Overall claim status vocabulary (brief s.12)
CLAIM_PASS, CLAIM_ATTENUATED, CLAIM_UNDECIDABLE, CLAIM_FAIL = "PASS", "ATTENUATED", "UNDECIDABLE", "FAIL"
CLAIM_NOT_EVALUATED = "NOT_EVALUATED"

# Resampling (house rule 2: protein- or component-clustered bootstrap, 5000 replicates, fixed seed)
BOOTSTRAP_REPS = 5000
DEFAULT_SEED = 20260922
CI_LEVEL = 0.95

# Minimum class sizes below which an interval-based test is UNDECIDABLE rather than computed
MIN_POSITIVES = 20
MIN_BACKGROUND = 20

# Minimal meaningful differences (equivalence margins) per statistic
MMD_LOG2_OR = 0.5          # |log2 odds ratio|; 0.5 = a 1.41-fold odds change
MMD_AUC = 0.05             # |AUC - 0.5|

# Cleavage-geometry bands, in residues from the cysteine (both sides). The distal band 6-12 is the
# attribute the manuscript's self-audit uses; the proximal band 1-3 matches the Phase 2 benchmark
# window. Two bands are tested, so each interval is Bonferroni-widened to 1 - 0.05/2.
CLEAVAGE_BANDS = {"proximal_1_3": (1, 3), "distal_6_12": (6, 12)}
CLEAVAGE_CI_LEVEL = 1 - (1 - CI_LEVEL) / len(CLEAVAGE_BANDS)

# Positive = detected coincidence thresholds (Artefact 4). Same cut points as the pre-registered
# K-HIGH / K-MID branches of the public Artefact 4 analysis; applied to the clustered-bootstrap bound.
OVERLAP_FAIL_LOWER = 0.90
OVERLAP_WARNING_LOWER = 0.50
OVERLAP_PASS_UPPER = 0.50

# Abundance test: rows with an abundance value must cover at least this share of the tested rows
ABUNDANCE_MIN_COVERAGE = 0.70

# Search-space test: share of positive peptides above the resolving mass that turns a missing
# alternative into FAIL (when observed peptide masses are supplied)
SEARCHSPACE_FAIL_SHARE = 0.01
DEFAULT_PEPTIDE_MASS_RANGE = (600.0, 4600.0)

# Claim retest
ATTENUATION_FLOOR = 0.5    # identical to scripts/ptm_detectability_diagnostics.py
ORPHANED_LIMIT = 0.05      # identical to scripts/ptm_detectability_diagnostics.py

# Monoisotopic masses (Da) used to build the isobaric-alternative table. Derived, not quoted.
MASS_S = 31.9720711744   # 32S, NIST atomic weights and isotopic compositions
MASS_O = 15.9949146196   # 16O, NIST
