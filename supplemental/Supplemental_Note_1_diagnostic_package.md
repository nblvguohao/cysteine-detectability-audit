# Supplemental Note 1. The three diagnostic checks, the label-semantics tiers and checks not shown in Table 1

The package comprises three checks, and each one is allowed to answer “cannot be computed” rather than return a number: the detectability share (undefined when the superset is itself barely above chance), the label-semantics tier (unclassifiable when the negative set is not defined), and the depth-confounding strength (not applicable when the depth covariate encodes the label).

Every threshold is written in the module docstring before any case is run, and the 55-test regression suite must pass before the package is used on anything. The suite mixes recomputation with replay of recorded precedent, and the audit marks which is which.

The tier is decided by the **negative set**, not by the chemistry. That is the whole point: the same enrichment chemistry can yield a tier that licenses a chemical claim or one that licenses only a ranking claim, depending on what the authors chose to call a negative.

**Label-semantics tiers.**

| tier | negative_semantics | claims_licensed | required_control | claims_in_this_corpus |
|---|---|---|---|---|
| T1_direct_adduct | same_run_mutually_exclusive | chemistry claims permitted | none beyond the usual grouping and interval discipline | 3 |
| T2_parallel_ambiguous | detected_background_proteome | chemistry claims permitted only with a detectability control reported | a detectability control must be reported with the claim | 3 |
| T3_capture_annotation | all_cys_in_identified_proteins; curated_negatives_other_db; random_or_simulated; unmodified_cys_same_protein; whole_proteome_annotation_background | detectability-ranking claims only; chemistry claims not supported | no chemistry claim; report detectability-only AUCs (global and within protein), not the share alone | 39 |
| unclassifiable | not_defined | no claim; the negative set must be defined first | define the negative set | 6 |

**Limitations.**

- The tier assignment is a reading of each paper's stated negative set. Where a paper does not state one, the claim is unclassifiable rather than assigned a default.
- The regression suite verifies the implementation against its declared expectations; the thresholds themselves are documented design choices.

**Source.** The tables above are printed from the stored analysis outputs.

**Checks not shown in Table 1.** Table 1 of the main text lists eight checks in decision order. The following checks complete the checklist.

| Check | Pre-analysis requirement | Reading rule if the check fails |
|---|---|---|
| Label swap (extends Table 1, item 1) | If needed, compute the bound 2d(1/p + 1/q) on how far an AUC can move when d instances are moved between a positive set of size p and a negative set of size q, then refit with the label "was this residue identified" and report the recovered share of above-chance discrimination with a cluster-bootstrap interval. | If the bound is smaller than the smallest AUC difference you would act on, the label swap cannot fail and carries no information; otherwise call the rules observability-confounded when the recovered-share interval has a lower bound of 0.5 or higher. |
| Co-residue depletion and protease reach | Compare the observed distribution of identified peptides by modifiable-residue count with the theoretical digest, then report the fraction observable individually under the protease used. | Residues failing the individual-observation criterion were not interrogated, not found unmodified. |
| Multiple testing | Correct across the full family of tests actually performed. | Report corrected values and preserve the complete tested family. |
