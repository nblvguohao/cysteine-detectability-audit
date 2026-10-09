# Supplemental Note 8. Structural site-preference claims, and what "undecidable" means here

Five published claims tie the modified cysteine to a structural attribute. Each was re-tested on the authors' own positive set, with the attribute recomputed from AlphaFold models (relative accessibility by the Shrake-Rupley algorithm, secondary structure by P-SEA, both as implemented in biotite) and with the same matched-background machinery used throughout.

**Each claim is read as its authors state it.** `SFE-001` states relative accessibility above 25%. `SNO-002` (Doulias et al. 2010) states that a relative accessibility of at most 10% denotes a buried cysteine, normalized to cysteine in the extended tripeptide Ala-Cys-Ala; the reference area itself is not printed, so two published tripeptide values, 140 A^2 (primary) and 135 A^2, are used. `SNO-012` (Marino and Gladyshev 2010) states a sulfur-atom criterion of 1.0 A^2 for burial and reports that about 35% of modified sulfur atoms are buried despite a slight overall enrichment for exposure; it is re-tested as a preference for exposure, which is what its authors report, and not as a null. An initial coding borrowed the 25% threshold for `SNO-002` in the belief that its authors gave none and read `SNO-012` as a null; both codings were corrected, and the rows produced by the initial coding are kept, marked superseded, in Supplemental Data 6.

**Verdicts.** `SFE-001` survives in both specifications. `SNO-002` is undecidable in every specification; its composition rebuilds to about a third of modified and unmodified cysteines exposed, against the 29% and 23% its authors report from experimental structures, so the recomputed attribute does not reconstruct the published one closely enough to settle the claim. `SNO-012` strengthens under the control in the primary specification (matched 1.75 [0.53, 3.22]), but its baseline interval covers zero, so it is undecidable there under the stated reading rule and the released tool; it survives with pLDDT >= 70, so its verdict depends on the specification. `SNO-001`'s categories do not reproduce under a different secondary-structure assignment and `SNO-009`'s per-position statistic does not reproduce; both are undecidable. Undecidable is not refuted: it says the claim could not be settled with this instrument on this material.

**All specifications.**

| claim_id | specification | verdict | n_observations | n_positive | baseline_log2_or | baseline_ci | matched_log2_or | matched_ci | attribute_rate_positive | attribute_rate_negative |
|---|---|---|---|---|---|---|---|---|---|---|
| SFE-001 | primary | survives | 8260 | 1030 | 2.5697 | [2.3428, 2.8151] | 2.2765 | [1.9773, 2.5863] | 0.6631 | 0.2488 |
| SFE-001 | pLDDT >= 70 | survives | 6929 | 718 | 2.7669 | [2.5005, 3.0347] | 2.6164 | [2.2435, 3.0140] | 0.5348 | 0.1444 |
| SNO-002 | primary | undecidable | 1468 | 238 | 0.0539 | [-0.3940, 0.4798] | 0.0269 | [-0.5445, 0.6326] | 0.3403 | 0.3325 |
| SNO-002 | reference area 135 A^2 | undecidable | 1468 | 238 | 0.0439 | [-0.4037, 0.4657] | 0.0 | [-0.5579, 0.6033] | 0.3445 | 0.3382 |
| SNO-002 | pLDDT >= 70 | undecidable | 1401 | 232 | 0.1986 | [-0.2502, 0.6536] | -0.1109 | [-0.7277, 0.4957] | 0.3276 | 0.2985 |
| SNO-012 | primary | undecidable | 366 | 55 | 0.5943 | [-0.2744, 1.5034] | 1.7482 | [0.5284, 3.2156] | 0.5636 | 0.4598 |
| SNO-012 | pLDDT >= 70 | survives | 343 | 54 | 0.8945 | [0.0478, 1.8060] | 1.5222 | [0.2667, 2.9092] | 0.5741 | 0.4187 |
| SNO-001 | primary | undecidable | 1468 | 238 | 0.2587 | [-0.1647, 0.6698] | 0.3625 | [-0.1772, 0.9152] | 0.3445 | 0.3057 |
| SNO-001 | pLDDT >= 70 | undecidable | 1401 | 232 | 0.2152 | [-0.2367, 0.6453] | 0.4601 | [-0.0819, 1.0031] | 0.3491 | 0.3165 |
| SNO-009 | primary | undecidable | 1468 | 238 | -0.0576 | [-0.4307, 0.3260] | -0.0242 | [-0.5195, 0.4550] | 0.521 | 0.5309 |
| SNO-009 | pLDDT >= 70 | undecidable | 1401 | 232 | -0.0417 | [-0.4135, 0.3475] | 0.0496 | [-0.4607, 0.5400] | 0.5129 | 0.5201 |

**Limitations.**

- The attribute is recomputed from predicted structures, not measured ones; the authors of `SNO-002` and `SNO-012` used experimental or homology structures, and a different secondary-structure assignment is itself one of the reasons a claim lands undecidable.
- `undecidable` is a statement about this instrument on this material. It licenses no conclusion about the underlying chemistry.

**Source.** The tables above are printed from the stored analysis outputs.
