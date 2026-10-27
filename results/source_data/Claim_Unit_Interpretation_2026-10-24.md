# Claim unit interpretation correction

24 October 2026. This note supersedes the description of the ten protein-level
claims as lacking a reference class or being unestimated. It does not overwrite
historical estimates or change any value in the 18-row current Claim_Ledger.csv.

Applying the retained final-selection rules yields 28 classified claims with
measured baselines: 18 with site-level and 10 with protein-level inference.
The latter are PERS-005, PERS-006, PERS-008, PERS-010, PERS-011, PERS-012,
SFE-011, SFI-001, SFI-002 and SNO-008. Their archived re-test records contain
baseline/matched estimates and classifications. The absence of those estimates
from a site-level display does not show that no reference existed.

The current ledger and Figure S1d summarize only the 18 site-level claims:
4 survives, 2 attenuated, 10 undecidable, 1 null broken, and 1 transfer-cohort
baseline contradiction. The ten protein-level claims remain outside this
denominator because their unit of inference differs. This correction changes
the explanation of their exclusion, not the site-level counts or outcomes.

Evidence: the historical Supplemental_Data_3/4/5_retest_round_*.csv tables,
their final-selection logic in claim_ledger.py, and the independently selected
28-row claim_selected_units.csv supplied with the final review records.
