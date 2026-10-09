# Supplemental Data 12. Analysis provenance

One row per analysis in the paper: whether it was registered, which artifact carries the registration,
whether its reading rule was fixed before the run, and which decisions were taken after a result had been
seen. This is the single place where the paper's pre-registration status is recorded, so that a reader can
check any analysis without searching the text.

**Reading the table.** `registration_status` is `registered`, `registered per round`, `partial` or `not stated`, and `not registered (post hoc)` for the analyses added after the primary analyses were complete (rows R14 to R22); "not stated" means
the manuscript does not claim a registration for that analysis, not that none exists in the authors'
records. `reading_rule_fixed_before_run` is about the decision rule, which can be fixed even where no
registration is claimed. `decision_after_result_seen` records every reading the manuscript labels post hoc,
including those where the registered analysis was replaced by a later reading; those are also stated in the
text where they affect a number. Registration artifacts named here are in the code repository (Data
Availability). In `where_reported`, R1–R7 denote the seven Results sections of the main text in order, and SR1–SR12 the sections of the Supplemental Results.

**Why this exists.** The paper asks readers to check other people's analyses. This table is the equivalent
audit of its own, and it is the reason the text can carry few provenance sentences: the per-analysis record
lives here rather than in the Results.
