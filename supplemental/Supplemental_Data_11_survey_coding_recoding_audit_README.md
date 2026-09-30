# Supplemental Data 11. Re-coding audit of the survey coding table

One row per record classified after full-text reading (74 rows), joined to the class assigned by each of
two re-coding passes that were blinded to the original code. Companion to
`Supplemental_Data_9_survey_coding_table.csv` (942 screened records; 74 classified).

**Procedure.** Two re-coding passes were run independently. Each received an extract of the 74 records
carrying only the title, the modification, the study type and the recorded evidence quotations, with
`background_class` and every other coded flag withheld, and each coded all 74 records under the same
rubric. Both passes were allowed to answer `unclassifiable from the retained evidence` (shown as `U`)
instead of assigning a class; pass 1 did so for 1 record and pass 2 for 28. Agreement was computed only
after both passes had been written to file.

**Rubric.** `a` all residues of that type, unrestricted; `b` matched on theoretical detectability or on
abundance; `c` residues detected in the same experiment without a modification assignment; `e` other, or
no statement found.

**Agreement with the original codes (74 records).**

| pass | exact agreement | Cohen's kappa |
|---|---|---|
| pass 1 | 46/74 (62.2%) | 0.31 |
| pass 2 | 37/74 (50.0%) | 0.31 |
| pass 1, crediting the 12 records whose background statement was found on re-grep but not retained | 57/74 (77.0%) | 0.58 |
| pass 2, same crediting | 48/74 (64.9%) | 0.46 |

Kappa is computed over the categories a/b/c/e/U with `U` as its own category.

**Where the disagreements are.**

| pass | disagreements | of which the retained quotation is cut off mid-clause | of which it is not |
|---|---|---|---|
| pass 1 | 28 of 74 | 22 | 6 |
| pass 2 | 37 of 74 | 32 | 5 |

Of the 74 retained background quotations, 51 are cut off mid-clause, in several cases immediately before
the clause that names the background. The disagreements that do not rest on truncation are the
consequential ones: in each, the retained evidence describes a negative set that is curated,
substrate-restricted or drawn from an external database rather than detected in the same experiment, and
the original code and the re-coding passes read that description differently. Disagreement is otherwise
concentrated on whether a background that the paper does not state should be read as the unrestricted
residue set or left unstated.

**Columns.** `original_code` is the code in Supplemental Data 9; `code_pass_1` / `code_pass_2` are the
re-coding passes; `pass_N_agrees_with_original` and `both_passes_agree_with_original` are string booleans;
`background_quotation_retained` is the evidence quotation the original coding rested on;
`quotation_truncated` is true when that quotation does not end in terminal punctuation;
`verification_note` is the original two-pass verification flag; `ambiguity_note` is where the original
coder recorded a doubt.

**Full-text re-check.** The five records on which the original count of matched backgrounds rested were then
re-read against the papers' full texts. The class of two (uid 441, uid 480) was confirmed; the other three
(uid 1111, uid 1537, uid 1553) were confirmed as not using a matched background, and their class in
Supplemental Data 9 was corrected with the original code retained in `background_class_original`. The
full-text re-check did not record verbatim quotations, so the quotations in this table remain the retained
extracts. `final_class` and `full_text_check` record the outcome.

**Limitations.** Both re-coding passes read the retained evidence quotations, not the papers' full texts,
so a background statement that exists in a paper but was not captured by the original extraction cannot be
recovered by either pass. The counts they give are therefore lower bounds, and the same limitation is what
makes the original `a` codes hard to reproduce. All passes used the same language model, so the agreement
figures bound how far the retained evidence determines the class; they are not a measurement of agreement
between human coders. The records that require a full-text decision are listed in the manual re-check table
accompanying the analysis records.
