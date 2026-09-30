# Supplemental Data 9. The survey coding table

One row per record retained at title and abstract screening (942 rows). The columns are as coded; `in_classified_sample` marks the 120 records drawn for full-text coding, and `background_class` is filled only for the 74 that proved to be site-level feature analyses.

**Funnel, recomputed from this table**

| step | records |
|---|---:|
| rows in the coding table (screened and included at title/abstract) | 942 |
| full text retrieved | 367 |
| sampled for classification | 120 |
| classified after full text | 74 |

**Background classes over the classified corpus**

`background_class` carries the class as released, after the blinded re-coding audit and the full-text
re-check of the records that changed; `background_class_original` carries the class as first coded, and
`recoding_note` records why a row changed. Three rows changed: uid 1111 and uid 1537 to `e`, uid 1553 to `a`.

| class | meaning | first coded | as released | rows that moved |
|---|---|---:|---:|---|
| a | all residues of that type, unrestricted | 34 | 35 | 1553 joins |
| b | matched on theoretical detectability or abundance | 1 | 0 | 1537 leaves |
| c | residues detected in the same experiment without a modification assignment | 4 | 2 | 1111 and 1553 leave |
| e | other, or no statement found | 35 | 37 | 1111 and 1537 join |

The released count of analyses using a matched background (classes b and c) is therefore **2 of 74 (2.7%)**,
against 5 as first coded. Both are S-palmitoylation studies. The per-record re-coding audit, with the class
assigned by each blinded pass, is Supplemental Data 11.

There is no member of class d (an equal-size shuffled set) in the classified corpus. `equal_size_random_background` is a separate, non-exclusive flag and is set for 4 papers.

**Modification composition of the classified corpus**

| modification | papers |
|---|---:|
| phosphorylation | 24 |
| acetylation | 12 |
| other-lysine-acyl | 12 |
| S-nitrosylation | 7 |
| multi-PTM | 5 |
| ubiquitination | 4 |
| other | 3 |
| palmitoylation | 3 |
| S-glutathionylation | 1 |
| S-sulfenylation | 1 |
| cysteine-oxidation-general | 1 |
| persulfidation | 1 |

**Coding procedure.** Titles and abstracts were screened by a large language model (Claude, Anthropic) against a written rubric; the false-negative rate of this screen was not estimated. The sampled full texts were coded against the same rubric by Claude from regex-extracted passages, with two verification passes over the full-text XML. Evidence quotations are recorded per coded field in the `evidence_quote_*` columns, and `ambiguity` records where the coder judged the paper's statement unclear.
