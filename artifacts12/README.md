# artifacts12/ - scripts and stored outputs behind Artifacts 1 and 2

The numbers printed for Artifact 1 (basic-residue enrichment of PlantPTMViewer cysteine sites, and its
reduction after matching the background on theoretical tryptic detectability), for Artifact 2 (the
neighbouring-cysteine statistic under asymmetric and symmetric single-cysteine filters) and the
4,000-permutation P of Experimental Procedures were computed in a separate working tree (the
"collaborator tree"). This folder deposits those scripts and outputs **unchanged** (byte for byte; two
CSV files keep the CRLF line endings Python's `csv` module wrote). They were copied on 2 October 2026.

`number_registry.csv` maps every printed number to its file, JSON path, stored value and producing script.

## Integrity

Twelve files were copied into a staging folder on 23 September 2026 and their SHA-256 recorded in
`MANIFEST.json` (which was itself recorded). All twelve deposited files match the staged SHA-256 on their
raw bytes (column "staged"). Every other file carries a SHA-256 computed when it was deposited. The
repository-wide `MANIFEST_SHA256.txt` lists the same values.

| file | bytes | SHA-256 | staged 2026-09-23 | what it is | manuscript numbers |
|---|---:|---|---|---|---|
| `MANIFEST.json` | 4,804 | `060996c3e6bc87514ad8544b2100faf7796ac4c594612751c76a0e9fbfc58ec5` | match | staging manifest: SHA-256 of the 11 staged files, PlantPTMViewer categories, Arabidopsis row counts, leading PubMed identifiers and the SHA-256 of the three exports (not redistributed) | PubMed identifiers of Artifact 1 |
| `scripts/09_map_large_scale.py` | 3,349 | `7e8ee590ca1cfc82be5a1eb36271c1f4ccd7585ae9f185748d95a73064aec8ec` | match | maps PlantPTMViewer Arabidopsis cysteine rows to protein sequences (step 1 of Artifact 1) | - (see Known gaps 1) |
| `scripts/10_large_scale_motif.py` | 4,438 | `e1873dfe5b49785d1852129f68d459a33a5f824eca1c2d95739bf45d7b704fd2` | match | within-protein residue-class enrichment; 4,000-permutation check | permutation P 0.7435 (Experimental Procedures) |
| `scripts/12_detectability_matched.py` | 4,333 | `78ec4e7d5c876c184b305db50105f988feecaf8dc2f7cdc617842f825f227096` | match | raw and detectability-matched KRH and DE z; rule 7-30 residues, 2 missed cleavages | Artifact 1 z values |
| `scripts/16_structural_matched_test.py` | 9,597 | `4cf678badecddb1b7e57fac4698e645d040e45658e97866496c49e91e6efd00e` | - | helper module imported by 21 and 36 (`detectable_positions`) | - |
| `scripts/17_disulfide_and_rank.py` | 6,273 | `3fcf23532689e586eabd5bbc4be312eed7957aad246c188dbcffdac46e503079` | - | helper module imported by 21 and 36 (`norm_rank`) | - |
| `scripts/21_crossspecies_test.py` | 13,305 | `f9c1f32d984e4644f28d3b85be9c1865e6ded9000bf9ff8b9d73cc3f0bf9fbc5` | match | cross-dataset structural features, asymmetric neighbouring-Cys test | rice -2.61, P 0.009, 756 sites |
| `scripts/24_build_persulfidation_table.py` | 4,794 | `80536b3ca89440cc37a9be4bd7dfb41dc7adbb942279c0c18bfc22f90201ef7f` | - | rice persulfidation site tables from PXD072089 `SS-all-peptides.tsv` | inputs of 21 and 36 |
| `scripts/36_ncys_decompose.py` | 8,909 | `f74133d4abf81544ecf854502721b16a99f2915c8439415e595c3f41348034e3` | match | near/far decomposition of the neighbouring-Cys statistic, both filters | rice -2.61 / +3.34, 689 sites (second store) |
| `results/09_mapped_sno.json` | 1,687,624 | `15f23a85cca1e0043601c6e409be40a0655cc706d7fd0302b271b93604483cc3` | - | mapped Arabidopsis S-nitrosylation site table written by 09, read by 10, 12, 16, 21, 36 | inputs of Artifact 1 |
| `results/09_mapped_so.json` | 1,052,532 | `a9d5ac16292a9052f961fa02a09cf429b7dd944252b30a19048a9042a634ff7d` | - | same, S-sulfenylation | inputs of Artifact 1 |
| `results/09_mapped_ox.json` | 916,011 | `a551398bd53f9540cff14f127ad477f1bf83788a7bdfe205b5685e35da1a2919` | - | same, reversible cysteine oxidation | inputs of Artifact 1 |
| `results/09_mapping_summary.json` | 361 | `a31da92f71be7a332d468cf262f4fefe041a65d250e6de323e64dac370acd0c2` | - | protein and site counts of the three mapped tables | - |
| `results/10_large_scale_motif.json` | 5,647 | `4126dcc66a6187d20e5bfa7af4bc79d566e7ec5b943fa06d4c37f6a01eeb60fe` | - | raw z (duplicate store) and the record `{type: sno, class: polar STNQ, p_perm_4000: 0.7435}` | 4,000; 0.7435 |
| `results/12_matched_sno.json` | 3,069 | `1de687bfd409fc22c2482b8b99776dd429c61c533a89b69b8e090177233c4b2d` | match | S-nitrosylation raw/matched z | +15.72, +8.22 |
| `results/12_matched_so.json` | 3,144 | `dd08e651f216be72f333e3291e0b4e6820f13cbf6db3ac79376afb783ce7c97e` | match | S-sulfenylation raw/matched z | +10.03, +4.71 |
| `results/12_matched_ox.json` | 3,206 | `716d1c8b35881d24a635c940832626d1fcc65e49eb89bd3a178598b669f1ea8e` | match | reversible oxidation raw/matched z, DE control, rule constants | +5.59, +0.67 (P 0.50), +6.04, +6.62, 7, 30, two |
| `results/14_accession_map.json` | 84,298 | `d60b9cd7a92c958baafacc055f5dcda885475d8bef51a24b471625b123a826d9` | - | Arabidopsis locus to UniProt accession map (input of 36) | - |
| `results/15_structural_features.json` | 4,256,659 | `2bca65b0326720ca9606686619c9c0e8b737f9c3550f68b66caba52fcefd2b4b` | - | per-cysteine AlphaFold features of the Arabidopsis proteins (input of 21) | - |
| `results/19_mapped_cre.json` | 618,403 | `71e217801c68669ee59dce229a00dd7516be95803f986d5186c3f64c4371e37e` | - | Chlamydomonas site tables (input of 21 and 36; not printed) | - |
| `results/19_mapped_dbptm.json` | 2,884,302 | `2e6fa070a58be0b2c39c6c6b6e2cb1b27c99d7c04f114950155c9041bdf1feef` | - | human + mouse dbPTM S-nitrosylation and S-glutathionylation sites mapped to UniProt sequences (derived table; input of 21, 36 and of the seven-dataset comparison) | inputs of Artifact 2 |
| `results/21_crossspecies_summary.json` | 42,221 | `34963cbdfd5b790aac9abdd59a7f54487ca846be424072e4db4cc529967aa6e0` | match | all cross-dataset tests; `rows[7]` is the rice neighbouring-Cys row | -2.61, P 0.008962, 756 |
| `results/21_crossspecies_summary.csv` | 12,697 | `09e99f87b1ae5fb04fc279bc730401d11f80fe14984b62eb3726ae2870076707` | match | same, CSV (CRLF) | same |
| `results/21_features_cross.json` | 9,075,284 | `a0535375b4d762a0cf9adcf2220ed230e117c215562c56d4a83d555ac7acc5aa` | - | per-cysteine AlphaFold feature cache of 21 (lets 21 run without the models) | - |
| `results/24_extra_sequences.json` | 30,958 | `dfd2ff86212384c9ff1b9d0d5f79b61f38fa9a29120f72c18e226b8888635fa5` | - | rice sequences fetched by accession outside the reference proteome (input of 24) | - (see Known gaps 3) |
| `results/24_mapped_osa_persulfidation.json` | 399,100 | `5c703885bab2df52d58737448213dbc04d52b138edca35006b7cffdd41b52252` | - | rice sites from single-cysteine peptides (input of 21 and 36); byte-identical to the rice self-audit cohort `PXD072089_mapped_rice.json` | inputs of Artifact 2 and of the self-audit |
| `results/24_mapped_osa_allpeptides.json` | 430,429 | `6e25f95aecd2087744b315a8ea990a3be7aed0637db029603a782da93232e11e` | - | rice sites from all cysteine-containing peptides (control) | - |
| `results/24_peptide_stats.json` | 168 | `8ca6c29097c797e8f07f1c775e30f234cdce116f6a8f7960fd6a2a6f98d8b391` | - | peptide counts of 24 | - |
| `results/27_mapped_palmitoylation.json` | 1,783,433 | `8db2602f0888bfb514c088628cfc10ca9e1a17b7390372f059b8705434d10b0b` | - | human + mouse dbPTM S-palmitoylation sites (input of 21 and 36; not printed) | - |
| `results/33_mapped_hsa_persulfidation.json` | 1,289,888 | `9e707e7e38693fff1d6e98d7eb6dec1e26764f18b6cb98e7431b1517b677d9bf` | - | human persulfidation sites from PXD044043 single-cysteine peptides (derived table); byte-identical to the human self-audit cohort `PXD044043_mapped_human.json` | inputs of Artifact 2 and of the self-audit |
| `results/35_symmetric_ncys_test.json` | 1,883 | `78e5b8657d6c058ade29241b9a15c367aa7ec08ea185b57fd9da34d84a525665` | match | seven-dataset asymmetric/symmetric comparison | +3.34 (P 0.00082), 689; -8.78, +2.67; seven; -0.98 to -11.09; +4.53 to -3.64; four; three; two |
| `results/36_ncys_decompose.json` | 10,744 | `cddbf4dadd18c2d07d40015253e481599088bdfa8cadee2419a0e718da75e966` | - | decomposition (second store of the rice pair) | -2.61 / +3.34, 689 |
| `results/36_ncys_decompose.csv` | 3,550 | `93ace78d0b177fe8caa73c469c6b1a5d951d90029a5c5ca436cef9a28a7e5f36` | - | same, CSV (CRLF) | same |
| `results/36_ncys_split.json` | 1,588,894 | `7236a532a52e34c33a9cbe530b364b79f038336005de94ab566c047f1b2864af` | - | per-protein near/far neighbour counts (cache of 36; lets 36 run without the models) | - |
| `positional_profile/scripts/positional_kr_profile_public_2026-09-17.py` | 30,885 | `aed19b7594202a6b52d04b30b54f1be897cc6582414a88fc32b62e0540f8d0dd` | - | positional K/R profile around the site (from the internal analysis tree) | 0.4576 [0.3681, 0.5487]; 1.4975 [1.3572, 1.6421]; 1.3831 [1.2417, 1.5235] |
| `positional_profile/results/positional_kr_profile_public_2026-09-17.csv` | 27,426 | `a76681d022846810bef9e2941e21a9bfa45f9b209c116d33e42bf55acd930ccb` | - | the profile values (CRLF) | same |
| `positional_profile/results/positional_kr_profile_public_audit.json` | 7,181 | `1978737971c8ecfaf294fa49d8f4cc4c76d6cdaae3afc06fa4ca49323e63cbde` | - | audit of the profile (input hashes) | - |

The three `positional_profile/` files match the SHA-256 of the internal analysis tree's artifact manifest.
Their dataset, `natcomm2023_ath_sno` (doi 10.1038/s41467-023-39078-0, PMID 37277371), is a separate file
from the three PlantPTMViewer sets but comes from the study listed first among the leading contributors to
the PlantPTMViewer S-nitrosylation set, so it is not independent of that set.

## Re-running

The scripts expect the layout of the tree they came from: `scripts/`, `results/` and `data/` side by side.
That is this folder (`ROOT = artifacts12/`); `data/` holds the third-party inputs listed below, which you
download yourself. Python 3.9 or later: 10, 12 and 24 use the standard library only; 21 and 36 import 16
and 17 from `scripts/` and need numpy (21 additionally imports biopython, pydssp and propka, but only for
proteins missing from its feature cache).

    cd artifacts12
    python scripts/10_large_scale_motif.py       # reads results/09_mapped_*.json
    python scripts/12_detectability_matched.py   # reads results/09_mapped_*.json
    python scripts/24_build_persulfidation_table.py   # needs data/persulfidation/SS-all-peptides.tsv, data/uniprot/osa_proteome.tsv
    python scripts/36_ncys_decompose.py          # uses the cache results/36_ncys_split.json
    python scripts/21_crossspecies_test.py       # uses the cache results/21_features_cross.json

Checked on 2 October 2026 in a copy of this folder with Python 3.9.6 and numpy 2.0.2, after deleting the
outputs being regenerated:

- 10 and 12: every value reproduced to a relative difference of at most 1.75e-14 (no printed value changes).
- 24 and 36: outputs byte-identical, including `36_ncys_split.json`.
- 21: the rice rows (z -2.61, P 0.008962, 756 sites) reproduce exactly. Six Arabidopsis rows (predicted
  disulfide and pKa without disulfides) differ, and with them the Benjamini-Hochberg q values of all rows,
  because 21 adds disulfide flags to the Arabidopsis proteins from the AlphaFold models in `data/af_large`,
  which are not redistributed (AlphaFold DB v6 models, fetched by accession). No printed number depends on
  those rows.
- 09 cannot be run as written (Known gaps 1).

## Inputs that are not redistributed

| input | used by | file name, date, SHA-256 |
|---|---|---|
| PlantPTMViewer cysteine exports (S-nitrosylation, S-sulfenylation, reversible cysteine oxidation) | 09 | `data/ptmviewer/sno_allSPECIES.csv` `d1ff195ea110c6220894f5060e7517c0b357aa8021c1c15aa44fcd6c5fb29800`; `so_allSPECIES.csv` `aad27d68bd8e8fbeffec40c0284916fe00e348448dbe2de5d385745d201d9633`; `ox_allSPECIES.csv` `3e15a85ca68785238240256c031f5c082102fc2e3cd502d998ef097dce81ce3b` (file dates 2026-09-10; the same values are in `MANIFEST.json`) |
| dbPTM experimentally supported site files (release not recorded) | producer of 19 (not retained) | `S-nitrosylation.gz` `ca74105632085a69b9e8580be7efa5617db973c564e3196b4aee1a9098682f29`; `Glutathionylation.gz` `f2859d86fec1c8d83a63abdec35c25c447b839e8eda83202122d5a4345751cb1` (both 2026-09-10); `S-palmitoylation.gz` `b7c1328ecbef2b9e0f9d8a217b95d10191f74741a54dab28d3e39b0220df0913` (2026-09-11) |
| PRIDE PXD072089 `SS-all-peptides.tsv` | 24 | `57c4ce64b5747eb591e879565916c04e0901ea228c20e511b2e1804fb384890c` (file date 2026-09-11; identical to the copy re-downloaded on 2026-09-30) |
| PRIDE PXD044043 MaxQuant archives | producer of 33 (not retained) | `txt_BTA.rar` `36381fb672bf75d892138d46e214f53b610473b1fbf19e5a907dda214110f149`; `txt_input.rar` `eccae06ebeee08898fe5d48f322bc0daa24b8a8aefc58336a5bb9cb6216a665f` |
| UniProt rice proteome table (columns accession, sequence; release not recorded) | 24 | `data/uniprot/osa_proteome.tsv` `8bf3d273236159460b28211351a3a484085a6719a6f55dfd83022a64f169b521` (2026-09-11) |
| AlphaFold DB models | 21 (Arabidopsis disulfide flags only, given the deposited caches) | `data/af_large/`, `data/af_cross/` |

## Known gaps

1. **09 cannot be re-run as written.** It reads an Arabidopsis proteome table from a temporary directory
   (`/private/tmp/.../scratchpad/ptmv/ath_proteome.tsv`, line 10) that no longer exists, and the
   UniProt/TAIR release of that table was not recorded. The mapped site tables 09 produced
   (`results/09_mapped_{sno,so,ox}.json`) are deposited, so 10, 12 and every downstream step re-run from
   them. Note added 2 October 2026: a UniProt table with the same columns,
   `data/uniprot/ath_proteome.tsv` (17,526,362 bytes, SHA-256
   `6e628e34c3e893e7efb7ba694ceb1bb4aab278e709b7c573cfbcc01558e19a5a`, file date 2026-09-10, release not
   recorded), was found in the collaborator tree; with line 10 pointed at it, 09 reproduces the three
   deposited mapped tables and the summary byte for byte. It is not redistributed here.
2. **The script that produced `results/35_symmetric_ncys_test.json` was not retained.** It is in neither
   the collaborator tree nor the internal analysis tree. The output is deposited with the SHA-256 recorded
   at staging, and the values it alone holds (human -8.78 and +2.67; the count of seven datasets; the ranges
   -0.98 to -11.09 and +4.53 to -3.64; four positive and three negative symmetric estimates) are reported as
   stored. The rice pair (+3.34, 689) is also stored in `36_ncys_decompose.json`, whose script re-runs
   byte-identically. The inputs such a script reads (the seven site tables `24_mapped_osa_persulfidation`,
   `33_mapped_hsa_persulfidation`, `09_mapped_{sno,so,ox}`, `19_mapped_dbptm`, and the neighbour cache
   `36_ncys_split`) are all deposited.
3. **No producing script for `results/24_extra_sequences.json`** (rice sequences fetched by accession for
   entries outside the reference proteome). The file is deposited; 24 reads it.
4. **No producing scripts for `results/19_mapped_dbptm.json`, `19_mapped_cre.json`,
   `27_mapped_palmitoylation.json`, `33_mapped_hsa_persulfidation.json`, `14_accession_map.json` and
   `15_structural_features.json`** (the collaborator tree has no scripts numbered 19, 27 or 33; 14 and 15
   exist there but need the AlphaFold models). The tables are deposited as the inputs of 21, 36 and the
   seven-dataset comparison.
5. An earlier draft quoted z = +2.54 (P = 0.011, 980 sites) for the unfiltered rice comparison; no stored
   table holds it and it is not used (recorded in `MANIFEST.json`).

## Licence

Code: MIT. Derived tables: CC BY 4.0. The dbPTM-derived tables (`19_mapped_dbptm.json`,
`27_mapped_palmitoylation.json`) and the PlantPTMViewer-derived tables (`09_mapped_*.json`) restate
site positions from those databases; the databases remain under their providers' terms.
