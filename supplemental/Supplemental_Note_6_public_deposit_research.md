# Supplemental Note 6. Re-searching PXD015307: what a public deposit can and cannot decide about sulfide versus dioxidation

The question, reading rule and outcome branches were registered before any search; search settings were amended after runs that failed on configuration (section 3). Sections 5 to 8 are post hoc and carry no verdict.

## 1. The deposit
PXD015307 [42] holds five raw files.

Fig-1D is the proteome-scale enrichment. Its persulfide sites come from a trapping-and-tagging chemistry read out with iodoTMT. The deposit's own search of this file (Fig-1D.msf: 25,014 peptides, 49,522 spectrum headers) placed:
- oxidation, on eleven residue types including cysteine;
- iodoTMT6plex, on C, D, E, H and K;
- one mass-shift entry.
It never placed sulfide or dioxidation.

The four Fig4-L files are named for NaHS, CTH with cysteine, heat-inactivated CTH and a control. The deposit's identifications in them are insulin chains only, consistent with an in vitro insulin experiment; the original paper was not re-read for this note. The deposit's search of these files returned 53 PSMs:
- all on human insulin chains: A chain GIVEQCCTSICSLYQLENYCN, 4 PSMs; intact B chain FVNQHLCGSHLVEALYLVCGERGFFYTPKT, 49 PSMs;
- 52 carrying N-ethylmaleimide;
- 14 in the CTH-with-cysteine file, 30 in the NaHS file, 9 in the heat-inactivated-CTH file and none in the control file;
- four carrying sulfide: two A-chain PSMs in the CTH-with-cysteine file and two B-chain PSMs in the NaHS file.

The consensus files carry no processing-node parameters. The modifications searched and the precursor tolerance therefore cannot be recovered; only the modifications placed can.

## 2. Search parameters (Table S6.1)
| item | setting | source |
|---|---|---|
| Raw files | Fig-1D.raw 498,563,222 B; Fig4-L_NaHS.raw 99,612,094; Fig4-L_CTH_Cys.raw 123,832,868; Fig4-L_heat_inactive_CTH.raw 150,754,884; Fig4-L_control.raw 102,688,780 (sizes match the PRIDE list); SHA-256 Fig-1D 1065817560a4768ddad8776750d6f41cb162b35a40ae2cbd1590d8ca9aaed14a, Fig4-L_NaHS 676abf25ac6ab0fc53769930eecd65b3fa09692c497071d87b46509bc42e8387, Fig4-L_CTH_Cys 6cf1bbb2ac14a2d588f4a091fc27b83baca8cbe1cb298097266e172565fcac28, Fig4-L_heat_inactive_CTH 41f3ead53dec8b4a521ce3b914c932502a843956c3c0931277f7f78f5c10b9aa, Fig4-L_control ddb40f5536280b900a84e994900aa2f344e8f06139dd9b85875a001075fd7da3 | plan; round report; checksums computed 2026-10-02 on the search host's copies |
| MS2 analyzer | Fig-1D FTMS MS2 (31,450 FTMS, 3 ITMS scans); Fig4-L files ion-trap MS2 (Fig4-L_NaHS 1,080 FTMS, 5,004 ITMS); survey scans FTMS in every file | scan headers of converted files |
| Conversion | ThermoRawFileParser 2.0.0.dev (bioconda build h9ee0642_1; written as 2.0.0.0 in the mzML software element) | conda environment record; mzML header |
| Engine | Comet 2026.01 rev. 1 (e4f767c), 24 threads, CPU | Comet output headers |
| Database | Fig-1D: UniProtKB/Swiss-Prot Mus musculus, reference proteome UP000000589, canonical + isoforms, 25,750 entries, release 2026_03 (downloaded 17 September 2026), SHA-256 070c674974b8b1b0da593248cf6d1565d532207241c390491f79d6f89f99d19f (with insulin appended: 1ef08af50c8c3b0eaa6f32b39e80aebc2fb500498939adac62a315c8c91d94ee). Fig4-L files: the same plus bovine insulin P01317 (25,751) | Comet output headers; plan amendment |
| Species | mouse: all 19,303 protein annotations of the deposit's Fig-1D.msf are mouse entries | plan amendment |
| Contaminants | none appended (the search database holds exactly the 25,750 reference-proteome entries, plus P01317 for the Fig4-L files) | database entry count |
| Decoys | Comet internal, reversed peptides with the C-terminal residue retained, concatenated; prefix DECOY_ | output sequences |
| Enzyme | trypsin, after K/R not before P, fully specific, up to 2 missed cleavages | plan; output sequences |
| Precursor | +/-10 ppm; isotope_error 0 | plan |
| Fragment | Fig-1D: fragment_bin_tol 0.02, fragment_bin_offset 0.0, theoretical_fragment_ions 0. Fig4-L: fragment_bin_tol 1.0005, fragment_bin_offset 0.4, theoretical_fragment_ions 1 | plan amendment; comet_hires.params (Fig-1D) |
| Static modifications | none (Comet's default add_C_cysteine disabled after run 1) | plan; amendment 6 |
| Variable modifications | C: sulfide +31.972071, dioxidation +31.989829, iodoTMT6plex +329.226595, carbamidomethyl +57.021464. M: oxidation +15.9949. At most 3 per peptide | report |
| Not searched | protein N-terminal acetylation; N-ethylmaleimide; iodoTMT6plex on D/E/H/K | plan amendment |
| Output | 5 candidates per spectrum | plan amendment |
| FDR | rank-1 PSMs, Comet E-value, decoys/targets at or below a threshold, largest threshold with ratio <= 0.01, per file, PSM level; no Percolator, no peptide- or protein-level control | reading script |
| Precursor error | (experimental - calculated)/calculated neutral mass x 10^6; baseline = 1%-FDR PSMs without a +32-class modification | reading script |
| Recovered from the Fig-1D parameter file (comet_hires.params, timestamp identical to the Fig-1D search start) | ion series b and y (no neutral losses); precursor charge from the spectrum, maximum 6; maximum fragment charge 3; peptide length 5-50; MH+ mass range 600-5,000 | parameter file on the search host |
| Outputs | Comet result tables, sha256: Fig-1D c07e7d329af2fb800d42dfc9e2dcac78caf1b6d96324e72c1ecf50f872ec7378; Fig4-L_CTH_Cys 088c995192ba864f89f1e5536fc81a2fe090b050f42f96806488b4061a0e03e9; Fig4-L_NaHS e0157a6c9919d3924e11564295d958dadd14f014955e396e76e939f36287333a; Fig4-L_control 5b2803fef02dba7b3e40c24497fe6a9a1aa069c95f72b80e1be790b4063bb701; Fig4-L_heat_inactive_CTH e5413073915c3cae218986c363c4ac194fe4cb5b5cdc0d09efea8465aaa26aad | reading audit |

## 3. Planning record
Before any search was run, the registration fixed:
- the question and the resolvability rule;
- the primary and secondary readouts;
- a minimum of twenty +32-class identifications for a rate, and four outcome branches;
- a search plan (bovine database, 0.02-Da fragment bins for all files, protein N-terminal acetylation among the variable modifications).

Seven amendments followed:
0. the compute host;
1. protein N-terminal acetylation dropped;
2. isotope_error 0 and at most three variable modifications;
3. five output lines;
4. the mouse instead of the registered bovine database, after the bovine run returned 174 PSMs at 1% FDR;
5. fragment settings per file;
6. four discarded runs, kept with their reasons:
   - Comet's default static carbamidomethyl term left on (a configuration error, not a change of plan);
   - the bovine database;
   - high-resolution fragment settings on the ion-trap files;
   - ion-trap settings on the high-resolution file.

Two further runs are recorded in the round report but not among the registered amendments: a fifth run of the insulin files without insulin in the database, after which bovine insulin was appended, and a wide-window diagnostic search stopped before completion. No amendment changed the question, the reading rule or the outcome branches.

The plan's timestamps are labeled UTC (00:40 to 02:00 on 17 September). They agree with the Comet run log (02:13 to 02:19 host time on 17 September) and with the reading audit (18:20 UTC on 16 September) when read as local times of the host, UTC+8.

The reading rule measures each error from the file's own baseline. The earlier consensus-file audit measured it from 0 ppm, so the two are not the same rule.

Two stored records still describe the discarded bovine run: the reading script's documentation and the audit field search_ran_on. The Comet headers in the same audit name the mouse databases that were actually read.

## 4. Results at 1% FDR (Table S6.2)
| file | MS2 | rank-1 PSMs | E-value threshold | PSMs at 1% FDR | +32-class Cys PSMs | baseline median (SD), ppm |
|---|---|---:|---:|---:|---:|---|
| Fig-1D | FTMS | 24,727 | 0.000932 | 245 | 0 | -3.01 (2.09) |
| Fig4-L_CTH_Cys | ion trap | 2,177 | 0.00292 | 14 | 0 | +5.55 (0.41) |
| Fig4-L_NaHS | ion trap | 1,896 | 0.00406 | 18 | 0 | +5.11 (2.04) |
| Fig4-L_control | ion trap | 726 | 0.00836 | 26 | 0 | +5.58 (0.59) |
| Fig4-L_heat_inactive_CTH | ion trap | 2,430 | 0.00549 | 18 | 0 | +4.75 (1.05) |

No file reaches twenty +32-class identifications, so by the registered rule the round is underpowered and no rate is quoted.

In the insulin files no decoy was accepted, so their counts are PSMs above the best-scoring decoy; with 14 to 26 targets the estimator cannot resolve 1%. The identities of these PSMs are not stored.

The Fig-1D yield is of the same order as the deposit's own 317 highest-confidence peptides.

## 5. Scope of the re-search
The deposit's four sulfide PSMs could not be re-matched, for three reasons:
- they are on human insulin chains, which differ from the bovine chains in the database at A8, A10 and B30;
- three also carry N-ethylmaleimide, which was not searched;
- the intact B chain is not a fully tryptic peptide of the precursor.

The planned re-examination of those four PSMs was therefore not achieved. Among spectra whose best-scoring candidates carry a +32-class modification, the only per-spectrum output stored, insulin appears only as eleven matches to the bovine A chain, all below the threshold (best E-value 0.046).

In Fig-1D the persulfide identity is inferred from the enrichment chemistry. A +32-class count of zero is therefore expected under the deposit's design and says nothing about its site calls.

## 6. The deposit's own four sulfide PSMs (post hoc; Table S6.3)
| file | chain | charge | placed | error (ppm) | file baseline, median (SD) ppm: reading | A- or B-chain baseline, median ppm (n): reading | across eight baselines |
|---|---|---|---|---:|---|---|---|
| CTH with Cys | A | 2 | NEM; sulfide | -6.50 | +5.77 (1.04): neither | -0.90 (2): neither | neither |
| CTH with Cys | A | 2 | sulfide | +4.60 | +5.77 (1.04): sulfide only | -0.90 (2): dioxidation only (file SD) or neither (chain SD) | depends on the baseline |
| NaHS | B | 4 | NEM; sulfide | +4.42 | +4.61 (2.19): sulfide only | +5.52 (47): sulfide only | sulfide only |
| NaHS | B | 4 | NEM; sulfide | +7.75 | +4.61 (2.19): either | +5.52 (47): either | either |

Readings are defined as follows:
- "sulfide only": the error lies within 2 SD of the sulfide expectation (the baseline) and not within 2 SD of the dioxidation expectation (baseline plus the separation: 7.0 ppm for the A chain, 5.0 ppm for the B chain);
- "dioxidation only": the converse;
- "either": within 2 SD of both;
- "neither": within 2 SD of neither.

The eight baselines are:
- the file's non-sulfide PSMs, as median with SD, as median with MAD-scaled SD, and as mean;
- the pooled median of all 49 non-sulfide PSMs;
- the same chain in the same file;
- the same chain in all files, with the chain's own SD or with the file SD;
- the same chain and charge in all files.

Not every baseline exists for every PSM: the CTH-with-cysteine file has no A-chain PSM without sulfide, and no A-chain PSM without sulfide has charge 2.

The deposit's two A-chain PSMs without sulfide (NaHS file, charge 3, m/z about 837) sit at -0.38 and -1.42 ppm. That is 6.4 ppm below the median of its 47 B-chain PSMs (+5.52 ppm; Mann-Whitney P = 0.019). B-chain ions at similar m/z (858 to 889; n = 21) sit at a median of +5.06 ppm, so the offset is not an m/z effect.

By the half-separation rule, the A-chain PSM at +4.60 ppm is resolvable as sulfide against the file baseline and as dioxidation against the A-chain baseline.

The deposit's tolerance is not recorded. With a window edge at +10 ppm, a dioxidized peptide reported as sulfide would have been expected at:
- +9.4 to +10.5 ppm for the B chain, inside or at the edge;
- +12.3 to +12.8 ppm (outside) or +6.1 ppm (inside) for the A chain, depending on the baseline.
The window can therefore remove part of the dioxidation-consistent error range, and the errors of reported sulfide PSMs are not independent evidence for sulfide.

The four PSMs illustrate why the check needs the run's own precursor-error distribution: their reading depends on the baseline chosen, and they are too few for a rate. Neutral masses, and hence separations in ppm, are recomputed approximations.

## 7. What precursor accuracy can decide (post hoc; Table S6.4)
| criterion | mass limit | theoretical mouse Cys peptides at or below (one +32; +32 with iodoTMT on the other Cys) |
|---|---|---|
| window, +/-20 ppm, unbiased | 888 Da | 2%; 2% |
| window, +/-10 ppm, unbiased | 1,776 Da | 34%; 30% |
| window, +/-4.5 ppm, unbiased | 3,946 Da | >99.9%; 98.7% |
| window, +/-10 ppm, Fig-1D offset -3.0 ppm: true sulfide; true dioxidation | 2,539 Da; 1,365 Da | 66%; 17% (one +32) |
| window, +/-10 ppm, insulin-file offsets +4.8 to +5.6 ppm: true sulfide; true dioxidation | 1,140-1,204 Da; 3,385-4,017 Da | 9-12%; 97-100% (one +32) |
| mass range, 2 SD each side, Fig-1D (SD 2.09 ppm, 245 PSMs) | 2,123 Da | 48%; 44% |
| mass range, 2 SD each side, insulin files (SD 0.41-2.04 ppm, 14-26 PSMs) | 2,181-10,863 Da | 51-100% |
| mass range, 3 SD each side, Fig-1D | 1,415 Da | 19%; 17% |

Window limits give the mass below which the alternative falls outside the window: 17,758/t Da for an unbiased measurement, and 17,758/(t + b) or 17,758/(t - b) Da with offset b when sulfide or dioxidation, respectively, is true. Mass-range limits are 17,758/(2z SD) Da.

Theoretical peptides: 516,143 unique tryptic cysteine peptides of the UniProt 2026_03 mouse reference proteome, 7 to 30 residues, up to 2 missed cleavages.

## 8. Fragment binning (post hoc; Tables S6.5 and S6.6)
A fragment that contains the modified cysteine moves by 0.017758/z in m/z between the two interpretations. The table gives the share of such b/y fragments whose Comet bin changes.

| binning | z = 1 | z = 2 | z = 3 |
|---|---:|---:|---:|
| 0.02 Da (Fig-1D), theoretical fragments of the 1,085 peptides of tied spectra | 0.89 | 0.45 | 0.30 |
| 1.0005 Da, offset 0.4 (insulin files), theoretical fragments of the 621 peptides of tied spectra | 0.0001 | 0.031 | 0.0001 |

The search output agrees. The sample is single-cysteine spectra whose best interpretation carries one +32-class modification. In these, the alternative chemistry was a scored candidate when its precursor error lay inside the +/-10-ppm window.

| | Fig-1D | insulin files |
|---|---:|---:|
| alternative outside the window | 977 | 453 |
| both scored, tied at the best score | 206 (30%) | 83 (83%) |
| both scored, alternative scored lower | 476 (70%) | 17 (17%) |
| within 0.3 ppm of the window edge (set aside) | 39 | 3 |
| tied among those with best E-value <= 1 | 6 of 60 | 21 of 26 |

With 0.02-Da bins, the fragment score is sensitive to the 0.0178-Da difference: the two candidates tie mainly when no fragment containing the modified cysteine is matched. This is more frequent among poorer matches (43% tied at best E-values above 100). With 1.0005-Da bins, the score is almost insensitive.

Most of these spectra lie below the identification threshold, so the comparison shows that the score responds to the shift, not which chemistry is correct; the contrast between the binnings holds among spectra with a best E-value of 1 or less (6 of 60 tied against 21 of 26).

## 9. Phosphorylation and sulfation: a precursor-mass check (Table S6.7)
| deposit | status | main-search tolerance | phosphopeptides | resolving mass M* | share at or above M* |
|---|---|---|---:|---:|---:|
| PXD071110 | assessed | 4.5 ppm | 33,137 | 2,115 Da | 0.479 |
| PXD071012 | assessed | 4.5 ppm | 27,372 | 2,115 Da | 0.427 |
| PXD068672 | assessed | 4.5 ppm | 21,646 | 2,115 Da | 0.246 |
| PXD068395 | assessed | 4.5 ppm | 23,153 | 2,115 Da | 0.213 |
| PXD071074, PXD071061 | not fetched (combined files 213 and 277 MB) | - | - | - | - |
| PXD068397, PXD064206 | phosphorylation not searched | - | - | - | - |

Configuration was read from each deposit's own mqpar.xml, and for PXD071110 also from its two Andromeda .apar files. The parameters.txt files state neither the modification list nor the tolerance.

None of the four assessed deposits searched sulfation. At 2 ppm, M* is 4,758 Da, above the heaviest phosphopeptide of each deposit (4,283 to 4,572 Da). The PXD071110 share was recomputed with the released tool (0.479343).

The check uses only the tolerance, the modification list and the peptide masses. It uses no fragment evidence and no measured mass errors, so it shows what the searches could not report, not that any identification is wrong.

## 10. Limitations
- The engine differs from the deposit's (Comet against Sequest HT), and no engine control was run.
- The database release differs from the authors'.
- The insulin files were searched with a bovine insulin sequence, without N-ethylmaleimide and with a fully tryptic enzyme, whereas the deposit identified intact human chains.
- No Percolator rescoring was applied.
- The deposit's A-chain baseline rests on two PSMs.
- Sections 5 to 8 are post hoc.

## 11. Files
Stored outputs:
- results/pxd015307_research_summary.csv
- results/pxd015307_research_psms.csv
- results/pxd015307_research_audit.json
- results/pxd015307_posthoc_score_ties.csv and its audit
- results/search_space_mass_accuracy.csv and its audit

Post hoc tables: posthoc_2026-09-30/results/H_artifact5_docs/. provenance.json lists every input and output with its sha256.
