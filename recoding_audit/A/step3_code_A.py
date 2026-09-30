"""Step 3 of the independent recoding protocol (coder A).

Coder A's independent codes for the 74 full-text-classified records.

Method notes (recorded before joining to the original table):

Decision order applied to every row, using ONLY the blinded extract
(title, modification, study_type, evidence_quote_background,
evidence_quote_detectability, ambiguity, note, verification):

  1. b  - the evidence explicitly matches the background on detectability
          (ionisation/MS-observability) or on protein abundance.
  2. c  - the negative set is residues *detected in the same experiment*
          that carry no modification assignment.
  3. a  - the evidence states the background is the complete, unrestricted
          set of residues of that type (whole proteome / all proteins /
          "all protein sequences" / all peptides of a reference database).
  4. e  - anything else, including "no background statement found".
  5. U  - the recorded extract is non-informative about any background and
          does not even permit the inference "no statement found".
          (Used for 1 of 74 rows; the truncation of a quotation is treated as
          "no statement retrievable" and therefore coded e, per the rubric's
          own wording for e.)

Strictness rule applied as instructed: "regulated vs all identified sites"
contrasts, externally-curated negative sets, and computational predictors
whose negatives were not generated in the same experiment are coded e (or a
if the evidence states an unrestricted whole-proteome/all-residues
background) - never c.

`coder_code` is coder A's single code per row; `code_A_strict` is the same
code re-derived under the strictness rule above (they coincide for every
row; the strict rule was applied throughout, not as a second pass).
"""

import csv
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

BLINDED = r"C:\Users\admin\Desktop\MCP\_recoding\A\blinded_74.csv"
OUT = r"C:\Users\admin\Desktop\MCP\_recoding\A\codes_A.csv"

# uid -> (code, one-sentence justification, confidence)
CODES = {
    1: ("e", "Computational predictor whose negatives are only called 'nonsulfhydration samples' with no stated provenance; not same-experiment residues and not stated as an unrestricted whole-proteome set.", "medium"),
    3: ("e", "Computational predictor; the negative description is truncated at 'while all', so an unrestricted all-cysteine proteome background cannot be established and the predictor default applies.", "low"),
    9: ("e", "Database resource in which only a foreground FASTA input is described; no background or negative set for the resource is named anywhere in the extract.", "medium"),
    36: ("e", "Predictor negatives are the remaining cysteines of the positive proteins only, i.e. a substrate-restricted complement, not an unrestricted proteome background and not same-experiment detection.", "medium"),
    43: ("e", "Predictor that labels 'all other samples' negative without stating their universe; provenance unstated so neither unrestricted-all-residues nor same-experiment status can be affirmed.", "low"),
    52: ("e", "Predictor; the negative sentence is truncated at 'all cysteine sites that', so the complement set cannot be confirmed as unrestricted proteome-wide.", "low"),
    102: ("a", "Negatives are randomly selected cysteines excluding known S-glutathionylation sites, i.e. an unrestricted sampling of all residues of that type with no detectability or abundance matching.", "medium"),
    128: ("e", "The analysis uploaded 586 reduced-cysteine peptides to pLogo but the extract names no background set for the motif comparison.", "medium"),
    196: ("a", "Explicit statement: a background set was constructed from 'all cysteine sites in human proteome'.", "high"),
    208: ("e", "Negatives are 810 'experimentally verified non-SNO sites' taken from an external curated database rather than generated in the same experiment.", "high"),
    211: ("e", "Predictor; the 2,581 negative sites are attributed only to '327 unique substrates' with no stated derivation rule.", "medium"),
    218: ("e", "Predictor; 478 negative samples are reported without any statement of how or from where they were drawn.", "medium"),
    441: ("c", "The negative group is the non-palmitoylated fraction of the same acyl-RAC dataset, i.e. proteins detected in the same experiment that carry no palmitoylation assignment.", "medium"),
    480: ("c", "Native MS with no protease digest means every cysteine of the analysed proteins is directly observed, so the comparison background is the same-experiment detected (unmodified) cysteine set.", "low"),
    713: ("e", "Motif-X is described with its +/-10 amino acid window but the sentence is cut off at 'with the', so no background set is retrievable.", "low"),
    736: ("e", "Negatives are annotated acetylation sites lacking low-throughput evidence, i.e. curated database annotations, not residues observed unmodified in the same experiment.", "high"),
    741: ("e", "Only flanking-sequence analysis of the identified Klac sites is described; no background set is named.", "medium"),
    769: ("e", "Soft MoMo is cited with a truncated template that never reaches a background specification.", "low"),
    800: ("e", "Only the MoMo tool and its version are cited; the extract names no background set for the motif analysis.", "medium"),
    802: ("e", "The extract states 'background matches' were referenced but the background set is never defined.", "high"),
    826: ("e", "MoMo usage is quoted only up to 'sequences constituted with amino', so the background specification is not retrievable.", "low"),
    828: ("U", "The recorded background evidence for this record is a trypsin reagent sentence, which is non-informative about any background and does not permit the inference that no statement exists.", "low"),
    836: ("e", "Motif-X motif searching is mentioned without any statement of the background set used.", "medium"),
    857: ("e", "Two independent regex passes over the full text returned no background statement.", "high"),
    859: ("e", "Two independent regex passes over the full text returned no background statement.", "high"),
    861: ("e", "rmotifx motif prediction is cited without any description of the background or reference set.", "medium"),
    877: ("e", "Soft motif-x is quoted with a truncated window description and no retrievable statement of the background.", "low"),
    885: ("e", "Context sequences around the acetyl sites were extracted but no background set for the comparison is named.", "medium"),
    892: ("e", "Motif analysis was run on differentially expressed phosphopeptides with no background set named.", "medium"),
    902: ("e", "Motif-x 21-mer analysis is described but the quoted sentence does not reach any background specification.", "low"),
    924: ("e", "Predictor; positives are experimentally identified succinylated lysines and the negative clause is truncated, so an unrestricted all-lysine background is not established.", "low"),
    982: ("e", "The study used a WebLogo frequency logo with no background set, i.e. no background was defined at all.", "medium"),
    989: ("a", "Foreground and background N-mers were both extracted from the protein-aligned sequences, i.e. the background is the complete set of N-mers of the analysed proteins with no detectability restriction.", "medium"),
    994: ("e", "Only Motif-X parameter settings are quoted; no background set is retrievable from the extract.", "low"),
    1005: ("e", "Two independent regex passes over the full text returned no background statement.", "high"),
    1022: ("e", "The MoMo/MEME function is cited without any statement of the background used for motif enrichment.", "medium"),
    1023: ("a", "Explicit statement that the Arabidopsis IPI database provided by Motif-X was used as the background, i.e. an unrestricted whole-proteome reference.", "high"),
    1067: ("a", "The background set was built by extracting all S/T/Y residues from the background protein set, i.e. every residue of that type with no detectability or abundance matching.", "medium"),
    1077: ("e", "Motifs were extracted from phosphopeptides but no background set is named.", "medium"),
    1111: ("e", "Foreground is regulated sites and background is all identified phosphosites from the same study, i.e. a contrast between two subsets of sites rather than unmodified residues, which the strictness rule places in e.", "high"),
    1114: ("e", "Motif-X is applied to enriched phosphopeptides with no retrievable statement of the background set.", "low"),
    1135: ("e", "Predictor; the extract states only that a window of 2n+1 was used for positives, leaving the negative-set source unstated.", "medium"),
    1253: ("e", "The 'background peptides' mentioned are non-target enrichment peptides and the motif-x background is unstated, so no residue-level background is described.", "medium"),
    1380: ("e", "Motif-X sequence logos are mentioned without any statement of the background set used.", "medium"),
    1537: ("e", "The abundance analysis in this study was applied to sites-per-protein, not to a sequence-level background, so it does not constitute an abundance-matched background.", "medium"),
    1553: ("e", "Predictor negatives were randomly drawn from identified human peptides held in an external database (MAPU), i.e. not from the same experiment.", "medium"),
    1558: ("e", "Soft MoMo is cited for lactylation motifs with no retrievable background specification.", "low"),
    1568: ("a", "Explicit statement that Motif-X analysed sequence models in all protein sequences, i.e. the unrestricted proteome.", "medium"),
    1593: ("a", "Explicit statement that modifier-21-mer motifs were analysed in all protein sequences.", "medium"),
    1596: ("e", "Predictor whose negatives are 'non-succinylated sites' that were not generated in the same experiment and are not stated to be an unrestricted all-lysine set.", "medium"),
    1608: ("e", "Motif-x is quoted with a truncated sentence that does not reach a background specification.", "low"),
    1609: ("e", "Two independent regex passes over the full text returned no background statement.", "high"),
    1614: ("e", "A complete sentence citing MoMo/motif-x over 'appraised proteins' that never names a background set.", "medium"),
    1620: ("e", "Soft motif-x is quoted with a truncated window template and no retrievable background statement.", "low"),
    1628: ("e", "The quoted passage concerns omission of ungroupable motifs and names no background set.", "low"),
    1636: ("e", "A complete sentence citing MoMo/motif-x on modification sites that never names a background set.", "medium"),
    1639: ("e", "Soft MoMo is quoted with a truncated template and no retrievable background statement.", "low"),
    1646: ("e", "MoMo Motif-X is cited for the acetylome but no background set is named.", "medium"),
    1674: ("e", "Soft motif-x is quoted with a truncated 21-mer template and no retrievable background statement.", "low"),
    1686: ("e", "The extract itself flags that the background set is not defined; 'motifs in the identified maize proteins' is not a background specification.", "medium"),
    1710: ("e", "Motif-X ubiquitination motif analysis is mentioned with no statement of the background set.", "medium"),
    1712: ("e", "Motif-X is named but the recorded sentence is cut off at 'which was' before any background set is specified, so no background statement is retrievable.", "low"),
    1723: ("e", "13-residue sequences centred on the phosphosite were extracted but no background set is named.", "medium"),
    1725: ("a", "The extract states the motif-x background was the SGD proteome, i.e. an unrestricted whole-proteome reference set.", "medium"),
    1726: ("e", "The quoted passage describes mutation enrichment over flanking sequences and never defines the reference region or background.", "low"),
    1729: ("e", "Only 'phosphorylation motifs predicted by Motif-X' is stated, with no background set named.", "medium"),
    1736: ("e", "IceLogo motif analysis is mentioned without any statement of the reference/background set.", "medium"),
    1742: ("e", "The quoted sentence concerns interpretation of ATR-regulated kinase action and names no background set.", "low"),
    1757: ("a", "Explicit statement that the Motif-X background was the IPI human proteome.", "high"),
    1789: ("e", "Computational predictor; the extract discusses decoy sites but never defines the comparison set used for conservation analysis.", "low"),
    1867: ("e", "Software tool paper that performs no enrichment of its own; it only recommends a background type for downstream tools, so no background used by this study is described.", "low"),
    1893: ("e", "Predictor whose negative-set source is explicitly noted as not stated in the extract; the quoted text is a decision formula.", "medium"),
    2038: ("e", "'Bioinformatics analyses' are invoked with no background specification of any kind.", "medium"),
    2113: ("a", "Background peptides were extracted from all length-13 peptide sequences of the UniProt human database, i.e. an unrestricted whole-proteome background.", "medium"),
}


def main():
    with open(BLINDED, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))

    missing = [r["uid"] for r in rows if int(r["uid"]) not in CODES]
    if missing:
        sys.exit("no code for uids: %s" % missing)
    extra = set(CODES) - {int(r["uid"]) for r in rows}
    if extra:
        sys.exit("codes for unknown uids: %s" % sorted(extra))

    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["uid", "doi", "coder_code", "code_A_strict", "justification", "confidence"])
        for r in rows:
            uid = int(r["uid"])
            code, just, conf = CODES[uid]
            w.writerow([r["uid"], r["doi"], code, code, just, conf])

    # distribution
    from collections import Counter
    dist = Counter(c for c, _, _ in CODES.values())
    print("codes_A.csv written: %d rows" % len(CODES))
    print("distribution coder A:", dict(sorted(dist.items())))
    print("strict (b or c):", dist.get("b", 0) + dist.get("c", 0))
    print("confidence:", dict(sorted(Counter(c for _, _, c in CODES.values()).items())))


if __name__ == "__main__":
    main()
