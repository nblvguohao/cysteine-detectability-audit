import csv, os, json

import os  # v3.1.1: paths made repository-relative (were absolute Windows paths of the authoring machine)
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
BLIND = os.path.join(_ROOT, 'recoding_audit', 'B', 'blinded_74.csv')
OUT = os.path.join(_ROOT, 'recoding_audit', 'B', 'codes_B.csv')

# Coder B (independent second pass). Codes assigned from blinded_74.csv ONLY.
# a = all residues of that type, unrestricted; b = matched on detectability/abundance;
# c = residues detected in the same experiment without a modification assignment;
# e = other, or no statement found; U = unclassifiable from the evidence retrieved.
CODES = {
1:    ("e", "Negative set described only as 'nonsulfhydration samples' with no stated source; not shown to be whole-proteome, same-experiment or detectability-matched.", "low"),
3:    ("U", "Retrieved quote is cut off mid-clause ('while all') before the negative set is described, so no class can be assigned.", "low"),
9:    ("U", "Retrieved quote is cut off inside the input specification, before any background descriptor appears.", "low"),
36:   ("e", "Negatives are the remaining cysteines of the proteins that carry the experimental SNO sites - a curated protein-restricted set, not the same experiment's unassigned detections.", "medium"),
43:   ("e", "Only 'all other samples were negative sites'; the origin of the negative samples is never given.", "low"),
52:   ("U", "Quote cut off at 'all cysteine sites that', so the qualifier that would decide between a and c is missing.", "low"),
102:  ("a", "'Randomly selected some cysteines except for known S-glutathionylation sites' describes an unrestricted cysteine pool drawn at random.", "low"),
128:  ("e", "Motif generation used 586 reduced-cysteine peptides; no background set is named anywhere in the retrieved evidence.", "medium"),
196:  ("a", "Explicit 'background set for all cysteine sites in human proteome' - an unrestricted whole-proteome cysteine background.", "high"),
208:  ("e", "Negatives are 810 'experimentally verified non-SNO sites', i.e. curated from prior experiments rather than generated in this study.", "medium"),
211:  ("e", "Training-set sizes are given but the origin of the 2,581 negative sites is never stated.", "medium"),
218:  ("e", "478 negative samples reported with no statement of where they came from.", "low"),
441:  ("c", "Negative group is the non-palmitoylated set detected within the same acyl-RAC dataset, i.e. residues detected in the same experiment without a modification assignment.", "medium"),
480:  ("c", "Native-MS experiment observes all cysteines of the analysed proteins and the recorded note states the background is that detected (unmodified) set.", "low"),
713:  ("U", "Quote cut off at 'with the' immediately before the Motif-X background specification.", "low"),
736:  ("e", "Negatives are acetylation sites annotated without low-throughput evidence - a curated annotation subset, not unmodified residues.", "medium"),
741:  ("U", "Quote concerns flanking-amino-acid analysis of identified Klac sites; no background set is described and the quote is cut off.", "low"),
769:  ("U", "Quote cut off inside the Motif-X parameter sentence before any background descriptor.", "low"),
800:  ("U", "Quote names the MoMo tool but is cut off inside its citation; no background set appears.", "low"),
802:  ("e", "The comparison is against undefined 'background matches' and the paper never defines the background set.", "medium"),
826:  ("U", "Motif tool named but quote cut off; the recorded note states a background statement exists in the full text but the statement itself was not recorded.", "low"),
828:  ("e", "Extracted 'background' sentence is a reagent purchase line (trypsin), and the recorded ambiguity confirms no per-residue background was defined.", "medium"),
836:  ("e", "Motif-X named but the quote is cut off; the only detectability evidence is an abundance remark recorded as discussion-only.", "medium"),
857:  ("e", "Recorded note states no background statement was found anywhere in the full text; motif analysis only.", "high"),
859:  ("e", "Recorded note states no background statement was found anywhere in the full text; motif analysis only.", "high"),
861:  ("e", "Only the motif-x/rmotifx software is named; no background set is stated.", "low"),
877:  ("U", "Recorded ambiguity states the motif-x background phrase was truncated in the extract, so the class cannot be recovered.", "low"),
885:  ("U", "Quote cut off mid-sentence ('extracted for') before any background description.", "low"),
892:  ("U", "Quote cut off at 'performed using' before the comparison and background are specified; the stated contrast is between two site subsets.", "low"),
902:  ("U", "Motif-x parameters described but quote cut off; the note records that an explicit background statement exists in the full text which was not recorded.", "low"),
924:  ("U", "Quote cut off at 'while the', exactly where the negative set would be described.", "low"),
982:  ("e", "Phosphosites in two groups visualised with WebLogo; the recorded ambiguity confirms no background set was used.", "medium"),
989:  ("a", "Background N-mers were extracted from the aligned protein sequences - an unrestricted set of all S/T/Y positions with no detectability restriction.", "medium"),
994:  ("U", "Only Motif-X default parameters are quoted; the note records that a background statement exists in the full text but it was not recorded.", "low"),
1005: ("e", "Recorded note states no background statement was found anywhere in the full text; only a general SP-motif enrichment remark.", "medium"),
1022: ("U", "MoMo/MEME named but quote cut off; the note records that an explicit background statement exists in the full text which was not recorded.", "low"),
1023: ("a", "Explicit statement that the Arabidopsis IPI (whole-proteome) database was used as the Motif-X background.", "high"),
1067: ("a", "Explicit statement that the background set was all S/T/Y residues extracted from the background protein set.", "high"),
1077: ("U", "Quote cut off mid-sentence while describing foreground motif extraction; no background set named.", "low"),
1111: ("e", "Background is all identified phosphosites while foreground is regulated sites - a regulation contrast between site subsets, which the strictness rule excludes from c.", "medium"),
1114: ("U", "Motif-X enrichment quote cut off mid-sentence; the note records that an explicit background statement exists in the full text which was not recorded.", "low"),
1135: ("e", "Only positive-set construction is quoted; the recorded ambiguity states the negative-set source is not stated.", "low"),
1253: ("e", "Quote concerns motif positions; the recorded ambiguity states the motif-x background is unstated (the recorded 'background peptides' means non-target peptides).", "medium"),
1380: ("e", "Motif-X version is named but no background set is stated.", "medium"),
1537: ("e", "Abundance was controlled for quantitation / sites-per-protein only; the recorded ambiguity confirms this is not a sequence-level background.", "medium"),
1553: ("a", "Negatives are 3,417 lysines randomly drawn from identified human peptides held in an external database (MAPU) - an unrestricted lysine pool, not same-experiment detections.", "low"),
1558: ("U", "Motif tool and 21-mers described but the quote is cut off before any background descriptor.", "low"),
1568: ("a", "Motif-X models were analysed over all protein sequences - an unrestricted all-proteins background.", "low"),
1593: ("a", "Motif-X models were analysed over all protein sequences - an unrestricted all-proteins background.", "low"),
1596: ("U", "Negatives described as 'non-succinylated sites which are those fragmented sequences' with the defining clause cut off.", "low"),
1608: ("U", "Motif-x parameters quoted but cut off; the note records that an explicit background statement exists in the full text which was not recorded.", "low"),
1609: ("e", "Recorded note states no background statement was found anywhere in the full text; only motif counts reported.", "medium"),
1614: ("e", "MoMo named and 'appraised proteins' mentioned, but no background set is defined.", "low"),
1620: ("U", "Motif-x 21-mer sentence cut off; the note records that an explicit background statement exists in the full text which was not recorded.", "low"),
1628: ("U", "Quote concerns excluding unconfident motifs and is cut off mid-clause; no background set is described.", "low"),
1636: ("e", "MoMo/motif-x named only; no background set is stated.", "medium"),
1639: ("U", "Motif-x sentence cut off mid-parameter list; the note records that an explicit background statement exists in the full text which was not recorded.", "low"),
1646: ("e", "Motif-X online tool named only; no background set is stated.", "medium"),
1674: ("U", "Motif-x 21-mer sentence cut off before any background descriptor appears.", "low"),
1686: ("e", "Quote cut off mid-sentence and the recorded ambiguity states the background set is not defined.", "medium"),
1710: ("U", "Motif-X named with no background given; the note records that an explicit background statement exists in the full text which was not recorded.", "low"),
1712: ("U", "Motif-X sentence cut off at 'which was'; the note records that an explicit background statement exists in the full text which was not recorded.", "low"),
1723: ("U", "Quote cut off mid-sentence describing extracted phosphosite sequences; no background set named.", "low"),
1725: ("a", "Recorded ambiguity states the motif-x background was the SGD (whole yeast) proteome.", "medium"),
1726: ("e", "Flanking-sequence analysis of variation datasets; the recorded ambiguity states the reference region is not defined.", "medium"),
1729: ("e", "Motif-X named only; no background set is stated.", "low"),
1736: ("e", "IceLogo described only in general terms; no background set is stated.", "medium"),
1742: ("U", "Quote cut off mid-clause and concerns ATR kinase action, not a background set.", "low"),
1757: ("a", "Explicit statement that the Motif-X background was the IPI human proteome - an unrestricted whole-proteome background.", "high"),
1789: ("e", "Quote concerns 'decoy' phosphosites and is cut off; the recorded ambiguity states the comparison set is not defined.", "low"),
1867: ("a", "Tool paper whose recorded ambiguity states it recommends a proteome-derived background and warns that species-mismatched backgrounds distort enrichment.", "low"),
1893: ("e", "Only the positive/negative assignment formula is quoted; the recorded ambiguity states the negative-set source is not stated.", "low"),
2038: ("e", "Only a generic statement that bioinformatics analyses were performed; the recorded ambiguity notes no background specification.", "medium"),
2113: ("a", "Explicit statement that background peptides were extracted from all length-13 peptides in the UniProt Human database - an unrestricted whole-proteome background.", "high"),
}

with open(BLIND, newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

assert len(rows) == 74, len(rows)
assert set(int(r["uid"]) for r in rows) == set(CODES), \
    set(int(r["uid"]) for r in rows) ^ set(CODES)

with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["uid", "doi", "coder_code", "code_B_strict", "justification", "confidence"])
    for r in rows:
        uid = int(r["uid"])
        code, just, conf = CODES[uid]
        # code_B_strict: the code after applying the strictness rule.
        # The rule was applied while assigning every code, so the two columns coincide;
        # no c-assignment required downgrading and no U row was resolved by guessing.
        w.writerow([uid, r["doi"], code, code, just, conf])

print("wrote", OUT, "rows:", len(rows))
from collections import Counter
c = Counter(v[0] for v in CODES.values())
print("coder_code distribution:", dict(c))
print("b or c (strict):", c["b"] + c["c"])
