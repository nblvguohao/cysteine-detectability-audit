"""Site-map the qPerS-SID supplementary tables. Records evidence; does not pick a label.

Source: Longen S, Richter F, Kohler Y, Wittig I, Beck KF, Pfeilschifter J.
Quantitative Persulfide Site Identification (qPerS-SID) Reveals Protein Targets of
H2S Releasing Donors in Mammalian Cells. Sci Rep 2016;6:29808, doi 10.1038/srep29808.
Gold open access; files and hashes in external/intake/scirep2016_qpers_sid/ and
results/qpers_sid_intake_audit.json.  Structural assessment:
reports/QPERS_SID_INTAKE_ASSESSMENT.md.

**Mapping is done by SEQUENCE, not by accession, and that is a decision not a
convenience.** The paper searched the full UniProt human set of 2015-04-08
(68,506 entries, isoforms and TrEMBL included); this project's reference is the
canonical Swiss-Prot human proteome in external/proteomes/hsa.fasta.gz (20,416
entries), which is what the v2 feature builder and the AlphaFold models are keyed
to.  Matching on accession loses 41.5% of the rows (1,086 of 1,855 accessions in
MOESM2/Beads are present).  So instead:

* tier A rows are located by their stripped peptide sequence,
* tier B rows by their 31-residue `Sequence window`, whose centre (index 15,
  0-based) is the modified cysteine.

A row is accepted only when its sequence occurs in EXACTLY ONE canonical protein
and the implied position really is a cysteine.  Rows that match nothing, match
several proteins, or land on a non-cysteine are counted and dropped, never
guessed.  The table's own accession is kept alongside as `accession_reported` so
the isoform-to-canonical collapse stays auditable.

**Adduct meaning depends on the fraction, so this script does not emit a
positive/negative label.** The enrichment works by eluting persulfides with TCEP
and then blocking the freed thiol with iodoacetamide, so the same adduct mass
means different things in different fractions:

    (Beads,       +414.19)  probe still on a free thiol          -> free thiol, same run
    (Beads,       +446.17)  probe on a persulfide                -> DIRECT adduct evidence
    (Total,       +446.17)  probe on a persulfide, whole lysate  -> DIRECT adduct evidence
    (Total,       +57.02)   iodoacetamide-blocked cysteine       -> free thiol, pre-enrichment
    (Total,       +414.19)  probe-labelled thiol, whole lysate   -> free thiol
    (Elution,     +57.02)   freed by TCEP, then blocked          -> ENRICHMENT-based persulfide evidence
    (Elution,     +414.19)  probe-labelled thiol in the eluate   -> carryover
    (Flowthrough, any)      not retained by streptavidin         -> unbound
    tier B (Elution, Carbamidomethyl C-Sites)          -> ENRICHMENT-based persulfide evidence
    tier B (Beads,   Iodoacetyl-PEG2-Biotin C-Sites)   -> free thiol, same run

Every row therefore carries `evidence_meaning` and the caller decides the label.
Cysteines present in a peptide but carrying no adduct are emitted as class
`C_unlabelled`; they were seen in the run but NOT seen by the probe, so they must
not be used as parallel-labelling negatives.

**Two evidence tiers, which reports/QPERS_SID_INTAKE_ASSESSMENT.md requires be
kept separate and never merged:**

* tier A, adduct-based: MOESM2 (main experiment, `Experiment` 1-2) and MOESM7
  (per-donor SILAC quadruplicate, `Experiment` L1/L2/H1/H2), both PEAKS 7
  searches that include +446.17.
* tier B, enrichment-based: MOESM3-6, one file per donor, MaxQuant quantification
  run with the SAME parameters but EXCLUDING +446.17, so a tier B positive means
  mean normalised H/L >= 1.30, significant in at least 3 of 4 experiments.

Outputs `results/qpers_sid_sites_normalised.csv` (long form, one row per
source/sheet/experiment/site/class), `results/qpers_sid_class_pairing.csv` (per
tier and donor: how many sites carry persulfide evidence, free-thiol evidence,
or both - the basis any negative set will be built from) and an audit.

Requires `xlrd` (the tables are legacy BIFF8 .xls, so the stdlib zip/XML route
used for .xlsx elsewhere in this tree does not apply). The installed version is
recorded in the audit; the project venv needs `pip install xlrd` to rerun this.
"""
from __future__ import annotations

import collections
import gzip
import json
import platform
import re
import sys
import time
from pathlib import Path

import xlrd

from common import ROOT, RESULTS, sha256, write_csv, write_json

SCRIPT = ROOT / "scripts" / "ingest_qpers_sid_supplementary.py"
INTAKE = ROOT / "external" / "intake" / "scirep2016_qpers_sid"
PROTEOME = ROOT / "external" / "proteomes" / "hsa.fasta.gz"
STEM = "41598_2016_BFsrep29808_MOESM{}_ESM.xls"

ADDUCTS = {"+57.02": "SH_iodoacetamide", "+414.19": "SH_probe", "+446.17": "SSH_probe"}
TIER_A = {"2": ("main_experiment", None), "7": ("silac_quadruplicate", None)}
TIER_B_DONORS = {"3": "GYY4137", "4": "Na2S4", "5": "Na2S", "6": "NaSH"}
MOESM7_DONOR_SHEETS = ("GYY4137", "Na2S4", "Na2S", "NaSH")
WINDOW_CENTRE = 15  # 0-based centre of the 31-residue Sequence window

MEANING = {
    ("Beads", "SH_probe"): "free thiol seen by the probe in the same run",
    ("Beads", "SSH_probe"): "direct adduct evidence of a persulfide",
    ("Beads", "SH_iodoacetamide"): "iodoacetamide-blocked cysteine on the beads",
    ("Total", "SH_probe"): "free thiol seen by the probe, whole lysate",
    ("Total", "SSH_probe"): "direct adduct evidence of a persulfide, whole lysate",
    ("Total", "SH_iodoacetamide"): "iodoacetamide-blocked cysteine, pre-enrichment",
    ("Elution", "SH_iodoacetamide"): "freed by TCEP then blocked: enrichment-based persulfide evidence",
    ("Elution", "SH_probe"): "probe-labelled thiol in the eluate, carryover",
    ("Elution", "SSH_probe"): "direct adduct evidence of a persulfide in the eluate",
    ("Flowthrough", "SH_probe"): "not retained by streptavidin",
    ("Flowthrough", "SSH_probe"): "not retained by streptavidin",
    ("Flowthrough", "SH_iodoacetamide"): "not retained by streptavidin",
}


def load_proteome():
    sequences = {}
    with gzip.open(PROTEOME, "rt") as handle:
        accession, chunks = None, []
        for line in handle:
            if line.startswith(">"):
                if accession:
                    sequences[accession] = "".join(chunks)
                match = re.match(r">(?:sp|tr)\|([^|]+)\|", line) or re.match(r">(\S+)", line)
                accession, chunks = match.group(1), []
            else:
                chunks.append(line.strip())
        if accession:
            sequences[accession] = "".join(chunks)
    return sequences


class Locator:
    """Exact lookup of a subsequence across the canonical proteome."""

    def __init__(self, sequences):
        self.order = sorted(sequences)
        self.blob = "|".join(sequences[a] for a in self.order)
        offsets, running = [], 0
        for accession in self.order:
            offsets.append(running)
            running += len(sequences[accession]) + 1
        self.offsets = offsets
        self.sequences = sequences
        self.cache = {}

    def _protein_at(self, index):
        low, high = 0, len(self.offsets) - 1
        while low < high:
            mid = (low + high + 1) // 2
            if self.offsets[mid] <= index:
                low = mid
            else:
                high = mid - 1
        return self.order[low], index - self.offsets[low] + 1

    def unique(self, subsequence):
        """Return (accession, 1-based start) when the subsequence occurs exactly once."""
        if subsequence in self.cache:
            return self.cache[subsequence]
        result, index, seen = None, self.blob.find(subsequence), 0
        while index != -1:
            seen += 1
            if seen > 1:
                result = "multi"
                break
            result = self._protein_at(index)
            index = self.blob.find(subsequence, index + 1)
        if result is None:
            result = "absent"
        self.cache[subsequence] = result
        return result


def strip_peptide(raw):
    """'R.TFC(+57.02)Q(+.98)LILDPIFK.V' -> ('TFCQLILDPIFK', {3: '+57.02', 4: '+0.98'})."""
    core = str(raw).strip()
    if len(core) > 2 and core[1] == ".":
        core = core[2:]
    if len(core) > 2 and core[-2] == ".":
        core = core[:-2]
    residues, mods, position = [], {}, 0
    for token in re.finditer(r"([A-Z])(\(([^)]*)\))?", core):
        position += 1
        residues.append(token.group(1))
        if token.group(3):
            mods[position] = "+" + token.group(3).lstrip("+")
    return "".join(residues), mods


def header_map(sheet, row):
    out = {}
    for column in range(sheet.ncols):
        key = re.sub(r"\s+", " ", str(sheet.cell_value(row, column))).strip()
        if key and key not in out:
            out[key] = column
    return out


def find_header_row(sheet, marker):
    for row in range(min(8, sheet.nrows)):
        if re.sub(r"\s+", " ", str(sheet.cell_value(row, 0))).strip() == marker:
            return row
    raise RuntimeError(f"no header row starting with {marker!r} in {sheet.name}")


def parse_accession(raw):
    raw = str(raw).strip()
    match (re.match(r"^(?:sp|tr)\|([^|]+)\|", raw) or re.match(r"^([A-Za-z0-9\-_]+)\|", raw)
           or re.match(r"^([A-Za-z0-9\-_]+)$", raw)):
        case None:
            return ""
        case found:
            return found.group(1)


def number(sheet, row, columns, *names):
    for name in names:
        if name in columns:
            value = sheet.cell_value(row, columns[name])
            if isinstance(value, float):
                return round(value, 6)
            text = str(value).strip()
            if text:
                try:
                    return round(float(text), 6)
                except ValueError:
                    return text
    return ""


def text(sheet, row, columns, *names):
    for name in names:
        if name in columns:
            return str(sheet.cell_value(row, columns[name])).strip()
    return ""


def main():
    started = time.time()
    sequences = load_proteome()
    locator = Locator(sequences)
    rows, stats = [], collections.Counter()

    # ---- tier A: peptide-level PEAKS searches that include +446.17 -------------
    for moesm, donor_sheets in (("2", (None,)), ("7", MOESM7_DONOR_SHEETS)):
        path = INTAKE / STEM.format(moesm)
        book = xlrd.open_workbook(path, on_demand=True)
        sheets = [s for s in book.sheet_names() if s != "Legend"]
        for sheet_name in sheets:
            sheet = book.sheet_by_name(sheet_name)
            head = find_header_row(sheet, "Experiment")
            columns = header_map(sheet, head)
            fraction = sheet_name if moesm == "2" else "Elution"
            donor = None if moesm == "2" else sheet_name
            for index in range(head + 1, sheet.nrows):
                peptide = text(sheet, index, columns, "Peptide")
                if not peptide:
                    continue
                stats[f"tierA_rows_{moesm}"] += 1
                stripped, mods = strip_peptide(peptide)
                located = locator.unique(stripped)
                if located == "absent":
                    stats["tierA_peptide_absent"] += 1
                    continue
                if located == "multi":
                    stats["tierA_peptide_multimaps"] += 1
                    continue
                accession, start = located
                protein = sequences[accession]
                for offset, residue in enumerate(stripped, start=1):
                    if residue != "C":
                        continue
                    site = start + offset - 1
                    if site > len(protein) or protein[site - 1] != "C":
                        stats["tierA_site_not_cysteine"] += 1
                        continue
                    adduct = mods.get(offset, "")
                    label = ADDUCTS.get(adduct, "C_unlabelled" if not adduct else "C_other_modification")
                    stats[f"tierA_class_{label}"] += 1
                    rows.append({
                        "evidence_tier": "A_adduct_based",
                        "source_table": f"MOESM{moesm}",
                        "sheet": sheet_name,
                        "fraction": fraction,
                        "donor": donor or "main_experiment_no_donor_split",
                        "experiment": text(sheet, index, columns, "Experiment"),
                        "species": "human",
                        "material": "HEK293 cells",
                        "accession": accession,
                        "site": site,
                        "accession_reported": parse_accession(text(sheet, index, columns, "Protein Accession")),
                        "gene": "",
                        "class_from_adduct": label,
                        "adduct": adduct,
                        "evidence_meaning": MEANING.get((fraction, label), "see the script docstring"),
                        "modified_sequence": peptide,
                        "mapping_route": "unique_peptide_in_canonical_proteome",
                        "ratio_average": "",
                        "ratio_sd": "",
                        "ratio_meaning": "",
                        "localization_score": "",
                        "significant": "",
                    })
        book.release_resources()

    # ---- tier B: MaxQuant site tables, +446.17 NOT searched -------------------
    for moesm, donor in TIER_B_DONORS.items():
        path = INTAKE / STEM.format(moesm)
        book = xlrd.open_workbook(path, on_demand=True)
        for sheet_name in (s for s in book.sheet_names() if s != "Legend"):
            sheet = book.sheet_by_name(sheet_name)
            head = find_header_row(sheet, "Protein")
            columns = header_map(sheet, head)
            if "Sequence window" not in columns:
                raise RuntimeError(f"MOESM{moesm}/{sheet_name} has no Sequence window column")
            site_column = ("Carbamidomethyl C-Sites" if "Carbamidomethyl C-Sites" in columns
                           else "Iodoacetyl-PEG2-Biotin C-Sites")
            label = "SSH_enriched" if sheet_name == "Elution" else "SH_probe_on_beads"
            meaning = MEANING.get(("Elution", "SH_iodoacetamide")) if sheet_name == "Elution" \
                else MEANING.get(("Beads", "SH_probe"))
            for index in range(head + 1, sheet.nrows):
                window = text(sheet, index, columns, "Sequence window")
                if not window:
                    continue
                stats[f"tierB_rows_{donor}_{sheet_name}"] += 1
                if len(window) != 31:
                    stats["tierB_window_wrong_length"] += 1
                    continue
                located = locator.unique(window)
                if located == "absent":
                    stats["tierB_window_absent"] += 1
                    continue
                if located == "multi":
                    stats["tierB_window_multimaps"] += 1
                    continue
                accession, start = located
                site = start + WINDOW_CENTRE
                if sequences[accession][site - 1] != "C":
                    stats["tierB_site_not_cysteine"] += 1
                    continue
                stats[f"tierB_class_{label}"] += 1
                rows.append({
                    "evidence_tier": "B_enrichment_based",
                    "source_table": f"MOESM{moesm}",
                    "sheet": sheet_name,
                    "fraction": sheet_name,
                    "donor": donor,
                    "experiment": "quadruplicate_SILAC_pooled",
                    "species": "human",
                    "material": "HEK293 cells",
                    "accession": accession,
                    "site": site,
                    "accession_reported": parse_accession(text(sheet, index, columns, "Protein")),
                    "gene": text(sheet, index, columns, "Gene names"),
                    "class_from_adduct": label,
                    "adduct": f"{site_column} (letter only, no position in the table)",
                    "evidence_meaning": meaning,
                    "modified_sequence": window,
                    "mapping_route": "unique_31mer_sequence_window_centre_is_the_site",
                    "ratio_average": number(sheet, index, columns, "Mean norm Ratio H/L",
                                            "Mean Norm Ratio H/L", "Mean  Norm Ratio H/L"),
                    "ratio_sd": number(sheet, index, columns, "STDEV"),
                    "ratio_meaning": ("MaxQuant normalised heavy-to-light ratio; a tier B positive means "
                                      "mean >= 1.30 and significant in at least 3 of 4 experiments"),
                    "localization_score": number(sheet, index, columns, "Score for localization"),
                    "significant": text(sheet, index, columns, "Significant"),
                })
        book.release_resources()

    write_csv(RESULTS / "qpers_sid_sites_normalised.csv", rows)

    # ---- pairing: what a negative set could be built from ---------------------
    direct_positive = {(r["accession"], r["site"]) for r in rows
                       if r["class_from_adduct"] == "SSH_probe"}
    probe_thiol = {(r["accession"], r["site"]) for r in rows
                   if r["class_from_adduct"] == "SH_probe" and r["fraction"] in ("Beads", "Total")}
    elution_enriched = {(r["accession"], r["site"]) for r in rows
                        if r["fraction"] == "Elution" and r["class_from_adduct"] == "SH_iodoacetamide"}
    pairing = [{
        "tier": "A_adduct_based", "scope": "all tier A tables pooled",
        "persulfide_evidence_sites": len(direct_positive),
        "free_thiol_probe_sites": len(probe_thiol),
        "both": len(direct_positive & probe_thiol),
        "free_thiol_only": len(probe_thiol - direct_positive),
        "persulfide_only": len(direct_positive - probe_thiol),
        "note": ("direct adduct evidence (+446.17) against probe-labelled free thiol (+414.19) in the "
                 "Beads and Total fractions; the free-thiol-only set is the candidate negative set"),
    }, {
        "tier": "A_adduct_based", "scope": "Elution fraction, enrichment route inside tier A tables",
        "persulfide_evidence_sites": len(elution_enriched),
        "free_thiol_probe_sites": len(probe_thiol),
        "both": len(elution_enriched & probe_thiol),
        "free_thiol_only": len(probe_thiol - elution_enriched),
        "persulfide_only": len(elution_enriched - probe_thiol),
        "note": ("iodoacetamide-blocked cysteines in the eluate are persulfides freed by TCEP, so this is "
                 "enrichment evidence appearing inside the tier A tables, not adduct evidence"),
    }]
    for donor in TIER_B_DONORS.values():
        positive = {(r["accession"], r["site"]) for r in rows
                    if r["donor"] == donor and r["class_from_adduct"] == "SSH_enriched"}
        negative = {(r["accession"], r["site"]) for r in rows
                    if r["donor"] == donor and r["class_from_adduct"] == "SH_probe_on_beads"}
        pairing.append({
            "tier": "B_enrichment_based", "scope": f"donor {donor}",
            "persulfide_evidence_sites": len(positive),
            "free_thiol_probe_sites": len(negative),
            "both": len(positive & negative),
            "free_thiol_only": len(negative - positive),
            "persulfide_only": len(positive - negative),
            "note": "Elution enriched sites against Beads probe-labelled thiols for the same donor",
        })
    union_b = {(r["accession"], r["site"]) for r in rows if r["class_from_adduct"] == "SSH_enriched"}
    arms = collections.Counter()
    for site in union_b:
        arms[sum(1 for d in TIER_B_DONORS.values()
                 if any(r["donor"] == d and r["class_from_adduct"] == "SSH_enriched"
                        and (r["accession"], r["site"]) == site for r in rows))] += 1
    pairing.append({
        "tier": "B_enrichment_based", "scope": "four donors pooled",
        "persulfide_evidence_sites": len(union_b),
        "free_thiol_probe_sites": len({(r["accession"], r["site"]) for r in rows
                                       if r["class_from_adduct"] == "SH_probe_on_beads"}),
        "both": "", "free_thiol_only": "", "persulfide_only": "",
        "note": "donor-arm depth distribution after mapping: " + json.dumps(dict(sorted(arms.items()))),
    })
    write_csv(RESULTS / "qpers_sid_class_pairing.csv", pairing)

    print(json.dumps({"rows": len(rows), "stats": dict(sorted(stats.items()))}, indent=2), flush=True)
    for entry in pairing:
        print(f"  [{entry['tier'][:1]}] {entry['scope'][:44]:46s} SSH {entry['persulfide_evidence_sites']:5} "
              f"SH {entry['free_thiol_probe_sites']:5} both {entry['both']:>5} SH_only {entry['free_thiol_only']:>5}",
              flush=True)

    write_json(RESULTS / "qpers_sid_ingest_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "purpose": ("map the qPerS-SID supplementary tables onto canonical human sites and record what "
                    "evidence each site carries; the positive/negative decision is deliberately left to "
                    "the scoring round"),
        "label_not_assigned_here": True,
        "sequence_reference": {
            "file": "external/proteomes/hsa.fasta.gz",
            "entries": len(sequences),
            "paper_used": "full UniProt human set of 2015-04-08, 68,506 entries, isoforms and TrEMBL included",
            "reconciliation": ("mapping is by sequence, not accession: accession matching would lose 41.5% of "
                              "rows (1,086 of 1,855 accessions in MOESM2/Beads are in the canonical set). "
                              "Isoform and TrEMBL accessions therefore collapse onto their canonical parent, "
                              "and the reported accession is kept as accession_reported for audit."),
            "acceptance_rule": ("a row is accepted only when its peptide (tier A) or 31-mer window (tier B) "
                               "occurs in exactly one canonical protein AND the implied position is a cysteine; "
                               "absent, multi-mapping and non-cysteine rows are counted and dropped"),
        },
        "evidence_tiers": {
            "A_adduct_based": "MOESM2 and MOESM7, PEAKS 7 searches including +446.17",
            "B_enrichment_based": ("MOESM3-6, MaxQuant with the same parameters but EXCLUDING +446.17, so a "
                                  "positive means mean normalised H/L >= 1.30 significant in at least 3 of 4 "
                                  "experiments after TCEP elution"),
            "must_not_be_merged": True,
        },
        "adduct_meaning_depends_on_fraction": {f"{k[0]}|{k[1]}": v for k, v in MEANING.items()},
        "counts": dict(sorted(stats.items())),
        "rows_written": len(rows),
        "outputs": ["results/qpers_sid_sites_normalised.csv", "results/qpers_sid_class_pairing.csv"],
        "limits": [
            "C_unlabelled cysteines were seen in the run but not seen by the probe, so they are not parallel-labelling negatives",
            "collapsing isoform accessions onto canonical parents can move a site's numbering relative to the paper's own table",
            "a peptide shared by several canonical proteins is dropped, which biases against paralogous families",
            "tier B positives are enrichment-based, not adduct-based, and must be reported as a reductive-elution enrichment label",
            "tier B's NaSH and Na2S arms contribute the most sites although the paper reports those donors as having hardly any detectable effect; that tension is unresolved",
            "no ProteomeXchange accession exists for this study, so nothing here can be re-searched from raw spectra",
        ],
        "input_hashes": {
            "external/proteomes/hsa.fasta.gz": sha256(PROTEOME),
            "scripts/ingest_qpers_sid_supplementary.py": sha256(SCRIPT),
            **{f"external/intake/scirep2016_qpers_sid/{STEM.format(m)}": sha256(INTAKE / STEM.format(m))
               for m in ("2", "3", "4", "5", "6", "7")},
        },
        "versions": {"python": sys.version, "xlrd": xlrd.__version__,
                     "note": ("xlrd is required because these are legacy BIFF8 .xls files; the stdlib zip/XML "
                              "route used for .xlsx elsewhere in this tree does not apply. The project venv "
                              "needs `pip install xlrd` to rerun this script.")},
    })


if __name__ == "__main__":
    main()
