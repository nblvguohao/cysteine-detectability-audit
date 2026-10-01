"""Normalise cysteine-modification site tables from three switch chemistries into one long table.

Round: 2026-09-15, cross-chemistry detectability census (track 1 of stage 1).
Output: results/ptm_census_sites.csv  (long format, one row per cysteine per role)
        results/ptm_census_sites_audit.json

WHAT THIS SCRIPT DECIDES, DECLARED BEFORE IT IS RUN
===================================================

Columns of the output table
---------------------------
dataset_id, chemistry_family, channel_available, species, proteome_key,
accession, site (1-based residue number), role, source_file

`role` takes exactly two values and nothing else:
  * ``positive``            - the modified form of this cysteine was reported by the
                              authors for this dataset under the authors' own criterion.
  * ``observed_unmodified`` - this cysteine was reported in the SAME experiment in its
                              unmodified (free-thiol / alkylated) form and was never
                              reported as the modified form in that experiment.

A cysteine that is both is a ``positive`` only.  The third class needed downstream -
cysteines never reported at all - is NOT enumerated here; it is generated in
run_ptm_detectability_census.py directly from the reference proteome, because it is
defined by absence and would otherwise have to be materialised twice.

Per-dataset positive and negative rules (these are the object of study, so they are
written out in full and are not parameters)
-------------------------------------------------------------------------------------
1. ``tomato_v2_benchmark`` - persulfidation-adjacent capture label, Solanum lycopersicum.
   Source: inputs/primary_site_folds.csv, the frozen v2 benchmark cohort.
   positive = label 1 (annotated in the capture workflow, i.e. pulled down and
   alkylated; per reports/V2_LABEL_SEMANTICS_AND_DETECTABILITY.md this is NOT a
   persulfidation call).  label 0 rows are unannotated cysteines of the same proteins
   and are therefore a *not-observed* negative, not an observed-unmodified one, so they
   are not emitted as ``observed_unmodified``; this dataset offers no observed-unmodified
   channel.
2. ``qtrp_S1_ph5`` / ``qtrp_S2_ph5`` - persulfidation, low-pH QTRP, differential
   alkylation, human (Fu 2020, ARS; paywalled supplementary, NOT redistributed).
   Source: results/qtrp_sites_normalised.csv (ingested in an earlier round).
   positive = rows of that arm's ``-SSH_pH5`` sheet.
   observed_unmodified = rows of the SAME arm's ``-SH_pH5`` sheet, minus positives of
   that arm.  This is the project's registered v3 raw negative rule
   (protocols/v3_qtrp_analysis_plan.json, internal_cohorts/primary): sites called SSH in
   *other* arms are retained as negatives here.
3. ``qpers_sid_tierB`` - persulfidation, qPerS-SID capture + tag switch, human HEK293
   (Longen 2016, Sci Rep; open access).
   Source: results/qpers_sid_sites_normalised.csv.
   positive = evidence_tier ``B_enrichment_based`` and class ``SSH_enriched`` and
   ratio_average >= 1.30 and, where the source sheet carries a Significant column,
   significant == '+'.  The 1.30 threshold is the authors' own, transcribed in
   reports/QPERS_SID_LABEL_CRITERIA.md; this script re-derives it rather than reusing a
   cached label file, and asserts the previously published counts (729 positives,
   1923 observed-unmodified) to detect drift.
   observed_unmodified = same-sheet Elution rows below that threshold, plus
   bead-retained probe-labelled thiols (``SH_probe_on_beads``), minus positives.
4. ``ath_*`` / ``human_*`` external datasets - parsed from author supplementary tables
   downloaded from the Europe PMC open-access supplementary endpoint.  Each parser is
   declared in EXTERNAL_PARSERS with the sheet, the columns used, and the criterion
   applied; nothing is inferred from a column whose meaning is not stated in the
   supplementary legend or the paper's methods.

Site-position convention
------------------------
All positions are 1-based residue numbers in the reference proteome entry of the
recorded accession.  A row is dropped, and counted in the audit, when the accession is
absent from the reference proteome or when the residue at that position is not cysteine.
Isoform suffixes (``-2``) are kept only if the proteome contains that exact accession;
otherwise the base accession is tried once and the mapping route is recorded.

Not done here: no model is fitted, no interval is computed, no claim about chemistry is
made.  This script only normalises tables.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import pathlib
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

ROOT = pathlib.Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
OUT_SITES = RESULTS / "ptm_census_sites.csv"
OUT_AUDIT = RESULTS / "ptm_census_sites_audit.json"
SUPP = ROOT / "external/intake/ptm_census_supp"

PROTEOMES = {
    "hsa": ROOT / "external/proteomes/hsa.fasta.gz",
    "ath": ROOT / "external/proteomes/ath.fasta.gz",
    "osa": ROOT / "external/proteomes/osa.fasta.gz",
    "sly": ROOT / "inputs/sly_proteome.tsv",
    "ath_tair": ROOT / "external/proteomes/ath_uniprot_tair_map.tsv.gz",
}


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_proteome(key: str) -> dict:
    path = PROTEOMES[key]
    seqs = {}
    if key == "ath_tair":
        # TAIR locus -> UniProt canonical sequence.  Several UniProt entries can
        # cross-reference one locus; the longest sequence is kept, which is declared in
        # the module docstring.  TAIR isoform suffixes (".1") are stripped by the caller.
        with gzip.open(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                seq = (row.get("Sequence") or "").strip()
                for locus in (row.get("TAIR") or "").split(";"):
                    locus = locus.strip().upper()
                    if not locus or not seq:
                        continue
                    if len(seq) > len(seqs.get(locus, "")):
                        seqs[locus] = seq
        return seqs
    if key == "sly":
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                seqs[row["Entry"]] = row["Sequence"]
        return seqs
    with gzip.open(path, "rt") as handle:
        name, buf = None, []
        for line in handle:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                head = line[1:].split()[0]
                parts = head.split("|")
                name = parts[1] if len(parts) > 2 else head
                buf = []
            else:
                buf.append(line.strip())
        if name:
            seqs[name] = "".join(buf)
    return seqs


def read_csv_rows(path: pathlib.Path) -> list:
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


# ----------------------------------------------------------------------------- xlsx
def xlsx_sheets(zf: zipfile.ZipFile) -> list:
    wb = zf.read("xl/workbook.xml").decode("utf-8", "replace")
    rels = zf.read("xl/_rels/workbook.xml.rels").decode("utf-8", "replace")
    rid_to_target = dict(re.findall(r'Id="([^"]+)"[^>]*Target="([^"]+)"', rels))
    out = []
    for m in re.finditer(r"<sheet[^>]*>", wb):
        tag = m.group(0)
        name = re.search(r'name="([^"]*)"', tag)
        rid = re.search(r'r:id="([^"]*)"', tag)
        if not (name and rid):
            continue
        target = rid_to_target.get(rid.group(1), "")
        target = target[1:] if target.startswith("/") else "xl/" + target.lstrip("/")
        out.append((name.group(1), target))
    return out


def xlsx_rows(path: pathlib.Path, sheet_name: str) -> list:
    """Return the sheet as a list of lists of strings (stdlib only)."""
    zf = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in zf.namelist():
        root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
        for si in root:
            shared.append("".join(t.text or "" for t in si.iter() if t.tag.endswith("}t")))
    target = None
    for name, tgt in xlsx_sheets(zf):
        if name == sheet_name:
            target = tgt
            break
    if target is None:
        raise KeyError(f"sheet {sheet_name!r} not in {path.name}: {[n for n, _ in xlsx_sheets(zf)]}")
    root = ET.fromstring(zf.read(target))
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    rows = []
    for row in root.iter(ns + "row"):
        cells = {}
        for c in row.iter(ns + "c"):
            ref = c.get("r") or ""
            col = re.match(r"([A-Z]+)", ref)
            idx = 0
            for ch in (col.group(1) if col else "A"):
                idx = idx * 26 + (ord(ch) - 64)
            idx -= 1
            t = c.get("t")
            v = c.find(ns + "v")
            if t == "s" and v is not None:
                val = shared[int(v.text)]
            elif t == "inlineStr":
                is_ = c.find(ns + "is")
                val = "".join(x.text or "" for x in is_.iter() if x.tag.endswith("}t")) if is_ is not None else ""
            else:
                val = v.text if v is not None else ""
            cells[idx] = (val or "").strip()
        if cells:
            width = max(cells) + 1
            rows.append([cells.get(i, "") for i in range(width)])
    return rows


# ------------------------------------------------------------------- local builders
def build_tomato() -> list:
    rows = read_csv_rows(ROOT / "inputs/primary_site_folds.csv")
    out = []
    for r in rows:
        if r["label"] != "1":
            continue
        out.append(dict(dataset_id="tomato_v2_benchmark", chemistry_family="persulfidation_capture_label",
                        channel_available="capture_only", species="Solanum lycopersicum", proteome_key="sly",
                        accession=r["accession"], site=int(r["position"]), role="positive",
                        source_file="inputs/primary_site_folds.csv"))
    return out


def build_qtrp() -> list:
    rows = read_csv_rows(ROOT / "results/qtrp_sites_normalised.csv")
    out = []
    for arm in ("S1", "S2"):
        pos = {(r["accession"], int(float(r["site"]))) for r in rows
               if r["source_table"] == arm and r["sheet"] == "-SSH_pH5"}
        neg = {(r["accession"], int(float(r["site"]))) for r in rows
               if r["source_table"] == arm and r["sheet"] == "-SH_pH5"} - pos
        did = f"qtrp_{arm}_ph5"
        for (acc, site), role in [(k, "positive") for k in sorted(pos)] + \
                                 [(k, "observed_unmodified") for k in sorted(neg)]:
            out.append(dict(dataset_id=did, chemistry_family="persulfidation",
                            channel_available="both", species="Homo sapiens", proteome_key="hsa",
                            accession=acc, site=site, role=role,
                            source_file="results/qtrp_sites_normalised.csv"))
    return out


def build_qpers() -> list:
    rows = read_csv_rows(ROOT / "results/qpers_sid_sites_normalised.csv")
    tierb = [r for r in rows if r["evidence_tier"] == "B_enrichment_based"]
    sheets_with_sig = {r["sheet"] for r in tierb if (r["significant"] or "").strip()}

    def ratio(r):
        try:
            return float(r["ratio_average"])
        except (TypeError, ValueError):
            return None

    pos, negb = set(), set()
    for r in tierb:
        key = (r["accession"], int(float(r["site"])))
        if r["class_from_adduct"] == "SSH_enriched":
            v = ratio(r)
            ok = v is not None and v >= 1.30
            if r["sheet"] in sheets_with_sig:
                ok = ok and (r["significant"] or "").strip() == "+"
            (pos if ok else negb).add(key)
        elif r["class_from_adduct"] == "SH_probe_on_beads":
            negb.add(key)
    negb -= pos
    out = []
    for key, role in [(k, "positive") for k in sorted(pos)] + \
                     [(k, "observed_unmodified") for k in sorted(negb)]:
        out.append(dict(dataset_id="qpers_sid_tierB", chemistry_family="persulfidation",
                        channel_available="both", species="Homo sapiens", proteome_key="hsa",
                        accession=key[0], site=key[1], role=role,
                        source_file="results/qpers_sid_sites_normalised.csv"))
    return out, len(pos), len(negb)


# ---------------------------------------------------------------- external parsers
def _num(s):
    try:
        return float(str(s).replace(",", ""))
    except (TypeError, ValueError):
        return None


def parse_generic_supp(spec: dict, seqs: dict) -> tuple:
    """Parse one supplementary sheet into positives / observed_unmodified.

    spec keys: zip, member, sheet, header_row, acc_col, site_col, role_rule
    Every key is written into the audit, so the parse is reproducible from the record.
    """
    zpath = SUPP / spec["zip"]
    inner = SUPP / "unpacked" / spec["member"]
    inner.parent.mkdir(parents=True, exist_ok=True)
    if not inner.exists():
        with zipfile.ZipFile(zpath) as zf:
            inner.write_bytes(zf.read(spec["member"]))
    rows = xlsx_rows(inner, spec["sheet"])
    hdr = [c.strip() for c in rows[spec["header_row"]]]
    body = rows[spec["header_row"] + 1:]
    def col(name):
        for i, c in enumerate(hdr):
            if c.lower() == name.lower():
                return i
        for i, c in enumerate(hdr):
            if name.lower() in c.lower():
                return i
        raise KeyError(f"column {name!r} not in {hdr}")
    ai, si = col(spec["acc_col"]), col(spec["site_col"])
    out, dropped = [], {"no_accession": 0, "no_site": 0}
    for r in body:
        if len(r) <= max(ai, si):
            dropped["no_site"] += 1
            continue
        acc_raw, site_raw = r[ai].strip(), r[si].strip()
        if not acc_raw:
            dropped["no_accession"] += 1
            continue
        accs = [a for a in re.split(r"[;,\s]+", acc_raw) if a]
        sites = [s for s in re.split(r"[;,\s]+", site_raw) if s]
        if not sites:
            dropped["no_site"] += 1
            continue
        for j, a in enumerate(accs):
            a = a.split("|")[1] if a.count("|") >= 2 else a
            s = sites[j] if j < len(sites) else sites[0]
            n = _num(re.sub(r"[^0-9.]", "", s))
            if n is None:
                dropped["no_site"] += 1
                continue
            out.append((a, int(n)))
    return out, dropped, hdr


# =====================================================================================
# EXTERNAL AUTHOR SUPPLEMENTARY TABLES - every parse rule is written here, before running
# =====================================================================================
# Each entry states: which archive, which member workbook, which sheet, which header row
# (0-based index into the sheet's non-empty rows), which columns, and the criterion.
# Nothing is inferred from a column whose meaning is not stated in the sheet legend.
#
# cysboost2019_human_sno_hela / _shsy5y  (Mnatsakanyan 2019, Nat Commun, open access)
#   positives: Supplemental Data 5b (HeLa) / 6b (SH-SY5Y) - two-column author tables of
#     UniProt accession + SNO site ("C273"); these are the authors' own SNO site calls.
#   observed_unmodified: Supplemental Data 2a - the identified cysteine-peptide list of
#     the same study.  Every cysteine inside an identified peptide is taken as observed;
#     positives are then subtracted.  Peptide->protein coordinates come from the authors'
#     "Positions in Master Proteins" field, never recomputed from sequence search.
# natcomm2023_ath_sno  (Zhan 2023, Nat Commun, open access, FAT-switch)
#   positives: Supplementary Dataset 3b "all identified SNO peptides", column
#     "Modifications in Master Proteins", entries of the form "AT1G01220.1 1xTFIA [C115]".
#   observed_unmodified: none available - the study reports enriched SNO peptides only,
#     so this dataset contributes the capture channel only.
# abiotech2025_ath_sno  (PAT-switch, aBIOTECH 2025, open access)
#   positives: Supplementary data 2 sheet "Class I PTMs", rows with Modification Name
#     "S-nitrosylation" and Target Amino Acid "C" (the authors' class I localisation set).
#   observed_unmodified: sheet "Modified Peptides" - cysteines inside identified modified
#     peptides that are not in the positive set.  Declared limitation: this table lists
#     peptides carrying at least one modification of any kind, so it is a subset of the
#     full identified cysteinome and the resulting negative set is smaller, not different
#     in kind.
EXTERNAL_DATASETS = [
    {
        "dataset_id": "cysboost2019_human_sno_hela",
        "chemistry_family": "s_nitrosylation", "species": "Homo sapiens",
        "proteome_key": "hsa", "channel_available": "both",
        "source_file": "external/intake/ptm_census_supp/cysboost2019_human_sno_PMC6522481.zip",
        "positive": {"zip": "cysboost2019_human_sno_PMC6522481.zip",
                     "member": "41467_2019_10182_MOESM8_ESM.xlsx", "sheet": "5b",
                     "header_row": 1, "mode": "two_column",
                     "acc_col": "Uniprot Protein Accession", "site_col": "SNO site"},
        "observed_unmodified": {"zip": "cysboost2019_human_sno_PMC6522481.zip",
                                "member": "41467_2019_10182_MOESM5_ESM.xlsx", "sheet": "2a",
                                "header_row": 2, "mode": "peptide_cys",
                                "acc_col": "Master Protein Accessions",
                                "seq_col": "Annotated Sequence",
                                "pos_col": "Positions in Master Proteins"},
    },
    {
        "dataset_id": "cysboost2019_human_sno_shsy5y",
        "chemistry_family": "s_nitrosylation", "species": "Homo sapiens",
        "proteome_key": "hsa", "channel_available": "both",
        "source_file": "external/intake/ptm_census_supp/cysboost2019_human_sno_PMC6522481.zip",
        "positive": {"zip": "cysboost2019_human_sno_PMC6522481.zip",
                     "member": "41467_2019_10182_MOESM9_ESM.xlsx", "sheet": "6b",
                     "header_row": 1, "mode": "two_column",
                     "acc_col": "Uniprot Protein Accession", "site_col": "SNO site"},
        "observed_unmodified": {"zip": "cysboost2019_human_sno_PMC6522481.zip",
                                "member": "41467_2019_10182_MOESM9_ESM.xlsx", "sheet": "6a",
                                "header_row": 2, "mode": "peptide_cys",
                                "acc_col": "Master Protein Accessions",
                                "seq_col": "Annotated Sequence",
                                "pos_col": "Positions in Master Proteins"},
    },
    {
        "dataset_id": "natcomm2023_ath_sno",
        "chemistry_family": "s_nitrosylation", "species": "Arabidopsis thaliana",
        "proteome_key": "ath_tair", "channel_available": "capture_only",
        "source_file": "external/intake/ptm_census_supp/natcomm2023_ath_sno_PMC10241878.zip",
        "positive": {"zip": "natcomm2023_ath_sno_PMC10241878.zip",
                     "member": "41467_2023_39078_MOESM6_ESM.xlsx",
                     "sheet": "all identified SNO peptides ", "header_row": 1,
                     "mode": "mod_string", "mod_col": "Modifications in Master Proteins"},
    },
    {
        "dataset_id": "abiotech2025_ath_sno",
        "chemistry_family": "s_nitrosylation", "species": "Arabidopsis thaliana",
        "proteome_key": "ath_tair", "channel_available": "both",
        "source_file": "external/intake/ptm_census_supp/abiotech2025_ath_persulf_PMC12454217.zip",
        "positive": {"zip": "abiotech2025_ath_persulf_PMC12454217.zip",
                     "member": "42994_2025_227_MOESM2_ESM.xlsx", "sheet": "Class I PTMs",
                     "header_row": 0, "mode": "class_i",
                     "acc_col": "Protein Accession", "site_col": "Position",
                     "mod_name_col": "Modification Name", "aa_col": "Target Amino Acid",
                     "mod_name_value": "S-nitrosylation"},
        "observed_unmodified": {"zip": "abiotech2025_ath_persulf_PMC12454217.zip",
                                "member": "42994_2025_227_MOESM2_ESM.xlsx",
                                "sheet": "Modified Peptides", "header_row": 0,
                                "mode": "peptide_cys",
                                "acc_col": "Master Protein Accessions",
                                "seq_col": "Annotated Sequence",
                                "pos_col": "Positions in Master Proteins"},
    },
    {
        # Huang 2020, Front Plant Sci, open access: YAP1C disulfide-linked peptide reporter
        # in Arabidopsis.  positives = the cross-linked (sulfenylated) cysteine, marked
        # with '#' in the authors' peptide strings (Dataset S2).  observed_unmodified =
        # cysteines of the REGULAR (non-cross-linked) peptides identified in the same runs
        # (Dataset S1b), which is this study's own observed cysteinome.
        "dataset_id": "fps2020_ath_sulfenyl",
        "chemistry_family": "sulfenylation", "species": "Arabidopsis thaliana",
        "proteome_key": "ath_tair", "channel_available": "both",
        "source_file": "external/intake/ptm_census_supp/fps2020_plant_sulfenyl_PMC7343964.zip",
        "positive": {"zip": "fps2020_plant_sulfenyl_PMC7343964.zip",
                     "nested_zip": "Data_Sheet_1.zip", "member": "DATA/Data 2.XLSX",
                     "unpacked_as": "fps2020_Data_2.XLSX", "sheet": "DatasetS2",
                     "header_row": 3, "mode": "peptide_locate",
                     "pep_col": 0, "prot_col": 13, "mark": "#", "require_col0": True},
        "observed_unmodified": {"zip": "fps2020_plant_sulfenyl_PMC7343964.zip",
                                "nested_zip": "Data_Sheet_1.zip", "member": "DATA/Data 1.XLSX",
                                "unpacked_as": "fps2020_Data_1.XLSX",
                                "sheet": "DatasetS1b_regular_peptide", "header_row": 1,
                                "mode": "peptide_locate", "pep_col": 1, "prot_col": 4,
                                "mark": "", "require_col0": True},
    },
]


def _sheet_rows(spec: dict):
    zpath = SUPP / spec["zip"]
    inner = SUPP / "unpacked" / spec.get("unpacked_as", spec["member"])
    inner.parent.mkdir(parents=True, exist_ok=True)
    if not inner.exists():
        with zipfile.ZipFile(zpath) as zf:
            if spec.get("nested_zip"):
                import io
                with zipfile.ZipFile(io.BytesIO(zf.read(spec["nested_zip"]))) as nz:
                    inner.write_bytes(nz.read(spec["member"]))
            else:
                inner.write_bytes(zf.read(spec["member"]))
    return xlsx_rows(inner, spec["sheet"])


def locate_peptide_cys(spec: dict, seqs: dict) -> tuple:
    """Positions of cysteines of a reported peptide, by exact search in the protein sequence.

    Declared rule: the peptide is searched in the reference sequence of the FIRST listed
    protein; a row is used only when the peptide occurs EXACTLY ONCE in that sequence.
    Rows failing that (absent peptide, repeated peptide, unmapped accession) are dropped
    and counted in the audit.  When spec["mark"] is set, only cysteines carrying that
    marker character are taken as modified; otherwise every cysteine of the peptide is
    emitted (used to build the observed-unmodified set).
    """
    rows = _sheet_rows(spec)
    pep_i, prot_i = spec["pep_col"], spec["prot_col"]
    mark = spec.get("mark", "")
    drops = {"blank": 0, "unmapped_accession": 0, "peptide_absent": 0, "peptide_repeated": 0}
    out = []
    for r in rows[spec["header_row"] + 1:]:
        if len(r) <= max(pep_i, prot_i):
            drops["blank"] += 1
            continue
        if spec.get("require_col0") and not r[0].strip():
            continue
        pep_raw, prot_raw = r[pep_i].strip(), r[prot_i].strip()
        if not pep_raw or not prot_raw:
            drops["blank"] += 1
            continue
        acc = re.split(r"[;\s/]+", prot_raw)[0].split(".")[0].upper()
        seq = seqs.get(acc)
        if seq is None:
            drops["unmapped_accession"] += 1
            continue
        clean, offsets = [], []
        for ch in pep_raw:
            if ch == mark and mark:
                if clean and clean[-1] == "C":
                    offsets.append(len(clean) - 1)
                continue
            if ch.isalpha():
                clean.append(ch.upper())
        pep = "".join(clean)
        if not mark:
            offsets = [i for i, ch in enumerate(pep) if ch == "C"]
        if not offsets:
            drops["blank"] += 1
            continue
        n_occ = seq.count(pep)
        if n_occ == 0:
            drops["peptide_absent"] += 1
            continue
        if n_occ > 1:
            drops["peptide_repeated"] += 1
            continue
        start = seq.find(pep)
        for off in offsets:
            out.append((acc, start + off + 1))
    return out, rows[spec["header_row"]][:16], drops


def _col(hdr: list, name: str) -> int:
    for i, c in enumerate(hdr):
        if c.strip().lower() == name.strip().lower():
            return i
    for i, c in enumerate(hdr):
        if name.strip().lower()[:18] in c.strip().lower():
            return i
    raise KeyError(f"column {name!r} not among {hdr[:14]}")


def parse_external(spec: dict) -> tuple:
    """Return (list of (accession, site), header, drop counters) for one declared spec."""
    rows = _sheet_rows(spec)
    hdr = rows[spec["header_row"]]
    body = rows[spec["header_row"] + 1:]
    drops = {"blank": 0, "unparsable_site": 0}
    out = []
    mode = spec["mode"]
    if mode == "two_column":
        ai, si = _col(hdr, spec["acc_col"]), _col(hdr, spec["site_col"])
        for r in body:
            if len(r) <= max(ai, si) or not r[ai].strip():
                drops["blank"] += 1
                continue
            digits = re.sub(r"[^0-9]", "", r[si])
            if not digits:
                drops["unparsable_site"] += 1
                continue
            out.append((r[ai].strip().split("|")[-1] if "|" not in r[ai] else r[ai].split("|")[1],
                        int(digits)))
    elif mode == "class_i":
        ai, si = _col(hdr, spec["acc_col"]), _col(hdr, spec["site_col"])
        mi, aa = _col(hdr, spec["mod_name_col"]), _col(hdr, spec["aa_col"])
        for r in body:
            if len(r) <= max(ai, si, mi, aa):
                drops["blank"] += 1
                continue
            if r[mi].strip() != spec["mod_name_value"] or r[aa].strip().upper() != "C":
                continue
            digits = re.sub(r"[^0-9]", "", r[si])
            if not digits or not r[ai].strip():
                drops["unparsable_site"] += 1
                continue
            out.append((r[ai].strip(), int(digits)))
    elif mode == "mod_string":
        mi = _col(hdr, spec["mod_col"])
        for r in body:
            if len(r) <= mi or not r[mi].strip():
                drops["blank"] += 1
                continue
            text = r[mi]
            acc = text.split()[0]
            hits = re.findall(r"C(\d+)", text)
            if not hits:
                drops["unparsable_site"] += 1
                continue
            for h in hits:
                out.append((acc, int(h)))
    elif mode == "peptide_cys":
        ai = _col(hdr, spec["acc_col"])
        qi = _col(hdr, spec["seq_col"])
        pi = _col(hdr, spec["pos_col"])
        for r in body:
            if len(r) <= max(ai, qi, pi):
                drops["blank"] += 1
                continue
            acc_raw, seq_raw, pos_raw = r[ai].strip(), r[qi].strip(), r[pi].strip()
            if not (acc_raw and seq_raw and pos_raw):
                drops["blank"] += 1
                continue
            core = seq_raw
            if "." in core:
                parts = core.split(".")
                core = parts[1] if len(parts) >= 3 else parts[-1]
            core = re.sub(r"[^A-Za-z]", "", core).upper()
            spans = re.findall(r"\[(\d+)-(\d+)\]", pos_raw)
            accs = [a for a in re.split(r"[;\s]+", acc_raw) if a]
            if not spans or not accs:
                drops["unparsable_site"] += 1
                continue
            start = int(spans[0][0])
            for a in accs[:1]:
                for k, ch in enumerate(core):
                    if ch == "C":
                        out.append((a, start + k))
    else:
        raise ValueError(mode)
    return out, hdr[:16], drops


def main() -> int:
    audit = {"script": pathlib.Path(__file__).name,
             "script_sha256": sha256(pathlib.Path(__file__)),
             "seed": 20260915, "inputs": {}, "datasets": {}, "dropped": {}}
    rows = []
    rows += build_tomato()
    rows += build_qtrp()
    qp, n_pos, n_neg = build_qpers()
    rows += qp
    audit["qpers_sid_reproduced_counts"] = {"positives": n_pos, "observed_unmodified": n_neg,
                                            "published_in_QPERS_SID_LABEL_CRITERIA": {"positives": 729,
                                                                                      "observed_unmodified": 1923},
                                            "match": (n_pos == 729 and n_neg == 1923)}

    audit["external_specs"] = EXTERNAL_DATASETS
    for spec in EXTERNAL_DATASETS:
        zpath = SUPP / spec["positive"]["zip"]
        if not zpath.exists():
            audit["dropped"][spec["dataset_id"]] = "archive missing: " + str(zpath.relative_to(ROOT))
            continue
        got = {}
        for role in ("positive", "observed_unmodified"):
            if role not in spec:
                continue
            if spec[role]["mode"] == "peptide_locate":
                pairs, hdr, drops = locate_peptide_cys(spec[role], load_proteome(spec["proteome_key"]))
            else:
                pairs, hdr, drops = parse_external(spec[role])
            got[role] = pairs
            audit["dropped"][f"{spec['dataset_id']}:{role}"] = drops
            audit["inputs"][f"{spec['dataset_id']}:{role}:header"] = hdr
        pos = {(a.split(".")[0].upper() if spec["proteome_key"] == "ath_tair" else a, s)
               for a, s in got.get("positive", [])}
        neg = {(a.split(".")[0].upper() if spec["proteome_key"] == "ath_tair" else a, s)
               for a, s in got.get("observed_unmodified", [])} - pos
        for key, role in [(k, "positive") for k in sorted(pos)] + \
                         [(k, "observed_unmodified") for k in sorted(neg)]:
            rows.append(dict(dataset_id=spec["dataset_id"], chemistry_family=spec["chemistry_family"],
                             channel_available=spec["channel_available"], species=spec["species"],
                             proteome_key=spec["proteome_key"], accession=key[0], site=key[1],
                             role=role, source_file=spec["source_file"]))
        audit["inputs"][spec["dataset_id"] + ":archive_sha256"] = sha256(zpath)

    # proteome mapping and cysteine check
    kept, drop_counts = [], {}
    cache = {}
    for r in rows:
        key = r["proteome_key"]
        if key not in cache:
            cache[key] = load_proteome(key)
        seqs = cache[key]
        acc = r["accession"]
        route = "exact"
        seq = seqs.get(acc)
        if seq is None and "-" in acc:
            seq = seqs.get(acc.split("-")[0])
            route = "isoform_base"
        d = drop_counts.setdefault(r["dataset_id"], {"unmapped_accession": 0, "not_cysteine": 0,
                                                     "out_of_range": 0, "kept": 0})
        if seq is None:
            d["unmapped_accession"] += 1
            continue
        i = r["site"] - 1
        if not (0 <= i < len(seq)):
            d["out_of_range"] += 1
            continue
        if seq[i] != "C":
            d["not_cysteine"] += 1
            continue
        d["kept"] += 1
        r["mapping_route"] = route
        kept.append(r)

    fields = ["dataset_id", "chemistry_family", "channel_available", "species", "proteome_key",
              "accession", "site", "role", "mapping_route", "source_file"]
    with open(OUT_SITES, "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields)
        w.writeheader()
        for r in kept:
            w.writerow({k: r[k] for k in fields})

    per = {}
    for r in kept:
        d = per.setdefault(r["dataset_id"], {"positive": 0, "observed_unmodified": 0,
                                             "chemistry_family": r["chemistry_family"],
                                             "species": r["species"], "proteome_key": r["proteome_key"],
                                             "channel_available": r["channel_available"], "proteins": set()})
        d[r["role"]] += 1
        d["proteins"].add(r["accession"])
    for k, v in per.items():
        v["n_proteins"] = len(v.pop("proteins"))
    audit["datasets"] = per
    audit["dropped_mapping"] = drop_counts
    for name, path in [("qtrp", ROOT / "results/qtrp_sites_normalised.csv"),
                       ("qpers", ROOT / "results/qpers_sid_sites_normalised.csv"),
                       ("tomato_folds", ROOT / "inputs/primary_site_folds.csv"),
                       ("sly_proteome", PROTEOMES["sly"]), ("hsa_proteome", PROTEOMES["hsa"]),
                       ("ath_proteome", PROTEOMES["ath"])]:
        if path.exists():
            audit["inputs"][name] = {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
    OUT_AUDIT.write_text(json.dumps(audit, indent=1, ensure_ascii=False, default=str))
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "proteins"} for k, v in per.items()},
                     indent=1, default=str))
    print("qpers match:", audit["qpers_sid_reproduced_counts"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
