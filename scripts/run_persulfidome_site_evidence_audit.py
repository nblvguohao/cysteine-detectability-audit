"""Does a deposited persulfidome actually carry cysteine-localised evidence?

The project's label provenance work concluded that tag-switch and dimedone-switch
enrichment is protein-level: the affinity step selects proteins, the digest then
yields peptides that need not carry the modified cysteine, and the recovered
cysteine carries a routine alkylation mass.  PXD005168 supported this indirectly
(`reports/V2_DIRECT_PERSULFIDE_EVIDENCE_CENSUS.md`).

This program tests it directly on PXD076106, "Persulfidome of A549 alveolar
epithelial cells" (2026, dimedone-switch, HPE-IAM alkylation, trypsin), using the
depositors' own Proteome Discoverer peptide-group export.  Two quantities decide
the question:

* how many peptide groups contain a cysteine at all, and
* where the Sulfide (+31.97) assignments are localised.

Counting rule for "contains a cysteine": Proteome Discoverer writes the annotated
sequence with flanking residues, as in ``[R].SAMPLEPEPTIDEK.[Y]``.  The flanks
belong to the protein, not the peptide, so they are stripped before the search;
counting them inflates the cysteine count several-fold.
"""
from __future__ import annotations

import collections
import json
import platform
import re
import statistics
import sys
import time
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from common import RESULTS, ROOT, sha256, write_csv, write_json

DATA = ROOT / "external/persulfidome_a549"
TABLE = DATA / "20210720_BignellE_AmichJ.peptide_groups.xlsx"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
FLANK = re.compile(r"^\[[A-Z\-]\]\.|\.\[[A-Z\-]\]$")
MOD_TOKEN = re.compile(r"(\d+)x([A-Za-z0-9\-\+/]+)\s*\[([^\]]+)\]")
ACCESSION = "PXD076106"


def cell_value(cell, shared):
    kind = cell.get("t")
    if kind == "inlineStr":
        inline = cell.find(f"{NS}is")
        if inline is None:
            return ""
        return "".join(node.text or "" for node in inline.iter() if node.tag.endswith("}t"))
    value = cell.find(f"{NS}v")
    if value is None:
        return ""
    if kind == "s" and shared:
        return shared[int(value.text)]
    return value.text


def read_sheet(path):
    archive = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in archive.namelist():
        for item in ET.fromstring(archive.read("xl/sharedStrings.xml")):
            shared.append(
                "".join(node.text or "" for node in item.iter() if node.tag.endswith("}t"))
            )
    sheet = next(n for n in archive.namelist() if n.startswith("xl/worksheets/sheet"))
    root = ET.fromstring(archive.read(sheet))
    rows = root.findall(f".//{NS}sheetData/{NS}row")
    header = [cell_value(c, shared) for c in rows[0].findall(f"{NS}c")]
    records = []
    for row in rows[1:]:
        values = [cell_value(c, shared) for c in row.findall(f"{NS}c")]
        if len(values) < len(header):
            values += [""] * (len(header) - len(values))
        records.append(dict(zip(header, values)))
    return header, records


def peptide_body(annotated):
    """Strip Proteome Discoverer flanking-residue notation."""
    return FLANK.sub("", annotated.strip()).replace(".", "").upper()


def main():
    started = time.time()
    if not TABLE.exists():
        raise SystemExit(f"Missing {TABLE}")
    header, records = read_sheet(TABLE)
    ratio_columns = [h for h in header if h.startswith("Abundance Ratio (log2)")]

    modification_residues = collections.defaultdict(collections.Counter)
    cysteine_rows = 0
    sulfide_rows = 0
    sulfide_on_cysteine = 0
    ratio_by_class = {"cysteine_containing": [], "no_cysteine": []}
    cysteine_site_rows = []

    for record in records:
        body = peptide_body(record.get("Annotated Sequence", ""))
        has_cysteine = "C" in body
        modifications = record.get("Modifications", "")
        for match in MOD_TOKEN.finditer(modifications):
            name = match.group(2)
            for part in match.group(3).split(";"):
                residue = part.strip()[:1]
                modification_residues[name][residue] += 1
        if "Sulfide" in modifications:
            sulfide_rows += 1
            if re.search(r"Sulfide \[[^\]]*C\d", modifications):
                sulfide_on_cysteine += 1
        if has_cysteine:
            cysteine_rows += 1
            cysteine_site_rows.append({
                "annotated_sequence": record.get("Annotated Sequence", ""),
                "peptide_body": body,
                "modifications": modifications,
                "master_protein_accessions": record.get("Master Protein Accessions", ""),
                "positions_in_master_proteins": record.get("Positions in Master Proteins", ""),
                "modifications_in_master_proteins": record.get(
                    "Modifications in Master Proteins", ""
                ),
                "n_psms": record.get("# PSMs", ""),
                "carries_sulfide_on_cysteine": bool(
                    re.search(r"Sulfide \[[^\]]*C\d", modifications)
                ),
                "carries_carbamidomethyl": "Carbamidomethyl" in modifications,
                **{column: record.get(column, "") for column in ratio_columns},
            })
        key = "cysteine_containing" if has_cysteine else "no_cysteine"
        first = record.get(ratio_columns[0], "") if ratio_columns else ""
        try:
            ratio_by_class[key].append(float(first))
        except (TypeError, ValueError):
            pass

    def describe(values):
        if not values:
            return {}
        values = sorted(values)
        return {
            "n": len(values),
            "median": round(statistics.median(values), 4),
            "q1": round(values[len(values) // 4], 4),
            "q3": round(values[3 * len(values) // 4], 4),
        }

    write_csv(RESULTS / "persulfidome_a549_cysteine_peptides.csv", cysteine_site_rows)
    audit = {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "deposition": ACCESSION,
        "title": "Persulfidome of A549 alveolar epithelial cells",
        "workflow": "dimedone-switch enrichment, DTT reduction, HPE-IAM alkylation, trypsin, Sequest HT",
        "source_table": str(TABLE.relative_to(ROOT)),
        "source_sha256": sha256(TABLE),
        "n_peptide_groups": len(records),
        "n_peptide_groups_containing_cysteine": cysteine_rows,
        "cysteine_peptide_fraction": round(cysteine_rows / len(records), 5),
        "n_rows_with_sulfide_assignment": sulfide_rows,
        "n_sulfide_assignments_localised_to_cysteine": sulfide_on_cysteine,
        "sulfide_assignment_residues": dict(modification_residues.get("Sulfide", {})),
        "carbamidomethyl_assignment_residues": dict(
            modification_residues.get("Carbamidomethyl", {})
        ),
        "all_modification_residues": {
            name: dict(counter) for name, counter in sorted(modification_residues.items())
        },
        "first_abundance_ratio_column": ratio_columns[0] if ratio_columns else None,
        "abundance_ratio_by_class": {
            name: describe(values) for name, values in ratio_by_class.items()
        },
        "counting_rule": (
            "flanking-residue notation stripped before testing for cysteine; counting the "
            "flanks inflates the cysteine peptide count several-fold"
        ),
        "interpretation": (
            "a persulfidome whose Sulfide assignments never localise to a cysteine carries no "
            "site-level persulfide evidence; the enrichment is protein-level, so most peptides "
            "of an enriched protein contain no cysteine at all"
        ),
        "what_this_is_not": [
            "a re-search of the raw files",
            "a claim that the depositors misreported anything",
            "evidence about whether the enriched proteins are persulfidated",
        ],
        "hardware": {"platform": platform.platform()},
        "versions": {"python": sys.version},
    }
    write_json(RESULTS / "persulfidome_a549_audit.json", audit)
    print(json.dumps({
        k: audit[k] for k in (
            "n_peptide_groups", "n_peptide_groups_containing_cysteine",
            "cysteine_peptide_fraction", "n_rows_with_sulfide_assignment",
            "n_sulfide_assignments_localised_to_cysteine", "sulfide_assignment_residues",
            "carbamidomethyl_assignment_residues", "abundance_ratio_by_class",
        )
    }, indent=2))


if __name__ == "__main__":
    main()
