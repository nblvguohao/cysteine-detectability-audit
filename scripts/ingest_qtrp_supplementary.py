"""Normalise the low-pH QTRP supplementary tables into one site-level table.

Source: Fu L, Liu K, He J, Tian C, Yu X, Yang J. *Direct Proteomic Mapping of
Cysteine Persulfidation.* Antioxid Redox Signal 2020;33(15):1061-1076,
doi 10.1089/ars.2019.7777, PMID 31411056. Supplementary Tables S1-S6, obtained
2026-09-15 through the reader's institutional subscription. The files stay under
`external/intake/` and are not redistributed in any deliverable.

Why this matters: the method labels free thiols and persulfides **in parallel in
the same experiment**, and the two classes are separated by the adduct on the
cysteine, `C(IPM)` for a free thiol against `C(IPMS)` for a persulfide. That
gives the detectability-matched negative set the project has lacked, where the
negative is a cysteine the same run detected rather than a cysteine nobody has
ever seen. PXD005168 failed this test because its channels overlapped at site
level and its intact-tag channel was majority decoy.

Two distinct labels come out:

* class, `SSH` against `SH`, a mutually exclusive identification per peptide;
* reactivity, only in Table S2, where R10:1 is the heavy(100 uM) to light(10 uM)
  probe ratio and the authors call below 3.0 hyperreactive, 3.0 to 5.0
  moderately reactive, and the rest weakly reactive.

In Tables S1 and S3 the light-to-heavy ratio is an identification quality filter
(kept between 0.67 and 1.5), not a biological quantity, so it is recorded but
must not be read as occupancy.
"""
from __future__ import annotations

import json
import platform
import re
import sys
import time
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

from common import ROOT, RESULTS, sha256, write_csv, write_json

INTAKE = ROOT / "external/intake/ars2020_qtrp"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

TABLE_CONTEXT = {
    "sj-xlsx-1": {"table": "S1", "material": "A549 lysate, NaHS-treated", "species": "human"},
    "sj-xlsx-2": {"table": "S2", "material": "A549 lysate, NaHS-treated, reactivity profiling", "species": "human"},
    "sj-xlsx-3": {"table": "S3", "material": "five cell lines, NaHS-treated", "species": "human"},
    "sj-xlsx-4": {"table": "S4", "material": "HEK293 expressing Flag-CSE", "species": "human"},
    "sj-xlsx-5": {"table": "S5", "material": "mouse tissues", "species": "mouse"},
    "sj-xlsx-6": {"table": "S6", "material": "pooled human serum", "species": "human"},
}


def sheet_values(archive, part, shared):
    def value(cell):
        kind = cell.get("t")
        if kind == "inlineStr":
            inline = cell.find(f"{NS}is")
            return "".join(n.text or "" for n in inline.iter() if n.tag.endswith("}t")) if inline is not None else ""
        node = cell.find(f"{NS}v")
        if node is None:
            return ""
        return shared[int(node.text)] if kind == "s" and shared else node.text

    root = ET.fromstring(archive.read(part))
    rows = []
    for row in root.findall(f".//{NS}sheetData/{NS}row"):
        cells = {}
        for cell in row.findall(f"{NS}c"):
            column = "".join(ch for ch in (cell.get("r") or "") if ch.isalpha())
            cells[column] = value(cell)
        rows.append(cells)
    return rows


def column_order(cells):
    def key(letters):
        total = 0
        for ch in letters:
            total = total * 26 + (ord(ch) - 64)
        return total
    return sorted(cells, key=key)


def find(header, *patterns):
    for column, name in header.items():
        lowered = (name or "").strip().lower()
        for pattern in patterns:
            if re.search(pattern, lowered):
                return column
    return None


def classify(sheet_name, modified_sequence):
    """Class from the adduct on the cysteine, cross-checked against the sheet name."""
    sequence = modified_sequence or ""
    has_persulfide = "(IPMS)" in sequence
    has_thiol = re.search(r"\(IPM\)", sequence) is not None
    from_sequence = None
    if has_persulfide and not has_thiol:
        from_sequence = "SSH"
    elif has_thiol and not has_persulfide:
        from_sequence = "SH"
    elif has_persulfide and has_thiol:
        from_sequence = "both_adducts_in_row"
    lowered = sheet_name.lower()
    from_sheet = None
    if "ssh" in lowered:
        from_sheet = "SSH"
    elif "-sh" in lowered or lowered.endswith("sh"):
        from_sheet = "SH"
    return from_sequence, from_sheet


def main():
    started = time.time()
    files = sorted(INTAKE.glob("sj-xlsx-*.xlsx"))
    if not files:
        raise SystemExit(f"No workbooks under {INTAKE}")

    rows = []
    disagreements = 0
    unparsed = 0
    per_sheet = []
    for path in files:
        stem = re.match(r"(sj-xlsx-\d+)", path.name).group(1)
        context = TABLE_CONTEXT.get(stem, {"table": stem, "material": "", "species": ""})
        archive = zipfile.ZipFile(path)
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            for item in ET.fromstring(archive.read("xl/sharedStrings.xml")):
                shared.append("".join(n.text or "" for n in item.iter() if n.tag.endswith("}t")))
        workbook = archive.read("xl/workbook.xml").decode("utf-8", errors="replace")
        names = re.findall(r'<sheet name="([^"]+)"', workbook)
        parts = sorted(
            [n for n in archive.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", n)],
            key=lambda n: int(re.search(r"sheet(\d+)", n).group(1)),
        )
        for sheet_name, part in zip(names, parts):
            if sheet_name.strip().lower() == "legend":
                continue
            values = sheet_values(archive, part, shared)
            if len(values) < 2:
                continue
            header = {c: values[0].get(c, "") for c in column_order(values[0])}
            accession_col = find(header, r"accession", r"^leading ac$")
            site_col = find(header, r"^site$", r"leading site")
            sequence_col = find(header, r"modified sequence")
            if accession_col is None or site_col is None:
                per_sheet.append({
                    "table": context["table"], "sheet": sheet_name,
                    "status": "skipped_no_site_columns", "rows": len(values) - 1,
                })
                continue
            gene_col = find(header, r"gene ?name", r"^genename$")
            spectra_col = find(header, r"spectral counts", r"spectrum count")
            average_col = find(header, r"^average", r"ratio ave")
            sd_col = find(header, r"^sd ")
            kept = 0
            for record in values[1:]:
                accession = (record.get(accession_col) or "").strip()
                site_raw = (record.get(site_col) or "").strip()
                if not accession or not site_raw:
                    unparsed += 1
                    continue
                try:
                    site = int(float(site_raw))
                except ValueError:
                    unparsed += 1
                    continue
                sequence = (record.get(sequence_col) or "") if sequence_col else ""
                from_sequence, from_sheet = classify(sheet_name, sequence)
                if from_sequence in ("SH", "SSH") and from_sheet in ("SH", "SSH") \
                        and from_sequence != from_sheet:
                    disagreements += 1
                rows.append({
                    "source_table": context["table"],
                    "sheet": sheet_name,
                    "material": context["material"],
                    "species": context["species"],
                    "ph": "5.0" if "5" in sheet_name and "7.6" not in sheet_name else
                          ("7.6" if "7.6" in sheet_name else ""),
                    "class_from_sheet": from_sheet or "",
                    "class_from_adduct": from_sequence or "",
                    "accession": accession,
                    "site": site,
                    "gene": (record.get(gene_col) or "").strip() if gene_col else "",
                    "modified_sequence": sequence,
                    "spectral_counts": (record.get(spectra_col) or "") if spectra_col else "",
                    "ratio_average": (record.get(average_col) or "") if average_col else "",
                    "ratio_sd": (record.get(sd_col) or "") if sd_col else "",
                    "ratio_meaning": (
                        "R10:1 heavy 100 uM over light 10 uM probe; below 3.0 hyperreactive, "
                        "3.0 to 5.0 moderate, above 5.0 weak"
                        if context["table"] == "S2" else
                        "light to heavy tag ratio used as an identification filter between 0.67 "
                        "and 1.5; not an occupancy measure"
                    ),
                })
                kept += 1
            per_sheet.append({
                "table": context["table"], "sheet": sheet_name,
                "status": "parsed", "rows": kept,
            })
            print(json.dumps({"table": context["table"], "sheet": sheet_name, "rows": kept}), flush=True)

    write_csv(RESULTS / "qtrp_sites_normalised.csv", rows)

    # overlap between the two classes inside each experiment
    grouped = defaultdict(lambda: {"SH": set(), "SSH": set()})
    for row in rows:
        label = row["class_from_adduct"] if row["class_from_adduct"] in ("SH", "SSH") else row["class_from_sheet"]
        if label not in ("SH", "SSH"):
            continue
        key = (row["source_table"], row["sheet"])
        grouped[key][label].add((row["accession"], row["site"]))
    overlap = []
    for (table, sheet), classes in sorted(grouped.items()):
        overlap.append({
            "table": table, "sheet": sheet,
            "n_sh_sites": len(classes["SH"]),
            "n_ssh_sites": len(classes["SSH"]),
            "n_shared_sites": len(classes["SH"] & classes["SSH"]),
        })
    # experiment-level pairing: SH and SSH sheets of the same table and pH
    pairs = defaultdict(lambda: {"SH": set(), "SSH": set()})
    for row in rows:
        label = row["class_from_adduct"] if row["class_from_adduct"] in ("SH", "SSH") else row["class_from_sheet"]
        if label not in ("SH", "SSH"):
            continue
        pairs[(row["source_table"], row["ph"])][label].add((row["accession"], row["site"]))
    pairing = []
    for (table, ph), classes in sorted(pairs.items()):
        pairing.append({
            "table": table, "ph": ph or "unspecified",
            "n_sh_sites": len(classes["SH"]),
            "n_ssh_sites": len(classes["SSH"]),
            "n_shared_sites": len(classes["SH"] & classes["SSH"]),
            "n_sh_only": len(classes["SH"] - classes["SSH"]),
            "n_ssh_only": len(classes["SSH"] - classes["SH"]),
        })
    write_csv(RESULTS / "qtrp_class_pairing.csv", pairing)

    write_json(RESULTS / "qtrp_ingest_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "host": platform.node(),
        "source": {
            "paper": "Fu et al., Antioxid Redox Signal 2020;33(15):1061-1076",
            "doi": "10.1089/ars.2019.7777",
            "pmid": "31411056",
            "obtained": "institutional subscription, 2026-09-15",
            "redistribution": "files stay under external/intake and are not included in any deliverable",
        },
        "n_rows": len(rows),
        "n_rows_unparsed": unparsed,
        "class_source_disagreements": disagreements,
        "class_rule": "adduct on the cysteine, C(IPM) means free thiol and C(IPMS) means persulfide; cross-checked against the sheet name",
        "per_sheet": per_sheet,
        "per_sheet_class_counts": overlap,
        "experiment_level_pairing": pairing,
        "label_semantics": {
            "class": "SSH is a persulfide identification, SH is a free-thiol identification in the same run",
            "S2_reactivity": "R10:1 heavy 100 uM over light 10 uM; authors' thresholds 3.0 and 5.0",
            "S1_S3_ratio": "light to heavy identification filter, not occupancy",
        },
        "file_hashes": {p.name: sha256(p) for p in files},
        "versions": {"python": sys.version},
    })
    print(json.dumps({"rows": len(rows), "unparsed": unparsed,
                      "class_disagreements": disagreements,
                      "pairing": pairing}, indent=2))


if __name__ == "__main__":
    main()
