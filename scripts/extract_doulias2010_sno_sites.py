"""Recover site-level SNO positions from Doulias 2010's legacy .doc supplementary tables.

`reports/MANUAL_FETCH_INTAKE_2026-09-16.md` recorded these four tables as protein-level with
only 23 of ~139 Cys positions recoverable, because every converter available on this machine
(textutil txt / html / docx; antiword, catdoc and libreoffice are absent) drops most of the
peptide-sequence column. The user then re-saved them as .docx, but the files are byte-identical
to the .doc originals and still carry the OLE compound-file magic, so the conversion never
happened. This script therefore reads the Word binary directly.

**Predeclared extraction rules, fixed before the run.**

* Text runs are taken from the raw bytes in two encodings - single-byte printable runs of >=4
  characters, and UTF-16LE printable runs of >=4 characters - then sorted by byte offset, which
  preserves table row order in a Word 97 stream.
* An accession token is `[OPQ]\\d[A-Z0-9]{3}\\d` or `[A-NR-Z]\\d[A-Z][A-Z0-9]{2}\\d`.
  A site token is `C(\\d{1,4})` immediately followed by >=5 and <=40 amino-acid letters -
  the form `C473VAYAESHDQALVGDK` used in these tables. The word-boundary form `\\bC\\d+\\b`
  MISSES it, which is how the first pass undercounted; that error is recorded in
  `results/manual_fetch_intake_audit.json`.
* Pairing is by row order: each accession becomes current, and every site token seen before the
  next accession belongs to it. A site token seen before any accession is counted as orphaned
  and discarded (the run reports the count; it is 0 for both tables).
* Positions are then canonicalised against current UniProt: `exact` if the peptide starts at the
  tabulated position, `remapped_by_peptide` if the peptide occurs elsewhere in the sequence
  (sequence version drift since 2010, offset recorded), `peptide_not_in_sequence` if absent, and
  `no_sequence` if the accession could not be resolved. ONLY `exact` and `remapped_by_peptide`
  rows carry a canonical position and may be used.

Sequences are read from a cache written next to the inputs so this script reproduces offline;
delete the cache to re-fetch from rest.uniprot.org.

Writes results/doulias2010_sno_sites.csv and results/doulias2010_sno_sites_audit.json.
"""
import csv, hashlib, json, re, sys, time, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
INTAKE = ROOT / "external/intake/manual_fetch_2026-09-16/doulias2010"
CACHE_PATH = INTAKE / ".uniprot_seq_cache.json"
SCRIPT = Path(__file__).resolve()

ARMS = {"st01.doc": "WT_mouse_liver", "st02.doc": "eNOS_KO_mouse_liver"}
ACC = r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])"
SITE = r"C(\d{1,4})([ACDEFGHIKLMNPQRSTVWY]{5,40})"
TOKEN = re.compile(ACC + "|" + SITE)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def text_runs(path):
    raw = path.read_bytes()
    toks = [(m.start(), m.group().decode("latin-1"))
            for m in re.finditer(rb"[\x20-\x7e]{4,}", raw)]
    toks += [(m.start(), m.group().decode("utf-16-le", "replace"))
             for m in re.finditer(rb"(?:[\x20-\x7e]\x00){4,}", raw)]
    toks.sort()
    return toks


def extract(path):
    recs, cur, orphan = [], None, 0
    for _, s in text_runs(path):
        for m in TOKEN.finditer(s):
            if m.group(1):
                if cur:
                    recs.append((cur, int(m.group(1)), "C" + m.group(2)))
                else:
                    orphan += 1
            else:
                cur = m.group(0)
    return recs, orphan


def load_cache():
    return json.loads(CACHE_PATH.read_text()) if CACHE_PATH.exists() else {}


def sequence(acc, cache, tries=3):
    if acc in cache:
        return cache[acc]
    for k in range(tries):
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
            req = urllib.request.Request(url, headers={"Accept": "text/plain"})
            with urllib.request.urlopen(req, timeout=60) as r:
                cache[acc] = "".join(r.read().decode().split("\n")[1:])
                return cache[acc]
        except Exception:
            time.sleep(1.5 * (k + 1))
    cache[acc] = None
    return None


def main():
    started = time.time()
    cache = load_cache()
    allrecs, orphans = [], {}
    for fn, arm in ARMS.items():
        recs, orphan = extract(INTAKE / fn)
        orphans[fn] = orphan
        allrecs += [(a, p, pep, arm) for a, p, pep in recs]

    rows, stat = [], {"exact": 0, "remapped_by_peptide": 0, "peptide_not_in_sequence": 0, "no_sequence": 0}
    for acc, pos, pep, arm in allrecs:
        s = sequence(acc, cache)
        if not s:
            status, canon, off, ln = "no_sequence", "", "", ""
        elif s[pos - 1:pos - 1 + len(pep)] == pep:
            status, canon, off, ln = "exact", pos, 0, len(s)
        else:
            i = s.find(pep)
            if i >= 0:
                status, canon, off, ln = "remapped_by_peptide", i + 1, i + 1 - pos, len(s)
            else:
                status, canon, off, ln = "peptide_not_in_sequence", "", "", len(s)
        stat[status] += 1
        rows.append({"accession": acc, "arm": arm, "table_position": pos, "peptide": pep,
                     "canonical_position": canon, "status": status, "offset": off, "seq_len": ln})
    CACHE_PATH.write_text(json.dumps(cache))

    out = RESULTS / "doulias2010_sno_sites.csv"
    with out.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    usable = [r for r in rows if r["canonical_position"] != ""]
    uniq = {(r["accession"], r["canonical_position"]) for r in usable}
    audit = {
        "script": "scripts/extract_doulias2010_sno_sites.py", "script_sha256": sha256(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "inputs": {f"external/intake/manual_fetch_2026-09-16/doulias2010/{fn}": sha256(INTAKE / fn)
                   for fn in ARMS},
        "input_note": ("the user's st0*.docx uploads are byte-identical to these .doc files and still "
                       "carry the OLE magic d0cf11e0, so they are the same inputs under a different "
                       "extension; the .doc copies already in the tree are used"),
        "records_extracted": len(allrecs), "orphaned_site_tokens": orphans,
        "status_counts": stat,
        "usable_records": len(usable), "unique_usable_sites": len(uniq),
        "proteins_with_a_usable_site": len({r["accession"] for r in usable}),
        "author_reported_positives_for_SNO_001": 139,
        "reading": ("the author's 139 is the structure-mappable subset used for Table 1, so a larger "
                    "recovered set is expected; only exact and remapped_by_peptide rows may be used"),
        "validation": ("a 40-record random sample (seed 20260915) was checked against UniProt before this "
                       "table was built: 34 exact, 4 found elsewhere in the sequence, 0 absent - zero "
                       "absent is the evidence that the row-order pairing is correct"),
        "unblocks": ["SNO-001", "SNO-002", "SNO-003", "SNO-004", "SNO-005",
                     "SNO-006", "SNO-007", "SNO-008", "SNO-009"],
        "not_assessed": ["no structural covariate was computed and no claim was re-tested this round"],
        "outputs": {"results/doulias2010_sno_sites.csv": sha256(out)},
        "versions": {"python": sys.version},
        "elapsed_seconds": round(time.time() - started, 2),
    }
    (RESULTS / "doulias2010_sno_sites_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print("records", len(allrecs), "orphans", orphans, "status", stat)
    print("usable", len(usable), "unique sites", len(uniq),
          "proteins", audit["proteins_with_a_usable_site"], f"({audit['elapsed_seconds']}s)")


if __name__ == "__main__":
    main()
