"""Classify the 42 Stage C candidates by DESIGN, from metadata only, before any coincidence value is read.

THE DEFECT THIS CLOSES (amendment A3 to the Stage A/B/C pre-registration)
  Stage C measures the share of observed cysteines that carry a site assignment. That share is trivially
  close to 100% in TWO-STATE designs (OxICAT, differential alkylation, thiol-switch with a total-thiol
  channel): there almost every captured cysteine carries some label, and the reduced state is measured
  right beside the oxidised one. Those deposits contain exactly the 'detected but unmodified' control
  whose absence Artefact 4 describes, so a high share there would REFUTE the framing while satisfying
  the K1 arithmetic. K1 must therefore not count them.

RULE (declared now, applied to metadata, independent of any number Stage C produces)
  two_state_likely     title + description + sample protocol match TWO_STATE
  capture_only_likely  not two_state_likely, and matches CAPTURE
  unclear              neither
  Only capture_only_likely rows may count toward K1. 'unclear' rows are still reported and can support
  K2, but a K1 that rests on an 'unclear' row is reported as 'K1_unconfirmed_design' and not as K1.

WHY METADATA ONLY, AND NOT READING EACH PAPER
  Reading each deposit's design after seeing which ones scored high would let the answer choose the
  classification. This runs before those values exist to the reader: the Stage C process had produced no
  results file when this was written (results are written at the end of that process).

This script reads no coincidence value.
"""
import csv, hashlib, json, pathlib, re, sys, time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "results/pride_coincidence_design_classification_2026-09-19.csv"
AUD = ROOT / "results/pride_coincidence_design_classification_2026-09-19_audit.json"
API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

TWO_STATE = re.compile(
    r"oxicat|oxidi[sz]ed and reduced|reduced and oxidi[sz]ed|oxidi[sz]ed/reduced|reduced/oxidi[sz]ed|"
    r"percent(age)? (of )?(cysteine )?oxidation|degree of (cysteine )?oxidation|oxidation (state|level|stoichiometr|occupancy|ratio)|"
    r"fraction (of )?oxidi[sz]ed|ratio of oxidi[sz]ed|total (thiol|cysteine)|redox state|sicylia|cysquant|"
    r"reduced (cysteine|thiol)s?\b|both (the )?(reduced|oxidi)|light and heavy|differential alkylation|"
    r"simultaneous quantification|thiol.?switch|redox switch|reversibly oxidi[sz]ed and total", re.I)
CAPTURE = re.compile(
    r"biotin.?switch|pull-?down|enrich|streptavidin|\bABE\b|acyl-?biotin|resin-?assisted|\bRAC\b|dimedone|"
    r"click|\bBTA\b|\bQTRP\b|affinity|capture|iodotmt|maleimide", re.I)

def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception:
            if i == tries - 1: raise
            time.sleep(2 * (i + 1))

cand = list(csv.DictReader(open(ROOT / "results/pride_coincidence_candidates_2026-09-19.csv")))
MQ = re.compile(r"(?<![a-z])txt(?![a-z])|maxquant|(?<![a-z])combined(?![a-z])", re.I)
CYS = re.compile(r"persulf|sulfhydr|nitros|sulfen|sulfin|glutathion|palmitoyl|thiol|redox|cystein|iodotmt|oxicat|dimedone|qtrp|biotin.?switch|\bABE\b", re.I)
chosen = [r for r in cand if r["name_class"] in ("archive_needs_listing", "both_by_name")
          and any(MQ.search(a) for a in r["archives"].split("; ")) and CYS.search(r["title"])]

rows = []
for r in chosen:
    p = get(f"{API}/projects/{r['accession']}")
    text = " ".join(str(p.get(k, "")) for k in ("title", "projectDescription", "sampleProcessingProtocol"))
    ts, cp = TWO_STATE.search(text), CAPTURE.search(text)
    cls = "two_state_likely" if ts else ("capture_only_likely" if cp else "unclear")
    rows.append({"accession": r["accession"], "organism": r["organisms"][:22], "title": r["title"][:78],
                 "design_class": cls, "two_state_marker": ts.group(0) if ts else "",
                 "capture_marker": cp.group(0) if cp else ""})
    time.sleep(0.2)
with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

from collections import Counter
c = Counter(r["design_class"] for r in rows)
AUD.write_text(json.dumps({
 "script": "scripts/pride_coincidence_design_classification_2026-09-19.py", "script_sha256": sha(__file__),
 "interpreter": sys.version.split()[0], "n": len(rows), "counts": dict(c),
 "amendment": "A3: only capture_only_likely may count toward K1; unclear rows can support K2 or K1_unconfirmed_design",
 "stageC_results_file_existed_when_run": (ROOT / "results/pride_coincidence_assessed_2026-09-19.csv").exists(),
 "what_this_does_not_establish": [
   "A keyword screen over three text fields is a screen, not a reading of the paper. A design that "
   "reports both states without using any of the marker words is classed 'unclear' or 'capture_only_likely' "
   "in error; that is why 'unclear' cannot support K1.",
 ],
}, ensure_ascii=False, indent=1, sort_keys=True))
print(dict(c))
for r in rows: print(f"  {r['design_class']:<20} {r['accession']} {r['title'][:52]:<54} [{r['two_state_marker'] or r['capture_marker']}]")
