"""Step 1b (POST HOC): why are the 400 candidate rice proteins missing from the stored cohort?

For each of the 400 accessions (s05_kept_vs_deleted_proteins.csv, kept_in_stored_cohort == False)
fetch the UniProtKB entry JSON (rest.uniprot.org only) and record entryType / inactiveReason, and
the UniSave last release in which the accession had an entry. Raw responses are saved.
Output: fetched/inactive_reasons_400.csv, fetched/inactive_reasons_400.json
"""
import time
import urllib.error
import urllib.request

import pandas as pd

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa

KV = pd.read_csv(RES / "s05_kept_vs_deleted_proteins.csv")
lost = sorted(KV.loc[~KV.kept_in_stored_cohort, "accession"])
out_dir = FETCHED / "inactive_reasons"
out_dir.mkdir(parents=True, exist_ok=True)
nbytes, rows = 0, []


def get(url):
    global nbytes
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "cys-revision-posthoc/2026-09-30"})
            with urllib.request.urlopen(req, timeout=60) as r:
                b = r.read(); nbytes += len(b); time.sleep(0.25)
                return b
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(3 * (attempt + 1)); continue
            return None
        except urllib.error.URLError:
            time.sleep(3 * (attempt + 1))
    return None


for a in lost:
    b = get(f"https://rest.uniprot.org/uniprotkb/{a}?format=json&fields=accession")
    rec = {"accession": a}
    if b is None:
        rec["status"] = "no_response"
    else:
        (out_dir / f"{a}.json").write_bytes(b)
        j = json.loads(b)
        rec["entryType"] = j.get("entryType")
        ir = j.get("inactiveReason") or {}
        rec["inactiveReasonType"] = ir.get("inactiveReasonType")
        rec["deletedReason"] = ir.get("deletedReason")
        rec["mergeDemergeTo"] = ";".join(ir.get("mergeDemergeTo", []) or [])
    rows.append(rec)
R = pd.DataFrame(rows)
R.to_csv(FETCHED / "inactive_reasons_400.csv", index=False)
summ = {"analysis_label": "POST HOC revision analysis 2026-09-30 (not registered)",
        "n": len(R), "bytes_downloaded": nbytes,
        "entryType": R["entryType"].value_counts(dropna=False).to_dict(),
        "inactiveReasonType": R["inactiveReasonType"].value_counts(dropna=False).to_dict(),
        "deletedReason": R["deletedReason"].value_counts(dropna=False).to_dict()}
write_json(FETCHED / "inactive_reasons_400.json", summ)
print(json.dumps(summ, indent=1))
