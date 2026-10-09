"""Step 1c (POST HOC): UniSave last release of 20 randomly chosen (seed 20260930) of the 400 dropped
rice accessions (rest.uniprot.org only). Output: fetched/unisave_last_release_sample.csv"""
import random, time, urllib.request
import pandas as pd
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
R = pd.read_csv(FETCHED / "inactive_reasons_400.csv")
acc = sorted(R.accession)
sample = sorted(random.Random(SEED).sample(acc, 20))
rows = []
for a in sample:
    req = urllib.request.Request(f"https://rest.uniprot.org/unisave/{a}?format=json",
                                 headers={"User-Agent": "cys-revision-posthoc/2026-09-30"})
    b = urllib.request.urlopen(req, timeout=60).read()
    (FETCHED / "inactive_reasons" / f"{a}_unisave.json").write_bytes(b)
    v = json.loads(b)["results"]
    top = max(v, key=lambda x: int(x["entryVersion"]))
    rows.append({"accession": a, "last_entry_version": int(top["entryVersion"]), "lastRelease": top["lastRelease"],
                 "lastReleaseDate": top["lastReleaseDate"], "sha256": hashlib.sha256(b).hexdigest()})
    time.sleep(0.5)
D = pd.DataFrame(rows); D.to_csv(FETCHED / "unisave_last_release_sample.csv", index=False)
print(D.to_string()); print(D.lastRelease.value_counts())
