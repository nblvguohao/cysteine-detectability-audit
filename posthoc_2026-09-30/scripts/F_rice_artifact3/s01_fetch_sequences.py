"""Step 1 (POST HOC): fetch sequences for rice accessions absent from UP000059680 (UniProt 2026_03).

Network: rest.uniprot.org only (allowed for this item). Total download capped at 19 MB.

  (a) active entries: GET /uniprotkb/accessions?accessions=...&format=tsv  (batches of 200)
  (b) everything not returned under its own accession by (a) -- entries DELETED in 2026_02
      ("Not part of a reference proteome") and merged/secondary accessions: UniParc, GET
      /uniparc/stream?query=dbid:A OR dbid:B ...&format=tsv&fields=upi,accession,length,sequence
      (batches of 100). UniParc lists every UniProtKB accession that ever pointed to a sequence as
      ACC.seqversion; for each requested accession the highest sequence version is kept.
      This deviates from the item brief's suggestion (UniSave per accession) because UniSave needs
      one version-listing call per accession (~10^4 calls); (c) checks UniParc against UniSave.
  (c) UniSave spot check: for 25 randomly chosen (seed 20260930) accessions resolved via UniParc,
      fetch the version listing and the last version's FASTA from /unisave and compare sequences.

Every raw response is saved under results/F_rice_artifact3/fetched/ with its sha256.
"""
import random
import time
import urllib.error
import urllib.parse
import urllib.request

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa

BASE = "https://rest.uniprot.org"
CAP = 19 * 1024 * 1024
state = {"bytes": 0, "requests": 0}
log = []


def get(url, tag):
    if state["bytes"] > CAP:
        raise SystemExit("download cap reached; stopping")
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "cys-revision-posthoc/2026-09-30"})
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
                rel = r.headers.get("X-UniProt-Release")
            state["bytes"] += len(body)
            state["requests"] += 1
            log.append({"tag": tag, "url": url[:300] + ("..." if len(url) > 300 else ""),
                        "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest(),
                        "uniprot_release": rel})
            time.sleep(0.6)
            return body.decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(5 * (attempt + 1)); continue
            if e.code == 404:
                log.append({"tag": tag, "url": url[:300], "bytes": 0, "http": 404})
                return ""
            raise
        except urllib.error.URLError:
            time.sleep(5 * (attempt + 1))
    raise SystemExit(f"failed after retries: {url[:200]}")


need = [a.strip() for a in (RES / "s00_accessions_to_fetch.txt").read_text(encoding="utf-8").split() if a.strip()]
FETCHED.mkdir(parents=True, exist_ok=True)

# ---------- (a) active entries ----------
active = {}
for i in range(0, len(need), 200):
    chunk = need[i:i + 200]
    url = (f"{BASE}/uniprotkb/accessions?accessions={','.join(chunk)}&format=tsv"
           f"&fields=accession,id,sequence_version,length,organism_id,sequence")
    txt = get(url, f"active_{i // 200:03d}")
    (FETCHED / f"active_batch_{i // 200:03d}.tsv").write_text(txt, encoding="utf-8")
    lines = txt.splitlines()
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) >= 6:
            active[f[0]] = {"entry_name": f[1], "seq_version": f[2], "length": int(f[3]),
                            "organism_id": f[4], "sequence": f[5]}
print("active returned:", len(active), "bytes so far", state["bytes"], flush=True)

rest = [a for a in need if a not in active]
print("to resolve via UniParc:", len(rest), flush=True)

# ---------- (b) UniParc for deleted / secondary accessions ----------
uniparc = {}   # acc -> (version, upi, sequence)
for i in range(0, len(rest), 100):
    chunk = rest[i:i + 100]
    q = " OR ".join(f"dbid:{a}" for a in chunk)
    url = (f"{BASE}/uniparc/stream?query={urllib.parse.quote(q)}&format=tsv"
           f"&fields=upi,accession,length,organism_id,sequence")
    txt = get(url, f"uniparc_{i // 100:03d}")
    (FETCHED / f"uniparc_batch_{i // 100:03d}.tsv").write_text(txt, encoding="utf-8")
    want = set(chunk)
    for ln in txt.splitlines()[1:]:
        f = ln.split("\t")
        if len(f) < 5:
            continue
        upi, accs, seq = f[0], f[1], f[4]
        for tok in accs.replace(",", ";").split(";"):
            tok = tok.strip()
            if not tok:
                continue
            base, _, ver = tok.partition(".")
            if base in want:
                # UniParc lists an ACTIVE cross-reference without a version suffix and an inactive
                # one as ACC.seqversion; an active xref (if any) wins, else the highest version
                v = int(ver) if ver.isdigit() else 10 ** 6
                if base not in uniparc or v > uniparc[base][0]:
                    uniparc[base] = (v, upi, seq, f[3])
print("resolved via UniParc:", len(uniparc), "bytes so far", state["bytes"], flush=True)

# ---------- (c) UniSave spot check ----------
rng = random.Random(SEED)
sample = sorted(rng.sample(sorted(uniparc), min(25, len(uniparc))))
spot = []
for a in sample:
    js = get(f"{BASE}/unisave/{a}?format=json", f"unisave_list_{a}")
    try:
        vers = json.loads(js)["results"]
        last = max(int(v["entryVersion"]) for v in vers)
    except Exception:
        spot.append({"accession": a, "status": "listing_failed"}); continue
    fa = get(f"{BASE}/unisave/{a}?format=fasta&versions={last}", f"unisave_fasta_{a}")
    seq = "".join(l.strip() for l in fa.splitlines() if l and not l.startswith(">"))
    spot.append({"accession": a, "unisave_last_entry_version": last,
                 "identical_to_uniparc_highest_seq_version": seq == uniparc[a][2],
                 "len_unisave": len(seq), "len_uniparc": len(uniparc[a][2])})
    (FETCHED / f"unisave_{a}_v{last}.fasta").write_text(fa, encoding="utf-8")

# ---------- write mapping ----------
rows = []
for a in need:
    if a in active:
        r = active[a]
        rows.append({"requested_accession": a, "source": "uniprot_rest_active", "upi": "",
                     "seq_version": r["seq_version"], "organism_id": r["organism_id"],
                     "length": len(r["sequence"]), "sequence": r["sequence"]})
    elif a in uniparc:
        v, upi, seq, org = uniparc[a]
        rows.append({"requested_accession": a, "source": "uniparc_inactive_or_secondary", "upi": upi,
                     "seq_version": v, "organism_id": org, "length": len(seq), "sequence": seq})
    else:
        rows.append({"requested_accession": a, "source": "not_found", "upi": "", "seq_version": "",
                     "organism_id": "", "length": 0, "sequence": ""})
with open(FETCHED / "fetched_sequences.tsv", "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)

summary = {
    "analysis_label": "POST HOC revision analysis 2026-09-30 (not registered)",
    "n_requested": len(need),
    "n_active_uniprotkb": sum(r["source"] == "uniprot_rest_active" for r in rows),
    "n_uniparc": sum(r["source"] == "uniparc_inactive_or_secondary" for r in rows),
    "n_not_found": sum(r["source"] == "not_found" for r in rows),
    "organism_ids": {k: sum(1 for r in rows if r["organism_id"] == k) for k in sorted({r["organism_id"] for r in rows})},
    "requests": state["requests"], "bytes_downloaded": state["bytes"],
    "unisave_spot_check": spot,
    "fetched_sequences_tsv_sha256": sha256_file(FETCHED / "fetched_sequences.tsv"),
    "request_log": log,
}
write_json(FETCHED / "fetch_summary.json", summary)
print(json.dumps({k: v for k, v in summary.items() if k not in ("request_log",)}, indent=1))
