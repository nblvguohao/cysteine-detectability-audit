"""Phase 1: screen public PRIDE deposits that could replicate each artefact on independent data.

All rules come from protocols/phase1_dataset_screening_preregistration_2026-09-22.json; gate A0 asserts
that file is byte-identical to what its pre-registration audit recorded. This script adds only the
mechanics below, fixed before the run:

MECHANICS (implementation choices, declared here, not rules)
* A pool-1 accession's record is fetched from the search endpoint with keyword=<accession>; the hit
  whose accession equals it is used. If absent, /projects/<acc> and /projects/<acc>/files are used and
  normalised to the same fields. Records of pool 2 come straight from the search results.
* Text for every regex = lower-cased join of title, projectDescription, sampleProcessingProtocol,
  dataProcessingProtocol and keywords. Protease regexes read the two protocols only (a protease named in
  the description is often a citation of someone else's work).
* Modification family = the first of: persulfidation, s_nitrosylation, sulfenylation, sulfinylation,
  glutathionylation, s_acylation, redox_general, cys_reactivity, phosphorylation whose pattern matches.
* Species and lab-PI last names of existing instances are read from their own records at run time.
  Existing instances named only by DOI (artefact 1's census) contribute DOIs, not species: the
  species-differs test for A1 therefore compares against accession instances only (a limit, reported).
* Shared-reference link: DOIs (10.\\d{4,}/\\S+) and PubMed ids (7-9 digit numbers following 'pmid' or
  'pubmed', or given as a reference's pubmedId) extracted from references + text; two projects sharing
  one are linked. Text links: PXD\\d{6} mentioned in a project's text.
* Every HTTP response is cached under external/phase1_screening_20260922/cache/ and every search call is
  logged with its URL, count and first accessions. Sleep 0.2 s between calls; 4 retries with back-off.
* File sizes (files endpoint) are fetched only for the projects that enter the download proposal.

AMENDMENT A1 (written before the first run; no screening output existed)
The registered independence rule compares lab-PI LAST names. Surnames such as Wang, Li, Zhang or Duan
are shared by many unrelated laboratories, so a surname match would mark unrelated deposits as
'not independent' - a false-negative bias against exactly the datasets Phase 1 is looking for.
Implemented instead: a lab PI is identified by the lower-cased 'first last' full name; only labPIs are
used (submitters are not lab PIs, which the registered text also implies). The registered JSON is not
edited; A0 still checks it byte for byte.

OUTPUTS (new names)
  results/phase1_dataset_screening_candidates_2026-09-22.csv   one row per topic-coherent project
  results/phase1_dataset_screening_by_artefact_2026-09-22.csv  one row per project x artefact
  results/phase1_download_proposal_2026-09-22.csv              files proposed for approval (not fetched)
  results/phase1_dataset_screening_2026-09-22_audit.json
  external/phase1_screening_20260922/query_log.json + cache/
Interpreter: standard library only.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PREREG = os.path.join(ROOT, "protocols", "phase1_dataset_screening_preregistration_2026-09-22.json")
PREREG_AUDIT = os.path.join(ROOT, "results", "phase1_dataset_screening_prereg_2026-09-22_audit.json")
POOL1 = os.path.join(ROOT, "results", "pride_coincidence_candidates_2026-09-19.csv")
EXT = os.path.join(ROOT, "external", "phase1_screening_20260922")
CACHE = os.path.join(EXT, "cache")
OUT_CAND = os.path.join(ROOT, "results", "phase1_dataset_screening_candidates_2026-09-22.csv")
OUT_ART = os.path.join(ROOT, "results", "phase1_dataset_screening_by_artefact_2026-09-22.csv")
OUT_DL = os.path.join(ROOT, "results", "phase1_download_proposal_2026-09-22.csv")
OUT_AUD = os.path.join(ROOT, "results", "phase1_dataset_screening_2026-09-22_audit.json")
API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
ARTEFACTS = ["A1_cleavage", "A2_multicys", "A3_abundance", "A4_overlap", "A5_search_space"]
FAMILIES = [("persulfidation", r"persulfid|sulfhydrat|polysulfid|sulfhydrome|persulfidome"),
            ("s_nitrosylation", r"s-?nitros|nitrosyl|nitrosothiol"),
            ("sulfenylation", r"sulfenyl|dimedone"),
            ("sulfinylation", r"sulfinyl|sulfinic"),
            ("glutathionylation", r"glutathionyl"),
            ("s_acylation", r"palmitoyl|s-acyl|acyl[- ]biotin|\babe\b"),
            ("redox_general", r"redox|oxidation|oxicat|iodotmt|thiol"),
            ("cys_reactivity", r"isotop-abpp|reactive cysteine|cysteine reactiv|cysteinome"),
            ("phosphorylation", r"phospho")]
LOG = []


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def get(url, tag):
    key = hashlib.sha256(url.encode()).hexdigest()[:24]
    cp = os.path.join(CACHE, f"{tag}_{key}.json")
    if os.path.exists(cp):
        with open(cp, encoding="utf-8") as fh:
            return json.load(fh)
    last = None
    for i in range(4):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                data = json.loads(r.read())
            with open(cp, "w", encoding="utf-8") as fh:
                json.dump(data, fh, ensure_ascii=False)
            time.sleep(0.2)
            return data
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"{url}: {last}")


def as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def names(x):
    out = []
    for v in as_list(x):
        if isinstance(v, dict):
            out.append(str(v.get("name") or v.get("value") or v.get("referenceLine") or v.get("lastName") or ""))
        else:
            out.append(str(v))
    return [o for o in out if o]


def pi_names(x):
    """full names, lower-cased and whitespace-normalised (amendment A1)"""
    out = set()
    for v in as_list(x):
        if isinstance(v, dict):
            n = " ".join(str(v.get(k) or "") for k in ("firstName", "lastName")).strip() or str(v.get("name") or "")
        else:
            n = str(v)
        n = re.sub(r"\s+", " ", n).strip().lower()
        if n:
            out.add(n)
    return out


def refs(rec):
    s = []
    for v in as_list(rec.get("references")):
        if isinstance(v, dict):
            s += [str(v.get("referenceLine") or ""), str(v.get("doi") or ""), str(v.get("pubmedId") or v.get("pubmedID") or "")]
        else:
            s.append(str(v))
    s.append(str(rec.get("doi") or ""))
    return " ".join(s)


def ids_from(text):
    dois = {d.rstrip(".,;)").lower() for d in re.findall(r"10\.\d{4,}/[^\s\"';,]+", text)}
    pmids = set(re.findall(r"(?:pmid|pubmed)[^0-9]{0,12}(\d{7,9})", text, re.I))
    return dois, pmids


def normalise(rec, files=None):
    return {
        "accession": rec.get("accession"),
        "title": rec.get("title") or "",
        "projectDescription": rec.get("projectDescription") or "",
        "sampleProcessingProtocol": rec.get("sampleProcessingProtocol") or "",
        "dataProcessingProtocol": rec.get("dataProcessingProtocol") or "",
        "keywords": " ".join(names(rec.get("keywords"))),
        "organisms": names(rec.get("organisms")),
        "labPIs": pi_names(rec.get("labPIs")),
        "references": refs(rec),
        "submissionType": rec.get("submissionType") or "",
        "publicationDate": rec.get("publicationDate") or "",
        "files": files if files is not None else list(as_list(rec.get("projectFileNames"))),
        "ptms": " ".join(names(rec.get("identifiedPTMStrings"))),
    }


def fetch_by_accession(acc):
    res = get(f"{API}/search/projects?keyword={urllib.parse.quote(acc)}&pageSize=10&page=0", "acc")
    LOG.append({"call": "search_by_accession", "accession": acc, "n": len(res) if isinstance(res, list) else None})
    for r in as_list(res):
        if isinstance(r, dict) and r.get("accession") == acc:
            return normalise(r), "search"
    p = get(f"{API}/projects/{acc}", "proj")
    f = get(f"{API}/projects/{acc}/files?pageSize=1000&page=0", "files")
    return normalise(p, [x.get("fileName") for x in as_list(f) if isinstance(x, dict)]), "project_endpoint"


def main():
    with open(PREREG_AUDIT, encoding="utf-8") as fh:
        pa = json.load(fh)
    if sha(PREREG) != pa["preregistration_sha256"]:
        sys.exit("A0 FAIL: pre-registration changed after it was recorded")
    for p in (OUT_CAND, OUT_ART, OUT_DL, OUT_AUD):
        if os.path.exists(p):
            sys.exit(f"REFUSE: {p} exists")
    with open(PREREG, encoding="utf-8") as fh:
        P = json.load(fh)
    os.makedirs(CACHE, exist_ok=True)
    RX = {k: re.compile(v, re.I) for k, v in P["regex"].items() if isinstance(v, str)}
    NT = {k: re.compile(v, re.I) for k, v in P["regex"]["NONTRYPTIC"].items()}
    t0 = time.time()

    # ---------------- retrieval ----------------
    records, how, failed = {}, {}, []
    with open(POOL1, newline="", encoding="utf-8") as fh:
        pool1 = [r["accession"] for r in csv.DictReader(fh)]
    pool2_hits = []
    for kw in P["retrieval"]["pool_2_keywords"]:
        for page in range(P["retrieval"]["pool_2_pages_per_keyword"]):
            url = f"{API}/search/projects?keyword={urllib.parse.quote(kw)}&pageSize=100&page={page}"
            res = as_list(get(url, "kw"))
            LOG.append({"call": "search_keyword", "url": url, "keyword": kw, "page": page, "n": len(res),
                        "first": [r.get("accession") for r in res[:5] if isinstance(r, dict)]})
            for r in res:
                if isinstance(r, dict) and r.get("accession"):
                    pool2_hits.append((kw, r["accession"]))
                    records.setdefault(r["accession"], normalise(r))
                    how.setdefault(r["accession"], f"pool2:{kw}")
            if len(res) < 100:
                break
    wanted = list(dict.fromkeys(pool1 + P["retrieval"]["pool_3_on_disk"]
                                + [a for v in P["existing_instances"].values() for a in v["accessions"]]))
    for acc in wanted:
        if acc in records:
            continue
        try:
            records[acc], src = fetch_by_accession(acc)
            how[acc] = ("pool1" if acc in pool1 else "pool3") + ":" + src
        except Exception as e:  # noqa: BLE001
            failed.append({"accession": acc, "error": str(e)[:200]})

    def text(r):
        return " ".join(r[f] for f in P["text_fields"]).lower()

    def proto(r):
        return (r["sampleProcessingProtocol"] + " " + r["dataProcessingProtocol"]).lower()

    def proteases(r):
        t = proto(r)
        nt = [k for k, rx in NT.items() if rx.search(t)]
        return nt, bool(RX["TRYPSIN"].search(t)), bool(RX["LYSC"].search(t))

    # topic coherence
    coherent, dropped_topic, fuzzy_fail = {}, Counter(), Counter()
    for acc, r in records.items():
        t = text(r)
        cys = bool(RX["CYS_MOD"].search(t))
        nt, _, _ = proteases(r)
        phos_orth = bool(RX["PHOSPHO"].search(t)) and bool(nt)
        src = how.get(acc, "")
        if src.startswith("pool2") and not nt:
            fuzzy_fail[src.split(":", 1)[1]] += 1
        if cys or phos_orth or acc in P["retrieval"]["pool_3_on_disk"]:
            coherent[acc] = dict(r, cys=cys, phos_orth=phos_orth)
        else:
            dropped_topic[src.split(":")[0]] += 1
    # pool 4: text links, fetched once
    linked_added = []
    for acc, r in list(coherent.items()):
        for m in set(re.findall(r"PXD\d{6}", text(r).upper())):
            if m == acc or m in coherent:
                continue
            if m not in records:
                try:
                    records[m], src = fetch_by_accession(m)
                    how[m] = "pool4:" + src
                except Exception as e:  # noqa: BLE001
                    failed.append({"accession": m, "error": str(e)[:200]})
                    continue
            rr = records[m]
            coherent[m] = dict(rr, cys=bool(RX["CYS_MOD"].search(text(rr))), phos_orth=False, linked_from=acc)
            linked_added.append((acc, m))

    # links (text mentions + shared DOI/PMID)
    idmap = defaultdict(set)
    for acc, r in coherent.items():
        d, p = ids_from(r["references"] + " " + text(r))
        for i in d | {f"pmid:{x}" for x in p}:
            idmap[i].add(acc)
    links = defaultdict(set)
    for acc, r in coherent.items():
        for m in set(re.findall(r"PXD\d{6}", text(r).upper())) - {acc}:
            links[acc].add(m)
            links[m].add(acc)
    for i, accs in idmap.items():
        for a in accs:
            links[a] |= accs - {a}

    # existing instances
    exist = {}
    for art, v in P["existing_instances"].items():
        sp, pis = set(), set()
        for a in v["accessions"]:
            r = records.get(a)
            if r:
                sp |= {o.lower() for o in r["organisms"]}
                pis |= r["labPIs"]
        exist[art] = {"acc": set(v["accessions"]), "dois": {d.lower() for d in v["dois"]}, "species": sp, "pis": pis}
    exist_fams = {"A1_cleavage": {"s_nitrosylation", "sulfenylation", "persulfidation"},
                  "A2_multicys": {"persulfidation"}, "A3_abundance": {"persulfidation"},
                  "A4_overlap": {"s_acylation"}, "A5_search_space": {"persulfidation"}}

    # ---------------- classification ----------------
    ext = P["file_ext"]

    def fclass(files):
        low = [str(f).lower() for f in files]
        raw = [f for f in low if any(f.endswith(e) for e in ext["raw"])]
        nonraw = [f for f in low if f not in raw]
        tables = [f for f in nonraw if any(f.endswith(e) for e in ext["table"])]
        archives = [f for f in nonraw if any(f.endswith(e) for e in ext["archive"])]
        site_t = [f for f in tables if RX["SITE_TABLE_NAME"].search(os.path.basename(f))]
        pep_t = [f for f in tables if RX["PEPTIDE_TABLE_NAME"].search(os.path.basename(f))]
        search_a = [f for f in archives if RX["SEARCH_ARCHIVE_NAME"].search(os.path.basename(f))]
        return {"n_files": len(low), "n_raw": len(raw), "n_nonraw": len(nonraw), "n_tables": len(tables),
                "n_archives": len(archives), "site_tables": site_t, "peptide_tables": pep_t, "search_archives": search_a}

    def family(r):
        t = text(r)
        for f, pat in FAMILIES:
            if re.search(pat, t, re.I):
                return f
        return "unclassified"

    cand_rows, art_rows = [], []
    for acc in sorted(coherent):
        r = coherent[acc]
        t = text(r)
        fc = fclass(r["files"])
        nt, tryp, lysc = proteases(r)
        fam = family(r)
        enrich = bool(RX["ENRICH"].search(t))
        block = bool(RX["BLOCK_ARM"].search(t))
        comp_text = bool(RX["COMPANION"].search(t))
        plus32 = bool(RX["PLUS32"].search(t))
        musts = {"M1": True, "M2": bool(re.fullmatch(r"PXD\d{6}", acc or "")),
                 "M3": r["cys"] or r["phos_orth"], "M4": fc["n_nonraw"] >= 1, "M5": bool(r["organisms"])}
        linked_nonenrich = sorted(m for m in links.get(acc, set())
                                  if m in coherent and not RX["ENRICH"].search(text(coherent[m])))
        orth = ("HIGH" if nt else "MEDIUM" if (lysc and not tryp) else "LOW" if tryp else "UNSTATED")
        cand_rows.append({
            "accession": acc, "retrieved_by": how.get(acc, ""), "title": r["title"][:160],
            "organisms": "; ".join(r["organisms"]), "family": fam, "cys_topic": r["cys"],
            "phospho_orthogonal": r["phos_orth"], "proteases_nontryptic": "; ".join(nt), "trypsin": tryp, "lysc": lysc,
            "orthogonal_level": orth, "enrich_design": enrich, "block_arm": block, "companion_text": comp_text,
            "plus32": plus32, "submissionType": r["submissionType"], "publicationDate": r["publicationDate"],
            "n_files": fc["n_files"], "n_raw": fc["n_raw"], "n_nonraw": fc["n_nonraw"],
            "site_table_files": "; ".join(fc["site_tables"][:5]), "peptide_table_files": "; ".join(fc["peptide_tables"][:5]),
            "search_archives": "; ".join(fc["search_archives"][:5]), "linked_accessions": "; ".join(sorted(links.get(acc, set()))),
            "linked_nonenriched": "; ".join(linked_nonenrich), "lab_pis": "; ".join(sorted(r["labPIs"])),
            "musts_failed": "; ".join(k for k, v in musts.items() if not v)})
        for art in ARTEFACTS:
            E = exist[art]
            indep = (acc not in E["acc"] and not (r["labPIs"] & E["pis"])
                     and not any(d in (r["references"] + " " + t).lower() for d in E["dois"]))
            sp_diff = bool(r["organisms"]) and not ({o.lower() for o in r["organisms"]} & E["species"])
            fam_diff = fam not in exist_fams[art]
            reason, core, pending = "", False, False
            if not all(musts.values()):
                reason = "must failed: " + ",".join(k for k, v in musts.items() if not v)
            elif art == "A1_cleavage":
                core = bool(fc["site_tables"] or fc["search_archives"] or fc["n_tables"])
                reason = "" if core else "no result table or archive"
            elif not r["cys"]:
                reason = "not a cysteine-modification deposit"
            elif art == "A2_multicys":
                core = bool(fc["peptide_tables"]) and not fc["site_tables"]
                pending = (not core) and bool(fc["n_archives"]) and not fc["site_tables"]
                reason = "" if core else ("archive contents unknown" if pending else "site table present or no peptide table")
            elif art == "A3_abundance":
                core = enrich and (comp_text or bool(linked_nonenrich))
                reason = "" if core else ("not an enrichment design" if not enrich else "no companion quantification found")
            elif art == "A4_overlap":
                core = enrich and bool(fc["site_tables"]) and bool(fc["peptide_tables"])
                pending = (not core) and enrich and bool(fc["search_archives"])
                reason = "" if core else ("search archive contents unknown" if pending else ("not an enrichment design" if not enrich else "no site+peptide tables"))
            elif art == "A5_search_space":
                core = plus32 and fc["n_raw"] >= 1
                reason = "" if core else ("no +32-class chemistry" if not plus32 else "no raw files")
            if core or pending:
                if not indep:
                    verdict = "NOT_INDEPENDENT"
                elif pending:
                    verdict = "PENDING_FILE_INSPECTION"
                elif (sp_diff or fam_diff) and (art != "A1_cleavage" or orth == "HIGH"):
                    verdict = "CANDIDATE_STRONG"
                else:
                    verdict = "CANDIDATE"
            else:
                verdict = "NOT_SUITABLE"
            art_rows.append({"accession": acc, "artefact": art, "verdict": verdict, "reason": reason,
                             "independent": indep, "species_differs": sp_diff, "family_differs": fam_diff,
                             "orthogonal_level": orth if art == "A1_cleavage" else "", "block_arm": block if art == "A4_overlap" else "",
                             "submissionType": r["submissionType"], "family": fam, "organisms": "; ".join(r["organisms"])})

    # ---------------- controls ----------------
    pc_ret = {a: (a in coherent) for a in P["controls"]["PC_retrieval"]["accessions"]}
    a, b = P["controls"]["PC_linking_A3"]["pair"]
    pc_link = (b in links.get(a, set())) or (a in links.get(b, set()))
    V = {(x["accession"], x["artefact"]): x for x in art_rows}
    pc_cls = {
        "PXD048216_A4": V.get(("PXD048216", "A4_overlap"), {}).get("verdict") in ("NOT_INDEPENDENT", "CANDIDATE", "CANDIDATE_STRONG", "PENDING_FILE_INSPECTION"),
        "PXD015307_A5": V.get(("PXD015307", "A5_search_space"), {}).get("verdict") in ("NOT_INDEPENDENT", "CANDIDATE", "CANDIDATE_STRONG"),
        "PXD072089_A2": V.get(("PXD072089", "A2_multicys"), {}).get("verdict") in ("NOT_INDEPENDENT", "CANDIDATE", "CANDIDATE_STRONG"),
    }
    a3_rule_downgrade = False
    if not pc_link:  # registered consequence
        a3_rule_downgrade = True
        for x in art_rows:
            if x["artefact"] == "A3_abundance" and x["reason"] == "no companion quantification found" and x["independent"]:
                x["verdict"], x["reason"] = "PENDING_FILE_INSPECTION", "no companion found; linking control failed (registered downgrade)"

    # ---------------- ranking and download proposal ----------------
    order = {"CANDIDATE_STRONG": 0, "CANDIDATE": 1, "PENDING_FILE_INSPECTION": 2}
    olev = {"HIGH": 0, "MEDIUM": 1, "LOW": 2, "UNSTATED": 3, "": 4}
    dl_rows, shortlist = [], {}
    size_cache = {}
    cand_by_acc = {c["accession"]: c for c in cand_rows}

    def files_of(acc):
        if acc not in size_cache:
            f = as_list(get(f"{API}/projects/{acc}/files?pageSize=1000&page=0", "files"))
            size_cache[acc] = [x for x in f if isinstance(x, dict)]
        return size_cache[acc]

    for art in ARTEFACTS:
        pool = [x for x in art_rows if x["artefact"] == art and x["verdict"] in order]
        pool.sort(key=lambda x: (order[x["verdict"]], olev[x["orthogonal_level"]] if art == "A1_cleavage" else 0,
                                 not x["species_differs"], not x["family_differs"], x["submissionType"] != "COMPLETE", x["accession"]))
        # size tie-break needs the files endpoint: fetch for the top 10 only, re-sort, keep top k
        top = pool[:10]
        for x in top:
            fl = files_of(x["accession"])
            need = [f for f in fl if ((f.get("fileCategory") or {}).get("value") == "RAW") == (art == "A5_search_space")]
            if art == "A3_abundance":  # the companion deposit's non-raw files are needed too (pre-registration)
                for m in [c for c in cand_by_acc[x["accession"]]["linked_nonenriched"].split("; ") if c][:1]:
                    need += [f for f in files_of(m) if (f.get("fileCategory") or {}).get("value") != "RAW"]
            x["_need_bytes"] = sum(int(f.get("fileSizeBytes") or 0) for f in need)
            x["_need"] = need
        top.sort(key=lambda x: (order[x["verdict"]], olev[x["orthogonal_level"]] if art == "A1_cleavage" else 0,
                                not x["species_differs"], not x["family_differs"], x["submissionType"] != "COMPLETE",
                                x["_need_bytes"], x["accession"]))
        k = P["download_proposal"]["top_k_per_artefact"]
        shortlist[art] = [x["accession"] for x in top[:k]]
        for rank, x in enumerate(top[:k], 1):
            for f in x["_need"]:
                locs = [l.get("value") for l in as_list(f.get("publicFileLocations")) if isinstance(l, dict)]
                url = next((u for u in locs if str(u).startswith(("ftp://", "http"))), locs[0] if locs else "")
                owner = f.get("projectAccessions") or [x["accession"]]
                dl_rows.append({"artefact": art, "rank": rank, "accession": x["accession"],
                                "file_owner": owner[0] if isinstance(owner, list) and owner else x["accession"], "verdict": x["verdict"],
                                "file_name": f.get("fileName"), "category": (f.get("fileCategory") or {}).get("value"),
                                "bytes": int(f.get("fileSizeBytes") or 0), "url": url, "approved": "no"})
        for rank, x in enumerate(pool, 1):
            x["rank_in_artefact"] = rank
    for x in art_rows:
        x.pop("_need", None)
        x.pop("_need_bytes", None)
        x.setdefault("rank_in_artefact", "")

    # ---------------- write ----------------
    for path, rows in ((OUT_CAND, cand_rows), (OUT_ART, art_rows), (OUT_DL, dl_rows)):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            fields = list(rows[0].keys()) if rows else ["empty"]
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
    with open(os.path.join(EXT, "query_log.json"), "w", encoding="utf-8") as fh:
        json.dump(LOG, fh, ensure_ascii=False, indent=0)
    vc = {art: dict(Counter(x["verdict"] for x in art_rows if x["artefact"] == art)) for art in ARTEFACTS}
    audit = {
        "script": "scripts/screen_phase1_datasets_2026-09-22.py", "script_sha256": sha(os.path.abspath(__file__)),
        "preregistration": os.path.relpath(PREREG, ROOT), "preregistration_sha256": sha(PREREG), "A0": "pass",
        "python": sys.version,
        "n_records_fetched": len(records), "n_topic_coherent": len(coherent), "fetch_failures": failed,
        "dropped_by_topic": dict(dropped_topic), "pool2_hits": len(pool2_hits),
        "NC_fuzzy_search_pool2_no_protease_in_protocol": dict(fuzzy_fail),
        "pool4_linked_added": linked_added,
        "controls": {"PC_retrieval": pc_ret, "PC_retrieval_recall": f"{sum(pc_ret.values())}/{len(pc_ret)}",
                     "PC_linking_A3": pc_link, "A3_registered_downgrade_applied": a3_rule_downgrade,
                     "PC_classifier": pc_cls,
                     "recall_undefined": P["controls"]["recall_undefined"]},
        "existing_instances_resolved": {k: {"species": sorted(v["species"]), "lab_pis": sorted(v["pis"])} for k, v in exist.items()},
        "verdict_counts_by_artefact": vc, "shortlist": shortlist,
        "download_proposal_bytes_by_artefact": {art: sum(d["bytes"] for d in dl_rows if d["artefact"] == art) for art in ARTEFACTS},
        "limits": ["A1 species-differs compares against accession instances only; the census datasets are DOI-only",
                   "classification reads metadata and file names; archive contents are not inspected",
                   "PRIDE only; MassIVE, jPOST, iProX not searched"],
        "outputs": {os.path.relpath(p, ROOT): sha(p) for p in (OUT_CAND, OUT_ART, OUT_DL)},
        "elapsed_s": round(time.time() - t0, 1),
    }
    with open(OUT_AUD, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=1, allow_nan=False)
    print(json.dumps({k: audit[k] for k in ("n_records_fetched", "n_topic_coherent", "controls", "verdict_counts_by_artefact",
                                            "shortlist", "download_proposal_bytes_by_artefact", "elapsed_s")}, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
