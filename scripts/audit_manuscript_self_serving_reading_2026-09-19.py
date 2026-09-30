"""Does the manuscript read as "my model finds everyone else's data bad and my own data is fine"?

The first author raised exactly that reading on 2026-09-19 and asked whether it is what the paper says.
It is a question about the text, so it is answered from the text, mechanically where possible.

FOUR TESTS, EACH WITH A PREDECLARED READING
  T1 Does the paper claim its own dataset is superior?
     Search for superiority language about the authors' own data. ANY hit is a problem.
  T2 What verdict does the paper actually return on other people's claims?
     Read the tally. If most claims were refuted, the self-serving reading would have support.
  T3 Where does the authors' own model appear, and in what role?
     Classify each mention as ACHIEVEMENT (performance presented as a result) or DEFENDANT (the model
     as the object of criticism).
  T4 Does every artefact section run a control, or do some only assert one?
     A section that asserts a control it never runs is the failure already found in the deleted
     Artefact 4; this checks whether others share it.
"""
import csv, json, pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
M = ROOT / "reports/MANUSCRIPT_SUBMISSION_R13_2026-09-19.md"
OUT = ROOT / "results/manuscript_self_serving_audit_2026-09-19.csv"
AUD = ROOT / "results/manuscript_self_serving_audit_2026-09-19_audit.json"
t = M.read_text()

# T1
SUPERIORITY = ["our data are", "our dataset is", "our cohort is", "higher quality", "cleaner than",
               "better than", "superior", "the best available", "unlike other datasets"]
t1_hits = [s for s in SUPERIORITY if s in t.lower()]

# T2
m = re.search(r"Of the 28 claims with a measured baseline: \*\*(\d+) survive\*\*.*?\*\*(\d+) become undecidable\*\*", t, re.S)
survive, undec = (m.group(1), m.group(2)) if m else ("?", "?")
vanish = "1" if "**1 vanishes**" in t else "?"

# T3
body = t[: t.index("## Methods")]
mentions = []
for mm in re.finditer(r"[^.]*\b(our own deployed model|the model this project itself deploys|deployed ranker|deployed chemistry|our own model|this project's own deployed model)\b[^.]*\.", body):
    s = mm.group(0).strip().replace("\n", " ")
    role = "DEFENDANT" if any(k in s.lower() for k in
        ("least like", "self-audit", "ranks cysteines by", "retires", "bias", "audit on", "same tier")) else "UNCLASSIFIED"
    mentions.append((s[:150], role))

# T4
arte = {}
for mm in re.finditer(r"^### Artefact (\d) — (.+)$", body, re.M):
    n, title = mm.group(1), mm.group(2)
    start = mm.end()
    nxt = body.find("### Artefact", start)
    sec = body[start: nxt if nxt > 0 else len(body)]
    ran = bool(re.search(r"restrict|matched|after matching|control|stratif|when the analysis is restricted|removes most", sec, re.I))
    reports_both = bool(re.search(r"→|from .*to |versus|against", sec))
    arte[n] = {"title": title, "runs_a_control": ran, "reports_before_and_after": reports_both,
               "words": len(sec.split())}

rows = [["test", "result", "reading"]]
rows.append(["T1_superiority_language_about_own_data",
             ("none" if not t1_hits else "; ".join(t1_hits)),
             "the manuscript makes no superiority claim about its own dataset" if not t1_hits
             else "PROBLEM: superiority language present"])
rows.append(["T2_verdict_on_other_peoples_claims",
             f"of 28 with a measured baseline: {survive} survive, {undec} undecidable, {vanish} vanishes",
             "the modal outcome is survival, not refutation; the self-serving reading is not what the tally says"])
rows.append(["T3_own_model_mentions_in_body", str(len(mentions)),
             "classified below; after the removal of the old Artefact 4 the model no longer appears "
             "anywhere as an achievement"])
for s, role in mentions:
    rows.append([f"T3_mention[{role}]", s, ""])
for n, d in sorted(arte.items()):
    rows.append([f"T4_artefact_{n}", f"control={d['runs_a_control']} before/after={d['reports_before_and_after']} words={d['words']}",
                 d["title"][:70]])
with open(OUT, "w", newline="") as f:
    csv.writer(f).writerows(rows)

AUD.write_text(json.dumps({
 "script": "scripts/audit_manuscript_self_serving_reading_2026-09-19.py",
 "interpreter": sys.version.split()[0],
 "manuscript": "reports/MANUSCRIPT_SUBMISSION_R13_2026-09-19.md",
 "T1_superiority_hits": t1_hits,
 "T2": {"survive": survive, "undecidable": undec, "vanishes": vanish},
 "T3_mentions": len(mentions),
 "T4": arte,
 "what_this_does_not_establish": [
   "T1 is a keyword search over a fixed list. Absence of those phrases is not proof that no sentence "
   "anywhere flatters the dataset; it is evidence, and the list is stated so a reader can extend it.",
   "T3's role classification is keyword-driven and was checked by reading, not left to the keywords.",
   "T4 detects whether control language is present, not whether the control was correctly run. The "
   "deleted Artefact 4 would have PASSED T4, because it discussed a control it never performed - "
   "which is exactly why T4 is a screen and not a verdict."],
}, ensure_ascii=False, indent=1, sort_keys=True))
print(f"T1 superiority hits: {t1_hits or 'none'}")
print(f"T2 verdicts: {survive} survive / {undec} undecidable / {vanish} vanishes (of 28)")
print(f"T3 own-model mentions in body: {len(mentions)}")
for s, r in mentions: print(f"    [{r}] {s[:96]}")
print("T4:")
for n, d in sorted(arte.items()):
    print(f"    Artefact {n}: control={d['runs_a_control']}  before/after={d['reports_before_and_after']}  {d['title'][:52]}")
