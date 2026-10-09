import csv, collections, os

SD9 = 'supplemental/Supplemental_Data_9_survey_coding_table.csv'
A = '_recoding/A/codes_A.csv'
B = '_recoding/B/codes_B.csv'

sd9 = list(csv.DictReader(open(SD9, encoding='utf-8-sig')))
cls = [r for r in sd9 if r['background_class'] in ('a', 'b', 'c', 'e')]
a = {r['uid']: r for r in csv.DictReader(open(A, encoding='utf-8-sig'))}
b = {r['uid']: r for r in csv.DictReader(open(B, encoding='utf-8-sig'))}

NOTE = 'explicit background statement found on full-text re-grep'
NOSTMT = 'no background statement anywhere in full text'
note_uids = [r['uid'] for r in cls if r['verification'].startswith(NOTE)]
nostmt_uids = [r['uid'] for r in cls if r['verification'].startswith(NOSTMT)]


def kappa(pairs):
    cats = sorted({x for p in pairs for x in p})
    n = len(pairs)
    po = sum(1 for x, y in pairs if x == y) / n
    ca = collections.Counter(x for x, _ in pairs)
    cb = collections.Counter(y for _, y in pairs)
    pe = sum(ca[c] * cb[c] for c in cats) / (n * n)
    return po, pe, (po - pe) / (1 - pe) if pe < 1 else float('nan')


pairs_a = [(r['background_class'], a[r['uid']]['coder_code']) for r in cls]
pairs_b = [(r['background_class'], b[r['uid']]['coder_code']) for r in cls]
print('=== agreement (all 74) ===')
for name, pr in [('A', pairs_a), ('B', pairs_b)]:
    po, pe, k = kappa(pr)
    print(f'  coder {name}: raw {sum(1 for x,y in pr if x==y)}/{len(pr)} = {po:.3f}  kappa {k:.3f}')

print('=== sensitivity: 12 note rows credited to the original code ===')
for name, src in [('A', a), ('B', b)]:
    pr = [('a', 'a') if r['uid'] in note_uids else (r['background_class'], src[r['uid']]['coder_code']) for r in cls]
    po, pe, k = kappa(pr)
    print(f'  coder {name}: raw {sum(1 for x,y in pr if x==y)}/{len(pr)} = {po:.3f}  kappa {k:.3f}')

# truncated quotations
def truncated(q):
    q = (q or '').strip()
    return bool(q) and q[-1] not in '.!?'
trunc = sum(1 for r in cls if truncated(r['evidence_quote_background']))
empty = sum(1 for r in cls if not (r['evidence_quote_background'] or '').strip())
print(f'\nretained background quotations: {trunc} of {len(cls)} cut off mid-clause, {empty} empty')

# ---- build supplemental audit table (74 rows) ----
os.makedirs('supplemental', exist_ok=True)
out = 'supplemental/Supplemental_Data_11_survey_coding_recoding_audit.csv'
with open(out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['uid', 'doi', 'pmid', 'year', 'journal', 'modification', 'study_type', 'title',
                'original_code', 'code_pass_1', 'code_pass_2',
                'pass_1_agrees_with_original', 'pass_2_agrees_with_original', 'both_passes_agree_with_original',
                'verification_note', 'ambiguity_note', 'background_quotation_retained',
                'quotation_truncated', 'justification_pass_1', 'justification_pass_2', 'confidence_pass_1', 'confidence_pass_2'])
    for r in sorted(cls, key=lambda x: int(x['uid'])):
        ca, cb = a[r['uid']], b[r['uid']]
        w.writerow([r['uid'], r['doi'], r['pmid'], r['year'], r['journal'], r['modification'],
                    r['study_type'], r['title'], r['background_class'], ca['coder_code'], cb['coder_code'],
                    str(ca['coder_code'] == r['background_class']).lower(),
                    str(cb['coder_code'] == r['background_class']).lower(),
                    str(ca['coder_code'] == r['background_class'] and cb['coder_code'] == r['background_class']).lower(),
                    r['verification'], r['ambiguity'], r['evidence_quote_background'],
                    str(truncated(r['evidence_quote_background'])).lower(),
                    ca['justification'], cb['justification'], ca['confidence'], cb['confidence']])
print('wrote', out)

# ---- build the human re-check checklist ----
bc = [r for r in cls if r['background_class'] in ('b', 'c')]
rows = []
for r in bc:
    rows.append(('1. headline count', r))
for r in cls:
    if r['uid'] in note_uids:
        rows.append(('2. statement not retained', r))
for r in cls:
    if r['uid'] in nostmt_uids:
        rows.append(('3. no statement found', r))

ACTION = {
    '1. headline count': 'Re-read the paper: find the sentence that says which residues the modification-preference comparison used as its reference. Confirmed as matched ONLY if it names (i) the same run\'s detected-but-unassigned residues, (ii) a detectability-matched set, or (iii) an abundance-matched set. A contrast between two subsets of sites (e.g. regulated vs all identified) is NOT a matched background; negatives drawn from an external database are NOT.',
    '2. statement not retained': 'Re-grep the full text for the background sentence (the note says one was found) and paste it verbatim into evidence_quote_background. Then assign a / b / c. If the recovered sentence turns out to be a / e, the original code stands.',
    '3. no statement found': 'Confirm that both regex passes really did miss nothing: search the methods and the motif-section text for background, reference set, negative control, random sequences. If nothing, e stands.',
}

md = ['# Manual re-check table — survey coding reliability',
      '',
      'Companion to `Supplemental_Data_11_survey_coding_recoding_audit.csv` (all 74 classified records).',
      'This file lists only the **21 records that a human must look at**. Everything else can be left as coded:',
      'the two blinded passes agreed with the original code, or disagreed only because the retained quotation is truncated.',
      '',
      'Groups:',
      f'- **Group 1 — the 5 records the headline count rests on** (originally coded b or c). Both blinded passes reject 3 of these 5.',
      f'- **Group 2 — {len(note_uids)} records** whose note says a background statement *was* found in the full text on re-grep, but the statement itself was not retained (all originally coded `a`).',
      f'- **Group 3 — {len(nostmt_uids)} records** whose note says no background statement exists anywhere (originally coded `e`).',
      '',
      'For each row below: read the paper, record the verbatim background sentence, assign one of',
      '`a` (all residues of that type, unrestricted) / `b` (matched on detectability or abundance) /',
      '`c` (residues detected in the same experiment without an assignment) / `e` (other, or none stated),',
      'and note whether it changes the count. Group 3 is a 5-minute confirmation; Group 1 decides whether the',
      'headline stays at 5 or becomes 2–3.',
      '']

for tag, r in rows:
    md.append(f'## {tag.split(".")[0]}.{r["uid"]} — `{r["background_class"]}` (original)')
    md.append(f'- **DOI**: {r["doi"] or "n/a"} | PMID {r["pmid"] or "n/a"} | {r["year"]} | {r["journal"]}')
    md.append(f'- **Title**: {r["title"]}')
    md.append(f'- **Modification / study type**: {r["modification"]} / {r["study_type"]}')
    md.append(f'- **Original code**: `{r["background_class"]}` — pass 1: `{a[r["uid"]]["coder_code"]}`, pass 2: `{b[r["uid"]]["coder_code"]}`')
    md.append(f'- **Retained quotation**: {(r["evidence_quote_background"] or "(empty)")[:300]}')
    if (r['ambiguity'] or '').strip():
        md.append(f'- **Original coder\'s own ambiguity note**: {r["ambiguity"]}')
    md.append(f'- **Pass 1 reading**: {a[r["uid"]]["coder_code"]} — {a[r["uid"]]["justification"]}')
    md.append(f'- **Pass 2 reading**: {b[r["uid"]]["coder_code"]} — {b[r["uid"]]["justification"]}')
    md.append(f'- **What to check**: {ACTION[tag]}')
    md.append('- **Your call**: `____` (a/b/c/e) | quotation: ____________________')
    md.append('')

md += ['## Roll-up (fill in when done)',
       '',
       '| group | records | changed away from original `b`/`c` | new count of matched-background papers |',
       '|---|---|---|---|',
       f'| 1. headline count | {len(bc)} | | |',
       f'| 2. statement not retained | {len(note_uids)} | n/a | n/a |',
       f'| 3. no statement found | {len(nostmt_uids)} | n/a | n/a |',
       '',
       'Headline numerator after the check: **2 / 3 / 4 / 5** of 74 — circle one.',
       '']
os.makedirs('_recoding', exist_ok=True)
open('_recoding/recheck_table.md', 'w', encoding='utf-8').write('\n'.join(md))
print('wrote _recoding/recheck_table.md with', len(rows), 'records')
