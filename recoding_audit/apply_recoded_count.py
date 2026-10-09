"""Apply the three full-text re-codings to Supplemental Data 9, Fig. 1c Source Data and Supplemental Data 11.

v3.1.1: paths are resolved from the repository root (run as `python recoding_audit/apply_recoded_count.py`
from anywhere); the Fig. 1 Source Data path is `source_data_submitted/` (the folder the repository ships;
the original pointed at the manuscript package's `source_data/`). The edit is idempotent: on the released
tables, which already carry the re-codings, it changes nothing and only prints the verification counts.
Files are written with LF line endings, as committed.
"""
import csv, os

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))

CORR = {
    '1111': ('c', 'e', 'Re-coded to e: foreground is regulated sites and background is all identified phosphosites, a contrast between two subsets of sites rather than against unassigned residues; both blinded re-coding passes read e.'),
    '1537': ('b', 'e', 'Re-coded to e: the abundance control was applied to sites-per-protein, not to a sequence-level background; both blinded re-coding passes read e.'),
    '1553': ('c', 'a', 'Re-coded to a: negatives were randomly drawn from identified human lysines held in an external database (MAPU), an unrestricted lysine pool rather than same-experiment detections; re-coding pass 2 read a and pass 1 read e, both rejecting c.'),
}

# ---------- 1. Supplemental Data 9 ----------
p = 'supplemental/Supplemental_Data_9_survey_coding_table.csv'
rows = list(csv.DictReader(open(p, encoding='utf-8-sig')))
cols = list(rows[0].keys())
if 'background_class_original' not in cols:
    cols.insert(cols.index('background_class_label') + 1, 'background_class_original')
if 'recoding_note' not in cols:
    cols.append('recoding_note')
changed = 0
for r in rows:
    r.setdefault('background_class_original', '')
    r.setdefault('recoding_note', '')
    if r['uid'] in CORR:
        old, new, note = CORR[r['uid']]
        if r['background_class'] == new and r['background_class_original'] == old:
            continue  # already applied (released table)
        assert r['background_class'] == old, (r['uid'], r['background_class'])
        r['background_class_original'] = old
        r['background_class'] = new
        r['recoding_note'] = note
        # keep the human-readable label in step with the class
        r['background_class_label'] = {
            'a': 'a: all residues of that type (proteome / all proteins), unrestricted',
            'b': 'b: matched on theoretical detectability or abundance',
            'c': 'c: residues detected in the same experiment without a modification assignment',
            'e': 'e: other, or no statement found',
        }[new]
        changed += 1
with open(p, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=cols, lineterminator='\n')
    w.writeheader()
    w.writerows(rows)
print(f'SD9: {changed} rows re-coded, columns now {len(cols)}')

# ---------- 2. Source data for Fig 1c ----------
p = 'source_data_submitted/Source_Data_Fig1_overview.csv'
rows = list(csv.DictReader(open(p, encoding='utf-8-sig')))
cols = list(rows[0].keys())
newvals = {'Other or not stated': '37', 'All residues, unrestricted': '35',
           'Detected but unmodified': '2', 'Detectability or abundance matched': '0'}
for r in rows:
    if r.get('figure') == 'Fig1' and r.get('panel') == 'c' and r.get('row') in newvals:
        if r['value'] != newvals[r['row']]:
            print(f"  Fig1c {r['row']}: {r['value']} -> {newvals[r['row']]}")
        r['value'] = newvals[r['row']]
with open(p, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=cols, lineterminator='\n')
    w.writeheader()
    w.writerows(rows)
print('Fig1 source data updated')

# ---------- 3. Supplemental Data 11 ----------
p = 'supplemental/Supplemental_Data_11_survey_coding_recoding_audit.csv'
rows = list(csv.DictReader(open(p, encoding='utf-8-sig')))
cols = list(rows[0].keys())
for c in ['final_class', 'full_text_check']:
    if c not in cols:
        cols.append(c)
for r in rows:
    r.setdefault('final_class', '')
    r.setdefault('full_text_check', '')
    if r['uid'] in CORR:
        r['final_class'] = CORR[r['uid']][1]
        r['full_text_check'] = 're-read in full text; class changed'
with open(p, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=cols, lineterminator='\n')
    w.writeheader()
    w.writerows(rows)
print('SD11 updated with final_class / full_text_check')

# ---------- verify ----------
import collections
rows = list(csv.DictReader(open('supplemental/Supplemental_Data_9_survey_coding_table.csv', encoding='utf-8-sig')))
cls = [r for r in rows if r['background_class'] in ('a', 'b', 'c', 'e')]
dist = collections.Counter(r['background_class'] for r in cls)
print('\nSD9 corrected distribution:', dict(dist), '| classified rows:', len(cls))
print('matched-background count (b + c):', dist['b'] + dist['c'])
print('overall rows:', len(rows))
