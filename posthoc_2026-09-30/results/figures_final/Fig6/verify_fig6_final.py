"""Independent read-back of the final Figure 6 PDF (run after fig6_final_2026-09-30.py). Reads only; writes nothing.

Nothing is taken from the plotting script: every position is measured from the vector content of
Fig6_claim_retests.pdf and mapped to data units through the drawn x tick marks, then compared with the Source Data.
  1  per claim: matched estimate (filled circle), baseline (open circle), matched 95% interval (thin line end points),
     shift segment (pale line) and the SFE-006 transfer diamond and interval
  2  colours: circle, interval, shift and strip cell of each row = released verdict colour of its final verdict
  3  row order: labels top to bottom = ranking by matched estimate within the site and protein blocks
  4  page, fonts and sizes
"""
import ast
import csv
import os
import re
import sys

import fitz

sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, 'Fig6_claim_retests.pdf')
SD = os.path.join(HERE, 'Source_Data_Fig6_claim_retests.csv')
FIG6_PY = '/path/to/local/巯基化/MCP/_revision/fig6_final.py'
TOL = 1e-4                         # log2 units; PDF coordinates carry ~1e-4 pt, i.e. ~4e-6 log2 units
results = []


def check(name, ok, detail=''):
    results.append(bool(ok))
    print(('PASS ' if ok else 'FAIL ') + name + (': ' + str(detail) if detail != '' else ''))


tree = ast.parse(open(FIG6_PY, encoding='utf-8').read())
VC = [ast.literal_eval(n.value) for n in ast.walk(tree) if isinstance(n, ast.Assign)
      and any(getattr(t, 'id', None) == 'VC' for t in n.targets)][0]
num, cat, panel = {}, {}, {}
for r in csv.DictReader(open(SD, encoding='utf-8', newline='')):
    if r['field'] == 'n_claims_in_block':
        continue
    panel[r['row']] = r['panel']
    try:
        num.setdefault(r['row'], {})[r['field']] = float(r['value'])
    except ValueError:
        cat.setdefault(r['row'], {})[r['field']] = r['value']

doc = fitz.open(PDF)
page = doc[0]
H = page.rect.height
spans = [s for b in page.get_text('dict')['blocks'] for ln in b.get('lines', []) for s in ln['spans']
         if s['text'].strip()]
dr = page.get_drawings()


def close(a, b, t=0.003):
    return a is not None and b is not None and len(a) == len(b) and max(abs(x - y) for x, y in zip(a, b)) < t


# ---- axes calibration from the drawn x tick marks (short vertical strokes) and their labels
xt = [d for d in dr if len(d['items']) == 1 and d['items'][0][0] == 'l' and abs(d['items'][0][1].x -
      d['items'][0][2].x) < 1e-6 and 1.5 < abs(d['items'][0][1].y - d['items'][0][2].y) < 2.5]
lab = {}
for s in spans:
    t = s['text'].replace('\u2212', '-')
    if re.fullmatch(r'-?\d+(\.\d+)?', t):
        lab[(round((s['bbox'][0] + s['bbox'][2]) / 2, 1), round(s['bbox'][1], 1))] = float(t)


def scale(ticks):
    """least-squares x = a + b * value over (tick x, label value) pairs; returns value(x)."""
    pairs = []
    for d in ticks:
        x = d['items'][0][1].x
        cand = [(abs(k[0] - x), v) for k, v in lab.items() if abs(k[0] - x) < 1.0 and k[1] > d['rect'].y1]
        pairs.append((x, min(cand)[1]))
    n = len(pairs)
    mx, mv = sum(p[0] for p in pairs) / n, sum(p[1] for p in pairs) / n
    b = sum((p[1] - mv) * (p[0] - mx) for p in pairs) / sum((p[1] - mv) ** 2 for p in pairs)
    a = mx - b * mv
    resid = max(abs(a + b * v - x) for x, v in pairs)
    return (lambda x: (x - a) / b), sorted(p[1] for p in pairs), resid


main_ticks = [d for d in xt if d['rect'].x0 < 400]
inset_ticks = [d for d in xt if d['rect'].x0 > 400]
val, main_vals, r1 = scale(main_ticks)
val_c, inset_vals, r2 = scale(inset_ticks)
check('x tick marks: main -2..6 and inset 0.8..1.6, linear to < 1e-3 pt', main_vals == [-2, 0, 2, 4, 6] and
      inset_vals == [0.8, 1.0, 1.2, 1.4, 1.6] and r1 < 1e-3 and r2 < 1e-3, (main_vals, inset_vals, r1, r2))

# ---- rows: label spans (claim id + suffix) with their vertical centre
rows = {}
for s in spans:
    m = re.match(r'^(PERS|SFE|SFI|SNO)-\d{3}', s['text'])
    if m and s['bbox'][0] < 100:
        rows[m.group(0)] = ((s['bbox'][1] + s['bbox'][3]) / 2, s['text'])
check('27 row labels in the hero panel', len(rows) == 27, len(rows))
lab_ok = {c: t for c, (y, t) in rows.items() if t != cat[c]['row_label']}
check('row labels = Source Data row_label', not lab_ok, lab_ok)


def near(y, kind, x_max=380):
    return [d for d in dr if kind(d) and abs((d['rect'].y0 + d['rect'].y1) / 2 - y) < 1.2 and d['rect'].x1 < x_max]


filled = lambda col: (lambda d: d['type'] == 'f' and d['items'][0][0] == 'c' and close(d['fill'], col))
opened = lambda d: d['type'] == 'fs' and d['items'][0][0] == 'c' and close(d['fill'], (1, 1, 1))
line = lambda w, col, op: (lambda d: d['type'] == 's' and len(d['items']) == 1 and d['items'][0][0] == 'l' and
                           abs(d['width'] - w) < 1e-3 and close(d['color'], col) and
                           abs((d.get('stroke_opacity') or 1) - op) < 0.01)
worst, bad, n_checked = 0.0, [], 0
for c, (y, _) in rows.items():
    col = tuple(VC[cat[c]['verdict']])
    fc = near(y, filled(col))
    oc = near(y, opened)
    ci = near(y, line(0.9, col, 1.0))
    sh = near(y, line(1.7, col, 0.45))
    if not (len(fc) == len(oc) == len(ci) == len(sh) == 1):
        bad.append((c, 'marks', len(fc), len(oc), len(ci), len(sh)))
        continue
    got = {'matched_log2_or': val((fc[0]['rect'].x0 + fc[0]['rect'].x1) / 2),
           'baseline_log2_or': val((oc[0]['rect'].x0 + oc[0]['rect'].x1) / 2),
           'matched_ci_low': val(min(ci[0]['items'][0][1].x, ci[0]['items'][0][2].x)),
           'matched_ci_high': val(max(ci[0]['items'][0][1].x, ci[0]['items'][0][2].x))}
    ends = sorted(val(p.x) for p in sh[0]['items'][0][1:3])
    got_shift = ends
    want_shift = sorted([num[c]['baseline_log2_or'], num[c]['matched_log2_or']])
    for k, v in got.items():
        e = abs(v - num[c][k])
        worst = max(worst, e)
        n_checked += 1
        if e > TOL:
            bad.append((c, k, v, num[c][k]))
    for a, b in zip(got_shift, want_shift):
        worst = max(worst, abs(a - b))
        n_checked += 1
        if abs(a - b) > TOL:
            bad.append((c, 'shift', a, b))
# SFE-006 transfer diamond (open, grey) and its interval, just below the row
y6 = rows['SFE-006'][0]
dia = [d for d in dr if d['type'] == 'fs' and len(d['items']) >= 3 and all(i[0] == 'l' for i in d['items'])
       and d['rect'].x1 < 380 and 0 < (d['rect'].y0 + d['rect'].y1) / 2 - y6 < 6]
tline = [d for d in dr if d['type'] == 's' and abs(d['width'] - 0.6) < 1e-3 and len(d['items']) == 1 and
         d['items'][0][0] == 'l' and d['rect'].x1 < 380 and 0 < d['rect'].y0 - y6 < 6]
if len(dia) == 1 and len(tline) == 1:
    got = {'transfer_matched_log2_or': val((dia[0]['rect'].x0 + dia[0]['rect'].x1) / 2),
           'transfer_matched_ci_low': val(min(tline[0]['items'][0][1].x, tline[0]['items'][0][2].x)),
           'transfer_matched_ci_high': val(max(tline[0]['items'][0][1].x, tline[0]['items'][0][2].x))}
    for k, v in got.items():
        e = abs(v - num['SFE-006'][k])
        worst = max(worst, e)
        n_checked += 1
        if e > TOL:
            bad.append(('SFE-006', k, v, num['SFE-006'][k]))
else:
    bad.append(('SFE-006 transfer marks', len(dia), len(tline)))
# PERS-010 inset
cc = tuple(VC[cat['PERS-010']['verdict']])
ins = lambda kind: [d for d in dr if kind(d) and d['rect'].x0 > 400 and d['rect'].y0 > H / 2]
fc, oc, ci = ins(filled(cc)), ins(opened), ins(line(0.9, cc, 1.0))
if len(fc) == len(oc) == len(ci) == 1:
    got = {'matched_log2_or': val_c((fc[0]['rect'].x0 + fc[0]['rect'].x1) / 2),
           'baseline_log2_or': val_c((oc[0]['rect'].x0 + oc[0]['rect'].x1) / 2),
           'matched_ci_low': val_c(min(ci[0]['items'][0][1].x, ci[0]['items'][0][2].x)),
           'matched_ci_high': val_c(max(ci[0]['items'][0][1].x, ci[0]['items'][0][2].x))}
    for k, v in got.items():
        e = abs(v - num['PERS-010'][k])
        worst = max(worst, e)
        n_checked += 1
        if e > TOL:
            bad.append(('PERS-010', k, v, num['PERS-010'][k]))
else:
    bad.append(('PERS-010 marks', len(fc), len(oc), len(ci)))
check('drawn positions = Source Data (matched, baseline, interval ends, shift ends; SFE-006 transfer; PERS-010)',
      not bad and n_checked == 27 * 6 + 3 + 4, '%d values, max |error| %.2e log2 units %s' % (n_checked, worst,
                                                                                         bad[:4] if bad else ''))

# ---- strip cells: colour per row = released colour of the final verdict
strip = [d for d in dr if d['type'] == 'f' and d['items'][0][0] == 're' and 375 < d['rect'].x0 < 392]
sbad = []
for c, (y, _) in rows.items():
    cells = [d for d in strip if abs((d['rect'].y0 + d['rect'].y1) / 2 - y) < 1.2]
    if len(cells) != 1 or not close(cells[0]['fill'], tuple(VC[cat[c]['verdict']])):
        sbad.append((c, [d['fill'] for d in cells]))
check('verdict strip: one cell per row, colour = released colour of the final verdict', not sbad and len(strip) == 27,
      sbad or '%d cells' % len(strip))

# ---- order: top to bottom = matched estimate descending within blocks
order = sorted(rows, key=lambda c: rows[c][0])
site = [c for c in order if panel[c] == 'a']
prot = [c for c in order if panel[c] == 'b']
check('row order: 18 site-level above 9 protein-level, each ranked by matched estimate',
      order == site + prot and len(site) == 18 and len(prot) == 9 and
      all(num[a]['matched_log2_or'] >= num[b]['matched_log2_or'] for blk in (site, prot)
          for a, b in zip(blk, blk[1:])), [c for c in order][:3])

# ---- the four verdict changes against the released figure are drawn as specified
spec = {'PERS-002': ('undecidable', 'PERS-002  2B  n'), 'SNO-012': ('undecidable', 'SNO-012  2E  S'),
        'SFE-006': ('attenuated', 'SFE-006  2C  t'), 'SNO-004': ('null_broken_by_control', 'SNO-004  2D  n'),
        'SNO-016': ('baseline_contradicts_claim', 'SNO-016  2B  T')}
check('PERS-002 undecidable (n), SNO-012 undecidable (S), SFE-006 attenuated own data (2C t), SNO-004, SNO-016',
      all(cat[c]['verdict'] == v and rows[c][1] == t for c, (v, t) in spec.items()),
      {c: (cat[c]['verdict'], rows[c][1]) for c in spec})
check('SFE-006 own-data values drawn: baseline 1.3840 -> matched 0.5950; transfer matched 0.1761',
      ('%.4f' % num['SFE-006']['baseline_log2_or'], '%.4f' % num['SFE-006']['matched_log2_or'],
       '%.4f' % num['SFE-006']['transfer_matched_log2_or']) == ('1.3840', '0.5950', '0.1761'))

# ---- page and type
box = doc.xref_get_key(page.xref, 'MediaBox')[1]
check('page 518.4 pt wide (7.2 in), <= 432 pt high', box.startswith('[0 0 518.4 ') and page.rect.height <= 432,
      box)
fonts = page.get_fonts()
check('fonts: Arial, embedded TrueType (fonttype 42)', all(f[1] == 'ttf' and 'Arial' in f[3] for f in fonts),
      sorted({f[3] for f in fonts}))
small = [(s['text'], round(s['size'], 2)) for s in spans if s['size'] < 5.7]
check('all text >= 5.7 pt except the two log2 subscripts', all(t == '2' for t, _ in small) and len(small) == 2,
      small)
print('ALL PASS' if all(results) else 'SOME CHECKS FAILED')
sys.exit(0 if all(results) else 1)
