"""Independent verification of the final Figure 3 PDF (run after fig3_final.py). Reads only; writes nothing.

Nothing is taken from the plotting script. Every drawn mark is measured in the vector content of
Fig3_search_space.pdf, mapped to data units through the drawn tick marks and their labels, and compared with
  (i)  the result files themselves (results/H_artifact5_docs/*, _cys_repo_work/repo/results/*), and
  (ii) the Source Data written next to the figure.
Checks: page and type; panel a crossings, strip bars and curve; panel b deposit PSMs, sulfide PSMs, expectation
markers and their fill (|z| < 2), readings; panel c bars and +32 counts; panel d tie shares and counts; printed labels;
text overlaps (text/text and text/marks); every Source Data value equals its named source cell.
"""
import csv
import json
import math
import pathlib
import re
import sys
from collections import Counter

import fitz

sys.stdout.reconfigure(encoding='utf-8')
HERE = pathlib.Path(__file__).resolve().parent
WORK = pathlib.Path('C:/Users/admin/Desktop/小论文/_cys_repo_work')
H5 = WORK / 'public/revision_2026-09-30/results/H_artifact5_docs'
REPO = WORK / 'repo/results'
PDF, SD = HERE / 'Fig3_search_space.pdf', HERE / 'Source_Data_Fig3_search_space.csv'
AXIS = (0.29, 0.31, 0.333)
results = []


def check(name, ok, detail=''):
    results.append(bool(ok))
    print(('PASS ' if ok else 'FAIL ') + name + (': ' + str(detail) if detail != '' else ''))


def rd(p):
    with open(p, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh))


def close(a, b, t=0.003):
    return a is not None and b is not None and len(a) == len(b) and max(abs(x - y) for x, y in zip(a, b)) < t


def ctr(r):
    return ((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2)


dec, sens, ties = rd(H5 / 'precursor_accuracy_decidability.csv'), rd(H5 / 'deposit_sulfide_psms_baseline_sensitivity.csv'), \
    rd(H5 / 'single_cys_tie_rate_when_both_scored.csv')
psms, summ, reass = rd(REPO / 'search_space_mass_accuracy.csv'), rd(REPO / 'pxd015307_research_summary.csv'), \
    rd(H5 / 'deposit_sulfide_psms_reassessed.csv')
frag = rd(H5 / 'fragment_bin_arithmetic.csv')
js = json.loads((H5 / 'deposit_plus32_summary.json').read_text(encoding='utf-8'))
WIN = 'window: alternative outside the precursor window (M < M*)'
ACC = 'accuracy: separation >= 4 SD (2 SD each side of the midpoint)'

doc = fitz.open(PDF)
page = doc[0]
spans = [s for b in page.get_text('dict')['blocks'] for ln in b.get('lines', []) for s in ln['spans']
         if s['text'].strip()]
dr = page.get_drawings()

# ---------------------------------------------------------------------------------------------------- 1 page, type
box = doc.xref_get_key(page.xref, 'MediaBox')[1]
check('page exactly 518.4 pt (7.2 in) wide and <= 432 pt (6 in) high', box.startswith('[0 0 518.4 ') and
      page.rect.height <= 432, '%s = %.2f x %.2f in' % (box, page.rect.width / 72, page.rect.height / 72))
fonts = page.get_fonts()
check('fonts: Arial only, embedded TrueType (pdf.fonttype 42)', fonts and all(f[1] == 'ttf' and 'Arial' in f[3]
                                                                            for f in fonts),
      sorted({(f[3], f[2]) for f in fonts}))
sizes = Counter(round(s['size'], 2) for s in spans)
check('all text >= 5.7 pt (the figure has no sub- or superscripts)', min(sizes) >= 5.7, dict(sorted(sizes.items())))
letters = [s for s in spans if s['text'] in 'abcd' and s['size'] == 8.0]
titles = [s for s in spans for L in letters if abs(s['origin'][1] - L['origin'][1]) < 0.01 and s is not L
          and 0 < s['origin'][0] - L['origin'][0] < 12]
check('panel letters 8 pt bold, titles 7 pt regular (house style of Figs 1-7)',
      len(letters) == 4 and all('Bold' in s['font'] for s in letters) and len(titles) == 4 and
      all(s['size'] == 7.0 and 'Bold' not in s['font'] for s in titles), [(s['text'], s['size']) for s in titles])


# ---------------------------------------------------------------------------------------------------- calibration
def num(t):
    t = t.replace('\u2212', '-').replace(',', '')
    return float(t) if re.fullmatch(r'-?\d+(\.\d+)?', t) else None


ticks = [d for d in dr if len(d['items']) == 1 and d['items'][0][0] == 'l' and abs((d.get('width') or 0) - 0.8) < 1e-3
         and close(d.get('color'), AXIS, 0.01)]
vt = [d for d in ticks if d['rect'].width < 0.01 and 2 < d['rect'].height < 3]
ht = [d for d in ticks if d['rect'].height < 0.01 and 2 < d['rect'].width < 3]
nums = [(s, num(s['text'])) for s in spans if num(s['text']) is not None]


def fit(pairs):
    n = len(pairs)
    mx, mv = sum(p[0] for p in pairs) / n, sum(p[1] for p in pairs) / n
    b = sum((p[1] - mv) * (p[0] - mx) for p in pairs) / sum((p[1] - mv) ** 2 for p in pairs)
    a = mx - b * mv
    return (lambda c: (c - a) / b), max(abs(a + b * v - c) for c, v in pairs)


def clusters(ds, along, key):
    """group tick marks that share a baseline ('key') and split them where the gap along the axis exceeds 60 pt
    (panels a and b share one x baseline, as do panels c and d)."""
    out = {}
    for k in sorted({round(key(d), 1) for d in ds}):
        row = sorted([d for d in ds if round(key(d), 1) == k], key=along)
        grp = [row[0]]
        for d in row[1:]:
            if along(d) - along(grp[-1]) > 60:
                out[(k, round(along(grp[0]), 1))] = grp
                grp = []
            grp.append(d)
        out[(k, round(along(grp[0]), 1))] = grp
    return out


xaxes = clusters(vt, lambda d: d['rect'].x0, lambda d: d['rect'].y0)
yaxes = clusters(ht, lambda d: d['rect'].y0, lambda d: d['rect'].x0)
XMAP, YMAP = {}, {}
for y0, ds in xaxes.items():
    pairs = []
    for d in ds:
        x = d['rect'].x0
        lab = [v for s, v in nums if abs((s['bbox'][0] + s['bbox'][2]) / 2 - x) < 1.5 and 0 < s['bbox'][1] - d['rect'].y1 < 8]
        if len(lab) == 1:
            pairs.append((x, lab[0]))
    if len(pairs) >= 3:
        XMAP[y0] = (fit(pairs), sorted(v for _, v in pairs), (min(d['rect'].x0 for d in ds), max(d['rect'].x0 for d in ds)))
for x0, ds in yaxes.items():
    pairs = []
    for d in ds:
        y = d['rect'].y0
        lab = [v for s, v in nums if abs((s['bbox'][1] + s['bbox'][3]) / 2 - y) < 2.0 and 0 < d['rect'].x0 - s['bbox'][2] < 8]
        if len(lab) == 1:
            pairs.append((y, lab[0]))
    if len(pairs) >= 3:
        YMAP[x0] = (fit(pairs), sorted(v for _, v in pairs), (min(d['rect'].y0 for d in ds), max(d['rect'].y0 for d in ds)))
xa = {tuple(v[1]): (k, v) for k, v in XMAP.items()}
ya = {tuple(v[1]): (k, v) for k, v in YMAP.items()}
want_x = {(0.0, 2000.0, 4000.0, 6000.0, 8000.0, 10000.0, 12000.0): 'a strip (shared by the curve)',
          (-5.0, 0.0, 5.0, 10.0): 'b', (0.0, 100.0, 200.0, 300.0): 'c'}
want_y = {(0.0, 5.0, 10.0, 15.0): 'a curve', (0.0, 0.25, 0.5, 0.75, 1.0): 'd'}
check('axes calibrated from drawn ticks and labels (linear to < 0.001 pt)',
      set(want_x) <= set(xa) and set(want_y) <= set(ya) and
      all(xa[k][1][0][1] < 1e-3 for k in want_x) and all(ya[k][1][0][1] < 1e-3 for k in want_y),
      {want_x.get(k, want_y.get(k)): round(v[1][0][1], 6) for k, v in {**xa, **ya}.items() if k in want_x or k in want_y})
xA, xB, xC = (xa[k][1][0][0] for k in want_x)
yA, yD = (ya[k][1][0][0] for k in want_y)
yA_axis_x = ya[(0.0, 5.0, 10.0, 15.0)][0]

errs = {}


def rec(panel, got, want, tol):
    e = abs(got - want)
    errs.setdefault(panel, []).append(e / tol)
    return e <= tol


# ---------------------------------------------------------------------------------------------------- 2 panel a
win = {r['parameter']: r for r in dec if r['criterion'] == WIN}
K = 10.0 * float(win['tolerance 10.0 ppm']['mass_limit_da'])
dark = (0.122, 0.137, 0.149)
dots = [d for d in dr if d['type'] == 'f' and d['items'][0][0] == 'c' and close(d.get('fill'), dark) and
        d['rect'].x1 < 250 and d['rect'].y1 < 135]
got = sorted((xA(ctr(d['rect'])[0]), yA(ctr(d['rect'])[1])) for d in dots)
want = sorted((float(win['tolerance %s ppm' % t]['mass_limit_da']), float(t)) for t in ('10.0', '4.5', '2.0'))
check('a: the three window crossings (M* in Da, tolerance in ppm) = precursor_accuracy_decidability.csv window rows',
      len(got) == 3 and all(rec('a', g[0], w[0], 0.5) and rec('a', g[1], w[1], 0.001) for g, w in zip(got, want)),
      [(round(g[0], 2), round(g[1], 4)) for g in got])
grey = (0.604, 0.631, 0.659)
sbars = sorted([d for d in dr if d['type'] == 'f' and d['items'][0][0] == 're' and close(d.get('fill'), grey) and
                d['rect'].x1 < 250 and 135 < d['rect'].y0 < 200], key=lambda d: d['rect'].y0)
acc = {}
for r in dec:
    m = re.match(r'PXD015307 re-search (\S+) \(', r['context'])
    if r['criterion'] == ACC and m:
        acc[m.group(1)] = float(r['mass_limit_da'])
arms = [r['arm'] for r in summ]
got = [xA(d['rect'].x1) for d in sbars]
check('a: strip bar ends (2-SD mass range per file, Da) = accuracy rows, in the order of panel c',
      len(got) == 5 and all(rec('a', g, acc[a], 0.5) for g, a in zip(got, arms)), [round(g, 2) for g in got])
curve = [d for d in dr if d['type'] == 's' and abs((d.get('width') or 0) - 1.4) < 1e-3]
pts = [p for d in curve for it in d['items'] if it[0] == 'l' for p in (it[1], it[2])]
dev = max(abs(yA(p.y) - K / xA(p.x)) for p in pts)
check('a: every vertex of the drawn curve lies on K / M ppm, K = 10 x M*(10 ppm) (%d vertices)' % len(pts),
      len(curve) == 1 and dev < 0.002, 'max |deviation| %.2e ppm' % dev)
txt = ' | '.join(s['text'] for s in spans)
lab_ok = True
for t in ('10.0', '4.5', '2.0'):
    r = win['tolerance %s ppm' % t]
    fr = float(r['fraction_theoretical_cys_peptides_at_or_below_plus32_only'])
    pct = '>99.9%' if 0.999 < fr < 1.0 else '%.0f%%' % (100 * fr)
    s = '±%s ppm: %s Da (%s)' % (t.rstrip('0').rstrip('.'), '{:,.0f}'.format(float(r['mass_limit_da'])), pct)
    lab_ok &= s in txt
lab_ok &= all('{:,.0f}'.format(acc[a]) in [s['text'] for s in spans] for a in arms)
lab_ok &= ('%s Da in ppm of M' % [r for r in frag if r['fragment_charge'] == '1'][0]['shift_mz']) in txt
check('a: printed crossing labels, strip values and the separation label = the result files', lab_ok)

# ---------------------------------------------------------------------------------------------------- 3 panel b
blue, orange = (0.122, 0.416, 0.647), (0.863, 0.604, 0.169)
leg = [s for s in spans if s['text'] in ('B chain', 'A chain', 'sulfide-assigned', 'expected if sulfide',
                                          'expected if dioxidation', 'filled: PSM within 2 SD')]
in_legend = lambda d: any(abs(ctr(d['rect'])[1] - (s['bbox'][1] + s['bbox'][3]) / 2) < 3 and
                          0 < s['bbox'][0] - ctr(d['rect'])[0] < 14 for s in leg)
rowlab = {s['text']: (s['bbox'][1] + s['bbox'][3]) / 2 for s in spans if s['text'] in
          ('CTH Cys', 'NaHS', 'control', 'heat inactive CTH') and s['bbox'][0] > 250 and s['bbox'][2] < 330}
chain_of = {r['sequence']: r['chain'] for r in reass}
arm_of_label = {'CTH Cys': 'CTH_Cys', 'NaHS': 'NaHS', 'control': 'control', 'heat inactive CTH': 'heat_inactive_CTH'}
circ = [d for d in dr if d['type'] == 'fs' and d['items'][0][0] == 'c' and abs(d['width'] - 0.3) < 1e-3
        and not in_legend(d)]
diam = [d for d in dr if d['type'] == 'fs' and len(d['items']) >= 4 and all(i[0] == 'l' for i in d['items'])
        and abs(d['width'] - 0.8) < 1e-3 and not in_legend(d)]
upper_y = max(rowlab.values()) + 14      # the swarm spreads up to ~10 pt; the lower section starts ~25 pt lower
ok_b, det_b = True, {}
for lab, yc in rowlab.items():
    arm = arm_of_label[lab]
    mine = [p for p in psms if p['arm'] == arm]
    # each mark of the upper section belongs to the nearest row label among the rows that have PSMs (the beeswarm
    # spreads a row by up to ~10 pt; the control row has none and carries the text 'no PSM in the deposit')
    has = {k for k in rowlab if any(p_['arm'] == arm_of_label[k] for p_ in psms)}
    nearest = lambda d: min(has, key=lambda k: abs(rowlab[k] - ctr(d['rect'])[1]))
    dots_ = [d for d in circ if ctr(d['rect'])[1] < upper_y and nearest(d) == lab]
    dias_ = [d for d in diam if ctr(d['rect'])[1] < upper_y and nearest(d) == lab]
    g = sorted((xB(ctr(d['rect'])[0]), 'A' if close(d['fill'], orange) else 'B') for d in dots_)
    w = sorted((float(p['observed_ppm']), chain_of[p['sequence']]) for p in mine if p['has_sulfide'] == 'False')
    gd = sorted((xB(ctr(d['rect'])[0]), 'A' if close(d['fill'], orange) else 'B') for d in dias_)
    wd = sorted((float(p['observed_ppm']), chain_of[p['sequence']]) for p in mine if p['has_sulfide'] == 'True')
    same = lambda u, v: len(u) == len(v) and all(rec('b', a[0], b[0], 0.001) and a[1] == b[1] for a, b in zip(u, v))
    ok_b &= same(g, w) and same(gd, wd)
    det_b[lab] = '%d dots, %d diamonds' % (len(g), len(gd))
n_ctrl = js['targeted_arms']['control']['n_psms']
check('b: every deposit PSM drawn at its observed_ppm (to 0.001 ppm) with its chain colour; sulfide PSMs as '
      'diamonds (search_space_mass_accuracy.csv; chain from deposit_sulfide_psms_reassessed.csv)',
      ok_b and sum(int(v.split()[0]) + int(v.split()[2]) for v in det_b.values()) == len(psms) == 53 and n_ctrl == 0
      and 'no PSM in the deposit' in txt, det_b)
SHOW = {('A', 'D1'), ('A', 'D6b'), ('B', 'D6a')}
rows = [r for r in sens if (r['chain'], r['baseline_definition'].split()[0]) in SHOW]
lower_labels = sorted([s for s in spans if re.match(r'^(B chain|A chain)', s['text']) and s['bbox'][2] < 330 and
                       (s['bbox'][1] + s['bbox'][3]) / 2 > upper_y], key=lambda s: s['bbox'][1])
readings = sorted([s for s in spans if s['text'] in ('neither', 'sulfide only', 'dioxidation only', 'either')],
                  key=lambda s: s['bbox'][1])
exp_marks = [d for d in dr if d['type'] == 'fs' and abs((d.get('width') or 0) - 1.0) < 1e-3 and d['rect'].x0 > 328
             and ctr(d['rect'])[1] > upper_y and d['rect'].width < 5 and not in_legend(d)]
ok_l, det_l, n_val = True, [], 0
order = []   # drawn order: for each sulfide PSM (arm, observed), B-chain baseline row then A-chain baseline row
for lab, rd_ in zip(lower_labels, readings):
    yc = (lab['bbox'][1] + lab['bbox'][3]) / 2
    cands = [r for r in rows if r['reading'] == rd_['text'] and
             lab['text'].startswith(('B chain, same file' if r['baseline_definition'].startswith('D1') else
                                     'A chain' if r['baseline_definition'].startswith('D6b') else 'B chain ('))]
    marks = [d for d in exp_marks if abs(ctr(d['rect'])[1] - yc) < 2]
    circle = [d for d in marks if d['items'][0][0] == 'c']
    square = [d for d in marks if d['items'][0][0] == 're']
    dia = [d for d in diam if abs(ctr(d['rect'])[1] - yc) < 2]
    if not (len(circle) == len(square) == len(dia) == 1):
        ok_l = False
        det_l.append((lab['text'], 'marks', len(circle), len(square), len(dia)))
        continue
    xo, xs, xd = xB(ctr(dia[0]['rect'])[0]), xB(ctr(circle[0]['rect'])[0]), xB(ctr(square[0]['rect'])[0])
    match = [r for r in cands if abs(float(r['observed_ppm']) - xo) < 0.01]
    if len(match) != 1:
        ok_l = False
        det_l.append((lab['text'], 'no unique source row', len(match)))
        continue
    r = match[0]
    zs, zo = float(r['z_if_sulfide']), float(r['z_if_dioxidation'])
    filled_s = not close(circle[0]['fill'], (1, 1, 1))
    filled_o = not close(square[0]['fill'], (1, 1, 1))
    n_base = int(re.search(r'\((\d+)\)', lab['text']).group(1))
    ok_l &= rec('b', xs, float(r['centre_ppm']), 0.002) and \
        rec('b', xd, float(r['expected_ppm_if_dioxidised_but_reported_as_sulfide']), 0.002) and \
        rec('b', xo, float(r['observed_ppm']), 0.002) and filled_s == (abs(zs) < 2) and filled_o == (abs(zo) < 2) \
        and n_base == int(r['n_baseline_psms'])
    n_val += 4
    det_l.append((lab['text'], r['arm'], round(xo, 3), r['baseline_definition'].split()[0], rd_['text']))
check('b: the six baseline rows: expectation circle = centre_ppm, square = expected_ppm_if_dioxidised, diamond = '
      'observed_ppm (to 0.002 ppm); filled exactly when |z| < 2; n and reading = deposit_sulfide_psms_baseline_'
      'sensitivity.csv', ok_l and len(det_l) == 6 and n_val == 24, det_l)

# ---------------------------------------------------------------------------------------------------- 4 panel c
cbars = sorted([d for d in dr if d['type'] == 'f' and d['items'][0][0] == 're' and close(d.get('fill'), grey) and
                d['rect'].x1 < 250 and d['rect'].y0 > 230], key=lambda d: d['rect'].y0)
got = [xC(d['rect'].x1) for d in cbars]
vals_c = [s['text'] for s in spans if s['bbox'][0] < 250 and s['bbox'][1] > 230 and s['bbox'][3] < 322 and
          re.fullmatch(r'\d+', s['text'])]                     # tick labels (below the axis at 324 pt) excluded
check('c: PSMs at 1% FDR per file (bar ends) = pxd015307_research_summary.csv; printed counts and the five +32 zeros',
      len(got) == 5 and all(rec('c', g, float(r['n_psms_at_1pct_fdr']), 0.01) for g, r in zip(got, summ)) and
      sorted(vals_c) == sorted([r['n_psms_at_1pct_fdr'] for r in summ] + [r['n_plus32_cys_psms'] for r in summ]) and
      all(r['n_plus32_cys_psms'] == '0' for r in summ), [round(g, 3) for g in got])

# ---------------------------------------------------------------------------------------------------- 5 panel d
red = (0.8, 0.353, 0.212)
dm = [d for d in dr if d['type'] == 'f' and d['rect'].x0 > 300 and d['rect'].y0 > 230 and
      (close(d.get('fill'), red) or close(d.get('fill'), dark))]
xt_d = sorted([d['rect'].x0 for d in vt if d['rect'].x0 > 300 and d['rect'].y0 > 230])
bins = ['(0.0, 1.0]', '(1.0, 10.0]', '(10.0, 100.0]', '(100.0, 1000.1]', 'all']
tie = {(r['group'].split()[0], r['evalue_bin']): r for r in ties}
ok_d, n_d = len(dm) == 10 and len(xt_d) == 5, 0
for d in dm:
    cx, cy = ctr(d['rect'])
    i = min(range(5), key=lambda k: abs(xt_d[k] - cx))
    grp = 'insulin' if close(d['fill'], red) else 'Fig-1D'
    r = tie[(grp, bins[i])]
    ok_d &= abs(xt_d[i] - cx) < 0.01 and rec('d', yD(cy), float(r['share_tied']), 0.0005)
    n_d += 1
printed = {s['text'] for s in spans if s['bbox'][0] > 300 and s['bbox'][1] > 230}
ok_d &= all('%s/%s' % (r['n_tied'], r['n_both_scored']) in printed for r in ties)
ok_d &= '%.0f%%' % (100 * float(tie[('insulin', 'all')]['share_tied'])) in printed and \
    '%.0f%%' % (100 * float(tie[('Fig-1D', 'all')]['share_tied'])) in printed
check('d: 10 tie shares (marker heights, to 0.0005) at their E-value bins, 10 tied/n counts and the two overall '
      'percentages = single_cys_tie_rate_when_both_scored.csv', ok_d and n_d == 10)

# ---------------------------------------------------------------------------------------------------- 6 overlaps
words = []
for b in page.get_text('rawdict')['blocks']:
    for ln in b.get('lines', []):
        horiz = abs(ln['dir'][0] - 1) < 1e-6
        for s in ln['spans']:
            cur = []
            for c in s['chars'] + [None]:
                if c is not None and c['c'].strip():
                    cur.append(c)
                    continue
                if cur:
                    t = ''.join(c_['c'] for c_ in cur)
                    x0, x1 = min(c_['bbox'][0] for c_ in cur), max(c_['bbox'][2] for c_ in cur)
                    if horiz:
                        oy, sz = s['origin'][1], s['size']
                        top = oy - (0.74 if any(ch in 'bdfhklt()[]{}|/%' for ch in t) else 0.716) * sz
                        bot = oy + (0.212 * sz if any(ch in 'gjpqy(),;[]{}|/_Q' for ch in t) else 0.0)
                        words.append((t, fitz.Rect(x0, top, x1, bot)))
                    else:
                        words.append((t, fitz.Rect(x0, min(c_['bbox'][1] for c_ in cur), x1,
                                                   max(c_['bbox'][3] for c_ in cur))))
                cur = []
tt = [(a[0], b[0]) for i, a in enumerate(words) for b in words[i + 1:]
      if (lambda r: not r.is_empty and r.width > 0.2 and r.height > 0.2)(fitz.Rect(a[1]) & b[1])]


def seg_hits(p, q, r, half):
    n = max(2, int(math.hypot(q.x - p.x, q.y - p.y) / 0.05))
    for k in range(n + 1):
        x, y = p.x + (q.x - p.x) * k / n, p.y + (q.y - p.y) * k / n
        if math.hypot(max(r.x0 - x, 0, x - r.x1), max(r.y0 - y, 0, y - r.y1)) < half - 0.05:
            return True
    return False


tm = []
for d in dr:
    big = d['rect'].width > 60 and d['rect'].height > 30
    for t, r in words:
        if d['type'] in ('f', 'fs') and not big:
            inter = fitz.Rect(d['rect']) & r
            if not inter.is_empty and inter.width > 0.2 and inter.height > 0.2:
                tm.append((t, 'fill', [round(v, 1) for v in d['rect']]))
                continue
        if d['type'] in ('s', 'fs') and (d.get('width') or 0) and d['rect'].intersects(r + (-2, -2, 2, 2)):
            for it in d['items']:
                if it[0] == 'l' and seg_hits(it[1], it[2], r, d['width'] / 2):
                    tm.append((t, 'line', [round(v, 1) for v in (it[1].x, it[1].y, it[2].x, it[2].y)]))
                    break
check('no text overlaps text (%d words, rotated axis labels included)' % len(words), not tt, tt[:5])
check('no text overlaps a stroked segment, marker or bar (knock-out boxes count as overlaps)', not tm, tm[:8])

# ---------------------------------------------------------------------------------------------------- 7 Source Data
sd = rd(SD)
src_ok, n_src = True, 0
for r in sd:
    p, row, f, v = r['panel'], r['row'], r['field'], r['value']
    if r['source_table'].startswith('derived'):
        if f == 'separation_ppm':
            m = float([x['value'] for x in sd if x['row'] == row and x['field'] == 'mass_da'][0])
            src_ok &= abs(float(v) - K / m) < 1e-9
        continue
    if p == 'a' and row.startswith('window'):
        t = row.split('+/-')[1].split()[0]
        want = t if f == 'tolerance_ppm' else win['tolerance %s ppm' % t][f]
    elif p == 'a' and row.startswith('2-SD'):
        want = [x for x in dec if x['criterion'] == ACC and x['context'].startswith(
            'PXD015307 re-search %s (' % row.split('|')[1])][0][f]
    elif p == 'a' and row == 'separation label':
        want = [x for x in frag if x['fragment_charge'] == '1'][0]['shift_mz']
    elif p == 'b' and row.startswith('deposit PSM '):
        want = psms[int(row.split('|')[0].split()[-1]) - 1][f]
    elif p == 'b' and row == 'deposit PSMs|control':
        want = str(n_ctrl)
    elif p == 'b' and row.startswith('sulfide PSM|'):
        parts = row.split('|')
        want = [x for x in sens if x['arm'] == parts[1] and x['chain'] == parts[2][0] and x['observed_ppm'] == parts[5]
                and x['baseline_definition'].split()[0] == parts[6].split()[1]][0][f]
    elif p == 'c':
        want = [x for x in summ if x['arm'] == row][0][f]
    elif p == 'd':
        g, b = row.split('|best E-value ')
        want = [x for x in ties if x['group'] == g and x['evalue_bin'] == b][0][f]
    else:
        want = None
    src_ok &= v == want
    n_src += 1
check('Source Data: every non-derived value is the verbatim cell of its named result file (%d values); the 186 '
      'curve points equal K / M' % n_src, src_ok and len(sd) == 523, '%d rows' % len(sd))
raw = SD.read_bytes()
check('Source Data CRLF and ASCII, as the released file', raw.count(b'\n') == raw.count(b'\r\n') and
      all(c < 128 for c in raw))

worst = {k: round(max(v), 3) for k, v in errs.items()}
print('largest drawn-value error as a fraction of its tolerance, per panel:', worst,
      '| values measured:', {k: len(v) for k, v in errs.items()})
print('ALL PASS' if all(results) else 'SOME CHECKS FAILED')
sys.exit(0 if all(results) else 1)
