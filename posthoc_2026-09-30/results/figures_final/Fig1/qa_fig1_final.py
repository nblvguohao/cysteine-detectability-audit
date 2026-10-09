"""Checks for the final Figure 1 (run after fig1_final.py). Reads only; writes nothing.

 1  page is 518.4 x 311.81 pt (7.2 in wide, <= 6 in high); every font is embedded Arial (Type0/TrueType, as
    pdf.fonttype 42 writes it); every text >= 5.7 pt; panel d uses 7-pt text, an 8-pt bold letter and 6-pt notes
 2  panels a-c are unchanged: text spans (string, origin, size, font, colour) identical to the released PDF, and the
    600-dpi rendering pixel-identical outside the panel d region
 3  nothing of the superseded panel d survives (no 'Vanishes' row, no 'Claims with a measured baseline', none of the
    superseded verdict classes); the panel d letter and title sit on the baseline of panels b and c
 4  every number printed in panel d, and every drawn segment (read back from the vector rectangles on the tick scale),
    equals a plotted Source Data row; bar fills (solid) and hatch lines equal the released Fig. 6 verdict colours
 5  Source Data: panels b and c rows identical to the released file; panel d values equal final_tally_counts.json
 6  recount (the only counting in this folder): final_tally.csv rows by level x final_verdict reproduce the stored
    counts drawn in panel d, and the stored verdicts equal Supplemental Data 2 for all 28 claims
 7  no two text boxes of panel d overlap, and no text box overlaps a bar
"""
import ast
import csv
import json
import pathlib
import re
import sys
from collections import Counter

import fitz

sys.stdout.reconfigure(encoding='utf-8')
HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
FINAL_CSV, FINAL_JSON = RESULTS / 'final' / 'final_tally.csv', RESULTS / 'final' / 'final_tally_counts.json'
MCP = pathlib.Path('/path/to/local/巯基化/MCP')
OLD_PDF, OLD_CSV = MCP / 'figures' / 'Fig1_overview.pdf', MCP / 'source_data' / 'Source_Data_Fig1_overview.csv'
SD2 = MCP / 'supplemental' / 'Supplemental_Data_2_claim_verdict_tally.csv'
FIG6_PY = MCP / '_revision' / 'fig6_final.py'
NEW_PDF, NEW_CSV, NEW_PNG = HERE / 'Fig1_overview.pdf', HERE / 'Source_Data_Fig1_overview.csv', HERE / 'Fig1_overview.png'
REGION = fitz.Rect(338.6, 178.0, 518.74, 311.81)
results = []


def check(name, ok, detail=''):
    results.append(bool(ok))
    print(('PASS ' if ok else 'FAIL ') + name + (': ' + str(detail) if detail != '' else ''))


def touches(r, R=REGION):
    """Closed-interval overlap; unlike fitz.Rect.intersects it also counts zero-width/-height rects (lines, ticks)."""
    return r.x1 >= R.x0 and r.x0 <= R.x1 and r.y1 >= R.y0 and r.y0 <= R.y1


def spans(page):
    return [s for b in page.get_text('dict')['blocks'] for ln in b.get('lines', []) for s in ln['spans']]


def released_palette():
    tree = ast.parse(FIG6_PY.read_text(encoding='utf-8'))
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(getattr(t, 'id', None) == 'VC' for t in node.targets):
            return {k: tuple(v) for k, v in ast.literal_eval(node.value).items()}


old_doc, new_doc = fitz.open(OLD_PDF), fitz.open(NEW_PDF)
old, new = old_doc[0], new_doc[0]
fc = json.loads(FINAL_JSON.read_text(encoding='utf-8'))
VC = released_palette()

# 1 ---------------------------------------------------------------------------------------------------------------
box = new_doc.xref_get_key(new.xref, 'MediaBox')[1]
check('page 518.4 pt (7.2 in) wide, 311.81 pt (<= 6 in) high', abs(new.rect.width - 518.4) < 0.001 and
      abs(new.rect.height - 311.81) < 0.01 and new.rect.height <= 6 * 72, 'MediaBox ' + box)
fonts = {(f[3], f[2], f[1]) for f in new.get_fonts(full=True)}
check('fonts: embedded Arial only (Type0 / TrueType)',
      all('Arial' in name and ftype == 'Type0' and ext == 'ttf' for name, ftype, ext in fonts), sorted(fonts))
sp_new = spans(new)
check('all text >= 5.7 pt', min(s['size'] for s in sp_new) >= 5.7, sorted({round(s['size'], 2) for s in sp_new}))
in_d = [s for s in sp_new if touches(fitz.Rect(s['bbox'])) and s['text'].strip()]
sizes_d = Counter(round(s['size'], 2) for s in in_d)
check('panel d text sizes {6 notes, 7 text, 8 letter}', set(sizes_d) == {6.0, 7.0, 8.0}, dict(sizes_d))
letter = [s for s in in_d if s['size'] == 8.0]
check('panel d letter is the only 8-pt text and is bold', len(letter) == 1 and letter[0]['text'] == 'd'
      and 'Bold' in letter[0]['font'], [(s['text'], s['font']) for s in letter])

# 2 ---------------------------------------------------------------------------------------------------------------
key = lambda s: (s['text'], round(s['origin'][0], 3), round(s['origin'][1], 3), round(s['size'], 3), s['font'],
                 s['color'])
outside = lambda sp: [key(s) for s in sp if not touches(fitz.Rect(s['bbox']))]
check('panels a-c: text spans identical to the released PDF', outside(spans(old)) == outside(sp_new),
      '%d spans' % len(outside(sp_new)))
DPI = 600
pa, pb = old.get_pixmap(dpi=DPI, alpha=False), new.get_pixmap(dpi=DPI, alpha=False)
s = DPI / 72
x_cut, y_cut = int(REGION.x0 * s) - 2, int(REGION.y0 * s) - 2               # 2-px margin inside the boundary
w_common = min(pa.width, pb.width) - 2


def rows_equal(x0, x1, y0, y1):
    n = pa.n
    for y in range(y0, y1):
        a = pa.samples_mv[(y * pa.width + x0) * n:(y * pa.width + x1) * n]
        b = pb.samples_mv[(y * pb.width + x0) * n:(y * pb.width + x1) * n]
        if a != b:
            return False, y
    return True, None


ok1, bad1 = rows_equal(0, x_cut, 0, pa.height)                              # everything left of panel d
ok2, bad2 = rows_equal(0, w_common, 0, y_cut)                                # everything above panel d
check('panels a-c: 600-dpi pixels identical outside the panel d region', ok1 and ok2,
      'left block %dx%d px, top block %dx%d px%s' % (x_cut, pa.height, w_common, y_cut,
                                                   '' if ok1 and ok2 else ' first bad row %s/%s' % (bad1, bad2)))

# 3 ---------------------------------------------------------------------------------------------------------------
all_text = ' '.join(s['text'] for s in sp_new)
gone = ('Vanishes', 'Claims with a measured baseline', 'Null contradicted', 'reading-dependent')
check('superseded panel d text absent', not any(t in all_text for t in gone), [t for t in gone if t in all_text])
b_letter = [s for s in sp_new if s['text'] == 'b' and s['size'] == 8.0]
d_title = [s for s in in_d if s['text'].startswith('Re-tested claims')]
old_d = [s for s in spans(old) if s['text'] == 'd' and s['size'] == 8.0]
one = len(letter) == 1 and len(b_letter) == 1 and len(d_title) == 1 and len(old_d) == 1   # else FAIL, not crash
check('letter and title on the baseline of panels b and c',
      one and abs(letter[0]['origin'][1] - b_letter[0]['origin'][1]) < 0.01
      and abs(d_title[0]['origin'][1] - b_letter[0]['origin'][1]) < 0.01 and d_title[0]['size'] == 7.0,
      (round(letter[0]['origin'][1], 3), round(b_letter[0]['origin'][1], 3)) if one else 'letter/title not unique')
check('letter d at the released x position', one and abs(letter[0]['origin'][0] - old_d[0]['origin'][0]) < 0.05,
      (round(letter[0]['origin'][0], 3), round(old_d[0]['origin'][0], 3)) if one else 'letter not unique')

# 4 ---------------------------------------------------------------------------------------------------------------
sd = list(csv.DictReader(NEW_CSV.open(encoding='utf-8', newline='')))
plotted = [r for r in sd if r['panel'] == 'd' and r['specification_label'].startswith('final tally')
           and 'plotted' in r['specification_label']]
plotted_values = {r['value'] for r in plotted}
ticks = {'0', '5', '10', '15'}
axis_y = max(s['origin'][1] for s in in_d if s['text'] in ticks)
printed = []
for s in in_d:
    if s['text'] in ticks and abs(s['origin'][1] - axis_y) < 0.5:     # x tick labels: the axis scale
        continue
    if 'Table 1' in s['text']:                                         # a cross-reference, not data
        continue
    printed += re.findall(r'\d+', s['text'])
check('every number printed in panel d is a plotted Source Data value (5 totals, title, 2 key counts)',
      len(printed) == 8 and all(p in plotted_values for p in printed), sorted(printed, key=int))
# read the bars back from the overlay's own content stream (the panel d Form XObject placed by show_pdf_page):
# every rectangle path with its paint (fill colour, hatch pattern), and the four tick paths for the scale.
frm = [x for x, name, *_ in new.get_xobjects() if name == 'fzFrm0']
full = int(re.search(r'/fullpage (\d+) 0 R', new_doc.xref_object(frm[0])).group(1))
ov = new_doc.xref_stream(full).decode('latin-1')
N = r'(-?\d+(?:\.\d*)?)'
Hpage = 311.81105                                   # overlay page height (= released page); PDF y is bottom-up
tick_x = sorted(float(m.group(1)) for m in re.finditer(r'%s %s m\s+\1 %s l\s+B' % (N, N, N), ov))
check('four x ticks in the overlay', len(tick_x) == 4, tick_x)
x_zero, pt_per_claim = tick_x[0], (tick_x[-1] - tick_x[0]) / 15
blocks = re.split(r'\bQ\b', ov)                      # matplotlib wraps each bar in q ... Q
rects = []
for blk in blocks:
    for m in re.finditer(r'%s %s m\s+%s \2 l\s+\3 %s l\s+\1 \4 l\s+h\s+(f|B)\b' % (N, N, N, N), blk):
        x0, y0, x1, y1, op = float(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4)), m.group(5)
        head = blk[:m.start()]
        # operands may be wrapped onto the next line by the PDF writer: separate them by any whitespace
        fills = [(f.start(), tuple(map(float, f.groups()))) for f in re.finditer(r'%s\s+%s\s+%s\s+rg\b' % (N, N, N),
                                                                                 head)]
        fills += [(f.start(), (float(f.group(1)),) * 3) for f in re.finditer(r'(?<![\d.])%s\s+g\b' % N, head)]
        stroke = re.findall(r'%s\s+%s\s+%s\s+RG' % (N, N, N), head)
        pattern = re.findall(r'/(H\d+) scn', head)
        rects.append(dict(x0=x0, x1=x1, top=Hpage - y1, bottom=Hpage - y0, op=op,
                          fill=max(fills)[1] if fills else None,
                          stroke=tuple(map(float, stroke[-1])) if stroke else None, hatched=bool(pattern)))
axes_bg = [r for r in rects if r['fill'] == (1.0, 1.0, 1.0) and r['op'] == 'f' and r['x0'] == x_zero]
check('axes background found', len(axes_bg) == 1, len(axes_bg))
top, bottom = axes_bg[0]['top'], axes_bg[0]['bottom']
labels = {s['text']: (s['bbox'][1] + s['bbox'][3]) / 2 for s in in_d if s['size'] == 7.0 and s['text'][0].isalpha()
          and s['bbox'][2] < x_zero and top < (s['bbox'][1] + s['bbox'][3]) / 2 < bottom}
CLASS = {'Survives': 'survives', 'Attenuated': 'attenuated', 'Undecidable': 'undecidable',
         'Null broken by control': 'null_broken_by_control', 'Baseline contradicts claim': 'baseline_contradicts_claim'}
check('the five category rows, in order, inside the axes',
      sorted(labels, key=labels.get) == list(CLASS), sorted(labels, key=labels.get))
bars = [r for r in rects if r is not axes_bg[0] and r['x0'] >= x_zero - 0.01 and top < r['top'] < bottom]
want = {(r['row'], r['field']): r['value'] for r in plotted}
seg_ok, detail, measured = True, {}, {}
for lab, yc in labels.items():
    segs = sorted([r for r in bars if r['top'] < yc < r['bottom']], key=lambda r: r['x0'])
    ends = [round((r['x1'] - x_zero) / pt_per_claim, 3) for r in segs]
    starts = [round((r['x0'] - x_zero) / pt_per_claim, 3) for r in segs]
    site, prot = int(want[(lab, 'site_level_claims')]), int(want[(lab, 'protein_level_claims')])
    expect = ([(0, site)] if site else []) + ([(site, site + prot)] if prot else [])
    kinds = ['hatched' if r['hatched'] else 'solid' for r in segs]
    seg_ok &= len(ends) == len(expect) and all(abs(a - e0) < 0.005 and abs(b - e1) < 0.005 for a, b, (e0, e1)
                                               in zip(starts, ends, expect)) \
        and kinds == ['solid'] * bool(site) + ['hatched'] * bool(prot)
    detail[lab] = list(zip(kinds, starts, ends))
    measured[lab] = segs
check('drawn segments (kind, start and end in claims on the tick scale) equal Source Data', seg_ok, detail)
col_ok, col_detail = True, {}
for lab, segs in measured.items():
    target = VC[CLASS[lab]]
    for r in segs:
        got = r['stroke'] if r['hatched'] else r['fill']        # hatch lines carry the verdict colour
        col_ok &= got is not None and max(abs(a - b) for a, b in zip(got, target)) < 0.003
        if r['hatched']:                                        # lighter shade: 40 % colour + 60 % white
            tinted = tuple(0.4 * c + 0.6 for c in target)
            col_ok &= max(abs(a - b) for a, b in zip(r['fill'], tinted)) < 0.003
        col_detail.setdefault(lab, []).append([round(v, 3) for v in got])
check('segment colours = released Fig. 6 verdict colours (fig6_final.py VC); protein fill = 40 % tint', col_ok,
      col_detail)

# 5 ---------------------------------------------------------------------------------------------------------------
old_rows = list(csv.DictReader(OLD_CSV.open(encoding='utf-8', newline='')))
bc_old = [(r['figure'], r['panel'], r['row'], r['field'], r['value'], r['source_table']) for r in old_rows
          if r['panel'] in 'bc']
bc_new = [(r['figure'], r['panel'], r['row'], r['field'], r['value'], r['source_table']) for r in sd
          if r['panel'] in 'bc']
check('Source Data panels b, c identical to the released file', bc_old == bc_new, '%d rows' % len(bc_new))
stored = {}
for lab, k in CLASS.items():
    stored[(lab, 'claims')] = fc['overall'][k]
    stored[(lab, 'site_level_claims')] = fc['by_level'].get('site|' + k, 0)
    stored[(lab, 'protein_level_claims')] = fc['by_level'].get('protein|' + k, 0)
stored[('All re-tested claims', 'claims')] = sum(fc['overall'].values())
stored[('All re-tested claims', 'site_level_claims')] = sum(v for k, v in fc['by_level'].items()
                                                            if k.startswith('site|'))
stored[('All re-tested claims', 'protein_level_claims')] = sum(v for k, v in fc['by_level'].items()
                                                               if k.startswith('protein|'))
mism = {k: (want.get(k), v) for k, v in stored.items() if str(v) != want.get(k)}
check('Source Data panel d plotted values = final_tally_counts.json', not mism and len(want) == len(stored),
      mism or '%d plotted values' % len(want))
header = list(sd[0].keys())
check('Source Data header', header == ['figure', 'panel', 'row', 'field', 'value', 'source_table',
                                       'specification_label'], header)
raw = NEW_CSV.read_bytes()
check('Source Data CRLF and ASCII, as the released file', b'\r\n' in raw and raw.count(b'\n') == raw.count(b'\r\n')
      and all(c < 128 for c in raw))

# 6 ---------------------------------------------------------------------------------------------------------------
tally = list(csv.DictReader(FINAL_CSV.open(encoding='utf-8', newline='')))
recount = Counter((r['level'], r['final_verdict']) for r in tally)
rc_ok = len(tally) == 28 and len({r['claim_id'] for r in tally}) == 28
rc_detail = {}
for lab, v in CLASS.items():
    got = (recount[('site', v)], recount[('protein', v)])
    exp = (int(want[(lab, 'site_level_claims')]), int(want[(lab, 'protein_level_claims')]))
    rc_ok &= got == exp
    rc_detail[lab] = got
rc_ok &= set(v for _, v in recount) == set(CLASS.values())
check('recount of final_tally.csv (level x final_verdict) = plotted', rc_ok, rc_detail)
sd2 = {r['claim_id']: r for r in csv.DictReader(SD2.open(encoding='utf-8-sig', newline=''))}
diff = {r['claim_id']: (r['stored_verdict'], sd2[r['claim_id']]['verdict'], r['final_verdict']) for r in tally
        if r['stored_verdict'] != sd2[r['claim_id']]['verdict'] or r['final_verdict'] != r['stored_verdict']}
check('stored verdicts = Supplemental Data 2 for all 28; the final differs only for SFE-006 (vanishes -> attenuated) '
      'and SNO-012 (attenuated -> undecidable)',
      diff == {'SFE-006': ('vanishes', 'vanishes', 'attenuated'),
               'SNO-012': ('attenuated', 'attenuated', 'undecidable')} and
      {k for k, r in sd2.items() if r['has_measured_baseline'] == '1'} == {r['claim_id'] for r in tally}, diff)

# 7 ---------------------------------------------------------------------------------------------------------------
boxes = [(s['text'], fitz.Rect(s['bbox'])) for s in in_d if s['text'].strip()]
clash = []
for i in range(len(boxes)):
    for j in range(i + 1, len(boxes)):
        inter = fitz.Rect(boxes[i][1]) & boxes[j][1]
        # the font boxes of adjacent lines touch by design (ascender/descender); count real overlaps only
        if not inter.is_empty and inter.width > 0.3 and inter.height > 1.0:
            clash.append((boxes[i][0], boxes[j][0], round(inter.width, 2), round(inter.height, 2)))
check('no overlapping text boxes in panel d', not clash, clash)
hit = []
for t, b in boxes:
    for r in bars:
        rb = fitz.Rect(r['x0'], r['top'], r['x1'], r['bottom'])
        inter = fitz.Rect(b) & rb
        if not inter.is_empty and inter.width > 0.3 and inter.height > 0.3:
            hit.append((t, [round(v, 1) for v in rb]))
check('no text box overlaps a bar', not hit, hit)

png = fitz.Pixmap(str(NEW_PNG))
check('PNG preview at 200 dpi', png.xres == 200 and png.yres == 200 and png.width == round(518.4 / 72 * 200),
      (png.width, png.height, png.xres))
print('ALL PASS' if all(results) else 'SOME CHECKS FAILED')
sys.exit(0 if all(results) else 1)
