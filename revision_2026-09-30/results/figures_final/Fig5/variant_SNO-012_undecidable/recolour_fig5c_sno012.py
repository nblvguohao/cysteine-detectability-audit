"""OPTIONAL VARIANT of Figure 5 (not the delivered copy): panel c's SNO-012 row recoloured from the 'attenuated' to
the 'undecidable' verdict colour, after the final tally was revised on 2026-09-30 (09:23): SNO-012 is undecidable.

Why: panels b and c of Fig. 5 colour each claim by its verdict (repo/scripts/plot_figs2to7_redesign_2026-09-23.py,
fig5(): c = S.VERDICT[r['verdict']]). The delivered Fig. 5 (../Fig5_three_axes.pdf, copied unchanged from
results/figures_rev/Fig5 as instructed) still draws SNO-012 in the attenuated light blue (0.435, 0.663, 0.839).
What: in the copied submitted page (the Form XObject that carries panels b and c), the four colour operators of the
SNO-012 row (connector line, filled marker fill and edge, open marker edge) are replaced by the undecidable grey
operator string that the same stream already uses for SNO-002, SNO-001 and SNO-009. Nothing else changes: no value,
no position, no text.
Source Data: the delivered file encodes no colour, but it lacks the two numbers printed in panel b's block headers
('n = 474', 'n = 54,855'; inherited from the submitted Source Data). This variant's Source Data is the delivered file
with those two rows added after the panel b rows (Supplemental Data 5, PERS-008, n_observations of the zf_background and
primary specifications, checked equal to n_proteins); every other byte is unchanged.
Checks: exactly four occurrences, all inside panel c's SNO-012 row; the 200-dpi rendering differs from the delivered
figure only inside that row; the recoloured marks now carry exactly the grey of the other undecidable rows; the added
Source Data rows equal Supplemental Data 5 and the rest of the file is byte-identical.

Run:  python recolour_fig5c_sno012.py      (writes Fig5_three_axes.pdf/.png and the Source Data into this folder only)
"""
import csv
import io
import hashlib
import pathlib
import re
import sys

import fitz

sys.stdout.reconfigure(encoding='utf-8')
HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE.parent / 'Fig5_three_axes.pdf'                 # the delivered, unchanged Fig. 5
OUT_PDF, OUT_PNG = HERE / 'Fig5_three_axes.pdf', HERE / 'Fig5_three_axes.png'
ATT = '.43529413 .6627451 .8392157'                        # released 'attenuated' colour as written in the stream
UND = '.6039216 .6313726 .65882357'                         # released 'undecidable' colour, same stream

doc = fitz.open(SRC)
page = doc[0]
hits = []
for xref in range(1, doc.xref_length()):
    try:
        raw = doc.xref_stream(xref)
    except Exception:
        continue
    if raw and ATT.encode() in raw:
        hits.append(xref)
assert len(hits) == 1, hits
xref = hits[0]
s = doc.xref_stream(xref).decode('latin-1')
occ = [m.start() for m in re.finditer(re.escape(ATT), s)]
assert len(occ) == 4, occ
# every occurrence is followed, before the next 'Q', by path coordinates of panel c's SNO-012 row: in this stream's
# user space (bottom-up) that row sits at y 101.2-105.4 pt, x 218-250 pt (0.4598 -> 0.5636 of the fraction axis)
for o in occ:
    seg = s[o:s.index('Q', o)]
    ys = [float(v) for v in re.findall(r'(-?\d+\.\d+|\d+) (-?\d+\.\d+|\d+) [mlc]\b', seg) for v in v[1:2]]
    xs = [float(v) for v in re.findall(r'(-?\d+\.\d+|\d+) (-?\d+\.\d+|\d+) [mlc]\b', seg) for v in v[0:1]]
    assert ys and all(100.5 < y < 106.0 for y in ys) and all(215 < x < 252 for x in xs), (o, ys, xs)
assert s.count(UND) > 0
doc.update_stream(xref, s.replace(ATT, UND).encode('latin-1'))
meta = dict(doc.metadata)
meta['subject'] = (meta.get('subject') or '') + ' | VARIANT: panel c SNO-012 recoloured attenuated -> undecidable ' \
                                               '(final tally revised 2026-09-30)'
doc.set_metadata(meta)
doc.save(OUT_PDF, garbage=3, deflate=True, no_new_id=True)
doc.close()

# ---- checks
a, b = fitz.open(SRC)[0], fitz.open(OUT_PDF)[0]
pa, pb = a.get_pixmap(dpi=200, alpha=False), b.get_pixmap(dpi=200, alpha=False)
assert (pa.width, pa.height) == (pb.width, pb.height)
n = pa.n
diff_rows = [y for y in range(pa.height) if pa.samples_mv[y * pa.width * n:(y + 1) * pa.width * n] !=
             pb.samples_mv[y * pb.width * n:(y + 1) * pb.width * n]]
diff_cols = set()
for y in diff_rows:
    ra, rb = pa.samples_mv[y * pa.width * n:(y + 1) * pa.width * n], pb.samples_mv[y * pb.width * n:(y + 1) * pb.width * n]
    diff_cols |= {x for x in range(pa.width) if ra[x * n:(x + 1) * n] != rb[x * n:(x + 1) * n]}
row_pt = (min(diff_rows) * 72 / 200, max(diff_rows) * 72 / 200)
col_pt = (min(diff_cols) * 72 / 200, max(diff_cols) * 72 / 200)
sno012 = [s_ for bl in b.get_text('dict')['blocks'] for ln in bl.get('lines', []) for s_ in ln['spans']
          if s_['text'] == 'SNO-012']
yc = (sno012[0]['bbox'][1] + sno012[0]['bbox'][3]) / 2
assert abs((row_pt[0] + row_pt[1]) / 2 - yc) < 3 and row_pt[1] - row_pt[0] < 8, (row_pt, yc)
marks = [d for d in b.get_drawings() if abs((d['rect'].y0 + d['rect'].y1) / 2 - yc) < 3 and d['rect'].x0 > 85 and
         d['rect'].x1 < 300 and d['type'] in ('s', 'fs')]
greys = {tuple(round(v, 3) for v in (d.get('color') or ())) for d in marks}
assert greys == {(0.604, 0.631, 0.659)}, greys
b.parent.load_page(0).get_pixmap(dpi=200).save(OUT_PNG)
print('recoloured 4 operators in xref %d; 200-dpi difference confined to x %.1f-%.1f pt, y %.1f-%.1f pt (SNO-012 row '
      'centre %.1f pt); SNO-012 marks now %s' % (xref, col_pt[0], col_pt[1], row_pt[0], row_pt[1], yc, greys))
for p in (OUT_PDF, OUT_PNG):
    print(hashlib.sha256(p.read_bytes()).hexdigest(), p.name)

# ---- Source Data: the delivered file plus the two numbers printed in panel b's block headers
SD_SRC, SD_OUT = HERE.parent / 'Source_Data_Fig5_three_axes.csv', HERE / 'Source_Data_Fig5_three_axes.csv'
SD5 = pathlib.Path('C:/Users/admin/Desktop/小论文/巯基化/MCP/supplemental/Supplemental_Data_5_retest_round_d.csv')
raw = SD_SRC.read_bytes()
assert raw.count(b'\r\n') == raw.count(b'\n')                      # CRLF, as the submitted file
lines = raw.decode('ascii').split('\r\n')
sd5 = {r['specification']: r for r in csv.DictReader(SD5.open(encoding='utf-8-sig', newline=''))
       if r['claim_id'] == 'PERS-008'}
printed = {s_['text'] for bl in fitz.open(SRC)[0].get_text('dict')['blocks'] for ln in bl.get('lines', [])
           for s_ in ln['spans']}
add = []
for spec, header in (('zf_background', "Authors\u2019 eligible list (n = 474)"),
                     ('primary', 'Whole proteome (n = 54,855)')):
    n = sd5[spec]['n_observations']
    assert n == sd5[spec]['n_proteins'] and header in printed and '{:,}'.format(int(n)) in header, (spec, n)
    label = {'zf_background': "zinc-finger-eligible background (the manuscript's primary verdict)",
             'primary': "whole-proteome background (the manuscript's secondary arm)"}[spec]
    buf = io.StringIO()
    csv.writer(buf, lineterminator='').writerow(['Fig5', 'b', spec, 'n_observations', n,
                                                 'Supplemental_Data_5_retest_round_d.csv',
                                                 label + '; printed in the block header (n = %s)' % '{:,}'.format(int(n))])
    add.append(buf.getvalue())
last_b = max(i for i, ln in enumerate(lines) if ln.startswith('Fig5,b,'))
out = lines[:last_b + 1] + add + lines[last_b + 1:]
SD_OUT.write_bytes('\r\n'.join(out).encode('ascii'))
new = SD_OUT.read_bytes().decode('ascii').split('\r\n')
assert [ln for ln in new if ln not in add] == lines and len(new) == len(lines) + 2
print('Source Data: %d rows + 2 added (panel b n = %s, %s)' % (len(lines) - 2, sd5['zf_background']['n_observations'],
                                                                 sd5['primary']['n_observations']))
print(hashlib.sha256(SD_OUT.read_bytes()).hexdigest(), SD_OUT.name)
