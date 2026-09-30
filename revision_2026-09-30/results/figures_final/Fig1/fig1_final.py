"""Figure 1, final revision (2026-09-30): panel d redrawn with the FINAL claim tally.

What changes
  Panel d only. Panels a-c are NOT redrawn: the released page (MCP/figures/Fig1_overview.pdf, read-only) is opened,
  the old panel d is deleted from its content stream (its 19 text objects and 12 vector units, so the superseded
  tally does not survive hidden under a white box; see 'assemble' for why PyMuPDF redaction is not used), and a newly
  drawn panel d (matplotlib, house style of panels b and c) is placed at the same position with fitz show_pdf_page at
  scale 1. The page is then trimmed on the right from 518.74 pt (183 mm) to 518.4 pt (7.2 in); panels a-c keep
  their coordinates and are unchanged to the pixel (checked by qa_fig1_final.py).
  Adapted from results/figures_rev/Fig1/fig1_rev.py (same page surgery, same house style); only the tally, the
  categories, the colours and the marks differ.

Tally (final verdict policy; results/final/final_tally.csv and final_tally_counts.json)
  Stored verdicts for 26 claims; two differ from the stored (released) verdict: SFE-006 is tallied on its authors'
  own data (attenuated; its transfer verdict, vanishes, is not tallied) and SNO-012 is undecidable (stored as
  attenuated by pipeline code that lacked the baseline-covers-zero branch; its baseline interval covers zero).
  11 survive, 2 attenuated, 13 undecidable, 1 null broken by control, 1 baseline contradicts claim, 0 vanish.
  Site level (18): 4 / 2 / 10 / 1 / 1; protein level (10): 7 / 0 / 3. (final_tally files as revised 2026-09-30 09:23)

Panel d design
  one bar per verdict class (Survives, Attenuated, Undecidable, Null broken by control, Baseline contradicts claim),
  split into a site-level segment (solid verdict colour) and a protein-level segment (lighter tint of the same colour,
  hatched in the verdict colour, so that a pale survives-blue segment cannot be read as the attenuated class; the
  released text reads protein-level commonality as not estimable, Table 1, item 4, which the key repeats);
  SNO-016's bar is annotated 'transfer cohort' (Supplemental Data 4: is_transfer). No 'Vanishes' bar: the class is
  empty in the final tally (a non-plotted Source Data row records the 0). Colours are the released Fig. 6 verdict
  colours (MCP/_revision/fig6_final.py, dict VC), which the final Fig. 6 reuses.

Values
  Read, never recomputed, from final_tally_counts.json (overall.*, by_level.'site|*' / 'protein|*'); the block sizes
  18 and 10 are the sums of the stored by_level cells. The Source Data names the publishable tables instead of the
  workspace files: Supplemental Data 2 (final_verdict, added 09:40) x Supplemental Data 10 (unit); the build refuses
  unless final_tally.csv equals them claim by claim and their recount equals final_tally_counts.json. A recount from
  final_tally.csv is also done in qa_fig1_final.py.

Run:  python fig1_final.py        (python 3.10, matplotlib 3.10, PyMuPDF 1.28; writes only into this folder)
Outputs: Fig1_overview.pdf, Fig1_overview.png (200 dpi), Source_Data_Fig1_overview.csv
"""
import ast
import csv
import hashlib
import io
import json
import pathlib
import re
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
import fitz  # noqa: E402

sys.stdout.reconfigure(encoding='utf-8')

# ----------------------------------------------------------------------------------------------- paths
HERE = pathlib.Path(__file__).resolve().parent                    # .../results/figures_final/Fig1
RESULTS = HERE.parents[1]                                          # .../revision_2026-09-30/results
FINAL_CSV = RESULTS / 'final' / 'final_tally.csv'
FINAL_JSON = RESULTS / 'final' / 'final_tally_counts.json'
MCP = pathlib.Path('C:/Users/admin/Desktop/小论文/巯基化/MCP')     # read-only
SRC_PDF = MCP / 'figures' / 'Fig1_overview.pdf'
SRC_CSV = MCP / 'source_data' / 'Source_Data_Fig1_overview.csv'
FIG6_PY = MCP / '_revision' / 'fig6_final.py'                     # released Fig. 6 verdict colours (dict VC)
SD4 = MCP / 'supplemental' / 'Supplemental_Data_4_retest_round_b.csv'
SD2 = MCP / 'supplemental' / 'Supplemental_Data_2_claim_verdict_tally.csv'      # carries final_verdict since 09:40
SD10 = MCP / 'supplemental' / 'Supplemental_Data_10_claim_independence.csv'    # unit = site / protein level
OUT_PDF, OUT_PNG, OUT_CSV = (HERE / 'Fig1_overview.pdf', HERE / 'Fig1_overview.png',
                             HERE / 'Source_Data_Fig1_overview.csv')
# source_table labels name the publishable tables: the values are READ from results/final (as specified) and gated
# equal to Supplemental Data 2 (final_verdict) x Supplemental Data 10 (unit) below
T2 = 'Supplemental_Data_2_claim_verdict_tally.csv'
T10 = 'Supplemental_Data_10_claim_independence.csv'

# ----------------------------------------------------------------------------------------------- house style
# rcParams of nc_figure_style_2026-09-21.py (V1) updated by nc_figure_style_v2_2026-09-23.py (the modules that drew
# panels b-d of the released figure), plus the mathtext setting of plot_fig1_biorender_r36_2026-09-23.py.
INK, AXIS, MUTE = '#1f2326', '#4a4f55', '#6f757c'
FS, FS_PANEL, FS_NOTE = 7, 8, 6          # FS_NOTE: in-panel note and key, as the 6-pt notes of Figs 2 and 6
LW = 1.0
RC = {
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'],
    'font.size': FS, 'axes.titlesize': FS, 'axes.labelsize': FS, 'xtick.labelsize': FS, 'ytick.labelsize': FS,
    'legend.fontsize': FS, 'axes.linewidth': LW, 'lines.linewidth': LW, 'xtick.major.width': LW,
    'ytick.major.width': LW, 'xtick.minor.width': LW, 'ytick.minor.width': LW, 'patch.linewidth': LW,
    'errorbar.capsize': 2, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.titlelocation': 'left', 'axes.titlepad': 4, 'legend.frameon': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'axes.edgecolor': AXIS, 'xtick.color': AXIS, 'ytick.color': AXIS, 'text.color': INK,
    'axes.labelcolor': INK, 'xtick.labelcolor': INK, 'ytick.labelcolor': INK,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.major.pad': 2.0, 'ytick.major.pad': 2.0,
    'axes.labelpad': 3.0, 'lines.markeredgewidth': LW, 'lines.solid_capstyle': 'butt',
    'axes.unicode_minus': True, 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
    'hatch.linewidth': 0.6,             # the hatch that marks protein-level segments
}
HATCH = '////'          # 3 pt spacing along x in the PDF backend (72-pt hatch cell, density 6 per '/')
TINT = 0.40             # protein-level fill = 40 % verdict colour + 60 % white ('lighter shade')
BAR_H = 0.62            # bar height in row units, as hbars() of the released script
XMAX, XTICKS = 15, [0, 5, 10, 15]
AX_RIGHT = 510.0        # pt from the left page edge; the released panel d ended at 505.77, 518.4 is the new edge
TRIM_W = 7.2 * 72       # 518.4 pt
SHOW_NOT_ESTIMABLE = False  # key line 'not estimable (Table 1, item 4)' under 'Protein level (10):'


def tint(c, f=TINT):
    r, g, b = to_rgb(c)
    return (f * r + 1 - f, f * g + 1 - f, f * b + 1 - f)


def released_palette():
    """The dict literal VC of the released Fig. 6 script, read with ast (no code is executed)."""
    tree = ast.parse(FIG6_PY.read_text(encoding='utf-8'))
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(getattr(t, 'id', None) == 'VC' for t in node.targets):
            return {k: tuple(v) for k, v in ast.literal_eval(node.value).items()}
    raise SystemExit('REFUSE: no VC dict in %s' % FIG6_PY)


# ----------------------------------------------------------------------------------------------- stored values
fc = json.loads(FINAL_JSON.read_text(encoding='utf-8'))
overall, by_level = fc['overall'], fc['by_level']
VC = released_palette()

# (row label as drawn, verdict key in final_tally*.*, colour key in the released VC)
CATS = [('Survives', 'survives'),
        ('Attenuated', 'attenuated'),
        ('Undecidable', 'undecidable'),
        ('Null broken by control', 'null_broken_by_control'),
        ('Baseline contradicts claim', 'baseline_contradicts_claim')]
# every verdict present in the final tally must have a row (nothing silently omitted)
assert set(overall) == {k for _, k in CATS}, sorted(overall)
assert {k.split('|')[1] for k in by_level} <= set(overall) and {k.split('|')[0] for k in by_level} == {'site',
                                                                                                        'protein'}
assert all(k in VC for _, k in CATS), sorted(VC)

ROWS = []
for label, k in CATS:
    n_site, n_prot, n_tot = by_level.get('site|' + k, 0), by_level.get('protein|' + k, 0), overall[k]
    assert n_site + n_prot == n_tot, label           # stored-versus-stored consistency; nothing is counted here
    ROWS.append(dict(label=label, key=k, colour=VC[k], site=n_site, prot=n_prot, total=n_tot,
                     src_site=f'{T2} (final_verdict = {k}) x {T10} (unit = site)',
                     src_prot=f'{T2} (final_verdict = {k}) x {T10} (unit = protein)',
                     src_tot=f'{T2} (final_verdict = {k})'))
N_SITE = sum(v for kk, v in by_level.items() if kk.startswith('site|'))
N_PROT = sum(v for kk, v in by_level.items() if kk.startswith('protein|'))
N_ALL = sum(overall.values())
assert N_SITE + N_PROT == N_ALL

# gate: the tally this figure is specified to show (task text of 2026-09-30, final policy); a changed input refuses,
# it does not silently redraw a different figure
EXPECTED = {'Survives': (4, 7), 'Attenuated': (2, 0), 'Undecidable': (10, 3), 'Null broken by control': (1, 0),
            'Baseline contradicts claim': (1, 0)}
got = {r['label']: (r['site'], r['prot']) for r in ROWS}
if got != EXPECTED or (N_ALL, N_SITE, N_PROT) != (28, 18, 10) or 'vanishes' in overall:
    sys.exit(f'REFUSE: stored tally differs from the specified one: {got} {(N_ALL, N_SITE, N_PROT)}')

# per-claim rows: membership of the segments (Source Data only) and the claims that carry a mark or a note
tally = list(csv.DictReader(FINAL_CSV.open(encoding='utf-8', newline='')))
members = {}
for t in tally:
    members.setdefault((t['final_verdict'], t['level']), []).append(t['claim_id'])
bcc = members[('baseline_contradicts_claim', 'site')]
sd4 = {r['claim_id']: r for r in csv.DictReader(SD4.open(encoding='utf-8-sig', newline=''))}
assert bcc == ['SNO-016'] and sd4['SNO-016']['is_transfer'] == 'True', bcc          # 'transfer cohort' mark
TRANSFER_COHORT = sd4['SNO-016']['transfer_cohort']
# the claims whose final verdict differs from the stored (released) one: noted in the Source Data, not plotted
CHANGED = [t for t in tally if t['final_verdict'] != t['stored_verdict']]
assert {t['claim_id']: (t['stored_verdict'], t['final_verdict']) for t in CHANGED} == {
    'SFE-006': ('vanishes', 'attenuated'), 'SNO-012': ('attenuated', 'undecidable')}, CHANGED
ROW_OF = {k: lab for lab, k in CATS}
ANNOT = {'Baseline contradicts claim': 'transfer cohort'}

# gate: the publishable tables named in the Source Data hold the same final tally, claim by claim
sd2 = {r['claim_id']: r for r in csv.DictReader(SD2.open(encoding='utf-8-sig', newline=''))}
sd10 = {r['claim_id']: r for r in csv.DictReader(SD10.open(encoding='utf-8-sig', newline=''))}
claim_ids = {t['claim_id'] for t in tally}
assert claim_ids == {c for c, r in sd2.items() if r['has_measured_baseline'] == '1'} == \
    {c for c, r in sd10.items() if r['has_measured_baseline'] == '1'}, 'claim sets differ'
bad = [t['claim_id'] for t in tally if (sd2[t['claim_id']]['final_verdict'], sd2[t['claim_id']]['verdict'],
                                        sd10[t['claim_id']]['unit']) != (t['final_verdict'], t['stored_verdict'],
                                                                         t['level'])
       or (t in CHANGED and sd2[t['claim_id']]['final_verdict_note'] != t['note'])]
assert not bad, ('final_tally.csv differs from Supplemental Data 2 / 10', bad)
recount = {}
for c in claim_ids:
    kk = '%s|%s' % (sd10[c]['unit'], sd2[c]['final_verdict'])
    recount[kk] = recount.get(kk, 0) + 1
assert recount == by_level, ('Supplemental Data 2 x 10 recount differs from final_tally_counts.json', recount)

# ----------------------------------------------------------------------------------------------- released geometry
src = fitz.open(SRC_PDF)
page = src[0]
W, H = page.rect.width, page.rect.height                             # 518.74 x 311.81 pt
spans = [s for b in page.get_text('dict')['blocks'] for ln in b.get('lines', []) for s in ln['spans']]
draws = page.get_drawings()
# panel b: bottom spine (stroked horizontal line) and axes background (white fill) -> panel d uses the same y range
b_spine = [d for d in draws if d['type'] == 's' and abs(d['rect'].x0 - 88.19) < 0.05 and d['rect'].height < 0.01]
b_axes = [d for d in draws if d['type'] == 'f' and abs(d['rect'].x0 - 88.19) < 0.05 and d['fill'] == (1.0, 1.0, 1.0)]
assert len(b_spine) == 1 and len(b_axes) == 1, (b_spine, b_axes)
Y_SPINE, Y_AXTOP = b_spine[0]['rect'].y0, b_axes[0]['rect'].y0          # 279.07, 199.56
old_d = [s for s in spans if s['text'] == 'd' and s['size'] == FS_PANEL]
assert len(old_d) == 1
X_LETTER = old_d[0]['bbox'][0]                                           # 339.71: the released panel d letter
REGION = fitz.Rect(338.6, 178.0, W, H)       # released panel d; panel c ends at 338.42 (spine cap), panel a at y<170
assert all(s['bbox'][0] >= REGION.x0 for s in spans if fitz.Rect(s['bbox']).intersects(REGION))

# ----------------------------------------------------------------------------------------------- panel d
plt.rcParams.update(RC)
fig = plt.figure(figsize=(W / 72, H / 72))


def frac(x, y):
    """PDF points (origin top-left, as fitz) -> figure fraction."""
    return x / W, 1 - y / H


def place(ax, left):
    ax.set_position([left / W, 1 - Y_SPINE / H, (AX_RIGHT - left) / W, (Y_SPINE - Y_AXTOP) / H])


ax = fig.add_axes([0, 0, 1, 1])
place(ax, 440.0)
n = len(ROWS)
ys = list(range(n))[::-1]
ax.set_xlim(0, XMAX)
ax.set_ylim(-0.6, n - 0.4)
ax.set_yticks(ys)
ax.set_yticklabels([r['label'] for r in ROWS])
ax.tick_params(axis='y', length=0, pad=3)                                # clean_y()
ax.spines['left'].set_visible(False)
ax.set_xticks(XTICKS)
ax.tick_params(axis='x', length=2)
ax.set_xlabel('Claims')

# solve the axes' left edge so that the category labels start where the released panel d's labels (and letter) did
fig.canvas.draw()
r_ = fig.canvas.get_renderer()
for _ in range(3):
    x0_pt = ax.get_tightbbox(r_).x0 * 72 / fig.dpi
    left = ax.get_position().x0 * W + (X_LETTER - x0_pt)
    place(ax, left)
    fig.canvas.draw()
AX_LEFT = ax.get_position().x0 * W
PT_X = (AX_RIGHT - AX_LEFT) / XMAX                       # pt per claim

for y, r in zip(ys, ROWS):
    c = r['colour']
    if r['site']:
        ax.barh(y, r['site'], height=BAR_H, color=c, lw=0, zorder=2)
    if r['prot']:                      # protein level: lighter shade, hatched in the verdict colour, no outline
        ax.barh(y, r['prot'], left=r['site'], height=BAR_H, facecolor=tint(c), edgecolor=c, lw=0, hatch=HATCH,
                zorder=2)
    t = ax.text(r['total'] + XMAX * 0.02, y, str(r['total']), va='center', ha='left', color=MUTE)
    if r['label'] in ANNOT:
        ax.annotate(ANNOT[r['label']], xy=(1, 0), xycoords=t, xytext=(2.5, 0), textcoords='offset points',
                    ha='left', va='bottom', fontsize=FS_NOTE, color=MUTE)


# panel letter and title: label_panels() of nc_figure_style_v2_2026-09-23.py, verbatim logic, one panel
def label_panels(fig, panels, gap_mm=1.4, above_mm=2.2):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    boxes = []
    for ax, letter, title in panels:
        ax.set_title("")
        bb = ax.get_tightbbox(r).transformed(inv)
        boxes.append((ax, letter, title, bb.x0, ax.get_position().x0))
    fw, fh = fig.get_figwidth() * 25.4, fig.get_figheight() * 25.4
    xs = {}
    for ax, _, _, _, frame_x0 in boxes:
        col = [b[3] for b in boxes if abs(b[4] - frame_x0) * fw < 2.0]
        xs[ax] = min(col)
    for ax, letter, title, _, _ in boxes:
        y = ax.get_position().y1 + above_mm / fh
        x = max(xs[ax], 0.002)
        t = fig.text(x, y, letter, fontsize=FS_PANEL, fontweight="bold", ha="left", va="baseline", color=INK)
        if title:
            fig.canvas.draw()
            w = t.get_window_extent(r).transformed(inv).width
            fig.text(x + w + gap_mm / fw, y, title, fontsize=FS, ha="left", va="baseline", color=INK)


label_panels(fig, [(ax, 'd', 'Re-tested claims (%d)' % N_ALL)])

# key, under the category labels (the panel's only free area; x tick labels and x label sit to its right). Its first
# line shares the baseline of the x tick labels; neutral grey swatches show the fill style, not a verdict.
tick_base = [s['origin'][1] for s in spans if s['text'] == '0' and abs(s['origin'][0] - 86.24) < 0.05]
assert len(tick_base) == 1                                               # panel b's '0' tick label: 288.09
KX, SW, SH, GAP, PITCH = X_LETTER, 9.0, 5.0, 2.0, 7.6   # swatch width / height, swatch-text gap, line pitch (pt)
CAP = 0.716 * FS_NOTE                                    # Arial cap height at 6 pt; swatches centred on it
KEY = [('solid', 'Site level (%d)' % N_SITE),
       ('hatched', 'Protein level (%d)%s' % (N_PROT, ':' if SHOW_NOT_ESTIMABLE else ''))]
if SHOW_NOT_ESTIMABLE:
    KEY.append((None, 'not estimable (Table 1, item 4)'))
for i, (kind, text) in enumerate(KEY):
    base = tick_base[0] + i * PITCH
    if kind:
        fc_, kw = (MUTE, {}) if kind == 'solid' else (tint(MUTE), dict(hatch=HATCH))
        y_bottom = base - CAP / 2 + SH / 2
        fig.add_artist(Rectangle(frac(KX, y_bottom), SW / W, SH / H, transform=fig.transFigure, facecolor=fc_,
                                 edgecolor=MUTE, lw=0, **kw))
    fig.text(*frac(KX + SW + GAP, base), text, fontsize=FS_NOTE, color=INK, ha='left', va='baseline')

buf = io.BytesIO()
fig.savefig(buf, format='pdf', facecolor='none', metadata={'CreationDate': None, 'ModDate': None})
plt.close(fig)

# ----------------------------------------------------------------------------------------------- assemble
# Remove the released panel d from the page's main content stream by deleting its drawing units whole: 19 text
# groups (q 1 0 0 1 x y cm ... BT ... ET Q), 6 clipped bar groups, 4 tick marks, the spine and the axes background.
# PyMuPDF redaction is NOT used: on this file its text filter re-wrote the next text object after each removed run
# with a spurious TJ offset, moving 'Enrichment' and the letter 'a' of panel a (seen in a first build of
# fig1_rev.py). Deleting self-contained units leaves every other operator of the stream byte-identical;
# qa_fig1_final.py checks panels a-c pixel by pixel against the released PDF.
NUM = r'(-?\d+(?:\.\d*)?|-?\.\d+)'
UNITS = {  # kind: (pattern, x group, y group, expected count inside REGION)
    'text': (r'q 1 0 0 1 %s %s cm [^qQ]*?BT.*?ET Q' % (NUM, NUM), 1, 2, 19),
    'bar': (r'q 1 G %s %s %s %s re W n 1 G/A2 gs (?:%s ){3}rg [-\d. ml]*?h f Q' % (NUM, NUM, NUM, NUM, NUM), 1, 2, 6),
    'tick': (r'(?<![-\d.])%s %s m \1 %s l B' % (NUM, NUM, NUM), 1, 2, 4),       # lookbehind: whole numbers only
    'spine': (r'(?<![-\d.])%s %s m %s \2 l S' % (NUM, NUM, NUM), 1, 2, 1),
    'axes_background': (r'(?<![-\d.])%s %s m %s \2 l \3 %s l \1 \4 l h f' % (NUM, NUM, NUM, NUM), 1, 4, 1),
}
main_xref = page.get_contents()[0]
stream0 = src.xref_stream(main_xref).decode('latin-1')
stream = stream0
removed = {}
for kind, (pat, gx, gy, want) in UNITS.items():
    hits = []

    def drop(m, gx=gx, gy=gy, hits=hits):
        x, y_top = float(m.group(gx)), H - float(m.group(gy))      # PDF user space is bottom-up
        if x >= REGION.x0 and y_top >= REGION.y0:
            hits.append((round(x, 2), round(y_top, 2)))
            return ''
        return m.group(0)

    stream = re.sub(pat, drop, stream)
    removed[kind] = hits
    assert len(hits) == want, (kind, hits)
assert stream.count('q') - stream.count('Q') == stream0.count('q') - stream0.count('Q')   # q/Q balance untouched
src.update_stream(main_xref, stream.encode('latin-1'))


def touches(r, R=REGION):
    """Closed-interval overlap; unlike fitz.Rect.intersects it also counts zero-width/-height rects (lines, ticks)."""
    return r.x1 >= R.x0 and r.x0 <= R.x1 and r.y1 >= R.y0 and r.y0 <= R.y1


left_text = [s['text'] for b in page.get_text('dict')['blocks'] for ln in b.get('lines', []) for s in ln['spans']
             if touches(fitz.Rect(s['bbox']))]
left_art = [d['rect'] for d in page.get_drawings() if touches(d['rect']) and d['rect'] != page.rect]
assert not left_text and not left_art, (left_text, left_art)             # the old panel d is gone (the page's own
#                                                                          white background is the only survivor)

ov = fitz.open('pdf', buf.getvalue())
assert abs(ov[0].rect.width - W) < 1e-3 and abs(ov[0].rect.height - H) < 1e-3
page.show_pdf_page(REGION, ov, 0, clip=REGION)                           # same size on both sides: scale 1
page.set_mediabox(fitz.Rect(0, 0, TRIM_W, H))
# byte-reproducible output: keep the document's permanent /ID and set the instance /ID (which MuPDF would otherwise
# draw at random on every save) from a hash of this build's inputs
id0 = re.search(r'<([0-9A-Fa-f]+)>', src.xref_get_key(-1, 'ID')[1]).group(1)
id1 = hashlib.sha256(SRC_PDF.read_bytes() + FINAL_JSON.read_bytes() + FINAL_CSV.read_bytes() +
                     pathlib.Path(__file__).read_bytes() + buf.getvalue()).hexdigest()[:32].upper()
src.xref_set_key(-1, 'ID', f'[<{id0}><{id1}>]')
meta = dict(src.metadata)                     # the released creation date stays: panels a-c were made then
meta.update(creator='%s (panels a-c, released figure); Matplotlib v%s (panel d)' % (
                meta['creator'].split(',')[0], mpl.__version__),
            producer='PyMuPDF %s (panel d placed on the released page)' % fitz.VersionBind,
            subject='Figure 1. Panel d revised 2026-09-30 with the final claim tally (results/final/'
                    'final_tally_counts.json); panels a-c unchanged from the released figure.')
src.set_metadata(meta)
src.save(OUT_PDF, garbage=3, deflate=True, no_new_id=True)
src.close()

out = fitz.open(OUT_PDF)
out[0].get_pixmap(dpi=200).save(OUT_PNG)
print(f'page {out[0].rect.width:.2f} x {out[0].rect.height:.2f} pt; panel d axes x {AX_LEFT:.2f}-{AX_RIGHT:.2f}, '
      f'y {Y_AXTOP:.2f}-{Y_SPINE:.2f}; {PT_X:.3f} pt per claim')
out.close()

# ----------------------------------------------------------------------------------------------- Source Data
old = list(csv.DictReader(SRC_CSV.open(encoding='utf-8', newline='')))
FIELDS = ['figure', 'panel', 'row', 'field', 'value', 'source_table', 'specification_label']
PLOT = 'final tally (final verdict policy, 2026-09-30); plotted'
rows = [dict(r, specification_label='') for r in old if r['panel'] in ('b', 'c')]
assert len(rows) == 8
for r in ROWS:
    rows += [dict(figure='Fig1', panel='d', row=r['label'], field='site_level_claims', value=r['site'],
                  source_table=r['src_site'], specification_label=PLOT + ' (solid segment)'),
             dict(figure='Fig1', panel='d', row=r['label'], field='protein_level_claims', value=r['prot'],
                  source_table=r['src_prot'], specification_label=PLOT + ' (lighter hatched segment)'),
             dict(figure='Fig1', panel='d', row=r['label'], field='claims', value=r['total'],
                  source_table=r['src_tot'], specification_label=PLOT + ' (bar total, printed at the bar end)')]
rows += [dict(figure='Fig1', panel='d', row='All re-tested claims', field='claims', value=N_ALL,
              source_table=f'{T2} (has_measured_baseline = 1)', specification_label=PLOT + ' (title)'),
         dict(figure='Fig1', panel='d', row='All re-tested claims', field='site_level_claims', value=N_SITE,
              source_table=f'{T2} (has_measured_baseline = 1) x {T10} (unit = site)',
              specification_label=PLOT + ' (key)'),
         dict(figure='Fig1', panel='d', row='All re-tested claims', field='protein_level_claims', value=N_PROT,
              source_table=f'{T2} (has_measured_baseline = 1) x {T10} (unit = protein)',
              specification_label=PLOT + ' (key)')]
MEMBERS = 'claim membership of the plotted segments (not plotted)'
for r in ROWS:
    for lvl, field in (('site', 'site_level_claim_ids'), ('protein', 'protein_level_claim_ids')):
        ids = members.get((r['key'], lvl), [])
        if ids:
            rows.append(dict(figure='Fig1', panel='d', row=r['label'], field=field, value='; '.join(ids),
                             source_table=f'{T2} (final_verdict = {r["key"]}) x {T10} (unit = {lvl})',
                             specification_label=MEMBERS))
rows += [dict(figure='Fig1', panel='d', row='Baseline contradicts claim', field='transfer_cohort_claim_id',
              value=bcc[0], source_table='Supplemental_Data_4_retest_round_b.csv (is_transfer = True; '
                                         'transfer_cohort = %s)' % TRANSFER_COHORT,
              specification_label='annotated transfer cohort (re-tested only on a transfer cohort)'),
         dict(figure='Fig1', panel='d', row='Vanishes', field='claims', value=0,
              source_table=f'{T2} (no final_verdict = vanishes)',
              specification_label='not plotted: no claim vanishes in the final tally')]
for t in CHANGED:          # final verdict differs from the stored (released) one: both recorded, with the reason
    cid, row = t['claim_id'], ROW_OF[t['final_verdict']]
    rows += [dict(figure='Fig1', panel='d', row=row, field=f'{cid}_final_verdict', value=t['final_verdict'],
                  source_table=f'{T2} (claim_id = {cid}; final_verdict, final_verdict_note)',
                  specification_label='not plotted: ' + t['note']),
             dict(figure='Fig1', panel='d', row=row, field=f'{cid}_stored_verdict', value=t['stored_verdict'],
                  source_table=f'{T2} (claim_id = {cid}; verdict)',
                  specification_label='not plotted: the stored (released) verdict; not tallied')]
with OUT_CSV.open('w', encoding='utf-8', newline='') as fh:    # CRLF and ASCII, as the released file
    w = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator='\r\n')
    w.writeheader()
    w.writerows(rows)
print(f'{OUT_CSV.name}: {len(rows)} rows ({sum(r["panel"] == "d" for r in rows)} in panel d)')
