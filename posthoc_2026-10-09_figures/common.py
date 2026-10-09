"""Shared style, inputs and output helpers for the MCP v3.3 figures (built 2026-10-09).

Rules carried by every build script:
  * values are only READ from existing files of this repository (paths below); nothing is re-estimated;
  * every plotted or printed value is written to that figure's Source Data CSV with the file it came from;
  * outputs go to OUT_DIR (default: posthoc_2026-10-09_figures/build_output/) and are refused if the target file
    already exists (no overwriting); a figure whose canvas QA reports an issue is refused;
  * Helvetica (macOS system face, extracted from Helvetica.ttc into a temporary folder at run time and embedded
    as a TrueType subset); never Arial unless Helvetica is missing (then the build stops).
Run with the system python3 (3.9; matplotlib 3.9.4, PyMuPDF 1.26), from any working directory, for example
    python3 posthoc_2026-10-09_figures/fig1_overview.py
"""
import csv
import io
import json
import os
import pathlib
import re
import tempfile

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.font_manager as fm  # noqa: E402
import matplotlib.text  # noqa: E402,F401
import matplotlib._mathtext as _mt  # noqa: E402

# sub/superscripts (log2) at 0.86 of the base size, so that on a 7-pt base they stay >= 6 pt (2.1 mm)
_mt.SHRINK_FACTOR = 0.86
import matplotlib.pyplot as plt  # noqa: E402

# ------------------------------------------------------------------------------------------------ paths
# every input is read from this repository, relative to its root
REPO = pathlib.Path(__file__).resolve().parent.parent
SD32 = REPO / 'source_data_submitted'      # Source Data of the earlier (v3.2) figures, whose values are re-used
SUP = REPO / 'supplemental'                # Supplemental Notes (v3.3; the tables read here are unchanged from v3.2)
FIG32 = REPO / 'figures_submitted'         # the earlier (v3.2) figures, re-used for Figs S1-S3, S5, S6 and Fig 1a
# provenance labels written into the Source Data (paths relative to the public repository)
SD32_LBL = 'cysteine-detectability-audit/source_data_submitted/'
FIG32_LBL = 'cysteine-detectability-audit/figures_submitted/'
OUT = pathlib.Path(os.environ.get('OUT_DIR', str(REPO / 'posthoc_2026-10-09_figures' / 'build_output'))).resolve()
for _protected in ('figures_final', 'source_data_final'):
    if OUT == REPO / _protected or (REPO / _protected) in OUT.parents:
        raise SystemExit('REFUSE: OUT_DIR inside the released folder %s/' % _protected)
OUT_FIG = OUT / 'figures'
OUT_SD = OUT / 'source_data'
QA_DIR = OUT / 'qa_logs'
for _d in (OUT_FIG, OUT_SD, QA_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------------------------------------ fonts
_FONT_DIR = pathlib.Path(tempfile.mkdtemp(prefix='mcp33_fonts_'))


def _register_helvetica():
    from fontTools.ttLib import TTCollection
    src = '/System/Library/Fonts/Helvetica.ttc'
    if not os.path.exists(src):
        raise SystemExit('Helvetica.ttc not found; refusing to fall back silently')
    col = TTCollection(src)
    want = {0: 'Helvetica-Regular', 1: 'Helvetica-Bold', 2: 'Helvetica-Oblique'}
    for i, name in want.items():
        got = col.fonts[i]['name'].getDebugName(4).strip()
        assert got.replace(' ', '-').startswith(name.replace('-Regular', '')), (i, got)
        p = _FONT_DIR / ('%s.ttf' % name)
        col.fonts[i].save(str(p))
        fm.fontManager.addfont(str(p))


_register_helvetica()

INK, AXIS, MUTE, REF, BAND = '#1f2326', '#4a4f55', '#6f757c', '#c3c8ce', '#f1f3f5'
FS, FS_S, FS_P = 7, 6, 8          # body, smallest, panel letter (pt). 6 pt = 2.12 mm em height
LW = 0.8
W_CM = 18.2
W_IN = W_CM / 2.54
# verdict colours of v3.2 Fig. 6 (read from Fig6_claim_retests.pdf fills) and artifact colours of v3.2 Fig. 1a
VC = {'survives': '#1f6aa5', 'attenuated': '#6fa9d6', 'undecidable': '#9aa1a8',
      'null_broken_by_control': '#c06c98', 'baseline_contradicts_claim': '#cc5a36'}
VLABEL = {'survives': 'survives', 'attenuated': 'attenuated', 'undecidable': 'undecidable',
          'null_broken_by_control': 'null broken by control',
          'baseline_contradicts_claim': 'baseline contradicts claim'}
AC = {1: '#1f6aa5', 2: '#dc9a2b', 3: '#23906f', 4: '#c06c98', 5: '#cc5a36'}
NEG_A_DEF = 'every other cysteine of a protein carrying a site'
NEG_B_DEF = 'cysteines the same experiment reported unmodified'
# bracket wording for Fig 1b and the graphical abstract (author decision 2026-10-09)
NEG_A_BRACKET = 'all other cysteines of the protein'
NEG_B_BRACKET = 'cysteines observed unmodified in the same run'

RC = {
    'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica'],
    'font.size': FS, 'axes.titlesize': FS, 'axes.labelsize': FS, 'xtick.labelsize': FS, 'ytick.labelsize': FS,
    'legend.fontsize': FS, 'axes.linewidth': LW, 'lines.linewidth': LW, 'xtick.major.width': LW,
    'ytick.major.width': LW, 'patch.linewidth': LW, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.titlelocation': 'left', 'axes.titlepad': 4, 'legend.frameon': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.edgecolor': AXIS, 'xtick.color': AXIS, 'ytick.color': AXIS,
    'text.color': INK, 'axes.labelcolor': INK, 'xtick.labelcolor': INK, 'ytick.labelcolor': INK,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.major.pad': 2.0, 'ytick.major.pad': 2.0,
    'axes.labelpad': 3.0, 'lines.solid_capstyle': 'butt', 'axes.unicode_minus': True,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Helvetica', 'mathtext.it': 'Helvetica:italic',
    'mathtext.bf': 'Helvetica:bold', 'pdf.compression': 6,
}
plt.rcParams.update(RC)


# ------------------------------------------------------------------------------------------------ Source Data
FIELDS = ['figure', 'panel', 'row', 'field', 'value', 'source_table', 'specification_label']


class SourceData:
    def __init__(self, fig_id, out_name):
        self.fig_id, self.rows = fig_id, []
        self.path = OUT_SD / out_name

    def add(self, panel, row, field, value, source_table, spec=''):
        self.rows.append(dict(figure=self.fig_id, panel=panel, row=row, field=field, value=value,
                              source_table=source_table, specification_label=spec))
        return value

    def write(self):
        refuse_existing(self.path)
        with open(self.path, 'w', encoding='utf-8', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator='\r\n')
            w.writeheader()
            w.writerows(self.rows)
        return self.path


def refuse_existing(p):
    p = pathlib.Path(p)
    if p.exists():
        raise SystemExit('REFUSE: %s exists (no overwriting)' % p)


def current_folder(source_table):
    """The earlier Source Data name some tables relative to posthoc_2026-09-30/results/; give the full path."""
    if source_table.startswith('H_artifact5_docs/'):
        return 'posthoc_2026-09-30/results/' + source_table
    return source_table


def read_csv(p):
    with open(p, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh))


def sd_lookup(rows, panel, row, field):
    hit = [r for r in rows if r['panel'] == panel and r['row'] == row and r['field'] == field]
    if len(hit) != 1:
        raise SystemExit('lookup failed (%d hits): %s %s %s' % (len(hit), panel, row, field))
    return hit[0]


def md_table(path, header_startswith):
    """Rows of the first markdown table whose header line starts with header_startswith."""
    lines = pathlib.Path(path).read_text(encoding='utf-8').splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith(header_startswith):
            hdr = [c.strip() for c in ln.strip().strip('|').split('|')]
            out = []
            for ln2 in lines[i + 2:]:
                if not ln2.startswith('|'):
                    break
                cells = [c.strip() for c in ln2.strip().strip('|').split('|')]
                out.append(dict(zip(hdr, cells)))
            return out
    raise SystemExit('table not found: %s in %s' % (header_startswith, path))


def num(s):
    s = s.replace('−', '-').replace(',', '').replace('`', '').strip()
    return float(s)


def pt_ci(s):
    """'0.815 [0.763, 0.858]' -> (0.815, 0.763, 0.858)"""
    m = re.match(r'\s*([-−+\d.]+)\s*\[\s*([-−+\d.]+)\s*,\s*([-−+\d.]+)\s*\]', s)
    if not m:
        raise SystemExit('cannot parse %r' % s)
    return tuple(num(g) for g in m.groups())


# ------------------------------------------------------------------------------------------------ drawing helpers
def mm(x):
    return x / 25.4


def panel_letter(fig, ax_or_xy, letter, title=None, dx_mm=0.0, dy_mm=2.2):
    """Bold letter (8 pt) + optional title (7 pt) placed above the left edge of the tight bbox of ax."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    if isinstance(ax_or_xy, tuple):
        x, y = ax_or_xy
    else:
        bb = ax_or_xy.get_tightbbox(r).transformed(inv)
        x = bb.x0 + mm(dx_mm) / fig.get_figwidth()
        y = ax_or_xy.get_position().y1 + mm(dy_mm) / fig.get_figheight()
    t = fig.text(x, y, letter, fontsize=FS_P, fontweight='bold', ha='left', va='baseline', color=INK)
    if title:
        fig.canvas.draw()
        w = t.get_window_extent(r).transformed(inv).width
        fig.text(x + w + mm(1.4) / fig.get_figwidth(), y, title, fontsize=FS, ha='left', va='baseline', color=INK)
    return t


def all_texts(fig):
    out = []
    for ax in fig.get_axes():
        if not ax.axison:
            out += [(t, ax) for t in list(ax.texts) + (list(ax.get_legend().get_texts()) if ax.get_legend() else [])]
            continue
        c = [ax.title, ax._left_title, ax._right_title, ax.xaxis.label, ax.yaxis.label]
        c += list(ax.get_xticklabels()) + list(ax.get_yticklabels()) + list(ax.texts)
        leg = ax.get_legend()
        if leg is not None:
            c += list(leg.get_texts())
        out += [(t, ax) for t in c]
    out += [(t, None) for t in fig.texts]
    return [(t, a) for t, a in out if t.get_visible() and t.get_text().strip()]


def check_figure(fig, name, number_pool, allow_numbers=(), tol_px=0.5):
    """Programmatic QA on the matplotlib canvas: text-text overlap, text inside the page, font size >= 6 pt,
    and every number printed in a non-tick text found in number_pool (Source Data values) or allow_numbers."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fb = fig.bbox
    items = all_texts(fig)
    boxes, issues = [], []
    ticks = set()
    for ax in fig.get_axes():
        ticks |= {id(t) for t in ax.get_xticklabels() + ax.get_yticklabels()}
    for t, ax in items:
        b = mpl.text.Text.get_window_extent(t, r)
        boxes.append((t, b))
        if t.get_fontsize() < FS_S - 1e-6:
            issues.append('font %.2f pt < 6: %r' % (t.get_fontsize(), t.get_text()))
        if b.x0 < fb.x0 - 0.5 or b.x1 > fb.x1 + 0.5 or b.y0 < fb.y0 - 0.5 or b.y1 > fb.y1 + 0.5:
            issues.append('outside page: %r' % t.get_text())
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            ix = min(a.x1, b.x1) - max(a.x0, b.x0) - tol_px
            iy = min(a.y1, b.y1) - max(a.y0, b.y0) - tol_px
            if ix > 0 and iy > 0:
                issues.append('overlap: %r <> %r' % (boxes[i][0].get_text(), boxes[j][0].get_text()))
    # text over data points (markers / line vertices) inside the same axes
    for ax in fig.get_axes():
        if not ax.axison:
            continue
        tb = [(t, mpl.text.Text.get_window_extent(t, r)) for t in ax.texts if t.get_visible() and t.get_text().strip()]
        for ln in ax.lines:
            if ln.get_label() == '_nocheck':
                continue
            pts = ax.transData.transform(ln.get_xydata())
            for t, b in tb:
                for x, y in pts:
                    if b.x0 + 1 < x < b.x1 - 1 and b.y0 + 1 < y < b.y1 - 1:
                        issues.append('text over data point: %r' % t.get_text())
                        break
    # numbers printed outside tick labels
    unmatched, matched = [], []
    for t, ax in items:
        if id(t) in ticks:
            continue
        s = t.get_text().replace('−', '-').replace(' ', '')
        s = re.sub(r'\$[^$]*\$', '', s)                     # mathtext (log_2 etc.)
        for tok in re.findall(r'(?<![A-Za-z\d.])-?\d[\d,]*\.?\d*', s):
            tok_c = tok.replace(',', '').rstrip('.')
            if tok_c in allow_numbers:
                matched.append((tok, 'allowed'))
                continue
            dec = len(tok_c.split('.')[1]) if '.' in tok_c else 0
            v = float(tok_c)
            ok = any(abs(round(p, dec) - v) < 10 ** (-dec) * 0.51 or abs(abs(round(p, dec)) - abs(v)) < 1e-12
                     for p in number_pool)
            (matched if ok else unmatched).append((tok, t.get_text()))
    return dict(name=name, n_texts=len(items), issues=issues, numbers_matched=len(matched),
                numbers_unmatched=unmatched)


def number_pool(sd):
    pool = []
    for r in sd.rows:
        try:
            pool.append(float(str(r['value']).replace(',', '')))
        except ValueError:
            pass
    return pool


def save(fig, stem, qa, png_dpi=300):
    """Write <stem>.pdf and <stem>.png (300 dpi) into OUT_DIR/figures; refuse if either exists."""
    pdf, png = OUT_FIG / (stem + '.pdf'), OUT_FIG / (stem + '.png')
    refuse_existing(pdf)
    refuse_existing(png)
    if qa['issues'] or qa['numbers_unmatched']:
        raise SystemExit('build refused: canvas QA issues %r' % (qa['issues'] + qa['numbers_unmatched'],))
    fig.savefig(pdf, metadata={'CreationDate': None, 'ModDate': None, 'Creator': 'MCP v3.3 figure_build'})
    fig.savefig(png, dpi=png_dpi)
    QA_DIR.mkdir(exist_ok=True)
    with open(QA_DIR / (stem + '_canvas_qa.json'), 'w', encoding='utf-8') as fh:
        json.dump(qa, fh, indent=1, ensure_ascii=False)
    w_cm, h_cm = fig.get_figwidth() * 2.54, fig.get_figheight() * 2.54
    print('%s: %.2f x %.2f cm; texts %d; issues %d; numbers matched %d, unmatched %d' % (
        stem, w_cm, h_cm, qa['n_texts'], len(qa['issues']), qa['numbers_matched'], len(qa['numbers_unmatched'])))
    for i in qa['issues']:
        print('   ISSUE', i)

    for u in qa['numbers_unmatched']:
        print('   UNMATCHED NUMBER', u)
    return pdf, png
