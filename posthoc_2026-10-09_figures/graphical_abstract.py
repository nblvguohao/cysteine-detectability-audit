"""Graphical abstract (v3.3), 1328 x 531 px canvas delivered at 2656 x 1062 px (PNG) and as PDF.

Figure plan (2026-10-09), 'Graphical abstract': title strip; left, the Fig 1a illustration (reused, author-confirmed not
AI-generated) with numbered badges 1-5 and a two-to-three-word label each; middle, the selection schematic with two
brackets; right, two text boxes linked to the brackets; three recommendations; no numbers, axes or data points.
Bracket wording per author decision 2026-10-09: 'all other cysteines of the protein' (NEG_A) and 'cysteines
observed unmodified in the same run' (NEG_B). Minimum text height 22 px at 1328 px width (16 pt at 100 px/in).
Deviation (reported): the three recommendations run as one row across the full width at the bottom, because three
22-px phrases do not fit in one row of the right zone (30% = ~400 px).
"""
import io

import fitz
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Wedge

import common as C

W, H = 1328, 531
DPI_BASE = 100
fig = plt.figure(figsize=(W / DPI_BASE, H / DPI_BASE))
P = fig.add_axes([0, 0, 1, 1])
P.set_xlim(0, W)
P.set_ylim(H, 0)
P.axis('off')
PT = 72 / DPI_BASE            # px -> pt at the base size
F_MIN = 22 * PT               # 15.84 pt = 22 px
F_TXT = 23 * PT
F_TITLE = 30 * PT

# title strip (12% of height)
P.add_patch(FancyBboxPatch((0, 0), W, 64, boxstyle='square,pad=0', fc='#f1f3f5', ec='none'))
P.text(W / 2, 33, 'What a cysteine site table can show depends on how its negatives were observed',
       fontsize=F_TITLE, fontweight='bold', ha='center', va='center', color=C.INK)

# ---------------------------------------------------------------------------------------------- left: workflow
doc = fitz.open(C.FIG32 / 'Fig1_overview.pdf')
info = doc.extract_image(doc[0].get_images(full=True)[0][0])
assert (info['width'], info['height']) == (4128, 754)
ILL = mpimg.imread(io.BytesIO(info['image']), format='png')
LX0, LX1 = 16, 606
LY0 = 112
LY1 = LY0 + (LX1 - LX0) * 754 / 4128
P.imshow(ILL, extent=(LX0, LX1, LY1, LY0), interpolation='none', zorder=1)
BBOX = (15.56, 510.92)
STAGE_X = [65.5, 164.5, 266.4, 369.4, 457.0]
LABELS = [(4, 'site/\nidentified\noverlap'), (1, 'cleavage\ngeometry'), (3, 'abundance'), (5, 'search\nspace'),
          (2, 'multi-Cys\nfiltering')]
for xs, (n, lab) in zip(STAGE_X, LABELS):
    x = LX0 + (xs - BBOX[0]) / (BBOX[1] - BBOX[0]) * (LX1 - LX0)
    P.add_patch(Circle((x, LY1 + 26), 15, color=C.AC[n], zorder=3))
    P.text(x, LY1 + 27, str(n), fontsize=F_TXT, fontweight='bold', color='white', ha='center', va='center', zorder=4)
    P.text(x, LY1 + 50, lab, fontsize=F_MIN, color=C.AC[n], ha='center', va='top', linespacing=1.0)
P.text(LX0, 88, 'Five workflow steps', fontsize=F_TXT, color=C.INK, ha='left', va='center', fontweight='bold')

# ---------------------------------------------------------------------------------------------- middle: schematic
MX0, MX1 = 652, 962
P.plot([628, 628], [80, 400], color='#d9dcdf', lw=1.2, label='_nocheck')
BY = 196
P.add_patch(FancyBboxPatch((MX0, BY - 7), MX1 - MX0, 14, boxstyle='round,pad=0,rounding_size=7', fc='#e4e7ea',
                           ec='#b4b9bf', lw=1.0, zorder=1))
STATES = ['M', 'D', 'U', 'D', 'M', 'U', 'D', 'U']
XS = [MX0 + 22 + i * (MX1 - MX0 - 44) / (len(STATES) - 1) for i in range(len(STATES))]
R, CY = 13, BY - 7 - 13 - 3
for x, s in zip(XS, STATES):
    if s == 'M':
        P.add_patch(Circle((x, CY), R, fc=C.AC[4], ec=C.AC[4], lw=1.6, zorder=3))
    elif s == 'D':
        P.add_patch(Circle((x, CY), R, fc='white', ec=C.INK, lw=1.8, zorder=3))
    else:
        P.add_patch(Circle((x, CY), R, fc='white', ec='#9aa1a8', lw=1.6, ls=(0, (3, 2)), zorder=3))
# key
KY = 104
for kx, s, lab in ((MX0, 'M', 'site'), (MX0 + 82, 'D', 'detected'), (MX0 + 200, 'U', 'not detected')):
    if s == 'M':
        P.add_patch(Circle((kx + 9, KY), 9, fc=C.AC[4], ec=C.AC[4], lw=1.2))
    elif s == 'D':
        P.add_patch(Circle((kx + 9, KY), 9, fc='white', ec=C.INK, lw=1.4))
    else:
        P.add_patch(Circle((kx + 9, KY), 9, fc='white', ec='#9aa1a8', lw=1.2, ls=(0, (3, 2))))
    P.text(kx + 23, KY + 1, lab, fontsize=F_MIN, ha='left', va='center', color=C.INK)


def comb(states, y, col, text, lw):
    xs = [x for x, s in zip(XS, STATES) if s in states]
    P.plot([min(xs), max(xs)], [y, y], color=col, lw=lw, solid_capstyle='butt', label='_nocheck')
    for x in xs:
        P.plot([x, x], [BY + 8, y], color=col, lw=1.0, ls=(0, (2, 2)), label='_nocheck')
    P.text(min(xs), y + 8, text, fontsize=F_MIN, color=col, ha='left', va='top', linespacing=1.0, zorder=5,
           bbox=dict(fc='white', ec='none', pad=1.0))
    return max(xs), y


ea = comb('DU', 240, '#4a4f55', C.NEG_A_BRACKET.replace('cysteines of', 'cysteines\nof'), 2.2)
eb = comb('D', 330, C.INK, C.NEG_B_BRACKET.replace('observed unmodified', 'observed\nunmodified'), 2.6)

# ---------------------------------------------------------------------------------------------- right: messages
RX0, RX1 = 1000, 1312
for (ex, ey), txt, fc in ((ea, 'Peptide detectability\nalone separates sites', '#e3edf6'),
                          (eb, 'Much of that\nseparation disappears', '#eceef0')):
    P.add_patch(FancyBboxPatch((RX0, ey - 34), RX1 - RX0, 68, boxstyle='round,pad=0,rounding_size=8', fc=fc,
                               ec='none', zorder=1))
    P.text((RX0 + RX1) / 2, ey, txt, fontsize=F_TXT, ha='center', va='center', color=C.INK, linespacing=1.05,
           zorder=2)
    P.add_patch(FancyArrowPatch((ex + 10, ey), (RX0 - 6, ey), arrowstyle='-|>', mutation_scale=16, lw=1.6,
                                color='#6f757c', zorder=2))
P.text(RX0, 160, 'Against these negatives:', fontsize=F_MIN, ha='left', va='center', color=C.MUTE)

# ---------------------------------------------------------------------------------------------- bottom row: advice
P.plot([22, W - 22], [418, 418], color='#d9dcdf', lw=1.2, label='_nocheck')
REC = [('report the identified-site share', 'share'), ('use observed negatives', 'neg'),
       ('report the search space', 'search')]
RY = 470
for i, (txt, icon) in enumerate(REC):
    x = 60 + i * 440
    if icon == 'share':
        P.add_patch(Circle((x, RY), 17, fc='#e4e7ea', ec='#4a4f55', lw=1.4))
        P.add_patch(Wedge((x, RY), 17, -90, 120, fc=C.AC[4], ec='none'))
    elif icon == 'neg':
        P.add_patch(Circle((x, RY), 15, fc='white', ec=C.INK, lw=2.4))
    else:
        P.add_patch(Circle((x - 3, RY - 3), 11, fc='white', ec=C.AC[5], lw=2.6))
        P.plot([x + 5, x + 15], [RY + 5, RY + 15], color=C.AC[5], lw=3.2, solid_capstyle='round', label='_nocheck')
    P.text(x + 30, RY + 1, txt, fontsize=F_TXT, ha='left', va='center', color=C.INK)

qa = C.check_figure(fig, 'Graphical_Abstract', [], allow_numbers=('1', '2', '3', '4', '5'))
for t, _ in C.all_texts(fig):
    if t.get_fontsize() < F_MIN - 1e-6:
        qa['issues'].append('GA text below 22 px: %r' % t.get_text())
qa['issues'] = [i for i in qa['issues'] if not i.startswith('font ')]   # 6-pt rule is for the main figures
C.save(fig, 'Graphical_Abstract', qa, png_dpi=2 * DPI_BASE)
sd = C.SourceData('Graphical_Abstract', 'Source_Data_Graphical_Abstract.csv')
sd.add('-', 'content', 'data', 'none: no numbers, axes or data points; workflow illustration reused from the original figure build',
       C.FIG32_LBL + 'Fig1_overview.pdf (embedded workflow image)', 'no data')
print(sd.write())
