"""Fig 1 (v3.3): where the five steps enter a site-mapping workflow, and the two negative sets.

Figure plan (2026-10-09), Fig 1:
  a  the workflow illustration of v3.2 Fig 1a, reused (author-confirmed, not AI-generated). In the v3.2 PDF it is one
     embedded 4128 x 754 px (600 dpi) image; it is extracted byte-for-byte with PyMuPDF and drawn unresampled
     (interpolation 'none'). Stage headers and all labels are redrawn in Helvetica: artifacts numbered 1-5 in text order
     (badges sit at the workflow positions, so they read 4-1-3-5-2 left to right), and under each a 'Diagnostic' and a
     'Control' line; artifacts 4 and 5 read 'requires a blocking arm or unenriched aliquot' and 'report search space'.
  b  selection schematic (no data): one protein, cysteines modified (filled), detected unmodified (open), not detected
     (grey dashed); two comb brackets NEG_A and NEG_B with the method definitions (author decision 2026-10-09);
     note 'A site table records modification only among detected cysteines'.
  c  background used by the classified analyses: v3.2 Fig 1c unchanged (Source_Data_Fig1_overview.csv panel c).
"""
import io

import fitz
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

import common as C

SRC_PDF = C.FIG32 / 'Fig1_overview.pdf'
F1 = C.SD32 / 'Source_Data_Fig1_overview.csv'
F1_LBL = 'cysteine-detectability-audit/source_data_submitted/Source_Data_Fig1_overview.csv'

sd = C.SourceData('Fig1', 'Source_Data_Fig1_overview.csv')

# ---------------------------------------------------------------------------------------------- illustration
doc = fitz.open(SRC_PDF)
imgs = doc[0].get_images(full=True)
assert len(imgs) == 1
info = doc.extract_image(imgs[0][0])
assert (info['width'], info['height']) == (4128, 754)
ILL = mpimg.imread(io.BytesIO(info['image']), format='png')
BBOX = doc[0].get_image_info()[0]['bbox']          # (15.56, 17.06, 510.92, 107.54) pt on the v3.2 page
sd.add('a', 'illustration', 'source', C.FIG32_LBL + 'Fig1_overview.pdf, embedded workflow image (4128 x 754 px)',
       C.FIG32_LBL + 'Fig1_overview.pdf', 'no data; reused unchanged (author-confirmed not AI-generated)')

# ---------------------------------------------------------------------------------------------- panel c data
f1 = C.read_csv(F1)
CROWS = ['Other or not stated', 'All residues, unrestricted', 'Detected but unmodified',
         'Detectability or abundance matched']
CV = []
for r in CROWS:
    n = int(C.sd_lookup(f1, 'c', r, 'classified_analyses')['value'])
    CV.append(n)
    sd.add('c', r, 'classified_analyses', n, '%s (panel c, row %s); from Supplemental_Data_9_survey_coding_table.csv' % (F1_LBL, r),
           'plotted; reused unchanged from the original figure build')
assert CV == [37, 35, 2, 0]

# ---------------------------------------------------------------------------------------------- page (pt coordinates)
W = C.W_CM / 2.54 * 72          # 515.9 pt
H = 342.0                       # 12.07 cm
fig = plt.figure(figsize=(W / 72, H / 72))
P = fig.add_axes([0, 0, 1, 1])
P.set_xlim(0, W)
P.set_ylim(H, 0)
P.axis('off')
S = W / 518.4                   # v3.2 page width -> 18.2 cm
x0, y0, x1, y1 = [v * S for v in BBOX]
y_off = 4.0
P.imshow(ILL, extent=(x0, x1, y1 + y_off, y0 + y_off), interpolation='none', zorder=1)

STAGES = [('Enrichment', 65.5), ('Digestion', 164.5), ('LC–MS/MS', 266.4), ('Database search', 369.4),
          ('Site table', 457.0)]
for name, xc in STAGES:
    P.text(xc * S, 10.5, name, fontsize=C.FS, fontweight='bold', ha='center', va='center')
P.text(3, 10.5, 'a', fontsize=C.FS_P, fontweight='bold', ha='left', va='center')

COLS = [  # (x_left v3.2 pt, artifact number, name, diagnostic, control)
    (36.5, 4, 'Positive set =\nidentified set', 'Diagnostic: share of\nidentified cysteines\ncarrying a site',
     'Needs: a blocking arm or\nunenriched aliquot'),
    (128.3, 1, 'Cleavage\ngeometry', 'Diagnostic: positional\nprofile', 'Control: background\nmatched on theoretical\ndetectability'),
    (218.1, 3, 'Protein\nabundance', 'Diagnostic: protein-level\nnegatives and abundance',
     'Control: background\nmatched on abundance\nand cysteine count'),
    (329.2, 5, 'Search\nspace', 'Diagnostic: separation in\nppm vs precursor error\nand fragment bins',
     'Needs: the search\nspace reported'),
    (423.2, 2, 'Multi-cysteine\npeptides', 'Diagnostic: site provenance\n(assigned or inferred)',
     'Control: filter both sets\nto single-cysteine peptides')]
YL = y1 + y_off + 8
for xl, n, name, diag, ctrl in COLS:
    x = xl * S
    col = C.AC[n]
    P.add_patch(Circle((x + 4.6, YL + 4.4), 4.6, color=col, zorder=3))
    P.text(x + 4.6, YL + 4.5, str(n), fontsize=C.FS, fontweight='bold', color='white', ha='center', va='center', zorder=4)
    P.text(x + 12.5, YL, name, fontsize=C.FS, fontweight='bold', color=col, ha='left', va='top', linespacing=1.05)
    P.text(x, YL + 21, diag, fontsize=C.FS_S, color=C.INK, ha='left', va='top', linespacing=1.1)
    P.text(x, YL + 45, ctrl, fontsize=C.FS_S, color=C.MUTE, ha='left', va='top', linespacing=1.1)
    sd.add('a', 'artifact %d' % n, 'badge_number', n, 'order of the five steps in Results (no data)', 'printed badge; no data')
YB = YL + 78   # bottom of panel a

# ---------------------------------------------------------------------------------------------- panel b (schematic)
TOP = YB + 14
P.text(3, TOP, 'b', fontsize=C.FS_P, fontweight='bold', ha='left', va='baseline')
P.text(13, TOP, 'Which cysteines are the negatives', fontsize=C.FS, ha='left', va='baseline')
BX0, BX1, BY = 40.0, 250.0, TOP + 62
P.add_patch(FancyBboxPatch((BX0, BY - 4), BX1 - BX0, 8, boxstyle='round,pad=0,rounding_size=4', fc='#e4e7ea',
                           ec='#b4b9bf', lw=0.6, zorder=1))
P.text(BX0 - 5, BY, 'protein', fontsize=C.FS_S, ha='right', va='center', color=C.MUTE)
STATES = ['M', 'D', 'U', 'D', 'M', 'U', 'D', 'U', 'D', 'U']
XS = [BX0 + 14 + i * (BX1 - BX0 - 28) / (len(STATES) - 1) for i in range(len(STATES))]
R = 5.0
CY = BY - 4 - R - 1.5
for x, s in zip(XS, STATES):
    P.plot([x, x], [BY - 4, CY + R], color='#8f969d', lw=0.6, zorder=2, label='_nocheck')
    if s == 'M':
        P.add_patch(Circle((x, CY), R, fc=C.AC[4], ec=C.AC[4], lw=0.9, zorder=3))
    elif s == 'D':
        P.add_patch(Circle((x, CY), R, fc='white', ec=C.INK, lw=0.9, zorder=3))
    else:
        P.add_patch(Circle((x, CY), R, fc='white', ec='#9aa1a8', lw=0.9, ls=(0, (2, 1.4)), zorder=3))
# legend of the three states (above)
LY = TOP + 16
for i, (s, lab) in enumerate((('M', 'modified (site)'), ('D', 'detected, reported unmodified'), ('U', 'not detected'))):
    lx = 40 + [0, 75, 205][i]
    if s == 'M':
        P.add_patch(Circle((lx, LY), 3.4, fc=C.AC[4], ec=C.AC[4], lw=0.8))
    elif s == 'D':
        P.add_patch(Circle((lx, LY), 3.4, fc='white', ec=C.INK, lw=0.8))
    else:
        P.add_patch(Circle((lx, LY), 3.4, fc='white', ec='#9aa1a8', lw=0.8, ls=(0, (2, 1.4))))
    P.text(lx + 6, LY, lab, fontsize=C.FS_S, ha='left', va='center')


def comb(xs, y, col, label, definition):
    P.plot([min(xs), max(xs)], [y, y], color=col, lw=1.0, solid_capstyle='butt', label='_nocheck')
    for x in xs:
        P.plot([x, x], [BY + 4 + 1.5, y], color=col, lw=0.6, ls=(0, (1.5, 1.2)), label='_nocheck')
    P.text(BX0 - 5, y, label, fontsize=C.FS, fontweight='bold', color=col, ha='right', va='center')
    P.text(min(xs), y + 4.0, definition, fontsize=C.FS_S, color=col, ha='left', va='top', zorder=5,
           bbox=dict(fc='white', ec='none', pad=0.6))


comb([x for x, s in zip(XS, STATES) if s in 'DU'], BY + 22, '#4a4f55', 'NEG_A', C.NEG_A_BRACKET)
comb([x for x, s in zip(XS, STATES) if s == 'D'], BY + 48, C.INK, 'NEG_B', C.NEG_B_BRACKET)
P.text(BX1 + 14, BY - 8, 'A site table\nrecords modification\nonly among\ndetected cysteines', fontsize=C.FS_S,
       color=C.INK, ha='left', va='center', linespacing=1.15)
P.text(XS[0], CY - R - 3, 'site', fontsize=C.FS_S, color=C.AC[4], ha='center', va='bottom')
P.text(XS[4], CY - R - 3, 'site', fontsize=C.FS_S, color=C.AC[4], ha='center', va='bottom')
sd.add('b', 'schematic', 'content', 'no data: illustrative protein with 10 cysteines (2 modified, 4 detected unmodified, 4 not detected)',
       'drawn schematic', 'no data')
sd.add('b', 'negative-set brackets', 'NEG_A', C.NEG_A_BRACKET, 'Figure 1b schematic wording (no data)', 'printed')
sd.add('b', 'negative-set brackets', 'NEG_B', C.NEG_B_BRACKET, 'Figure 1b schematic wording (no data)', 'printed')

# ---------------------------------------------------------------------------------------------- panel c
axc_l, axc_r = 405.0, 500.0
ax = fig.add_axes([axc_l / W, 1 - (TOP + 92) / H, (axc_r - axc_l) / W, 74 / H])
cols = ['#9aa1a8', '#9aa1a8', C.AC[1], C.AC[1]]
for i, (n, c) in enumerate(zip(CV, cols)):
    ax.barh(i, n, height=0.62, color=c, lw=0)
    ax.text(n + 0.8, i, str(n), va='center', ha='left', color=C.MUTE, fontsize=C.FS)
ax.set_yticks(range(4))
ax.set_yticklabels(['Other or not stated', 'All residues, unrestricted', 'Detected but unmodified',
                    'Detectability or\nabundance matched'], fontsize=C.FS)
ax.invert_yaxis()
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.set_xlim(0, 42)
ax.set_xticks([0, 20, 40])
ax.set_xlabel('Classified analyses')
P.text(320, TOP, 'c', fontsize=C.FS_P, fontweight='bold', ha='left', va='baseline')
P.text(330, TOP, 'Background used in the literature survey', fontsize=C.FS, ha='left', va='baseline')

qa = C.check_figure(fig, 'Fig1', C.number_pool(sd), allow_numbers=('1', '2', '3', '4', '5'))
C.save(fig, 'Fig1_overview', qa)
print(sd.write())
