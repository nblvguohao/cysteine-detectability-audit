"""Supplementary figures S1, S2, S3, S5, S6, S7 (v3.3). Fig S4 is drawn by fig5_claim_retests.py.

Figure plan (2026-10-09), 'Figures moved to the Supplemental':
  S1 = v3.2 Fig 2 whole (four proteases)        - v3.2 PDF reused byte-for-byte
  S2 = v3.2 Fig 3 whole (PXD015307)             - v3.2 PDF reused byte-for-byte
  S3 = v3.2 Fig 7 whole (self-audit)            - v3.2 PDF reused byte-for-byte
  S5 = v3.2 Fig 5b (PERS-008 background) and 5c (five structural claims) - the two panels are placed 1:1 from the v3.2
       PDF (vector, PyMuPDF show_pdf_page with clips that leave out the old letters and titles; v3.2 Fig 5b and 5c were
       themselves carried over from the submitted PDF, see fig5_rebuild.py); letters a, b and titles redrawn
  S6 = v3.2 Fig 4 whole (synthetic benchmark)   - v3.2 PDF reused byte-for-byte
  S7 = v3.2 Fig 1b (survey funnel)              - redrawn from Source_Data_Fig1_overview.csv panel b (4 values)
PNG previews are rendered from the PDFs at 300 dpi. Reused figures keep their Arial text (unchanged).
"""
import shutil

import fitz
import matplotlib.pyplot as plt

import common as C

REUSE = [('FigS1_public_four_protease', 'Fig2_public_four_protease', 'Source_Data_Fig2_public_four_protease.csv', 'Figure 2'),
         ('FigS2_search_space', 'Fig3_search_space', 'Source_Data_Fig3_search_space.csv', 'Figure 3'),
         ('FigS3_self_audit_public', 'Fig7_self_audit_public', 'Source_Data_Fig7_self_audit_public.csv', 'Figure 7'),
         ('FigS6_synthetic_ground_truth', 'Fig4_synthetic_ground_truth', 'Source_Data_Fig4_synthetic_ground_truth.csv', 'Figure 4')]


REUSED = 'reused unchanged from the original figure build'


def render_png(pdf, png):
    C.refuse_existing(png)
    d = fitz.open(pdf)
    assert len(d) == 1
    d[0].get_pixmap(dpi=300).save(str(png))


for stem, old, sdname, label in REUSE:
    src = C.FIG32 / (old + '.pdf')
    pdf = C.OUT_FIG / (stem + '.pdf')
    C.refuse_existing(pdf)
    shutil.copyfile(src, pdf)
    render_png(pdf, C.OUT_FIG / (stem + '.png'))
    rows = C.read_csv(C.SD32 / sdname)
    fid = stem.split('_')[0]
    sd = C.SourceData(fid, 'Source_Data_%s.csv' % stem)
    origin = C.SD32_LBL + sdname
    for r in rows:
        if 'benchmark' in r:      # v3.2 Fig 4 layout: figure,panel,benchmark,quantity,strength,value
            sd.add(r['panel'], '%s|strength %s' % (r['benchmark'], r['strength']), r['quantity'], r['value'], origin,
                   REUSED)
        else:
            st = C.current_folder(r.get('source_table', ''))
            sd.add(r['panel'], r['row'], r['field'], r['value'], '%s (originally: %s)' % (origin, st) if st else origin,
                   REUSED + ('; ' + r['specification_label'] if r.get('specification_label') else ''))
    print(stem, len(rows), 'rows ->', sd.write())

# ---------------------------------------------------------------------------------------------- S5
SRC5 = C.FIG32 / 'Fig5_three_axes.pdf'
d5 = fitz.open(SRC5)
p5 = d5[0]
assert abs(p5.rect.width - 518.4) < 0.01
spans = [s for b in p5.get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans']]
# clips (v3.2 page points): panel b body, panel c header line (right column), panel c body
CLIP_B = fitz.Rect(310.0, 17.0, 518.4, 236.5)
CLIP_CH = fitz.Rect(290.0, 277.0, 518.4, 288.6)
CLIP_C = fitz.Rect(0.0, 288.6, 518.4, 418.0)
for s in spans:   # the old letters and titles lie outside every clip
    if s['text'].strip() in ('b', 'c', 'Background choice', 'Attribute definition'):
        r = fitz.Rect(s['bbox'])
        assert not any(r.intersects(c) for c in (CLIP_B, CLIP_C)), s['text']
        assert not (r.intersects(CLIP_CH) and s['text'].strip() != ''), s['text']
# Remove everything outside the three clips from an in-memory copy of the v3.2 page, so that the old panel a, letters
# and titles do not survive as hidden content of the embedded page; then check the kept regions render unchanged.
src_bytes = d5.tobytes()
d5 = fitz.open('pdf', src_bytes)
p5 = d5[0]
pr = p5.rect
KEEP = [CLIP_B, CLIP_CH, CLIP_C]
cut = [fitz.Rect(0, 0, pr.x1, CLIP_B.y0), fitz.Rect(0, CLIP_B.y0, CLIP_B.x0, CLIP_CH.y0),
       fitz.Rect(0, CLIP_B.y1, pr.x1, CLIP_CH.y0), fitz.Rect(0, CLIP_CH.y0, CLIP_CH.x0, CLIP_C.y0),
       fitz.Rect(0, CLIP_C.y1, pr.x1, pr.y1)]
for r in cut:
    p5.add_redact_annot(r)
p5.apply_redactions(images=fitz.PDF_REDACT_IMAGE_REMOVE, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                    text=fitz.PDF_REDACT_TEXT_REMOVE)
orig = fitz.open(SRC5)[0]
for c in KEEP:
    a = orig.get_pixmap(dpi=200, clip=c).samples
    b = p5.get_pixmap(dpi=200, clip=c).samples
    diff = sum(1 for x, y in zip(a, b) if x != y)
    print('kept region', c, 'differing bytes after redaction:', diff)
    assert diff == 0, c
left = p5.get_text()
for w in ('NEG_A', 'QTRP', 'Background choice', 'Attribute definition', 'Global AUC'):
    assert w not in left, w
WP = C.W_IN * 72
SC = WP / 518.4
out = fitz.open()
HP = 16 + (CLIP_B.height + 26 + CLIP_CH.height + CLIP_C.height) * SC + 4
page = out.new_page(width=WP, height=HP)
yb = 16
page.show_pdf_page(fitz.Rect(0, yb, CLIP_B.width * SC, yb + CLIP_B.height * SC), d5, 0, clip=CLIP_B)
yc = yb + CLIP_B.height * SC + 26
page.show_pdf_page(fitz.Rect(CLIP_CH.x0 * SC, yc, CLIP_CH.x1 * SC, yc + CLIP_CH.height * SC), d5, 0, clip=CLIP_CH)
page.show_pdf_page(fitz.Rect(0, yc + CLIP_CH.height * SC, WP, yc + (CLIP_CH.height + CLIP_C.height) * SC), d5, 0, clip=CLIP_C)
page.insert_font(fontname='HelvB', fontfile=str(C._FONT_DIR / 'Helvetica-Bold.ttf'))
page.insert_font(fontname='HelvR', fontfile=str(C._FONT_DIR / 'Helvetica-Regular.ttf'))
for x, y, letter, title in ((3, 12, 'a', 'Background choice (PERS-008)'),
                            (3, yc + 8.5, 'b', 'Attribute definition (five structural claims)')):
    page.insert_text((x, y), letter, fontname='HelvB', fontsize=8)
    page.insert_text((x + 7.5, y), title, fontname='HelvR', fontsize=7)
pdf = C.OUT_FIG / 'FigS5_background_and_attribute.pdf'
C.refuse_existing(pdf)
out.save(str(pdf), garbage=3, deflate=True)
render_png(pdf, C.OUT_FIG / 'FigS5_background_and_attribute.png')
sd = C.SourceData('FigS5', 'Source_Data_FigS5_background_and_attribute.csv')
for r in C.read_csv(C.SD32 / 'Source_Data_Fig5_three_axes.csv'):
    if r['panel'] in ('b', 'c'):
        new = {'b': 'a', 'c': 'b'}[r['panel']]
        sd.add(new, r['row'], r['field'], r['value'],
               'cysteine-detectability-audit/source_data_submitted/Source_Data_Fig5_three_axes.csv (panel %s; originally: %s)' % (r['panel'], C.current_folder(r['source_table'])),
               REUSED + ('; ' + r['specification_label'] if r['specification_label'] else ''))
print('FigS5 ->', sd.write())

# ---------------------------------------------------------------------------------------------- S7
f1 = C.read_csv(C.SD32 / 'Source_Data_Fig1_overview.csv')
ROWS = ['Passed screening', 'Full text retrieved', 'Sampled for coding', 'Classified']
V = [int(C.sd_lookup(f1, 'b', r, 'records')['value']) for r in ROWS]
assert V == [942, 367, 120, 74]
sd = C.SourceData('FigS7', 'Source_Data_FigS7_literature_survey.csv')
for r, v in zip(ROWS, V):
    sd.add('-', r, 'records', v, 'cysteine-detectability-audit/source_data_submitted/Source_Data_Fig1_overview.csv (panel b, row %s); from '
           'Supplemental_Data_9_survey_coding_table.csv' % r, 'plotted; redrawn from the original figure build (same values and colours)')
COLS = ['#c9d9e8', '#95b6d3', '#5f90bb', '#1f6aa5']     # v3.2 Fig 1b fills (read from the PDF)
fig = plt.figure(figsize=(C.mm(90), C.mm(45)))
ax = fig.add_axes([0.36, 0.25, 0.56, 0.62])
for i, (v, c) in enumerate(zip(V, COLS)):
    ax.barh(i, v, height=0.62, color=c, lw=0)
    ax.text(v + 15, i, format(v, ','), va='center', ha='left', color=C.MUTE)
ax.set_yticks(range(4))
ax.set_yticklabels(ROWS)
ax.invert_yaxis()
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.set_xlim(0, 1100)
ax.set_xticks([0, 500, 1000])
ax.set_xticklabels(['0', '500', '1,000'])
ax.set_xlabel('Records')
fig.text(0.02, 0.93, 'Literature survey', fontsize=C.FS, va='baseline')
qa = C.check_figure(fig, 'FigS7', C.number_pool(sd))
C.save(fig, 'FigS7_literature_survey', qa)
print(sd.write())
