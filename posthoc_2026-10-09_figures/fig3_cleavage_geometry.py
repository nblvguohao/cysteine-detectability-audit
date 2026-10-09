"""Fig 3 (v3.3): tryptic cleavage geometry creates a basic-residue signature that follows the protease.

Figure plan (2026-10-09), Fig 3:
  a  positional K/R profile, dataset_id == natcomm2023_ath_sno only (the file's other rows are dropped on read and never
     written anywhere), offsets -7..+7 without 0, ratio_vs_proteome with ci95_lo/hi_vs_proteome, dashed y = 1,
     offsets -1, +5 and +7 labelled;
  b  KRH enrichment z before and after theoretical-detectability matching (three Arabidopsis sets) plus the acidic
     (D/E) control of the reversible-oxidation set; artifacts12/results/12_matched_{sno,so,ox}.json;
  c  four-protease deposit PXD063463, distal band (6-12) log2 OR against the proteome background, own rule (filled)
     and trypsin rule (open), 97.5% (Bonferroni) intervals, identified cysteines per arm;
  d  the same statistic against the observed background.
"""
import json

import matplotlib.pyplot as plt
import pandas as pd

import common as C

POS = C.REPO / 'artifacts12/positional_profile/results/positional_kr_profile_public_2026-09-17.csv'
POS_LBL = 'cysteine-detectability-audit/artifacts12/positional_profile/results/positional_kr_profile_public_2026-09-17.csv'
MJ = {k: C.REPO / ('artifacts12/results/12_matched_%s.json' % k) for k in ('sno', 'so', 'ox')}
MJ_LBL = 'cysteine-detectability-audit/artifacts12/results/12_matched_%s.json'
FP = C.SD32 / 'Source_Data_Fig2_public_four_protease.csv'
FP_LBL = 'cysteine-detectability-audit/source_data_submitted/Source_Data_Fig2_public_four_protease.csv'
N4 = C.SUP / 'Supplemental_Note_4_four_protease_public_deposit.md'
N4_LBL = 'Supplemental Note 4'
DS = 'natcomm2023_ath_sno'

sd = C.SourceData('Fig3', 'Source_Data_Fig3_cleavage_geometry.csv')

# ---------------------------------------------------------------------------------------------- panel a data
pos = pd.read_csv(POS)
pos = pos[pos['dataset_id'] == DS].copy()             # filter first; nothing else from the file is kept
assert len(pos) == 14 and sorted(pos['offset']) == [-7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7]
PA = {}
for _, r in pos.sort_values('offset').iterrows():
    o = int(r['offset'])
    PA[o] = (float(r['ratio_vs_proteome']), float(r['ci95_lo_vs_proteome']), float(r['ci95_hi_vs_proteome']))
    for f in ('ratio_vs_proteome', 'ci95_lo_vs_proteome', 'ci95_hi_vs_proteome'):
        sd.add('a', 'offset %+d' % o, f, float(r[f]), '%s (dataset_id %s, offset %d, %s)' % (POS_LBL, DS, o, f),
               'plotted; K/R fraction at the offset relative to the proteome background, 95% interval')
assert PA[-1][0] == 0.4576 and PA[5][0] == 1.4975 and PA[7][0] == 1.3831
for o in (-1, 5, 7):   # printed at the manuscript's rounding (revision 2026-10-09); the point itself is unrounded
    sd.add('a', 'offset %+d' % o, 'ratio_vs_proteome_printed', '%.2f' % PA[o][0], '%s (dataset_id %s, offset %d, ratio_vs_proteome)' % (POS_LBL, DS, o),
           'printed label, rounded to 2 decimals as in the text')
sd.add('a', 'dataset', 'dataset_id', DS, POS_LBL,
       'Arabidopsis S-nitrosylation site table (FAT-switch cohort, Supplemental Data 7); label drawn in panel title')

# ---------------------------------------------------------------------------------------------- panel b data
NAMES = {'sno': 'S-nitrosylation', 'so': 'S-sulfenylation', 'ox': 'Reversible oxidation'}
KEYS = {'KRH': '碱性 KRH', 'DE': '酸性 DE'}
PB = []
for k, cls in (('sno', 'KRH'), ('so', 'KRH'), ('ox', 'KRH'), ('ox', 'DE')):
    j = json.loads(MJ[k].read_text(encoding='utf-8'))
    c = j['classes'][KEYS[cls]]
    raw, mat = c['raw']['z'], c['matched']['z']
    lab = '%s, %s' % (NAMES[k], 'K/R/H' if cls == 'KRH' else 'D/E (control)')
    PB.append((lab, raw, mat))
    for stage, v in (('raw', raw), ('matched', mat)):
        sd.add('b', lab, 'z_%s' % stage, v, '%s (classes["%s"].%s.z)' % (MJ_LBL % k, 'basic KRH' if cls == 'KRH' else 'acidic DE', stage),
               'plotted; open = all other cysteines of the proteins, filled = theoretically detectable cysteines')
assert [(round(a, 2), round(b, 2)) for _, a, b in PB] == [(15.72, 8.22), (10.03, 4.71), (5.59, 0.67), (6.04, 6.62)]

# ---------------------------------------------------------------------------------------------- panels c, d data
fp = C.read_csv(FP)
n4 = C.md_table(N4, '| arm | identified cysteines (+hydroxylamine)')
NID = {r['arm']: int(r['identified cysteines (+hydroxylamine)']) for r in n4}
assert NID == {'trypsin': 5328, 'AspN': 187, 'chymotrypsin': 728, 'GluC': 884}
ARMS = [('AspN', 'aspn', 'AspN', 'AspN'), ('CT', 'chymotrypsin', 'Chymotrypsin', 'chymotrypsin'),
        ('GluC', 'gluc', 'GluC', 'GluC'), ('Trypsin', 'trypsin', 'Trypsin', 'trypsin')]
PC, PD = {}, {}


def trio(panel, row):
    v = [float(C.sd_lookup(fp, panel, row, f)['value']) for f in ('estimate', 'ci_low', 'ci_high')]
    return tuple(v)


for code, own, disp, n4name in ARMS:
    rows_c = [('own', '%s|%s|proteome|distal_6_12' % (code, own))]
    if code != 'Trypsin':
        rows_c.append(('trypsin', '%s|trypsin|proteome|distal_6_12' % code))
    for rule, row in rows_c:
        PC[(code, rule)] = trio('a', row)
        for f, v in zip(('estimate', 'ci_low', 'ci_high'), PC[(code, rule)]):
            sd.add('c', row, f, v, '%s (panel a, row %s, field %s)' % (FP_LBL, row, f),
                   'plotted; distal band 6-12, proteome background, %s rule; 97.5%% (Bonferroni) interval' % rule)
    row = '%s|%s|observed|distal_6_12' % (code, own)
    PD[code] = trio('b', row)
    for f, v in zip(('estimate', 'ci_low', 'ci_high'), PD[code]):
        sd.add('d', row, f, v, '%s (panel b, row %s, field %s)' % (FP_LBL, row, f),
               'plotted; distal band 6-12, observed background (identified cysteines without a site), own rule; 97.5% interval')
    for pnl in ('c', 'd'):
        sd.add(pnl, disp, 'identified_cysteines', NID[n4name],
               '%s (input table, identified cysteines (+hydroxylamine), arm %s)' % (N4_LBL, n4name), 'printed next to the arm')
# cross-check against Note 4's cleavage-geometry table
n4g = C.md_table(N4, '| arm | rule | background |')
for r in n4g:
    code = {'trypsin': 'Trypsin', 'AspN': 'AspN', 'chymotrypsin': 'CT', 'GluC': 'GluC'}[r['arm']]
    est, lo, hi = C.pt_ci(r['distal 6-12'])
    if r['background'] == 'proteome':
        got = PC[(code, 'trypsin' if (r['rule'] == 'trypsin' and code != 'Trypsin') else 'own')]
    else:
        got = PD[code]
    assert all(abs(a - b) < 1e-6 for a, b in zip(got, (est, lo, hi))), (r, got)

# ---------------------------------------------------------------------------------------------- draw
fig = plt.figure(figsize=(C.W_IN, C.mm(118)))
BLUE, PINK, GREY = C.AC[1], '#c06c98', '#6f757c'

# a
ax = fig.add_axes([0.075, 0.60, 0.44, 0.32])
ax.axhline(1, color=C.REF, lw=0.8, ls=(0, (3, 2)), zorder=0)
for o, (v, lo, hi) in PA.items():
    hl = o in (-1, 5, 7)
    col = C.INK if hl else '#8f969d'
    ax.plot([o, o], [lo, hi], color=col, lw=0.8)
    ax.plot(o, v, 'o', ms=3.4, color=col)
ax.annotate('−1: %.2f' % PA[-1][0], xy=(-1, PA[-1][0]), xytext=(-3.0, 0.48), fontsize=C.FS_S, ha='right', va='center',
            arrowprops=dict(arrowstyle='-', lw=0.5, color=C.MUTE, shrinkA=0, shrinkB=2.5))
ax.annotate('+5: %.2f' % PA[5][0], xy=(5, PA[5][0]), xytext=(3.0, 1.72), fontsize=C.FS_S, ha='center', va='bottom',
            arrowprops=dict(arrowstyle='-', lw=0.5, color=C.MUTE, shrinkA=0, shrinkB=2.5))
ax.annotate('+7: %.2f' % PA[7][0], xy=(7, PA[7][0]), xytext=(7.0, 1.72), fontsize=C.FS_S, ha='center', va='bottom',
            arrowprops=dict(arrowstyle='-', lw=0.5, color=C.MUTE, shrinkA=0, shrinkB=2.5))
ax.set_xlim(-7.6, 7.9)
ax.set_xticks([-7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7])
ax.set_xticklabels(['−7', '', '−5', '', '−3', '', '−1', '+1', '', '+3', '', '+5', '', '+7'])
ax.set_ylim(0.3, 1.95)
ax.set_yticks([0.5, 1.0, 1.5])
ax.set_xlabel('Offset from the modified cysteine (residues)')
ax.set_ylabel('K/R frequency / proteome')
C.panel_letter(fig, ax, 'a', 'K/R by offset, FAT-switch S-nitrosylation sites (Arabidopsis)')

# b
ax = fig.add_axes([0.715, 0.60, 0.27, 0.32])
ys = [3, 2, 1, -0.2]
for y, (lab, raw, mat) in zip(ys, PB):
    col = '#8f969d' if 'D/E' in lab else BLUE
    ax.plot([raw, mat], [y, y], color=col, lw=1.4, alpha=0.45, solid_capstyle='butt')
    ax.plot(raw, y, 'o', ms=4, mfc='white', mec=col, mew=0.9)
    ax.plot(mat, y, 'o', ms=4, mfc=col, mec=col)
    if abs(raw - mat) > 3:
        for v in (raw, mat):
            ax.text(max(v, 1.7), y + 0.22, '%+.2f' % v, fontsize=C.FS_S, va='bottom', ha='center', color=C.MUTE)
    else:
        for v in (raw, mat):
            right = (v >= max(raw, mat))
            ax.text(v + (0.6 if right else -0.6), y, '%+.2f' % v, fontsize=C.FS_S, va='center',
                    ha='left' if right else 'right', color=C.MUTE)
ax.axvline(0, color=C.AXIS, lw=0.8, zorder=0)
ax.axhline(0.4, color=C.REF, lw=0.6, ls=(0, (3, 2)))
ax.set_yticks(ys)
ax.set_yticklabels(['S-nitrosylation', 'S-sulfenylation', 'Rev. oxidation', 'Rev. oxidation'])
ax.text(19.5, 3.62, 'K/R/H', fontsize=C.FS_S, ha='right', va='bottom', color=BLUE)
ax.text(19.5, -0.2, 'D/E (control)', fontsize=C.FS_S, ha='right', va='center', color='#6f757c')
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.set_xlim(-2.5, 19.8)
ax.set_ylim(-0.75, 4.0)
ax.set_xticks([0, 5, 10, 15])
ax.set_xlabel('Enrichment z vs other cysteines')
C.panel_letter(fig, ax, 'b', 'Before (open) and after (filled) matching')

# c
LABEL = {'AspN': 'AspN', 'CT': 'Chymotrypsin', 'GluC': 'GluC', 'Trypsin': 'Trypsin'}
NN = {'AspN': NID['AspN'], 'CT': NID['chymotrypsin'], 'GluC': NID['GluC'], 'Trypsin': NID['trypsin']}
ax = fig.add_axes([0.20, 0.10, 0.36, 0.33])
yy = {'AspN': 3, 'CT': 2, 'GluC': 1, 'Trypsin': 0}
for code, y in yy.items():
    v = PC[(code, 'own')]
    ax.plot([v[1], v[2]], [y + 0.12, y + 0.12], color=BLUE, lw=0.9)
    ax.plot(v[0], y + 0.12, 'o', ms=4, color=BLUE)
    if (code, 'trypsin') in PC:
        t = PC[(code, 'trypsin')]
        ax.plot([t[1], t[2]], [y - 0.14, y - 0.14], color=GREY, lw=0.9)
        ax.plot(t[0], y - 0.14, 'o', ms=4, mfc='white', mec=GREY, mew=0.9)
ax.axvline(0, color=C.AXIS, lw=0.8, zorder=0)
ax.set_yticks(list(yy.values()))
ax.set_yticklabels(['%s (n = %s)' % (LABEL[c], format(NN[c], ',')) for c in yy])
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.set_ylim(-0.6, 3.6)
ax.set_xlim(-0.75, 2.65)
ax.set_xlabel('Distal-band (6–12) log$_2$ odds ratio')
ax.plot([1.55], [0.42], 'o', ms=4, color=BLUE, label='_nocheck')
ax.text(1.65, 0.42, "arm's own rule", fontsize=C.FS_S, va='center')
ax.plot([1.55], [0.02], 'o', ms=4, mfc='white', mec=GREY, mew=0.9, label='_nocheck')
ax.text(1.65, 0.02, 'trypsin rule', fontsize=C.FS_S, va='center')
C.panel_letter(fig, ax, 'c', 'Proteome background, PXD063463')

# d
ax = fig.add_axes([0.715, 0.10, 0.27, 0.33])
for code, y in yy.items():
    v = PD[code]
    ax.plot([v[1], v[2]], [y, y], color=PINK, lw=0.9)
    ax.plot(v[0], y, 'o', ms=4, color=PINK)
ax.axvline(0, color=C.AXIS, lw=0.8, zorder=0)
ax.set_yticks(list(yy.values()))
ax.set_yticklabels(['%s (n = %s)' % (LABEL[c], format(NN[c], ',')) for c in yy])
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.set_ylim(-0.6, 3.6)
ax.set_xlim(-1.5, 1.5)
ax.set_xlabel('Distal-band (6–12) log$_2$ odds ratio')
C.panel_letter(fig, ax, 'd', 'Observed background, own rule')
fig.text(0.985, 0.005, 'Intervals: a, 95%; c, d, 97.5% (Bonferroni, two bands)', fontsize=C.FS_S, ha='right',
         va='bottom', color=C.MUTE)
for f, v, s in (('interval_level_percent_a', 95, POS_LBL + ' (ci95_* columns)'),
                ('interval_level_percent_cd', 97.5, N4_LBL + ' (Cleavage geometry: ci_level 0.975)'),
                ('distal_band_from', 6, N4_LBL), ('distal_band_to', 12, N4_LBL),
                ('reference_line_a', 1, 'ratio 1 = proteome background (no data)')):
    sd.add('a-d', 'printed constants', f, v, s, 'printed label')

qa = C.check_figure(fig, 'Fig3', C.number_pool(sd), allow_numbers=('2',))
C.save(fig, 'Fig3_cleavage_geometry', qa)
print(sd.write())
