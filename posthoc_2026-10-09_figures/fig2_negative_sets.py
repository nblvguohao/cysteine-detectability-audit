"""Fig 2 (v3.3): the negative set decides what detectability-only features can separate (v3.2 Fig 5a redrawn).

Figure plan (2026-10-09), Fig 2: three paired panels (a global VIS10, b within-protein VIS10, c within-protein DIG25);
x = NEG_A / NEG_B; one line per cohort with 95% intervals at both ends; FAT-switch NEG_A only; YAP1C and QTRP S2
coloured and labelled directly; NEG_B end filled when comparable (comparable, comparable_threshold), open when
not_comparable (Source_Data_text_artefact4_public_cohorts.csv); corner note '6 of 7 cohorts fall by 0.13-0.28'
(Supplemental Note 20, Reading, bullet 2).
Deviation (reported): y axis 0.25-1.0 instead of 0.4-0.95, because intervals reach 0.2836 (QTRP S2, NEG_B, within
protein) and 0.9557 (QTRP S1, NEG_A, DIG25 within protein); a 0.4-0.95 axis would clip them.
Values are read from source_data_submitted/Source_Data_Fig5_three_axes.csv (panel a) and cross-checked against
Supplemental Note 20 Table S20.1 (3 decimals).
"""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import common as C

SRC = C.SD32 / 'Source_Data_Fig5_three_axes.csv'
SRC_COMP = C.SD32 / 'Source_Data_text_artefact4_public_cohorts.csv'
N20 = C.SUP / 'Supplemental_Note_20_fig5a_absolute_discrimination.md'
SRC_LBL = 'cysteine-detectability-audit/source_data_submitted/Source_Data_Fig5_three_axes.csv'
COMP_LBL = 'cysteine-detectability-audit/source_data_submitted/Source_Data_text_artefact4_public_cohorts.csv'
N20_LBL = 'Supplemental Note 20'

COHORTS = [('qtrp_S1_ph5', 'QTRP S1 pH5'), ('qtrp_S2_ph5', 'QTRP S2 pH5'), ('qpers_sid_tierB', 'qPerS-SID tier B'),
           ('cysboost2019_human_sno_hela', 'Cys-BOOST HeLa'), ('cysboost2019_human_sno_shsy5y', 'Cys-BOOST SH-SY5Y'),
           ('natcomm2023_ath_sno', 'FAT-switch'), ('abiotech2025_ath_sno', 'PAT-switch'),
           ('fps2020_ath_sulfenyl', 'YAP1C reporter')]
PANELS = [('a', 'auc_global_VIS10', 'All sites, VIS10', 'VIS10 global'),
          ('b', 'auc_within_protein_VIS10', 'Within proteins, VIS10', 'VIS10 within protein'),
          ('c', 'auc_within_protein_DIG25', 'Within proteins, DIG25', 'DIG25 within protein')]
NEUTRAL = '#8f969d'
HIGHLIGHT = {'qtrp_S2_ph5': '#1f6aa5', 'fps2020_ath_sulfenyl': '#dc9a2b'}
SHORT = {'qtrp_S2_ph5': 'QTRP S2', 'fps2020_ath_sulfenyl': 'YAP1C'}

rows = C.read_csv(SRC)
comp = {r['dataset_id']: r for r in C.read_csv(SRC_COMP)}
s201 = C.md_table(N20, '| Cohort | Neg. set |')
sd = C.SourceData('Fig2', 'Source_Data_Fig2_negative_sets.csv')

# ---------------------------------------------------------------------------------------------- read + cross-check
VAL = {}
for cid, name in COHORTS:
    for neg in ('NEG_A', 'NEG_B'):
        n20 = [r for r in s201 if r['Cohort'] == name and r['Neg. set'] == neg]
        assert len(n20) == 1, (name, neg)
        n20 = n20[0]
        for pnl, field, _, n20col in PANELS:
            key = '%s|%s' % (cid, neg)
            hits = [r for r in rows if r['panel'] == 'a' and r['row'] == key and r['field'] == field]
            if cid == 'natcomm2023_ath_sno' and neg == 'NEG_B':
                assert not hits or hits[0]['value'] in ('', 'nan', 'NaN'), hits
                assert 'no observed-unmodified class' in n20['sites / negatives']
                continue
            est = float(C.sd_lookup(rows, 'a', key, field)['value'])
            lo = float(C.sd_lookup(rows, 'a', key, field + '_lo')['value'])
            hi = float(C.sd_lookup(rows, 'a', key, field + '_hi')['value'])
            p, l, h = C.pt_ci(n20[n20col])
            assert all(abs(x - y) <= 0.00051 for x, y in zip((est, lo, hi), (p, l, h))), (name, neg, field, est, lo, hi, p, l, h)
            VAL[(cid, neg, field)] = (est, lo, hi)
            spec = 'plotted; %s; %s; point and 95%% cluster-bootstrap interval (cross-checked with Note 20 Table S20.1)' % (name, neg)
            for f, v in ((field, est), (field + '_lo', lo), (field + '_hi', hi)):
                sd.add(pnl, key, f, v, '%s (panel a, row %s, field %s)' % (SRC_LBL, key, f), spec)

COMP = {}
for cid, name in COHORTS:
    c = comp[cid]['comparability']
    COMP[cid] = c
    sd.add('a-c', cid, 'NEG_B_comparability', c, '%s (dataset_id %s, comparability)' % (COMP_LBL, cid),
           'NEG_B end filled if comparable or comparable_threshold, open if not_comparable, absent if not_computable')
    sd.add('a-c', cid, 'display_name', name, '%s, Table S20.1 (Cohort)' % N20_LBL, 'cohort display name')
assert {k for k, v in COMP.items() if v in ('comparable', 'comparable_threshold')} == {
    'qtrp_S1_ph5', 'qtrp_S2_ph5', 'qpers_sid_tierB', 'fps2020_ath_sulfenyl'}
assert {k for k, v in COMP.items() if v == 'not_comparable'} == {
    'abiotech2025_ath_sno', 'cysboost2019_human_sno_hela', 'cysboost2019_human_sno_shsy5y'}

# corner note: quoted from Note 20 (Reading, bullet 2) - checked to be present verbatim
n20txt = N20.read_text(encoding='utf-8')
assert 'VIS10 global AUC fell by 0.13–0.28 in six cohorts' in n20txt
for f, v in (('cohorts_falling', 6), ('cohorts_with_NEG_B', 7), ('fall_min', 0.13), ('fall_max', 0.28)):
    sd.add('a', 'corner note', f, v, '%s (Reading, bullet 2: "VIS10 global AUC fell by 0.13-0.28 in six cohorts")' % N20_LBL,
           'printed: 6 of 7 cohorts fall by 0.13-0.28')
# consistency of the quoted note with the plotted values (no new statistic is reported; assertion only)
drops = [VAL[(c, 'NEG_A', 'auc_global_VIS10')][0] - VAL[(c, 'NEG_B', 'auc_global_VIS10')][0]
         for c, _ in COHORTS if (c, 'NEG_B', 'auc_global_VIS10') in VAL]
falls = sorted(d for d in drops if d > 0)
assert len(drops) == 7 and len(falls) == 6 and round(falls[0], 2) == 0.13 and round(falls[-1], 2) == 0.28, drops
sd.add('a-c', 'negative-set definitions', 'NEG_A', C.NEG_A_DEF, 'Experimental Procedures (Detectability-only models: negative-set definitions)', 'printed key')
sd.add('a-c', 'negative-set definitions', 'NEG_B', C.NEG_B_DEF, 'Experimental Procedures (Detectability-only models: negative-set definitions)', 'printed key')

# ---------------------------------------------------------------------------------------------- draw
fig = plt.figure(figsize=(C.W_IN, C.mm(80)))
L, R, B, T, GAP = 0.075, 0.985, 0.30, 0.90, 0.055
pw = (R - L - 2 * GAP) / 3
DX = dict(zip(['natcomm2023_ath_sno'] + [c for c, _ in COHORTS if c != 'natcomm2023_ath_sno'], [-0.15, -0.107, -0.064, -0.021, 0.021, 0.064, 0.107, 0.15]))
axes = []
for i, (pnl, field, title, _) in enumerate(PANELS):
    ax = fig.add_axes([L + i * (pw + GAP), B, pw, T - B])
    axes.append(ax)
    ax.axhline(0.5, color=C.REF, lw=0.8, ls=(0, (3, 2)), zorder=0)
    for cid, name in COHORTS:
        col = HIGHLIGHT.get(cid, NEUTRAL)
        z = 3 if cid in HIGHLIGHT else 2
        a = VAL[(cid, 'NEG_A', field)]
        xa = 0 + DX[cid]
        ax.plot([xa, xa], [a[1], a[2]], color=col, lw=0.7, zorder=z)
        ax.plot(xa, a[0], 'o', ms=3.4, mfc=col, mec=col, mew=0.8, zorder=z + 1)
        if (cid, 'NEG_B', field) in VAL:
            b = VAL[(cid, 'NEG_B', field)]
            xb = 1 + DX[cid]
            ax.plot([xa, xb], [a[0], b[0]], color=col, lw=0.9 if cid in HIGHLIGHT else 0.7, zorder=z)
            ax.plot([xb, xb], [b[1], b[2]], color=col, lw=0.7, zorder=z)
            filled = COMP[cid] in ('comparable', 'comparable_threshold')
            ax.plot(xb, b[0], 'o', ms=3.4, mfc=col if filled else 'white', mec=col, mew=0.8, zorder=z + 1)
            if cid in SHORT:
                ax.text(1.24, b[0], SHORT[cid], color=col, fontsize=C.FS_S, va='center', ha='left')
        else:   # FAT-switch: NEG_A only
            ax.plot(xa, a[0], 'D', ms=3.2, mfc=NEUTRAL, mec=NEUTRAL, zorder=z + 1)
            ax.text(-0.21, a[0], 'FAT-switch', fontsize=C.FS_S, color=C.MUTE, ha='right', va='center')
    ax.set_xlim(-0.75, 1.62)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['NEG_A', 'NEG_B'])
    ax.set_ylim(0.25, 1.0)
    ax.set_yticks([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.spines['bottom'].set_bounds(-0.3, 1.3)
    if i == 0:
        ax.set_ylabel('Out-of-fold AUC')
    C.panel_letter(fig, ax, pnl, title)
axes[0].text(1.60, 0.27, '6 of 7 cohorts fall\nby 0.13–0.28', fontsize=C.FS_S, color=C.INK, ha='right', va='bottom')

# key under the panels
kax = fig.add_axes([L, 0.02, R - L, 0.17])
kax.axis('off')
h = [Line2D([], [], color=NEUTRAL, marker='o', ms=3.4, lw=0.7, mfc=NEUTRAL, mec=NEUTRAL),
     Line2D([], [], color=NEUTRAL, marker='o', ms=3.4, lw=0, mfc='white', mec=NEUTRAL),
     Line2D([], [], color=NEUTRAL, marker='D', ms=3.2, lw=0, mfc=NEUTRAL, mec=NEUTRAL)]
kax.legend(h, ['one cohort; NEG_B end filled: same-experiment readout (QTRP S1, S2, qPerS-SID, YAP1C)',
               'NEG_B end open: observed class not comparable (PAT-switch, Cys-BOOST HeLa, SH-SY5Y)',
               'FAT-switch: NEG_A only (no observed-unmodified class)'], loc='upper left', bbox_to_anchor=(0, 1.08), fontsize=C.FS_S, ncol=1,
           handlelength=1.6, handletextpad=0.5, labelspacing=0.25, borderaxespad=0)
kax.text(0.62, 0.86, 'NEG_A: ' + C.NEG_A_DEF, fontsize=C.FS_S, va='center', transform=kax.transAxes)
kax.text(0.62, 0.62, 'NEG_B: ' + C.NEG_B_DEF, fontsize=C.FS_S, va='center', transform=kax.transAxes)
kax.text(0.62, 0.38, 'Bars: 95% cluster-bootstrap intervals; dashed line: AUC 0.5', fontsize=C.FS_S,
         va='center', transform=kax.transAxes)

qa = C.check_figure(fig, 'Fig2', C.number_pool(sd), allow_numbers=('0.5', '95', '10', '25'))
sd.add('a-c', 'printed constants', 'reference_line', 0.5, 'chance level (no data)', 'dashed line')
sd.add('a-c', 'printed constants', 'interval_level_percent', 95,
       '%s (panel a, specification_label)' % SRC_LBL, 'printed in key')
C.save(fig, 'Fig2_negative_sets', qa)
print(sd.write())
