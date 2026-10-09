"""Fig 4 (v3.3): two steps remove the comparison - coincident positive and identified sets, undecidable search spaces.

Figure plan (2026-10-09), Fig 4:
  a  share of identified cysteines carrying a site, point + 95% interval, sorted by value: PXD048216 (Note 7, PRIMARY,
     CAM_only, Wilson), PXD063463 four arms (v3.2 Source Data Fig 2 panel c, protein-clustered bootstrap, Note 4),
     the four published cohorts with used_in_main_text == yes (Source_Data_text_artefact4_public_cohorts.csv, Wilson);
     no 0.50 / 0.90 decision lines;
  b  v3.2 Fig 3a without the per-file precision strip: 0.017758-Da separation in ppm of M, +/-10, 4.5, 2 ppm windows;
  c  v3.2 Fig 3d: share of single-cysteine spectra tied at the best score by fragment binning.
"""
import matplotlib.pyplot as plt

import common as C

N7 = C.SUP / 'Supplemental_Note_7_artefact4_public_deposit.md'
N7_LBL = 'Supplemental Note 7'
FP = C.SD32 / 'Source_Data_Fig2_public_four_protease.csv'
FP_LBL = 'cysteine-detectability-audit/source_data_submitted/Source_Data_Fig2_public_four_protease.csv'
COH = C.SD32 / 'Source_Data_text_artefact4_public_cohorts.csv'
COH_LBL = 'cysteine-detectability-audit/source_data_submitted/Source_Data_text_artefact4_public_cohorts.csv'
SS = C.SD32 / 'Source_Data_Fig3_search_space.csv'
SS_LBL = 'cysteine-detectability-audit/source_data_submitted/Source_Data_Fig3_search_space.csv'
N4 = C.SUP / 'Supplemental_Note_4_four_protease_public_deposit.md'

sd = C.SourceData('Fig4', 'Source_Data_Fig4_removed_comparisons.csv')
PINK, BLUE, GREY, VERMIL = C.AC[4], C.AC[1], '#4a4f55', C.AC[5]

# ---------------------------------------------------------------------------------------------- a
A = []   # (label, est, lo, hi, group)
n7 = [r for r in C.md_table(N7, '| definition | reading |')
      if r['definition'] == 'PRIMARY_leading_both_sides' and r['reading'] == 'CAM_only']
assert len(n7) == 1
r = n7[0]
est, lo, hi = float(r['p']), float(r['wilson_low']), float(r['wilson_high'])
assert (round(est, 4), round(lo, 4), round(hi, 4)) == (0.9418, 0.9336, 0.949)
A.append(('PXD048216 (ABE)', est, lo, hi, 'deposit'))
for f, v in (('share', est), ('wilson_low', lo), ('wilson_high', hi), ('identified', int(r['observed_cysteines'])),
             ('with_site', int(r['matched']))):
    sd.add('a', 'PXD048216', f, v, '%s (coincidence table, PRIMARY_leading_both_sides, CAM_only)' % N7_LBL,
           'plotted (share, Wilson 95% interval)' if f in ('share', 'wilson_low', 'wilson_high') else 'not plotted; count')
fp = C.read_csv(FP)
for code, disp in (('Trypsin', 'trypsin'), ('CT', 'chymotrypsin'), ('GluC', 'GluC'), ('AspN', 'AspN')):
    v = [float(C.sd_lookup(fp, 'c', code, f)['value']) for f in ('estimate', 'ci_low', 'ci_high')]
    A.append(('PXD063463 %s' % disp, v[0], v[1], v[2], 'four'))
    for f, x in zip(('estimate', 'ci_low', 'ci_high'), v):
        sd.add('a', 'PXD063463 %s' % disp, f, x, '%s (panel c, row %s, field %s)' % (FP_LBL, code, f),
               'plotted; protein-clustered 95% interval (Supplemental Note 4, Coincidence)')
assert '**Coincidence.** The share of identified cysteines that carry a site, with a protein-clustered 95 per cent interval' \
    in N4.read_text(encoding='utf-8')
coh = [r for r in C.read_csv(COH) if r['used_in_main_text'] == 'yes']
assert {r['dataset_id'] for r in coh} == {'fps2020_ath_sulfenyl', 'qpers_sid_tierB', 'qtrp_S1_ph5', 'qtrp_S2_ph5'}
NAME = {'fps2020_ath_sulfenyl': 'YAP1C reporter', 'qpers_sid_tierB': 'qPerS-SID tier B', 'qtrp_S1_ph5': 'QTRP S1 pH5',
        'qtrp_S2_ph5': 'QTRP S2 pH5'}
for r in coh:
    v = float(r['share']), float(r['wilson_low']), float(r['wilson_high'])
    A.append((NAME[r['dataset_id']], v[0], v[1], v[2], 'published'))
    for f, x in zip(('share', 'wilson_low', 'wilson_high'), v):
        sd.add('a', r['dataset_id'], f, x, '%s (dataset_id %s, %s)' % (COH_LBL, r['dataset_id'], f),
               'plotted; %s; Wilson 95%% interval; used_in_main_text = yes' % NAME[r['dataset_id']])
A.sort(key=lambda t: t[1])

# ---------------------------------------------------------------------------------------------- b
ss = C.read_csv(SS)
curve = sorted((float(C.sd_lookup(ss, 'a', r['row'], 'mass_da')['value']),
                float(C.sd_lookup(ss, 'a', r['row'], 'separation_ppm')['value']), r['row'])
               for r in ss if r['panel'] == 'a' and r['field'] == 'mass_da' and r['row'].startswith('separation curve'))
for m, p, row in curve:
    sd.add('b', row, 'mass_da', m, '%s (panel a, row %s)' % (SS_LBL, row), 'drawn curve grid point (reused unchanged from the original figure build)')
    sd.add('b', row, 'separation_ppm', p, '%s (panel a, row %s)' % (SS_LBL, row), 'drawn curve grid point (reused unchanged from the original figure build)')
SEP = C.sd_lookup(ss, 'a', 'separation label', 'separation_da')['value']
sd.add('b', 'separation label', 'separation_da', float(SEP), '%s (panel a, separation label)' % SS_LBL, 'printed')
WIN = []
for t in ('10.0', '4.5', '2.0'):
    row = 'window +/-%s ppm' % t
    tol = float(C.sd_lookup(ss, 'a', row, 'tolerance_ppm')['value'])
    ms = float(C.sd_lookup(ss, 'a', row, 'mass_limit_da')['value'])
    fr = float(C.sd_lookup(ss, 'a', row, 'fraction_theoretical_cys_peptides_at_or_below_plus32_only')['value'])
    WIN.append((tol, ms, fr))
    for f, v in (('tolerance_ppm', tol), ('mass_limit_da', ms),
                 ('fraction_theoretical_cys_peptides_at_or_below_plus32_only', fr)):
        sd.add('b', row, f, v, '%s (panel a, row %s, field %s)' % (SS_LBL, row, f), 'plotted / printed')
    sd.add('b', row, 'mass_limit_da_printed', round(ms), SS_LBL, 'mass_limit_da rounded to 1 Da for printing')
    pct = 99.9 if 0.999 < fr < 1.0 else round(fr * 100)
    sd.add('b', row, 'percent_printed', pct, SS_LBL,
           'fraction x 100 rounded for printing (">99.9%" printed when 0.999 < fraction < 1)')

# ---------------------------------------------------------------------------------------------- c
BINS = ['(0.0, 1.0]', '(1.0, 10.0]', '(10.0, 100.0]', '(100.0, 1000.1]', 'all']
BLAB = ['≤1', '1–10', '10–100', '>100']
GROUPS = [('insulin files (ion-trap MS2, 1.0005-Da bins)', VERMIL, 's', '1.0005-Da bins (ion trap)'),
          ('Fig-1D (FTMS MS2, 0.02-Da bins)', C.INK, 'o', '0.02-Da bins (FTMS)')]
TIE = {}
for g, *_ in GROUPS:
    for b in BINS:
        row = '%s|best E-value %s' % (g, b)
        vals = {f: float(C.sd_lookup(ss, 'd', row, f)['value']) for f in ('share_tied', 'n_tied', 'n_both_scored')}
        TIE[(g, b)] = vals
        for f, v in vals.items():
            sd.add('c', row, f, v, '%s (panel d, row %s, field %s)' % (SS_LBL, row, f),
                   'plotted (share) / printed under the axis (tied/n)')
    sd.add('c', '%s|best E-value all' % g, 'percent_printed', round(TIE[(g, 'all')]['share_tied'] * 100), SS_LBL,
           'share_tied x 100 rounded for printing')

# ---------------------------------------------------------------------------------------------- draw
fig = plt.figure(figsize=(C.W_IN, C.mm(78)))
ax = fig.add_axes([0.18, 0.17, 0.185, 0.71])
for i, (lab, e, lo, hi, grp) in enumerate(A):
    col = {'deposit': PINK, 'four': BLUE, 'published': GREY}[grp]
    ax.plot([lo, hi], [i, i], color=col, lw=0.9)
    ax.plot(e, i, 'o', ms=3.6, color=col)
ax.set_yticks(range(len(A)))
ax.set_yticklabels([t[0] for t in A])
for tl, t in zip(ax.get_yticklabels(), A):
    tl.set_color({'deposit': PINK, 'four': BLUE, 'published': GREY}[t[4]])
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.set_xlim(0, 1)
ax.set_ylim(-0.6, len(A) - 0.4)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_xticklabels(['0', '0.25', '0.50', '0.75', '1'])
ax.set_xlabel('Identified cysteines carrying a site')
C.panel_letter(fig, ax, 'a', 'Positive set vs identified set')

# b
ax = fig.add_axes([0.435, 0.17, 0.235, 0.71])
mx, sp = [c[0] for c in curve], [c[1] for c in curve]
YMAX = 16.0
ax.fill_between(mx, sp, YMAX, color=VERMIL, alpha=0.08, lw=0, zorder=0)
ax.plot(mx, sp, color=VERMIL, lw=1.3, zorder=3, label='_nocheck')
ax.text(1700, 14.3, '%s Da in ppm of M' % SEP, color=VERMIL, fontsize=C.FS_S, ha='left', va='center')
ax.text(13900, 12.3, 'both candidates\ninside the window', color=C.MUTE, fontsize=C.FS_S, ha='right', va='center')
for tol, ms, fr in WIN:
    ax.axhline(tol, color='#8f969d', lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.plot([ms, ms], [0, tol], color='#8f969d', lw=0.6, ls=(0, (1, 1.6)), zorder=1, label='_nocheck')
    ax.plot([ms], [tol], 'o', ms=3.2, color=C.INK, zorder=5)
    pct = '>99.9%' if 0.999 < fr < 1.0 else '%d%%' % round(fr * 100)
    lab = '±%s ppm: %s Da (%s)' % (('%g' % tol), format(round(ms), ','), pct)
    if tol == 2.0:
        ax.text(13900, 2.3, lab.replace(' Da (', ' Da\n('), color=C.INK, fontsize=C.FS_S, ha='right', va='bottom',
                multialignment='right')
    else:
        ax.text(ms + 250, tol + 0.35, lab, color=C.INK, fontsize=C.FS_S, ha='left', va='bottom')
ax.text(13900, 7.4, '(%): theoretical Cys\npeptides ≤ that mass', color=C.MUTE, fontsize=C.FS_S, ha='right', va='center')
ax.set_xlim(0, 14000)
ax.set_ylim(0, YMAX)
ax.set_yticks([0, 5, 10, 15])
ax.set_xticks([0, 4000, 8000, 12000])
ax.set_xticklabels(['0', '4,000', '8,000', '12,000'])
ax.set_xlabel('Peptide neutral mass (Da)')
ax.set_ylabel('Separation (ppm)')
C.panel_letter(fig, ax, 'b', 'Sulfide vs dioxidation, precursor window')

# c
ax = fig.add_axes([0.755, 0.33, 0.225, 0.55])
XALL = 4.8
for gi, (g, col, mk, lab) in enumerate(GROUPS):
    ys = [TIE[(g, b)]['share_tied'] for b in BINS[:4]]
    ax.plot(range(4), ys, color=col, lw=0.9, zorder=2)
    ax.plot(range(4), ys, mk, color=col, ms=3.4, mew=0, zorder=3)
    ya = TIE[(g, 'all')]['share_tied']
    ax.plot([XALL], [ya], mk, color=col, ms=4.6, mew=0, zorder=3)
    ax.text(XALL + 0.3, ya, '%d%%' % round(ya * 100), color=col, fontsize=C.FS_S, va='center')
    for xi, b in enumerate(BINS):
        v = TIE[(g, b)]
        ax.annotate('%d/%d' % (v['n_tied'], v['n_both_scored']), xy=(XALL if b == 'all' else xi, 0),
                    xycoords=ax.get_xaxis_transform(), xytext=(0, -12.5 - 7.0 * gi), textcoords='offset points',
                    color=col, fontsize=C.FS_S, ha='center', va='top')
ax.annotate('tied/n', xy=(-0.45, 0), xycoords=ax.get_xaxis_transform(), xytext=(-3, -12.5),
            textcoords='offset points', color=C.MUTE, fontsize=C.FS_S, ha='right', va='top')
ax.text(-0.3, 0.68, GROUPS[0][3], color=VERMIL, fontsize=C.FS_S, ha='left', va='center')
ax.text(0.6, 0.06, GROUPS[1][3], color=C.INK, fontsize=C.FS_S, ha='left', va='center')
ax.axvline(4.2, color=C.REF, lw=0.6, zorder=0)
ax.set_xlim(-0.45, XALL + 0.6)
ax.set_ylim(0, 1.05)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_yticklabels(['0', '0.25', '0.50', '0.75', '1'])
ax.set_xticks([0, 1, 2, 3, XALL])
ax.set_xticklabels(BLAB + ['all'], fontsize=C.FS_S)
ax.set_ylabel('Share tied at best score')
ax.text(0.5, -0.42, 'Best E-value of the spectrum', transform=ax.transAxes, ha='center', va='top', fontsize=C.FS)
C.panel_letter(fig, ax, 'c', 'Fragment binning')

# key for panel a
fig.text(0.012, 0.035, 'a: ', fontsize=C.FS_S, fontweight='bold', va='bottom')
fig.text(0.03, 0.035, 'PXD048216 and published cohorts, Wilson 95% interval; PXD063463 arms, protein-clustered '
         'bootstrap 95% interval', fontsize=C.FS_S, va='bottom', color=C.MUTE)
sd.add('a', 'printed constants', 'interval_level_percent', 95, N7_LBL + '; ' + COH_LBL, 'printed in key')
sd.add('b', 'printed constants', 'ppm_label_tolerances', '10; 4.5; 2', SS_LBL, 'printed window labels')
for g, *_ in GROUPS:
    pass
for f, v in (('fragment_bin_ion_trap_da', 1.0005), ('fragment_bin_ftms_da', 0.02)):
    sd.add('c', 'printed constants', f, v, '%s (panel d, row names)' % SS_LBL, 'printed group label')
for f, v in (('evalue_bin_edge_1', 1), ('evalue_bin_edge_10', 10), ('evalue_bin_edge_100', 100)):
    sd.add('c', 'printed constants', f, v, '%s (panel d, row names: best E-value bins)' % SS_LBL, 'x tick label')

qa = C.check_figure(fig, 'Fig4', C.number_pool(sd), allow_numbers=('048216', '063463', '10', '4.5', '2', '1'))
C.save(fig, 'Fig4_removed_comparisons', qa)
print(sd.write())
