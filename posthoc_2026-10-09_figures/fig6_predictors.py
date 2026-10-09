"""Fig 6 (v3.3): detectability-only features against a published predictor and our own ranker.

Figure plan (2026-10-09), Fig 6:
  a  pLMSNOSite benchmark, two small panels (AUROC over all test sites; within-protein AUROC, mean of per-protein
     AUROCs): pLMSNOSite (probabilities), VIS10, DIG25 with 95% intervals; the all-sites panel adds an open marker for
     the tabulated 0.754 (thresholded calls); paired differences listed under the panels:
     all sites VIS10 0.106 [0.051, 0.153] resolved, DIG25 0.039 [-0.021, 0.105] unresolved; within protein VIS10
     0.013 [-0.041, 0.065] unresolved. Values: posthoc_2026-09-30/results/G_plmsnosite_rescore/ (metrics_main.csv,
     paired_differences.csv, withinprot_mean_summary.csv), cross-checked with Note 22 Table S22.1.
  b  own ranker, human and rice: within-protein AUC of the detectability-only set (25 columns) and of the full model
     (1,107 columns), with the recovered share and its interval (v3.2 Source Data Fig 7 panel c; Note 12).
"""
import matplotlib.pyplot as plt

import common as C

G = C.REPO / 'posthoc_2026-09-30/results/G_plmsnosite_rescore'
G_LBL = 'cysteine-detectability-audit/posthoc_2026-09-30/results/G_plmsnosite_rescore/'
N22 = C.SUP / 'Supplemental_Note_22_plmsnosite_rescored.md'
N22_LBL = 'Supplemental Note 22'
F7 = C.SD32 / 'Source_Data_Fig7_self_audit_public.csv'
F7_LBL = 'cysteine-detectability-audit/source_data_submitted/Source_Data_Fig7_self_audit_public.csv'
N12 = C.SUP / 'Supplemental_Note_12_public_detectability_only.md'
N12_LBL = 'Supplemental Note 12'

sd = C.SourceData('Fig6', 'Source_Data_Fig6_predictors.csv')
t22 = {r['Model or comparison']: r for r in C.md_table(N22, '| Model or comparison |')}


def chk(v, cell):
    ref = C.pt_ci(cell)
    assert all(abs(a - b) <= 0.00051 for a, b in zip(v, ref)), (v, ref)


# ---------------------------------------------------------------------------------------------- a: all sites
mm_ = {r['score']: r for r in C.read_csv(G / 'metrics_main.csv')}
ALL = {}
for key, row22 in (('pLMSNOSite', 'pLMSNOSite, released models (probabilities)'), ('VIS10', 'VIS10'),
                   ('DIG25', 'DIG25'), ('pLMSNOSite_call_0.5', 'pLMSNOSite, calls thresholded at 0.5ᶜ')):
    r = mm_[key]
    v = (float(r['auroc']), float(r['auroc_ci_low']), float(r['auroc_ci_high']))
    chk(v, t22[row22]['AUROC, all test sites'].replace('ᵈ', ''))
    ALL[key] = v
    for f, x in zip(('auroc', 'auroc_ci_low', 'auroc_ci_high'), v):
        sd.add('a', key, 'all_sites_' + f, x, G_LBL + 'metrics_main.csv (score %s, %s); = %s, Table S22.1' % (key, f, N22_LBL),
               'plotted; 95% percentile interval, 5,000 protein resamples (seed 20260930)'
               + ('; open marker: tabulated value of thresholded calls (post hoc re-scoring)' if 'call' in key else ''))
# ---------------------------------------------------------------------------------------------- a: within protein
wp = C.read_csv(G / 'withinprot_mean_summary.csv')
PRIM = 'all 278 test proteins resampled, seed 20260930, resamples shared with step 30 (primary)'
EST = 'mean of per-protein AUROCs (manuscript estimand; primary)'


def wrow(quantity, model):
    h = [r for r in wp if r['estimand'] == EST and r['scheme'] == PRIM and r['quantity'] == quantity and r['model'] == model]
    assert len(h) == 1, (quantity, model)
    return h[0]


WIN = {}
for key, row22 in (('pLMSNOSite', 'pLMSNOSite, released models (probabilities)'), ('VIS10', 'VIS10'), ('DIG25', 'DIG25')):
    r = wrow('within-protein AUROC', key)
    v = (float(r['value']), float(r['ci_low']), float(r['ci_high']))
    chk(v, t22[row22]['Within-protein AUROC, mean of per-protein AUROCsᵃ'])
    WIN[key] = v
    for f, x in zip(('value', 'ci_low', 'ci_high'), v):
        sd.add('a', key, 'within_protein_auroc_' + f, x,
               G_LBL + 'withinprot_mean_summary.csv (estimand mean of per-protein AUROCs, primary scheme, model %s, %s); = %s, Table S22.1'
               % (key, f, N22_LBL), 'plotted; 220 test proteins with both labels; 95% percentile interval')
# ---------------------------------------------------------------------------------------------- a: paired differences
pdiff = C.read_csv(G / 'paired_differences.csv')
DIFF = []
for comp, row22, col in (('pLMSNOSite - VIS10', 'pLMSNOSite − VIS10, paired difference', 'AUROC, all test sites'),
                         ('pLMSNOSite - DIG25', 'pLMSNOSite − DIG25, paired difference', 'AUROC, all test sites')):
    h = [r for r in pdiff if r['comparison'] == comp and r['metric'] == 'auroc']
    assert len(h) == 1
    v = (float(h[0]['difference']), float(h[0]['ci_low']), float(h[0]['ci_high']))
    chk(v, t22[row22][col])
    DIFF.append(('All sites', comp, v))
r = wrow('within-protein AUROC difference', 'pLMSNOSite - VIS10')
v = (float(r['value']), float(r['ci_low']), float(r['ci_high']))
chk(v, t22['pLMSNOSite − VIS10, paired difference']['Within-protein AUROC, mean of per-protein AUROCsᵃ'])
DIFF.append(('Within proteins', 'pLMSNOSite - VIS10', v))
for scope, comp, v in DIFF:
    res = 'resolved' if (v[1] > 0 or v[2] < 0) else 'unresolved'
    src = (G_LBL + 'paired_differences.csv (comparison %s, metric auroc)' % comp if scope == 'All sites' else
           G_LBL + 'withinprot_mean_summary.csv (within-protein AUROC difference, %s, primary)' % comp)
    for f, x in zip(('difference', 'ci_low', 'ci_high'), v):
        sd.add('a', '%s | %s' % (scope, comp), f, x, src + '; = %s, Table S22.1' % N22_LBL,
               'printed under the panels (3 decimals); resolved = interval excludes 0')
    sd.add('a', '%s | %s' % (scope, comp), 'reading', res, N22_LBL + ' (Statistics: resolved when the interval excludes zero)',
           'printed')
assert [d[2][0] > 0 and d[2][1] > 0 for d in DIFF] == [True, False, False]

# ---------------------------------------------------------------------------------------------- b
f7 = C.read_csv(F7)
n12 = C.md_table(N12, '| cohort |')
B = {}
for coh, disp in (('human_PXD044043', 'Human'), ('rice_PXD072089', 'Rice')):
    det = float(C.sd_lookup(f7, 'c', coh + '|detect_only', 'within_protein_auc')['value'])
    full = float(C.sd_lookup(f7, 'c', coh + '|full', 'within_protein_auc')['value'])
    row = [r for r in n12 if r['cohort'] == coh]
    assert len(row) == 1
    vals = list(row[0].values())
    B[coh] = (disp, det, full)
    sd.add('b', coh, 'within_protein_auc_detect_only', det, '%s (panel c, row %s|detect_only)' % (F7_LBL, coh),
           'plotted; detectability-only set, 25 columns')
    sd.add('b', coh, 'within_protein_auc_full', full, '%s (panel c, row %s|full)' % (F7_LBL, coh),
           'plotted; full model, 1,107 columns')
    B[coh] += (row[0],)
hdr = list(n12[0].keys())
print('Note 12 header:', hdr)
for coh, (disp, det, full, row) in B.items():
    cells = list(row.values())
    # columns: cohort | n | full | detect-only | ratio | lo | hi  (asserted against Fig 7 Source Data)
    assert abs(float(cells[2]) - full) < 1e-9 and abs(float(cells[3]) - det) < 1e-9, cells
    ratio, lo, hi = float(cells[4]), float(cells[5]), float(cells[6])
    for f, col, x in (('recovery_share', 'recovery_ratio', ratio), ('recovery_share_ci_low', 'ci_low', lo),
                      ('recovery_share_ci_high', 'ci_high', hi)):
        sd.add('b', coh, f, x, '%s (LightGBM ranking member (primary) table, row %s, column %s)' % (N12_LBL, coh, col),
               'printed (2 decimals); paired protein-clustered bootstrap interval')
    B[coh] = (disp, det, full, ratio, lo, hi)
for f, v in (('n_columns_detect_only', 25), ('n_columns_full', 1107)):
    sd.add('b', 'model sizes', f, v, '%s (first paragraph)' % N12_LBL, 'printed in key')

# ---------------------------------------------------------------------------------------------- draw
fig = plt.figure(figsize=(C.W_IN, C.mm(78)))
MODELS = [('pLMSNOSite', 'pLMSNOSite', '#1f2326'), ('VIS10', 'VIS10 (10 features)', '#1f6aa5'),
          ('DIG25', 'DIG25 (25 features)', '#6fa9d6')]
YS = [2, 1, 0]
axs = [fig.add_axes([0.155, 0.47, 0.2, 0.41]), fig.add_axes([0.40, 0.47, 0.2, 0.41])]
for ax, data, title in ((axs[0], ALL, 'All test sites'), (axs[1], WIN, 'Within proteins')):
    for y, (k, lab, col) in zip(YS, MODELS):
        v = data[k]
        ax.plot([v[1], v[2]], [y, y], color=col, lw=1.0)
        ax.plot(v[0], y, 'o', ms=4, color=col)
    ax.axvline(0.5, color=C.REF, lw=0.6, ls=(0, (3, 2)), zorder=0)
    ax.set_xlim(0.45, 0.9)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9])
    ax.set_ylim(-0.6, 2.95)
    ax.set_yticks(YS)
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('AUROC')
    ax.set_title(title, fontsize=C.FS, loc='center', pad=2)
axs[0].set_yticklabels([m[1] for m in MODELS])
axs[1].set_yticklabels([])
v = ALL['pLMSNOSite_call_0.5']
axs[0].plot([v[1], v[2]], [2.5, 2.5], color='#8f969d', lw=0.7)
axs[0].plot(v[0], 2.5, 'o', ms=4, mfc='white', mec='#1f2326', mew=0.8)
axs[0].text(0.455, 2.5, 'tabulated %.3f\n(thresholded calls)' % v[0], fontsize=C.FS_S, ha='left', va='center',
            color=C.MUTE)
C.panel_letter(fig, (0.012, 0.95), 'a', 'Published predictor (pLMSNOSite benchmark, test set)')
lines = ['Paired difference, pLMSNOSite minus:']
for scope, comp, v in DIFF:
    res = 'resolved' if (v[1] > 0 or v[2] < 0) else 'unresolved'
    lines.append('%s, %s: %.3f [%.3f, %.3f], %s' % (scope.lower(), comp.split(' - ')[1], v[0], v[1], v[2], res))
for i, ln in enumerate(lines):
    fig.text(0.155, 0.285 - i * 0.05, ln.replace('-', '−'), fontsize=C.FS_S,
             color=C.INK if i else C.MUTE, va='top')
fig.text(0.155, 0.285 - 4 * 0.05, 'Within proteins: mean of per-protein AUROCs; 95% intervals', fontsize=C.FS_S,
         color=C.MUTE, va='top')

# b
ax = fig.add_axes([0.72, 0.22, 0.26, 0.66])
for x, coh in ((0, 'human_PXD044043'), (1, 'rice_PXD072089')):
    disp, det, full, ratio, lo, hi = B[coh]
    ax.plot([x - 0.12, x + 0.12], [det, full], color='#b4b9bf', lw=0.8, zorder=1)
    ax.plot(x - 0.12, det, 'o', ms=4.4, mfc='white', mec='#1f6aa5', mew=0.9, zorder=3)
    ax.plot(x + 0.12, full, 'o', ms=4.4, color='#1f2326', zorder=3)
    ax.text(x, 0.66, 'share\n%.2f [%.2f, %.2f]' % (ratio, lo, hi), fontsize=C.FS_S, ha='center', va='center',
            color=C.INK)
ax.axhline(0.5, color=C.REF, lw=0.6, ls=(0, (3, 2)), zorder=0)
ax.set_xlim(-0.55, 1.55)
ax.set_ylim(0.5, 0.95)
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9])
ax.set_xticks([0, 1])
ax.set_xticklabels(['Human\nPXD044043', 'Rice\nPXD072089'])
ax.set_ylabel('Within-protein AUC')
C.panel_letter(fig, (0.64, 0.95), 'b', 'Our ranker, refit within each cohort')
kb = fig.add_axes([0.72, 0.0, 0.27, 0.08])
kb.axis('off')
from matplotlib.lines import Line2D  # noqa: E402
kb.legend([Line2D([], [], ls='', marker='o', ms=4.4, mfc='white', mec='#1f6aa5', mew=0.9),
           Line2D([], [], ls='', marker='o', ms=4.4, color='#1f2326')],
          ['detectability only (25 columns)', 'full model (1,107 columns)'], loc='lower left', fontsize=C.FS_S,
          borderaxespad=0, handletextpad=0.3, bbox_to_anchor=(-0.05, 0.0))
sd.add('a', 'printed constants', 'interval_level_percent', 95, N22_LBL, 'printed')
qa = C.check_figure(fig, 'Fig6', C.number_pool(sd), allow_numbers=('044043', '072089', '0.5'))
C.save(fig, 'Fig6_predictors', qa)
print(sd.write())
