"""Figure 6 of the manuscript, rebuilt: detectability control shifts each re-tested claim by a
visible amount. Run from the scripts/ directory; reads the released Source Data and the two
supplemental tables, writes figures_r31/Fig6_claim_retests.pdf and .png.
"""
"""Figure 6, finalised.

Contract
  core conclusion : under the matched detectability control almost every re-tested claim keeps its
                    effect; the few that move, and the few that sit on the wrong side of zero,
                    are readable at a glance. The figure shows the shift, not a list of intervals.
  evidence chain  : per claim, the baseline estimate under the authors' own definition, the matched
                    estimate with its interval, and the verdict class; PERS-010 kept separate
                    because its estimand is an observed-to-expected ratio.
  archetype       : quantitative grid - one hero panel plus one small separately-scaled inset.
  backend         : Python / matplotlib (exclusive).
  export          : PDF + PNG, 7 pt Arial, editable text, 183 mm wide (the paper's figure width).
"""
import csv, os, sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.stdout.reconfigure(encoding='utf-8')

mpl.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'font.size': 7, 'axes.linewidth': 0.8,
    'axes.spines.right': False, 'axes.spines.top': False, 'legend.frameon': False,
})
DARK, GREY = (31 / 255, 35 / 255, 38 / 255), (111 / 255, 117 / 255, 124 / 255)
VC = {'survives': (0.122, 0.416, 0.647), 'undecidable': (0.604, 0.631, 0.659),
      'attenuated': (0.435, 0.663, 0.839), 'null_broken_by_control': (0.753, 0.424, 0.596),
      'vanishes': (0.863, 0.604, 0.169), 'baseline_contradicts_claim': (0.8, 0.353, 0.212)}
LAB = {'survives': 'survives', 'undecidable': 'undecidable', 'attenuated': 'attenuated',
       'null_broken_by_control': 'null broken by control', 'vanishes': 'vanishes',
       'baseline_contradicts_claim': 'baseline contradicts claim'}
# the round and flag suffixes the released figure carries, carried over verbatim
SUF = {'SFE-001': '2E', 'SFE-002': '2D', 'SNO-006': '2D', 'SNO-021': '2', 'PERS-009': '2B',
       'SNO-012': '2E  S', 'SNO-004': '2D  n', 'SFE-006': '2  T', 'SNO-016': '2B  T',
       'PERS-002': '2B  n', 'PERS-003': '2B', 'SFE-007': '2B  T', 'SFE-008': '2B  T',
       'SNO-001': '2E', 'SNO-002': '2E', 'SNO-005': '2D', 'SNO-009': '2E', 'SNO-014': '2D',
       'PERS-005': '2', 'PERS-011': '2', 'PERS-012': '2', 'PERS-006': '2',
       'SFI-001': '2D', 'SFI-002': '2D', 'SNO-008': '2D', 'PERS-008': '2D  SR', 'SFE-011': '2B  T'}

SRC = '../source_data_submitted/Source_Data_Fig6_claim_retests.csv'
vals = {}
for r in csv.DictReader(open(SRC, encoding='utf-8-sig')):
    try:
        vals.setdefault(r['row'], {}).setdefault(r['panel'], {})[r['field']] = float(r['value'])
    except ValueError:
        pass
verdict = {r['claim_id']: r['verdict'] for r in csv.DictReader(
    open('../supplemental/Supplemental_Data_2_claim_verdict_tally.csv', encoding='utf-8-sig'))}
unit = {r['claim_id']: r['unit'] for r in csv.DictReader(
    open('../supplemental/Supplemental_Data_10_claim_independence.csv', encoding='utf-8-sig'))
    if r['has_measured_baseline'] == '1'}

def panel_of(c):
    for k in ('a', 'b'):
        if k in vals.get(c, {}):
            return k
    return None
main = [c for c in SUF if panel_of(c)]
site = sorted([c for c in main if unit.get(c) == 'site'], key=lambda c: -vals[c][panel_of(c)].get('matched_log2_or', 0))
prot = sorted([c for c in main if unit.get(c) == 'protein'], key=lambda c: -vals[c][panel_of(c)].get('matched_log2_or', 0))
order = site + prot
n = len(order)
assert n == 27, n
print(f'figure 6: {n} claims in the hero panel (site {len(site)} + protein {len(prot)}) + PERS-010 inset')
os.makedirs('../figures_r31', exist_ok=True)
W, H = 518.7 / 72, 330 / 72

fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0.125, 0.115, 0.585, 0.865])
ys = np.arange(n)[::-1]
for y, c in zip(ys, order):
    v = vals[c][panel_of(c)]; b, m = v.get('baseline_log2_or', 0), v.get('matched_log2_or', 0)
    col = VC.get(verdict[c], GREY)
    ax.plot([b, m], [y, y], color=col, lw=1.7, alpha=0.45, solid_capstyle='butt', zorder=2)
    ax.plot([v.get('matched_ci_low', m), v.get('matched_ci_high', m)], [y, y], color=col, lw=0.9, zorder=3)
    ax.scatter([m], [y], s=24, color=col, zorder=4, linewidths=0)
    ax.scatter([b], [y], s=15, facecolor='white', edgecolor=DARK, linewidths=0.7, zorder=4)
ax.axvline(0, color=DARK, lw=0.9, zorder=1)
ax.set_yticks(ys)
ax.set_yticklabels([f'{c}  {SUF[c]}' for c in order], fontsize=6.3, color=DARK)
ax.set_xlabel('log$_2$ odds ratio', color=DARK, fontsize=7)
ax.tick_params(axis='x', colors=GREY, labelsize=6.5, length=2)
ax.tick_params(axis='y', length=2, colors=DARK)
ax.spines['left'].set_visible(False)
ax.set_ylim(-0.9, n - 0.1)
ax.set_xlim(-2.6, max(6.6, ax.get_xlim()[1]))
split = len(prot) - 0.5
ax.axhline(split, color=GREY, lw=0.6, ls=(0, (3, 2)), xmin=0.0, xmax=1.05)
ax.text(0.012, split, 'site-level claims (%d)' % len(site), transform=ax.get_yaxis_transform(),
        fontsize=6, color=GREY, va='bottom', ha='left')
ax.text(0.012, split, 'protein-level claims (%d)' % len(prot), transform=ax.get_yaxis_transform(),
        fontsize=6, color=GREY, va='top', ha='left')
bx, bw = 1.035, 0.020
for y, c in zip(ys, order):
    ax.add_patch(Rectangle((bx, y - 0.34), bw, 0.68, transform=ax.get_yaxis_transform(),
                           color=VC.get(verdict[c], GREY), clip_on=False))
ax.text(bx + bw / 2, n - 0.55, 'verdict', transform=ax.get_yaxis_transform(), fontsize=5.9,
        color=GREY, ha='center', va='bottom')

# ---- small separately scaled panel for PERS-010 ----
axc = fig.add_axes([0.795, 0.135, 0.185, 0.115])
p = vals['PERS-010']['c']
b, m = p['baseline_log2_or'], p['matched_log2_or']
col = VC.get(verdict['PERS-010'], GREY)
axc.plot([b, m], [0, 0], color=col, lw=1.7, alpha=0.45, solid_capstyle='butt')
axc.plot([p['matched_ci_low'], p['matched_ci_high']], [0, 0], color=col, lw=0.9)
axc.scatter([m], [0], s=24, color=col, zorder=4, linewidths=0)
axc.scatter([b], [0], s=15, facecolor='white', edgecolor=DARK, linewidths=0.7, zorder=4)
axc.set_yticks([]); axc.set_ylim(-0.6, 0.6)
axc.set_xlabel('log$_2$ (observed / expected shared)', fontsize=6, color=DARK)
axc.tick_params(axis='x', colors=GREY, labelsize=6, length=2)
axc.spines['left'].set_visible(False)
axc.set_title('PERS-010, own estimand', fontsize=6.3, color=DARK, pad=3, loc='left')

handles = [plt.Line2D([], [], marker='o', ls='', ms=4, color=GREY,
                      label='after matching, with its 95% interval'),
           plt.Line2D([], [], marker='o', ls='', ms=3.4, mfc='white', mec=DARK, mew=0.7,
                      label='baseline (authors’ definition)')]
handles += [plt.Line2D([], [], marker='s', ls='', ms=4, color=v, label=LAB[k])
            for k, v in VC.items() if any(verdict[x] == k for x in order)]
ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.075, 1.0), fontsize=5.9,
          labelcolor=DARK, handletextpad=0.45, borderpad=0.2, labelspacing=0.30)
# the T/S/R/n flags are defined in the figure legend, as in the released figure
META = {'CreationDate': None, 'ModDate': None}
fig.savefig('../figures_r31/Fig6_claim_retests.pdf', metadata=META); fig.savefig('../figures_r31/Fig6_claim_retests.png', dpi=300)
plt.close(fig)
print('../figures_r31/Fig6_claim_retests.pdf written')
