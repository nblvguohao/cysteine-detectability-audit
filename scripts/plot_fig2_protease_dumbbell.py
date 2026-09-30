"""Figure 2 of the manuscript, rebuilt: panel a is a dumbbell so that what the trypsin rule
loses is a visible length. Run from the scripts/ directory; reads the released Source Data and
writes figures_r31/Fig2_public_four_protease.pdf and .png.
"""
"""Figure 2 rebuilt so that its three panels carry three different readings of the same deposit:
a, the signature follows the protease used (dumbbell); b, it is absent against the observed
background (dot plot); c, how much of the negative class is left (lollipop with the registered
decision levels). Same data, same dimensions as the released figure (183 x 104 mm).

backend: Python / matplotlib (exclusive).
"""
import csv, os, sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')

mpl.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none', 'font.size': 7,
                     'axes.linewidth': 0.8, 'axes.spines.right': False, 'axes.spines.top': False,
                     'legend.frameon': False})
DARK, GREY = (31 / 255, 35 / 255, 38 / 255), (111 / 255, 117 / 255, 124 / 255)
BLUE, PINK = (0.122, 0.416, 0.647), (0.753, 0.424, 0.596)

D = {}
for r in csv.DictReader(open('../source_data_submitted/Source_Data_Fig2_public_four_protease.csv',
                             encoding='utf-8-sig')):
    D.setdefault(r['panel'], {}).setdefault(r['row'], {})[r['field']] = float(r['value'])
A, B, C = D['a'], D['b'], D['c']
ARMS = ['Trypsin', 'AspN', 'CT', 'GluC']
NICE = {'Trypsin': 'Trypsin', 'AspN': 'AspN', 'CT': 'Chymotrypsin', 'GluC': 'GluC'}
RULE = {'Trypsin': 'trypsin', 'AspN': 'aspn', 'CT': 'chymotrypsin', 'GluC': 'gluc'}
BANDS = [('proximal_1_3', 'proximal, 1–3 residues'), ('distal_6_12', 'distal, 6–12 residues')]

os.makedirs('../figures_r31', exist_ok=True)
fig = plt.figure(figsize=(518.7 / 72, 294.8 / 72))

# ---------------- a: dumbbell ----------------
ax = fig.add_axes([0.125, 0.47, 0.60, 0.48])
ylab, ys, y = [], [], 0.0
for band, _ in BANDS:
    for arm in ARMS:
        own = A.get(f'{arm}|{RULE[arm]}|proteome|{band}')
        tryp = A.get(f'{arm}|trypsin|proteome|{band}')
        if own:
            ax.plot([own['ci_low'], own['ci_high']], [y, y], color=BLUE, lw=0.9, zorder=3)
            ax.scatter([own['estimate']], [y], s=26, color=BLUE, zorder=5, linewidths=0)
        if tryp and (not own or abs(tryp['estimate'] - own['estimate']) > 1e-9):
            ax.plot([tryp['ci_low'], tryp['ci_high']], [y, y], color=GREY, lw=0.9, zorder=3)
            ax.scatter([tryp['estimate']], [y], s=20, facecolor='white', edgecolor=GREY,
                       linewidths=0.8, zorder=5)
            if own:
                ax.plot([own['estimate'], tryp['estimate']], [y, y], color=GREY, lw=1.5,
                        alpha=0.45, zorder=2, solid_capstyle='butt')
        ylab.append(NICE[arm]); ys.append(y); y += 1
    y += 0.9
ax.axvline(0, color=DARK, lw=0.9, zorder=1)
ax.set_yticks(ys); ax.set_yticklabels(ylab, fontsize=6.4, color=DARK)
ax.tick_params(axis='y', length=2, colors=DARK)
ax.tick_params(axis='x', colors=GREY, labelsize=6.4, length=2)
ax.set_xlabel('cleavage-geometry log$_2$ odds ratio', color=DARK, fontsize=6.8)
ax.spines['left'].set_visible(False)
ax.set_ylim(-0.9, y - 1.8)
split = ys[len(ARMS) - 1]
ax.axhline(split + 0.45, color=GREY, lw=0.6, ls=(0, (3, 2)))
ax.text(0.010, split + 0.45, BANDS[0][1], transform=ax.get_yaxis_transform(), fontsize=6,
        color=GREY, va='top', ha='left')
ax.text(0.010, split + 0.45, BANDS[1][1], transform=ax.get_yaxis_transform(), fontsize=6,
        color=GREY, va='bottom', ha='left')
h = [plt.Line2D([], [], marker='o', ls='', ms=4, color=BLUE, label='protease’s own cleavage rule'),
     plt.Line2D([], [], marker='o', ls='', ms=3.6, mfc='white', mec=GREY, mew=0.8,
                label='scored with the trypsin rule'),
     plt.Line2D([], [], color=GREY, lw=1.5, alpha=0.45, label='what the trypsin rule loses')]
ax.legend(handles=h, loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize=6, labelcolor=DARK,
          handletextpad=0.5, borderpad=0.2, labelspacing=0.35)
fig.text(0.055, 0.952, 'a', fontsize=8, fontweight='bold', color=DARK)
fig.text(0.088, 0.956, 'Proteome background: the signature follows the protease that was used',
         fontsize=6.8, color=DARK, va='bottom')

# ---------------- b: observed background, dot plot ----------------
axb = fig.add_axes([0.125, 0.115, 0.35, 0.24])
for i, arm in enumerate(ARMS):
    v = B.get(f'{arm}|{RULE[arm]}|observed|distal_6_12')
    if not v:
        continue
    axb.plot([v['ci_low'], v['ci_high']], [i, i], color=PINK, lw=0.9, zorder=3)
    axb.scatter([v['estimate']], [i], s=22, color=PINK, zorder=5, linewidths=0)
axb.axvline(0, color=DARK, lw=0.9, zorder=1)
axb.set_yticks(range(len(ARMS))); axb.set_yticklabels([NICE[a] for a in ARMS], fontsize=6.2, color=DARK)
axb.tick_params(axis='y', length=2, colors=DARK)
axb.tick_params(axis='x', colors=GREY, labelsize=6.2, length=2)
axb.set_xlabel('cleavage-geometry log$_2$ odds ratio', color=DARK, fontsize=6.6)
axb.spines['left'].set_visible(False)
axb.set_ylim(-0.7, len(ARMS) - 0.3)
fig.text(0.055, 0.368, 'b', fontsize=8, fontweight='bold', color=DARK)
fig.text(0.088, 0.372, 'Observed background: the distal signature is absent',
         fontsize=6.6, color=DARK, va='bottom')

# ---------------- c: coincidence, lollipop with the decision levels ----------------
axc = fig.add_axes([0.63, 0.115, 0.35, 0.24])
for i, arm in enumerate(ARMS):
    v = C.get(arm)
    if not v:
        continue
    axc.plot([0, v['estimate']], [i, i], color=GREY, lw=1.0, zorder=2)
    axc.plot([v['ci_low'], v['ci_high']], [i, i], color=GREY, lw=2.6, alpha=0.30, zorder=1)
    axc.scatter([v['estimate']], [i], s=22, color=GREY, zorder=5, linewidths=0)
    axc.text(v['estimate'] + 0.03, i, f"{v['estimate']:.2f}", fontsize=5.8, color=DARK, va='center')
for level in (0.50, 0.90):
    axc.axvline(level, color=GREY, lw=0.6, ls=(0, (3, 2)), zorder=0)
axc.set_yticks(range(len(ARMS))); axc.set_yticklabels([NICE[a] for a in ARMS], fontsize=6.2, color=DARK)
axc.tick_params(axis='y', length=2, colors=DARK)
axc.tick_params(axis='x', colors=GREY, labelsize=6.2, length=2)
axc.set_xlim(0, 1.05)
axc.set_xlabel('identified cysteines carrying a site assignment', color=DARK, fontsize=6.6)
axc.spines['left'].set_visible(False)
axc.set_ylim(-0.7, len(ARMS) - 0.3)
fig.text(0.560, 0.368, 'c', fontsize=8, fontweight='bold', color=DARK)
fig.text(0.593, 0.372, 'How much of the negative class is left',
         fontsize=6.6, color=DARK, va='bottom')

META = {'CreationDate': None, 'ModDate': None}
fig.savefig('../figures_r31/Fig2_public_four_protease.pdf', metadata=META); fig.savefig('../figures_r31/Fig2_public_four_protease.png', dpi=300)
plt.close(fig)
print('../figures_r31/Fig2_public_four_protease.pdf written')
