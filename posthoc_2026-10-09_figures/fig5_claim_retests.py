"""Fig 5 (v3.3): re-testing published claims under a matched control (v3.2 Fig 6 redrawn) + Fig S4 (protein level).

Figure plan (2026-10-09), Fig 5:
  a  18 site-level claims, blocked by attribute class (basic, acidic, structural, other; order as in the figure plan);
     row label 'ID . attribute'; per row baseline (open circle + interval), matched (filled, verdict colour + interval),
     same-size random control (grey diamond + interval); verdict text column; T = re-tested only in a transfer cohort.
     Values: source_data_submitted/Source_Data_Fig6_claim_retests.csv panel a (SFE-006 = own-data re-run, tallied).
  b  retained fraction matched/baseline, 6 claims (Note 14 Table S14.2), grey ticks = random/baseline of the same
     table (together 0.91-1.02); no intervals; note 'point estimates; 3 basic vs 2 acidic; post hoc'.
  c  verdict counts: site level 18 stacked (4, 2, 10, 1, 1); protein level 10 as one grey bar,
     'not estimable (Table 1, item 6)' (author decision 2026-10-09); source_data_submitted/Source_Data_Fig1_overview.csv d.
Fig S4: the protein-level block (9 claims, v3.2 Fig 6 panel b rows) and PERS-010 on its own estimand (panel c row),
     every mark grey, labelled 'not estimable (Table 1, item 6)'.
"""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import common as C

F6 = C.SD32 / 'Source_Data_Fig6_claim_retests.csv'
F6_LBL = 'cysteine-detectability-audit/source_data_submitted/Source_Data_Fig6_claim_retests.csv'
F1 = C.SD32 / 'Source_Data_Fig1_overview.csv'
F1_LBL = 'cysteine-detectability-audit/source_data_submitted/Source_Data_Fig1_overview.csv'
N14 = C.SUP / 'Supplemental_Note_14_claim_retest_tally.md'
N14_LBL = 'Supplemental Note 14'

BLOCKS = [('Basic residue', ['SFE-006', 'PERS-009', 'SNO-021', 'SFE-008', 'SNO-016']),
          ('Acidic residue', ['SFE-002', 'SNO-006', 'SFE-007', 'SNO-014']),
          ('Structural', ['SFE-001', 'SNO-001', 'SNO-002', 'SNO-009', 'SNO-012']),
          ('Other', ['PERS-002', 'PERS-003', 'SNO-004', 'SNO-005'])]
SHORT = {  # abbreviated from the Attribute column of Note 14 Table S14.1
    'SFE-006': 'K/R flanking', 'PERS-009': '≥3 A/K/R/V within ±5', 'SNO-021': 'K at ±10, +6, +5 or +1',
    'SFE-008': 'K/R at SFE-006 offsets', 'SNO-016': '≥3 K/R/H within ±5',
    'SFE-002': 'E at −4, −3, +1, +3 to +5', 'SNO-006': 'D/E at +3', 'SFE-007': 'E at −3, +1, +3 or +4',
    'SNO-014': 'D at −1 or +1',
    'SFE-001': 'accessibility > 25%', 'SNO-001': 'α-helix', 'SNO-002': 'accessibility > 10%',
    'SNO-009': 'coil at 0 to +3', 'SNO-012': 'sulfur exposed',
    'PERS-002': 'C at ±2 or ±3 (null)', 'PERS-003': '≥3 hydrophobic within ±5', 'SNO-004': '≥30% hydrophobic (null)',
    'SNO-005': 'G at −1'}
VF = ['baseline_log2_or', 'baseline_ci_low', 'baseline_ci_high', 'matched_log2_or', 'matched_ci_low',
      'matched_ci_high', 'random_control_log2_or', 'random_control_ci_low', 'random_control_ci_high']

f6 = C.read_csv(F6)
s141 = {r['Claim'].strip('`'): r for r in C.md_table(N14, '| Claim | Level | Attribute |')}
sd = C.SourceData('Fig5', 'Source_Data_Fig5_claim_retests.csv')


def claim_vals(panel, cid):
    v = {f: float(C.sd_lookup(f6, panel, cid, f)['value']) for f in VF}
    v['verdict'] = C.sd_lookup(f6, panel, cid, 'verdict')['value']
    v['flags'] = C.sd_lookup(f6, panel, cid, 'flags')['value'] if panel != 'c' else ''
    v['level'] = C.sd_lookup(f6, panel, cid, 'level')['value']
    # cross-check with Note 14 Table S14.1 (4 decimals)
    n = s141[cid]
    for col, keys in (('Baseline [95% CI]', VF[0:3]), ('Matched [95% CI]', VF[3:6]),
                      ('Random control [95% CI]', VF[6:9])):
        ref = C.pt_ci(n[col])
        assert all(abs(v[k] - x) < 6e-5 for k, x in zip(keys, ref)), (cid, col, ref, [v[k] for k in keys])
    return v


# ---------------------------------------------------------------------------------------------- a
A = {}
site_ids = [c for _, ids in BLOCKS for c in ids]
assert sorted(site_ids) == sorted(r['row'] for r in f6 if r['panel'] == 'a' and r['field'] == 'verdict')
for cid in site_ids:
    v = claim_vals('a', cid)
    assert v['level'] == 'site'
    nverd = s141[cid]['Verdict'].replace(' ', '_')
    assert nverd == v['verdict'], (cid, nverd, v['verdict'])
    A[cid] = v
    for f in VF:
        src = C.sd_lookup(f6, 'a', cid, f)
        sd.add('a', cid, f, v[f], '%s (panel a, row %s, field %s; originally %s)' % (F6_LBL, cid, f, src['source_table']),
               'plotted; log2 odds ratio, 95% cluster-bootstrap interval; ' + (src['specification_label'] or 'primary'))
    sd.add('a', cid, 'verdict', v['verdict'], '%s (panel a, row %s, verdict = Supplemental Data 2 final_verdict)' % (F6_LBL, cid),
           'plotted: matched-estimate colour and verdict column')
    sd.add('a', cid, 'transfer_only', 'T' if 'T' in v['flags'] else '', '%s (panel a, row %s, flags)' % (F6_LBL, cid),
           'printed T: re-tested only in a transfer cohort')
    sd.add('a', cid, 'attribute_label', SHORT[cid], '%s, Table S14.1 (Attribute), abbreviated' % N14_LBL, 'row label')
assert {c for c in A if 'T' in A[c]['flags']} == {'SFE-007', 'SFE-008', 'SNO-016'}

# ---------------------------------------------------------------------------------------------- b
s142 = C.md_table(N14, '| Claim | Attribute class | Baseline |')
B = []
for r in s142:
    cid = r['Claim'].strip('`')
    B.append((cid, r['Attribute class'], C.num(r['Matched / baseline']), C.num(r['Random / baseline'])))
    for f, col in (('baseline_log2_or', 'Baseline'), ('matched_log2_or', 'Matched'), ('matched_over_baseline', 'Matched / baseline'),
                   ('random_control_log2_or', 'Random control'), ('random_over_baseline', 'Random / baseline')):
        sd.add('b', cid, f, C.num(r[col]), '%s, Table S14.2 (%s; row %s)' % (N14_LBL, col, cid),
               'plotted' if 'over' in f else 'not plotted; context')
    sd.add('b', cid, 'attribute_class', r['Attribute class'], '%s, Table S14.2' % N14_LBL, 'x grouping')
assert [b[0] for b in B] == ['SFE-006', 'PERS-009', 'SNO-021', 'SFE-002', 'SNO-006', 'SFE-001']
assert [b[2] for b in B] == [0.430, 0.371, 0.552, 0.829, 1.020, 0.886]
assert min(b[3] for b in B) == 0.911 and max(b[3] for b in B) == 1.020
sd.add('b', 'random control range', 'min_random_over_baseline', 0.911, '%s, Table S14.2 (Random / baseline, min)' % N14_LBL,
       'printed 0.91 (rounded)')
sd.add('b', 'random control range', 'max_random_over_baseline', 1.020, '%s, Table S14.2 (Random / baseline, max)' % N14_LBL,
       'printed 1.02 (rounded)')
for f, v in (('basic_claims', 3), ('acidic_claims', 2)):
    sd.add('b', 'note', f, v, '%s, Table S14.2 (Attribute class counts)' % N14_LBL, 'printed: 3 basic vs 2 acidic')

# ---------------------------------------------------------------------------------------------- c
f1 = C.read_csv(F1)
CATS = [('Survives', 'survives'), ('Attenuated', 'attenuated'), ('Undecidable', 'undecidable'),
        ('Null broken by control', 'null_broken_by_control'), ('Baseline contradicts claim', 'baseline_contradicts_claim')]
CC = []
for lab, k in CATS:
    n = int(C.sd_lookup(f1, 'd', lab, 'site_level_claims')['value'])
    CC.append((lab, k, n))
    sd.add('c', lab, 'site_level_claims', n, '%s (panel d, row %s, site_level_claims)' % (F1_LBL, lab), 'plotted segment')
assert [c[2] for c in CC] == [4, 2, 10, 1, 1]
assert sum(1 for c in A.values() if c['verdict'] == 'undecidable') == 10
NSITE = int(C.sd_lookup(f1, 'd', 'All re-tested claims', 'site_level_claims')['value'])
NPROT = int(C.sd_lookup(f1, 'd', 'All re-tested claims', 'protein_level_claims')['value'])
assert (NSITE, NPROT) == (18, 10)
sd.add('c', 'All re-tested claims', 'site_level_claims', NSITE, F1_LBL + ' (panel d, All re-tested claims)', 'printed bar label')
sd.add('c', 'All re-tested claims', 'protein_level_claims', NPROT, F1_LBL + ' (panel d, All re-tested claims)',
       'printed bar label; grey bar, not estimable (Table 1, item 6); no verdict colour')
sd.add('c', 'protein level', 'table1_item', 6, 'Table 1, item 6',
       'printed: not estimable (Table 1, item 6)')

# ---------------------------------------------------------------------------------------------- draw Fig 5
fig = plt.figure(figsize=(C.W_IN, C.mm(128)))
ax = fig.add_axes([0.215, 0.115, 0.285, 0.825])
y = 0
YS, HEAD = {}, []
for bi, (bname, ids) in enumerate(BLOCKS):
    HEAD.append((bname, y))
    y -= 0.85
    for cid in ids:
        YS[cid] = y
        y -= 1
    y -= 0.35
YMIN = y + 0.35
OFF = 0.24
for cid, yy in YS.items():
    v = A[cid]
    col = C.VC[v['verdict']]
    ax.plot([v['baseline_ci_low'], v['baseline_ci_high']], [yy + OFF] * 2, color='#4a4f55', lw=0.6)
    ax.plot(v['baseline_log2_or'], yy + OFF, 'o', ms=3.2, mfc='white', mec=C.INK, mew=0.7, zorder=4)
    ax.plot([v['matched_ci_low'], v['matched_ci_high']], [yy] * 2, color=col, lw=1.0)
    ax.plot(v['matched_log2_or'], yy, 'o', ms=3.8, color=col, zorder=5)
    ax.plot([v['random_control_ci_low'], v['random_control_ci_high']], [yy - OFF] * 2, color='#b4b9bf', lw=0.6)
    ax.plot(v['random_control_log2_or'], yy - OFF, 'D', ms=2.6, color='#9aa1a8', zorder=4)
    ax.text(1.075, yy, C.VLABEL[v['verdict']], color=col if v['verdict'] != 'undecidable' else '#6f757c',
            fontsize=C.FS_S, va='center', ha='left', clip_on=False, transform=ax.get_yaxis_transform())
    if 'T' in v['flags']:
        ax.text(1.035, yy, 'T', fontsize=C.FS_S, va='center', ha='center', color=C.MUTE, clip_on=False,
                transform=ax.get_yaxis_transform())
ax.axvline(0, color=C.AXIS, lw=0.8, zorder=0)
ax.set_ylim(YMIN - 0.3, 0.45)
ax.set_xlim(-2.8, 3.3)
ax.set_xticks([-2, -1, 0, 1, 2, 3])
ax.set_yticks(list(YS.values()))
ax.set_yticklabels(['%s · %s' % (c, SHORT[c]) for c in YS], fontsize=C.FS_S)
ax.tick_params(axis='y', length=0, pad=2)
ax.spines['left'].set_visible(False)
ax.set_xlabel('log$_2$ odds ratio')
for bname, yy in HEAD:
    ax.text(-0.02, yy - 0.3, bname, transform=ax.get_yaxis_transform(), fontsize=C.FS, fontweight='bold',
            ha='right', va='center')
    if yy < 0:
        ax.axhline(yy + 0.05, color=C.REF, lw=0.5, xmin=-0.72, xmax=1.5, clip_on=False)
C.panel_letter(fig, (0.012, 0.955), 'a', 'Site-level claims (18) under a matched control')
# key for a
kx = 0.40
fig.text(0.205, 0.012, '', fontsize=C.FS_S)
handles = [Line2D([], [], color='#4a4f55', lw=0.6, marker='o', ms=3.2, mfc='white', mec=C.INK),
           Line2D([], [], color=C.VC['survives'], lw=1.0, marker='o', ms=3.8),
           Line2D([], [], color='#b4b9bf', lw=0.6, marker='D', ms=2.6, mfc='#9aa1a8', mec='#9aa1a8')]
kax = fig.add_axes([0.012, 0.0, 0.66, 0.035])
kax.axis('off')
kax.legend(handles, ['reconstructed baseline', 'matched (colour = verdict)', 'same-size random control'],
           ncol=3, loc='lower left', fontsize=C.FS_S, handlelength=1.8, columnspacing=1.0, borderaxespad=0,
           bbox_to_anchor=(0.0, 0.0))
fig.text(0.55, 0.008, 'T, transfer cohort only; 95% intervals', fontsize=C.FS_S, color=C.MUTE, ha='left', va='bottom')

# b
axb = fig.add_axes([0.775, 0.60, 0.205, 0.335])
xs = [0, 1, 2, 3.4, 4.4, 5.8]
for x, (cid, cls, mr, rr) in zip(xs, B):
    axb.plot([x - 0.28, x + 0.28], [rr, rr], color='#9aa1a8', lw=1.6, solid_capstyle='butt')
    axb.plot(x, mr, 'o', ms=4.2, color=C.INK)
axb.axhline(1, color=C.REF, lw=0.6, ls=(0, (3, 2)), zorder=0)
axb.set_xticks(xs)
axb.set_xticklabels([b[0] for b in B], rotation=90, fontsize=C.FS_S)
axb.set_ylim(0, 1.25)
axb.set_yticks([0, 0.5, 1.0])
axb.set_xlim(-0.6, 6.4)
axb.set_ylabel('Matched / baseline')
for xc, lab in ((1, 'basic'), (3.9, 'acidic'), (5.8, 'struct.')):
    axb.text(xc, 1.12, lab, ha='center', va='center', fontsize=C.FS_S, color=C.MUTE)
axb.text(6.35, 0.06, 'point estimates;\n3 basic vs 2 acidic;\npost hoc', fontsize=C.FS_S, ha='right', va='bottom',
         color=C.MUTE)
axb.text(-0.5, 0.06, 'grey: random\ncontrol, 0.91–1.02', fontsize=C.FS_S, ha='left', va='bottom', color='#6f757c')
C.panel_letter(fig, axb, 'b', 'Retained fraction')

# c
axc = fig.add_axes([0.815, 0.235, 0.165, 0.16])
left = 0
for lab, k, n in CC:
    axc.barh(1, n, left=left, height=0.55, color=C.VC[k], lw=0)
    left += n
axc.barh(0, NPROT, height=0.55, color='#d9dcdf', lw=0)
axc.text(NPROT + 0.4, 0, 'not estimable\n(Table 1,\nitem 6)', fontsize=C.FS_S, va='center', ha='left', color=C.MUTE)
axc.set_yticks([1, 0])
axc.set_yticklabels(['Site (%d)' % NSITE, 'Protein (%d)' % NPROT])
axc.tick_params(axis='y', length=0)
axc.spines['left'].set_visible(False)
axc.set_xlim(0, 18.6)
axc.set_ylim(-0.5, 1.5)
axc.set_xticks([0, 6, 12, 18])
axc.set_xlabel('Claims')
# per-verdict key with counts, under the axis
kax2 = fig.add_axes([0.735, 0.03, 0.25, 0.10])
kax2.axis('off')
for i, (lab, k, n) in enumerate(CC):
    cx, cy = (0.0 if i < 3 else 0.0), 0.92 - i * 0.2
    kax2.add_patch(Rectangle((cx, cy - 0.07), 0.035, 0.14, color=C.VC[k], transform=kax2.transAxes, clip_on=False))
    kax2.text(cx + 0.055, cy, '%s (%d)' % (C.VLABEL[k], n), fontsize=C.FS_S, va='center', transform=kax2.transAxes)
C.panel_letter(fig, axc, 'c', 'Verdicts')

qa = C.check_figure(fig, 'Fig5', C.number_pool(sd), allow_numbers=('95', '0', '1', '2', '3', '4', '5', '10', '25', '30', '6', '7'))
C.save(fig, 'Fig5_claim_retests', qa)
sd.add('a', 'printed constants', 'interval_level_percent', 95, F6_LBL + ' (95% cluster-bootstrap, Note 14)', 'printed in key')
print(sd.write())

# ---------------------------------------------------------------------------------------------- Fig S4
sd4 = C.SourceData('FigS4', 'Source_Data_FigS4_protein_level_claims.csv')
prot_ids = [r['row'] for r in f6 if r['panel'] == 'b' and r['field'] == 'verdict']
P = {}
for cid in prot_ids:
    v = claim_vals('b', cid)
    assert v['level'] == 'protein'
    P[cid] = v
    for f in VF:
        src = C.sd_lookup(f6, 'b', cid, f)
        sd4.add('a', cid, f, v[f], '%s (panel b, row %s, field %s; originally %s)' % (F6_LBL, cid, f, src['source_table']),
                'plotted grey; log2 odds ratio, 95% cluster-bootstrap interval')
    sd4.add('a', cid, 'stored_verdict', v['verdict'], '%s (panel b, row %s, verdict)' % (F6_LBL, cid),
            'not plotted: protein-level commonality not estimable (Table 1, item 6)')
    sd4.add('a', cid, 'attribute_label', s141[cid]['Attribute'], '%s, Table S14.1 (Attribute)' % N14_LBL, 'row label')
assert len(P) == 9
pc = claim_vals('c', 'PERS-010')
for f in VF:
    sd4.add('b', 'PERS-010', f, pc[f], '%s (panel c, row PERS-010, field %s)' % (F6_LBL, f),
            "plotted grey; log2(observed / expected shared), the claim's own estimand")
sd4.add('b', 'PERS-010', 'stored_verdict', pc['verdict'], '%s (panel c, verdict)' % F6_LBL,
        'not plotted: not estimable (Table 1, item 6)')
sd4.add('a-b', 'protein level', 'table1_item', 6, 'Table 1, item 6',
        'printed: not estimable (Table 1, item 6)')
sd4.add('a-b', 'protein level', 'claims', 10, F1_LBL + ' (panel d, All re-tested claims, protein_level_claims)', 'printed')

PSHORT = {'PERS-005': 'glycolysis', 'PERS-006': 'primary metabolism', 'PERS-008': 'ubiquitin-dependent catabolism',
          'PERS-011': 'thylakoid', 'PERS-012': 'carboxylic acid metabolism', 'SFE-011': 'photosynthesis',
          'SFI-001': 'extracellular exosome', 'SFI-002': 'oxidoreductase activity', 'SNO-008': 'oxidoreductase activity'}
order = sorted(P, key=lambda c: -P[c]['matched_log2_or'])
fig = plt.figure(figsize=(C.W_IN, C.mm(80)))
ax = fig.add_axes([0.25, 0.17, 0.43, 0.72])
GR, GL = '#6f757c', '#b4b9bf'
for i, cid in enumerate(order):
    v, yy = P[cid], -i
    ax.plot([v['baseline_ci_low'], v['baseline_ci_high']], [yy + OFF] * 2, color=GL, lw=0.6)
    ax.plot(v['baseline_log2_or'], yy + OFF, 'o', ms=3.2, mfc='white', mec=GR, mew=0.7, zorder=4)
    ax.plot([v['matched_ci_low'], v['matched_ci_high']], [yy] * 2, color=GR, lw=1.0)
    ax.plot(v['matched_log2_or'], yy, 'o', ms=3.8, color=GR, zorder=5)
    ax.plot([v['random_control_ci_low'], v['random_control_ci_high']], [yy - OFF] * 2, color=GL, lw=0.6)
    ax.plot(v['random_control_log2_or'], yy - OFF, 'D', ms=2.6, color='#9aa1a8', zorder=4)
    if 'T' in C.sd_lookup(f6, 'b', cid, 'flags')['value']:
        ax.text(8.1, yy, 'T', fontsize=C.FS_S, va='center', ha='center', color=C.MUTE)
ax.axvline(0, color=C.AXIS, lw=0.8, zorder=0)
ax.set_yticks([-i for i in range(len(order))])
ax.set_yticklabels(['%s · %s' % (c, PSHORT[c]) for c in order], fontsize=C.FS_S)
ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.set_xlim(-1.5, 8.4)
ax.set_ylim(-len(order) + 0.4, 0.6)
ax.set_xlabel('log$_2$ odds ratio')
C.panel_letter(fig, (0.012, 0.94), 'a', 'Protein-level claims (10): not estimable (Table 1, item 6)')
axi = fig.add_axes([0.78, 0.45, 0.2, 0.22])
for f, yy, mk, col in ((('baseline_log2_or', 'baseline_ci_low', 'baseline_ci_high'), 0.3, 'o', 'white'),
                       (('matched_log2_or', 'matched_ci_low', 'matched_ci_high'), 0.0, 'o', GR),
                       (('random_control_log2_or', 'random_control_ci_low', 'random_control_ci_high'), -0.3, 'D', '#9aa1a8')):
    axi.plot([pc[f[1]], pc[f[2]]], [yy, yy], color=GL if col != GR else GR, lw=0.8)
    axi.plot(pc[f[0]], yy, mk, ms=3.4 if mk == 'o' else 2.6, mfc=col, mec=GR if col == 'white' else col, mew=0.7)
axi.set_ylim(-0.6, 0.6)
axi.set_yticks([])
axi.spines['left'].set_visible(False)
axi.set_xlim(0.7, 1.9)
axi.set_xticks([0.8, 1.2, 1.6])
axi.set_xlabel('log$_2$ (observed / expected shared)\nmodified in all four organs', fontsize=C.FS)
C.panel_letter(fig, (0.75, 0.74), 'b', 'PERS-010, own estimand')
kax = fig.add_axes([0.25, 0.0, 0.7, 0.05])
kax.axis('off')
kax.legend([Line2D([], [], color=GL, lw=0.6, marker='o', ms=3.2, mfc='white', mec=GR),
            Line2D([], [], color=GR, lw=1.0, marker='o', ms=3.8),
            Line2D([], [], color=GL, lw=0.6, marker='D', ms=2.6, mfc='#9aa1a8', mec='#9aa1a8')],
           ['reconstructed baseline', 'matched', 'same-size random control'], ncol=3, loc='lower left',
           fontsize=C.FS_S, handlelength=1.8, borderaxespad=0, bbox_to_anchor=(-0.02, 0))
fig.text(0.78, 0.012, 'T, transfer cohort only; 95% intervals', fontsize=C.FS_S, color=C.MUTE, va='bottom')
sd4.add('a-b', 'printed constants', 'interval_level_percent', 95, F6_LBL, 'printed in key')
sd4.add('b', 'PERS-010', 'organs', 4, '%s, Table S14.1 (Attribute: Modified in all four organs)' % N14_LBL, 'printed title')
qa = C.check_figure(fig, 'FigS4', C.number_pool(sd4), allow_numbers=('9', '0'))
C.save(fig, 'FigS4_protein_level_claims', qa)
print(sd4.write())
