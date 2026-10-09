"""Figure 6, final revision (2026-09-30): the 28 re-tested claims under the FINAL verdict policy.

Design kept from MCP/_revision/fig6_final.py (the released figure) and from the candidate rebuild
results/figures_rev/Fig6/fig6_rev_2026-09-30.py: one hero panel with the 18 site-level claims above a dashed line and
the nine protein-level claims below, each ranked by its matched estimate; open circle = reconstructed baseline,
filled circle = matched estimate, pale segment = the shift, thin line = the matched 95% interval, colour = verdict,
repeated in a strip at the right; PERS-010 in an inset on its own estimand.

Final verdict policy (results/final/final_tally.csv, final_tally_counts.json; as revised 2026-09-30 09:23)
  * every claim keeps its stored verdict (= Supplemental Data 2) except two: SFE-006 is tallied on its authors' own
    data (Yang 2014 component; attenuated; its transfer verdict, vanishes, is not tallied) and SNO-012 is undecidable
    (its baseline interval covers zero; stored as attenuated by pipeline code that lacked that branch; it survives
    under the pLDDT >= 70 specification, hence its flag S); 11 survive, 2 attenuated, 13 undecidable, 1 null broken
    by control (SNO-004), 1 baseline contradicts claim (SNO-016, transfer), 0 vanish;
  * colour classes and colours are the released ones (fig6_final.py, dict VC); 'vanishes' is no longer in use;
  * row labels are the released ones (fig6_final.py, dict SUF: round code and flags T, S, R, n), so PERS-002 keeps
    'n' and SNO-012 keeps 'S'; SFE-006 alone changes, to '2C  t': round C is its own-data re-run
    (Source_Data_text_phase2c_sfe006_reproduction.json) and t marks the transfer re-test drawn under its row.

Kept from the candidate rebuild (not verdict-related):
  * SFE-006 is drawn from its own-data re-run and labelled 'own data'; its transfer re-test (the released row, matched
    0.1761 [-0.0814, 0.4379]) is a grey open diamond with its interval, labelled 'transfer', just below the row;
  * the block labels sit in a gap between the blocks, so they no longer collide with the PERS-002 interval;
  * 'reconstructed baseline' replaces 'baseline (authors' definition)', which is wrong for transfer rows;
  * width exactly 7.2 in; mathtext in Arial.

No statistic is computed here. The script
  1. builds Source_Data_Fig6_claim_retests.csv by copying stored cells verbatim (numbers as stored strings) from the
     released Source Data (27 claims), the SFE-006 own-data JSON, Supplemental Data 4 (SFE-006 transfer row), the
     final tally, and the released figure script (labels, colours), cross-checking every number against its primary
     table (Supplemental Data 4/5/6) and every verdict against Supplemental Data 2 (gates B0-B8);
  2. re-reads that CSV and draws from it only (sorting by the matched estimate is the only operation on values);
  3. gates the drawing and the written PDF (G0-G7): drawn values present in Source Data, text collisions with text
     and with data marks, canvas, font sizes and families, page size.
Outputs are written to temporary names and renamed only after every gate has passed; a failed gate stops the run,
removes the temporary files and leaves any earlier outputs untouched.

Run:  python fig6_final_2026-09-30.py      (system python 3.10, matplotlib 3.10, PyMuPDF)
"""
import ast
import csv
import hashlib
import json
import os
import sys
from collections import Counter

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.text import Text  # noqa: E402

sys.stdout.reconfigure(encoding='utf-8')

# ------------------------------------------------------------------------------------------------ paths (read-only)
HERE = os.path.dirname(os.path.abspath(__file__))
MCP = '/path/to/local/巯基化/MCP'
REV = '/path/to/local/_cys_repo_work/public/revision_2026-09-30'
IN = {
    'final_tally': REV + '/results/final/final_tally.csv',
    'final_counts': REV + '/results/final/final_tally_counts.json',
    'fig6_released_script': MCP + '/_revision/fig6_final.py',
    'sfe006': MCP + '/source_data/Source_Data_text_phase2c_sfe006_reproduction.json',
    'sd_released': MCP + '/source_data/Source_Data_Fig6_claim_retests.csv',
    'SD2': MCP + '/supplemental/Supplemental_Data_2_claim_verdict_tally.csv',     # carries final_verdict (09:40)
    'SD4': MCP + '/supplemental/Supplemental_Data_4_retest_round_b.csv',
    'SD5': MCP + '/supplemental/Supplemental_Data_5_retest_round_d.csv',
    'SD6': MCP + '/supplemental/Supplemental_Data_6_retest_round_e.csv',
    'SD10': MCP + '/supplemental/Supplemental_Data_10_claim_independence.csv',   # unit = site / protein level
}
# names written into the source_table column (basenames, as in the released Source Data). The verdicts are READ from
# results/final/final_tally.csv (as specified) and gated equal, claim by claim, to the publishable Supplemental Data 2
# (final_verdict, verdict, final_verdict_note) and Supplemental Data 10 (unit), which the Source Data therefore names.
TAB = {'SD4': 'Supplemental_Data_4_retest_round_b.csv', 'SD5': 'Supplemental_Data_5_retest_round_d.csv',
       'SD6': 'Supplemental_Data_6_retest_round_e.csv',
       'SFE006_OWN_JSON': 'Source_Data_text_phase2c_sfe006_reproduction.json',
       'SD2': 'Supplemental_Data_2_claim_verdict_tally.csv',
       'SD10': 'Supplemental_Data_10_claim_independence.csv',
       'labels': 'row labels of the submitted Figure 6 (MCP/_revision/fig6_final.py, SUF)',
       'colours': 'verdict colours of the submitted Figure 6 (MCP/_revision/fig6_final.py, VC)'}
OUT = {'sd': os.path.join(HERE, 'Source_Data_Fig6_claim_retests.csv'),
       'pdf': os.path.join(HERE, 'Fig6_claim_retests.pdf'),
       'png': os.path.join(HERE, 'Fig6_claim_retests.png')}
TMP = {k: os.path.join(HERE, '_tmp_' + os.path.basename(v)) for k, v in OUT.items()}
OUT_LOG = os.path.join(HERE, 'build_log.txt')

VALUE_FIELDS = ['baseline_log2_or', 'baseline_ci_low', 'baseline_ci_high',
                'matched_log2_or', 'matched_ci_low', 'matched_ci_high',
                'random_control_log2_or', 'random_control_ci_low', 'random_control_ci_high']
TRANSFER_FIELDS = ['transfer_' + f for f in VALUE_FIELDS]
NUMERIC = set(VALUE_FIELDS) | set(TRANSFER_FIELDS)
CLASS_ORDER = ['survives', 'attenuated', 'undecidable', 'null_broken_by_control', 'baseline_contradicts_claim']
CLASS_LABEL = {'survives': 'survives', 'attenuated': 'attenuated', 'undecidable': 'undecidable',
               'null_broken_by_control': 'null broken by control',
               'baseline_contradicts_claim': 'baseline contradicts claim'}     # as the released legend
SFE006_SUFFIX = '2C  t'           # round C = own-data re-run; t = transfer re-test drawn under the row
LOG = []


def log(msg):
    LOG.append(msg)
    print(msg)


QUIET = {}                        # per-value gates: counted, logged once as a total


def gate(ok, name, detail='', quiet=False):
    if not ok:
        raise SystemExit('GATE FAILED  %s  %s' % (name, detail))
    if quiet:
        QUIET[name] = QUIET.get(name, 0) + 1
    else:
        log('gate ok   %s%s' % (name, ('  ' + detail) if detail else ''))


def flush_quiet():
    for name, n in QUIET.items():
        log('gate ok   %s  (%d checks)' % (name, n))
    QUIET.clear()


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def read_csv(p):
    return list(csv.DictReader(open(p, encoding='utf-8-sig', newline='')))


def released_dicts():
    """The dict literals VC (verdict colours) and SUF (row-label suffixes) of the released Fig. 6 script, read with
    ast; nothing is executed."""
    tree = ast.parse(open(IN['fig6_released_script'], encoding='utf-8').read())
    out = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if getattr(t, 'id', None) in ('VC', 'SUF'):
                    out[t.id] = ast.literal_eval(node.value)
    return {k: tuple(v) for k, v in out['VC'].items()}, out['SUF']


# ================================================================================ 1. build Source Data (copy only)
def build_source_data():
    final = read_csv(IN['final_tally'])
    fin = {r['claim_id']: r for r in final}
    counts = json.load(open(IN['final_counts'], encoding='utf-8'))
    sd2 = {r['claim_id']: r for r in read_csv(IN['SD2'])}
    own = json.load(open(IN['sfe006'], encoding='utf-8'), parse_float=str)     # numbers kept as stored text
    VC, SUF = released_dicts()
    released = {}
    for r in read_csv(IN['sd_released']):
        released.setdefault(r['row'], {'panel': r['panel']})[r['field']] = (r['value'], r['source_table'])
    sd = {'SD4': [], 'SD5': [], 'SD6': []}
    for k in sd:
        sd[k] = read_csv(IN[k])
    sd4 = {}
    for r in sd['SD4']:
        gate(r['claim_id'] not in sd4, 'B0 SD4 one row per claim', r['claim_id'], quiet=True)
        sd4[r['claim_id']] = r
    by_tab = {TAB[k]: k for k in sd}

    # ---- the claim set and the verdicts
    gate(len(final) == len(fin) == 28, 'B1 final_tally.csv lists 28 distinct claims')
    gate(set(fin) == set(released) == {c for c, r in sd2.items() if r['has_measured_baseline'] == '1'},
         'B1 same 28 claims in final_tally.csv, the released Source Data and SD2 (has_measured_baseline = 1)')
    for key, sel in (('overall', lambda r: r['final_verdict']),
                     ('by_level', lambda r: '%s|%s' % (r['level'], r['final_verdict']))):
        got = dict(Counter(sel(r) for r in final))
        gate(got == counts[key], 'B2 recount of final_tally.csv = final_tally_counts.json %s' % key,
             str(dict(sorted(got.items()))))
    changed = {c: (r['stored_verdict'], r['final_verdict']) for c, r in fin.items()
               if r['final_verdict'] != r['stored_verdict']}
    gate(all(r['stored_verdict'] == sd2[c]['verdict'] for c, r in fin.items()),
         'B3 stored_verdict = Supplemental Data 2 verdict for all 28 claims')
    sd10 = {r['claim_id']: r for r in read_csv(IN['SD10'])}
    gate(all(sd2[c]['final_verdict'] == r['final_verdict'] and sd10[c]['unit'] == r['level'] and
             (r['final_verdict'] == r['stored_verdict'] or sd2[c]['final_verdict_note'] == r['note'])
             for c, r in fin.items()),
         'B3 final_verdict (and its note) = Supplemental Data 2, level = Supplemental Data 10 unit, for all 28 claims')
    gate(changed == {'SFE-006': ('vanishes', 'attenuated'), 'SNO-012': ('attenuated', 'undecidable')},
         'B3 the final verdict differs from the stored one only for SFE-006 (vanishes -> attenuated) and SNO-012 '
         '(attenuated -> undecidable)', str(changed))
    gate('covers zero' in fin['SNO-012']['note'] and float(released['SNO-012']['baseline_ci_low'][0]) < 0 <
         float(released['SNO-012']['baseline_ci_high'][0]),
         'B3 SNO-012: the baseline interval drawn from the Source Data covers zero, as the final tally note states',
         '[%s, %s]' % (released['SNO-012']['baseline_ci_low'][0], released['SNO-012']['baseline_ci_high'][0]))
    gate(set(fin[c]['final_verdict'] for c in fin) == set(CLASS_ORDER) and all(k in VC for k in CLASS_ORDER),
         'B3 the five final classes all have a released colour', ', '.join(CLASS_ORDER))
    gate(set(SUF) == set(fin) - {'PERS-010'}, 'B3 released row labels exist for the 27 claims of the hero panel')

    # ---- blocks: level from the final tally; PERS-010 (own estimand) in the inset, as released
    panel_of = {c: ('c' if c == 'PERS-010' else 'a' if r['level'] == 'site' else 'b') for c, r in fin.items()}
    gate(all(panel_of[c] == released[c]['panel'] for c in fin), 'B4 block of every claim = released panel')

    rows = []
    n_same_as_released = 0
    for cid in sorted(fin):
        f, panel = fin[cid], panel_of[cid]
        if cid == 'SFE-006':
            src_key, spec_label = 'SFE006_OWN_JSON', ("own_data_yang2014_component: re-run on the Yang 2014 "
                                                      "component of the authors' cohort (tallied)")
            for field in VALUE_FIELDS:
                if field.endswith('_ci_low') or field.endswith('_ci_high'):
                    raw = own[field.rsplit('_ci_', 1)[0] + '_interval'][0 if field.endswith('_low') else 1]
                else:
                    raw = own[field]
                rows.append(dict(figure='Fig6', panel=panel, row=cid, field=field, value=raw,
                                 source_table=TAB[src_key], specification_label=spec_label))
            gate(own['claim_id'] == 'SFE-006' and own['verdict'] == f['final_verdict'] == 'attenuated',
                 'B5 own-data JSON is SFE-006 and its stored verdict is the final one', own['verdict'])
            gate(('%.4f' % float(own['baseline_log2_or']), '%.4f' % float(own['matched_log2_or'])) ==
                 ('1.3840', '0.5950') and '1.3840 -> 0.5950' in f['note'],
                 'B5 own-data baseline -> matched = 1.3840 -> 0.5950, as in the final tally note')
            # the transfer re-test: the released SFE-006 row = Supplemental Data 4, drawn as the grey overlay
            tr = sd4['SFE-006']
            gate(tr['is_transfer'] == 'True' and tr['verdict'] == f['stored_verdict'] == 'vanishes',
                 'B5 SD4 SFE-006 is the transfer re-test, verdict vanishes = stored verdict', tr['transfer_cohort'])
            tspec = 'primary; transfer cohort %s (drawn under the row, not tallied)' % tr['transfer_cohort']
            for field in VALUE_FIELDS:
                raw = tr[field]
                gate(released['SFE-006'][field] == (raw, TAB['SD4']),
                     'B5 transfer value identical to the released SFE-006 row', field, quiet=True)
                rows.append(dict(figure='Fig6', panel=panel, row=cid, field='transfer_' + field, value=raw,
                                 source_table=TAB['SD4'], specification_label=tspec))
            for field, value in (('transfer_verdict', tr['verdict']), ('transfer_cohort', tr['transfer_cohort']),
                                 ('transfer_baseline_reproduction', tr['baseline_reproduction'])):
                rows.append(dict(figure='Fig6', panel=panel, row=cid, field=field, value=value,
                                 source_table=TAB['SD4'], specification_label=tspec))
            rows.append(dict(figure='Fig6', panel=panel, row=cid, field='baseline_reproduction',
                             value='own_data_component', source_table=TAB[src_key],
                             specification_label="assigned in revision from reproduction_notes: positives = the Yang "
                                                 "2014 component of the authors' cohort, negatives by the authors' "
                                                 "rule (no released baseline_reproduction term covers a component "
                                                 "re-run)"))
            direction_row = tr
        else:
            rel = released[cid]
            tab = {rel[fl][1] for fl in VALUE_FIELDS} | {rel['baseline_reproduction'][1]}
            gate(len(tab) == 1 and tab <= set(by_tab), 'B4 one primary table per claim', cid, quiet=True)
            tab = tab.pop()
            match = [r for r in sd[by_tab[tab]] if r['claim_id'] == cid and
                     all(r[fl] == rel[fl][0] for fl in VALUE_FIELDS)]
            gate(len(match) == 1, 'B4 the 9 released numbers = exactly one row of their primary table', cid,
                 quiet=True)
            m = match[0]
            spec = m.get('specification', 'primary')
            spec_label = ('zf_background: zinc-finger complement, the authors\' primary background (tallied)'
                          if spec == 'zf_background' else spec)
            gate(m['baseline_reproduction'] == rel['baseline_reproduction'][0],
                 'B4 baseline_reproduction = primary table', cid, quiet=True)
            for field in VALUE_FIELDS + ['baseline_reproduction']:
                rows.append(dict(figure='Fig6', panel=panel, row=cid, field=field, value=rel[field][0],
                                 source_table=rel[field][1], specification_label=spec_label))
                n_same_as_released += field != 'baseline_reproduction'
            direction_row = m
        # ---- categories: what the colours, the strip and the row labels show
        cat = [('verdict', f['final_verdict'], TAB['SD2'] + ' (final_verdict)'),
               ('verdict_class', CLASS_LABEL[f['final_verdict']], TAB['colours']),
               ('stored_verdict', f['stored_verdict'], TAB['SD2'] + ' (verdict)'),
               ('level', f['level'], TAB['SD10'] + ' (unit)'),
               ('claim_direction', direction_row['claim_direction'],
                TAB['SD4'] if cid == 'SFE-006' else rel['matched_log2_or'][1])]
        if cid == 'PERS-010':
            cat.append(('inset_title', 'PERS-010, own estimand', 'MCP/_revision/fig6_final.py (inset title)'))
        else:
            suffix = SFE006_SUFFIX if cid == 'SFE-006' else SUF[cid]
            parts = suffix.split('  ')
            round_code, flags = parts[0], ' '.join(''.join(parts[1:]))
            src = ('revision 2026-09-30 (round 2C = the own-data re-run, %s; flag t = transfer re-test drawn '
                   'under the row; released suffix: %s)' % (TAB['SFE006_OWN_JSON'], SUF[cid])
                   if cid == 'SFE-006' else TAB['labels'])
            cat += [('row_label', cid + '  ' + suffix, src), ('round_code', round_code, src), ('flags', flags, src)]
        if f['final_verdict'] != f['stored_verdict']:          # SFE-006 and SNO-012: the reason, verbatim
            cat.append(('final_verdict_note', f['note'], TAB['SD2'] + ' (final_verdict_note)'))
        for field, value, s in cat:
            rows.append(dict(figure='Fig6', panel=panel, row=cid, field=field, value=value, source_table=s,
                             specification_label=spec_label))
    flush_quiet()
    log('B4 detail: %d numbers of the 27 claims other than SFE-006 are copied verbatim (text and source table) from '
        'the released Source Data, each = one row of Supplemental Data 4/5/6' % n_same_as_released)

    # ---- the two block labels print claim counts; they are stored values too
    n_site = sum(1 for c in fin if panel_of[c] == 'a')
    n_prot = sum(1 for c in fin if panel_of[c] == 'b')
    gate((n_site, n_prot) == (sum(v for k, v in counts['by_level'].items() if k.startswith('site|')),
                              sum(v for k, v in counts['by_level'].items() if k.startswith('protein|')) - 1),
         'B6 block sizes = final_tally_counts by_level sums (protein minus the PERS-010 inset)', '%d, %d' % (n_site,
                                                                                                           n_prot))
    block_rows = {
        'a': dict(figure='Fig6', panel='a', row='site-level claims', field='n_claims_in_block', value=str(n_site),
                  source_table=TAB['SD2'] + ' (has_measured_baseline = 1) x ' + TAB['SD10'] + ' (unit = site)',
                  specification_label=''),
        'b': dict(figure='Fig6', panel='b', row='protein-level claims', field='n_claims_in_block', value=str(n_prot),
                  source_table=TAB['SD2'] + ' (has_measured_baseline = 1) x ' + TAB['SD10'] +
                  ' (unit = protein), minus PERS-010 (drawn in the inset)', specification_label='')}

    # ---- the tally written into the Source Data equals the final counts (overall, site, protein)
    ver = {r['row']: r['value'] for r in rows if r['field'] == 'verdict'}
    lev = {r['row']: r['value'] for r in rows if r['field'] == 'level'}
    for lvl in ('site', 'protein'):
        got = Counter('%s|%s' % (lvl, v) for c, v in ver.items() if lev[c] == lvl)
        want = {k: v for k, v in counts['by_level'].items() if k.startswith(lvl + '|')}
        gate(dict(got) == want, 'B7 Source Data verdicts = final_tally_counts by_level %s' % lvl,
             str(dict(sorted(got.items()))))
    gate('vanishes' not in ver.values(), 'B7 no vanishes class in the tally')

    # ---- order rows as they are drawn: panel a top to bottom, b top to bottom, c; field order within a claim kept
    matched = {r['row']: float(r['value']) for r in rows if r['field'] == 'matched_log2_or'}
    rank = {c: ('abc'.index(panel_of[c]), -matched[c]) for c in fin}
    rows.sort(key=lambda r: rank[r['row']])
    out = []
    for p in 'abc':
        if p in block_rows:
            out.append(block_rows[p])
        out += [r for r in rows if r['panel'] == p]
    rows = out
    gate(all(isinstance(r['value'], str) and all(ord(ch) < 128 for ch in r['value'] + r['source_table'] +
                                                   r['specification_label'])
             for r in rows), 'B8 Source Data is ASCII text', '%d rows' % len(rows))
    with open(TMP['sd'], 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['figure', 'panel', 'row', 'field', 'value', 'source_table',
                                           'specification_label'], lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    log('Source Data built: %d rows' % len(rows))
    return VC


# ================================================================================ 2. draw from the Source Data only
DARK, GREY = (31 / 255, 35 / 255, 38 / 255), (111 / 255, 117 / 255, 124 / 255)
W_IN = 7.2                        # exactly 518.4 pt
ROW_PT = 10.3                     # vertical pitch of one claim, as in the released figure
GAP = 1.4                         # extra rows between the blocks: room for the two block labels
FS_ROW, FS_TICK, FS_LAB, FS_LEG, FS_NOTE = 6.3, 6.5, 7.0, 5.9, 6.0

mpl.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'],
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold', 'mathtext.sf': 'Arial', 'mathtext.default': 'regular',
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'font.size': 7, 'axes.linewidth': 0.8,
    'axes.spines.right': False, 'axes.spines.top': False, 'legend.frameon': False, 'axes.unicode_minus': True,
})


def draw(VC):
    sdrows = read_csv(TMP['sd'])
    num, cat, panel, nblock = {}, {}, {}, {}
    for r in sdrows:
        if r['field'] == 'n_claims_in_block':
            nblock[r['panel']] = int(r['value'])
            continue
        panel[r['row']] = r['panel']
        if r['field'] in NUMERIC:
            num.setdefault(r['row'], {})[r['field']] = float(r['value'])
        else:
            cat.setdefault(r['row'], {})[r['field']] = r['value']
    drawn = []                                   # (claim, field, value) for every number put on the page

    def V(cid, field):
        v = num[cid][field]
        drawn.append((cid, field, v))
        return v

    col = {cid: VC[cat[cid]['verdict']] for cid in cat}

    site = sorted([c for c in num if panel[c] == 'a'], key=lambda c: -num[c]['matched_log2_or'])
    prot = sorted([c for c in num if panel[c] == 'b'], key=lambda c: -num[c]['matched_log2_or'])
    inset = [c for c in num if panel[c] == 'c']
    gate(len(site) == nblock['a'] == 18 and len(prot) == nblock['b'] == 9 and inset == ['PERS-010'],
         'G0 blocks 18 + 9 + inset PERS-010, equal to the stored block counts')
    y = {}
    for i, c in enumerate(reversed(prot)):
        y[c] = float(i)
    for i, c in enumerate(reversed(site)):
        y[c] = len(prot) - 1 + 1 + GAP + i
    top = max(y.values())
    split = (len(prot) - 1 + (len(prot) + GAP)) / 2.0
    ylim = (-0.9, top + 0.9)

    # ---- geometry in points
    probe = plt.figure(figsize=(W_IN, 3))
    r0 = probe.canvas.get_renderer()
    lab_w = max(Text(0, 0, cat[c]['row_label'], fontsize=FS_ROW, figure=probe).get_window_extent(r0).width
                for c in site + prot) * 72.0 / probe.dpi
    plt.close(probe)
    W_PT = W_IN * 72.0
    L = max(66.0, lab_w + 2.0 + 3.5 + 5.0)          # tick length + pad + margin
    R = 368.0                                       # right edge of the main axes (released figure: 368.1 pt)
    B, T = 38.0, 9.0
    HM = (ylim[1] - ylim[0]) * ROW_PT
    H_PT = B + HM + T
    fig = plt.figure(figsize=(W_IN, H_PT / 72.0))
    data_artists = []

    ax = fig.add_axes([L / W_PT, B / H_PT, (R - L) / W_PT, HM / H_PT])
    for c in site + prot:
        yy, cc = y[c], col[c]
        b, m = V(c, 'baseline_log2_or'), V(c, 'matched_log2_or')
        lo, hi = V(c, 'matched_ci_low'), V(c, 'matched_ci_high')
        data_artists += ax.plot([b, m], [yy, yy], color=cc, lw=1.7, alpha=0.45, solid_capstyle='butt', zorder=2)
        data_artists += ax.plot([lo, hi], [yy, yy], color=cc, lw=0.9, solid_capstyle='butt', zorder=3)
        data_artists.append(ax.scatter([m], [yy], s=24, color=cc, zorder=4, linewidths=0))
        data_artists.append(ax.scatter([b], [yy], s=15, facecolor='white', edgecolor=DARK, linewidths=0.7, zorder=5))
        if 'transfer_matched_log2_or' in num[c]:
            # SFE-006: tallied on its own data (above); its transfer re-test, a grey open diamond, just below the row
            yt = yy - 0.31
            tl, th = V(c, 'transfer_matched_ci_low'), V(c, 'transfer_matched_ci_high')
            tm = V(c, 'transfer_matched_log2_or')
            data_artists += ax.plot([tl, th], [yt, yt], color=GREY, lw=0.6, solid_capstyle='butt', zorder=3)
            data_artists.append(ax.scatter([tm], [yt], s=13, marker='D', facecolor='white', edgecolor=GREY,
                                           linewidths=0.7, zorder=5))
            ax.annotate('own data', xy=(max(b, m, hi), yy), xytext=(4.5, 0), textcoords='offset points',
                        fontsize=FS_LEG, color=GREY, va='center', ha='left')
            ax.annotate('transfer', xy=(tl, yt), xytext=(-3.0, 0), textcoords='offset points', fontsize=FS_LEG,
                        color=GREY, va='center', ha='right')
    data_artists.append(ax.axvline(0, color=DARK, lw=0.9, zorder=1))
    ticks = [y[c] for c in site + prot]
    ax.set_yticks(ticks)
    ax.set_yticklabels([cat[c]['row_label'] for c in site + prot], fontsize=FS_ROW, color=DARK)
    ax.set_xlabel('log$_2$ odds ratio', color=DARK, fontsize=FS_LAB)
    ax.set_xlim(-2.6, 6.6)
    ax.set_xticks([-2, 0, 2, 4, 6])
    ax.tick_params(axis='x', colors=GREY, labelsize=FS_TICK, length=2)
    ax.tick_params(axis='y', length=2, colors=DARK)
    ax.spines['left'].set_visible(False)
    ax.set_ylim(*ylim)

    # block divider and labels, in the gap between the blocks (clear of every row)
    bx, bw = 1.035, 0.020
    data_artists += ax.plot([0, bx + bw], [split, split], transform=ax.get_yaxis_transform(), color=GREY, lw=0.6,
                            ls=(0, (3, 2)), clip_on=False, zorder=1)
    ax.text(0.012, split + 0.14, 'site-level claims (%d)' % nblock['a'], transform=ax.get_yaxis_transform(),
            fontsize=FS_NOTE, color=GREY, va='bottom', ha='left')
    ax.text(0.012, split - 0.14, 'protein-level claims (%d)' % nblock['b'], transform=ax.get_yaxis_transform(),
            fontsize=FS_NOTE, color=GREY, va='top', ha='left')

    # verdict strip
    for c in site + prot:
        p = Rectangle((bx, y[c] - 0.34), bw, 0.68, transform=ax.get_yaxis_transform(), facecolor=col[c],
                      edgecolor='none', lw=0, clip_on=False)
        ax.add_patch(p)
        data_artists.append(p)
    ax.text(bx + bw / 2, top + 0.52, 'verdict', transform=ax.get_yaxis_transform(), fontsize=FS_LEG, color=GREY,
            ha='center', va='bottom')

    # ---- inset: PERS-010 on its own estimand (observed / expected shared), not comparable with the odds ratios
    c = inset[0]
    axc = fig.add_axes([412.0 / W_PT, (B + 7.0) / H_PT, 96.0 / W_PT, 38.0 / H_PT])
    b, m = V(c, 'baseline_log2_or'), V(c, 'matched_log2_or')
    lo, hi = V(c, 'matched_ci_low'), V(c, 'matched_ci_high')
    data_artists += axc.plot([b, m], [0, 0], color=col[c], lw=1.7, alpha=0.45, solid_capstyle='butt')
    data_artists += axc.plot([lo, hi], [0, 0], color=col[c], lw=0.9, solid_capstyle='butt')
    data_artists.append(axc.scatter([m], [0], s=24, color=col[c], zorder=4, linewidths=0))
    data_artists.append(axc.scatter([b], [0], s=15, facecolor='white', edgecolor=DARK, linewidths=0.7, zorder=5))
    axc.set_yticks([])
    axc.set_ylim(-0.6, 0.6)
    axc.set_xlim(0.76, 1.72)
    axc.set_xticks([0.8, 1.0, 1.2, 1.4, 1.6])
    axc.set_xlabel('log$_2$ (observed / expected shared)', fontsize=FS_NOTE, color=DARK)
    axc.tick_params(axis='x', colors=GREY, labelsize=FS_NOTE, length=2)
    axc.spines['left'].set_visible(False)
    axc.set_title(cat[c]['inset_title'], fontsize=FS_ROW, color=DARK, pad=3, loc='left')

    # ---- legend (the row flags are defined in the figure legend text, as in the released figure)
    used = {cat[c]['verdict'] for c in cat}
    gate(used == set(CLASS_ORDER), 'G0 colour classes in use = the five final classes', ', '.join(CLASS_ORDER))
    handles = [Line2D([], [], color=GREY, lw=0.9, marker='o', ms=4, mfc=GREY, mec=GREY, mew=0,
                      label='after matching, with its 95% interval'),
               Line2D([], [], ls='', marker='o', ms=3.4, mfc='white', mec=DARK, mew=0.7,
                      label='reconstructed baseline'),
               Line2D([], [], color=GREY, lw=0.6, marker='D', ms=3.0, mfc='white', mec=GREY, mew=0.7,
                      label='SFE-006 transfer re-test,\nafter matching (t)')]
    handles += [Line2D([], [], ls='', marker='s', ms=4, color=VC[k], label=CLASS_LABEL[k]) for k in CLASS_ORDER]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.075, 1.0), fontsize=FS_LEG, labelcolor=DARK,
              handletextpad=0.45, borderpad=0.2, labelspacing=0.30, handlelength=1.5)
    return fig, drawn, data_artists, num


# ================================================================================ 3. gates on the drawing
def check_drawing(fig, drawn, data_artists, num):
    for cid, field, v in drawn:
        if not (field in num[cid] and num[cid][field] == v):
            raise SystemExit('GATE FAILED  G1 drawn value not in Source Data: %s %s' % (cid, field))
    want = {c: {'baseline_log2_or', 'matched_log2_or', 'matched_ci_low', 'matched_ci_high'} for c in num}
    want['SFE-006'] |= {'transfer_matched_log2_or', 'transfer_matched_ci_low', 'transfer_matched_ci_high'}
    got = {}
    for cid, field, v in drawn:
        got.setdefault(cid, set()).add(field)
    gate(got == want, 'G1 every drawn number read from Source Data (%d numbers, 28 claims)' % len(drawn))

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    Wpx, Hpx = fig.bbox.width, fig.bbox.height
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    boxes = [(t, t.get_window_extent(r)) for t in texts]
    tol = 0.3 * fig.dpi / 72.0
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if (min(a.x1, b.x1) - max(a.x0, b.x0) > tol) and (min(a.y1, b.y1) - max(a.y0, b.y0) > tol):
                bad.append((boxes[i][0].get_text(), boxes[j][0].get_text()))
    gate(not bad, 'G2 no text overlaps text (%d text items)' % len(boxes), str(bad[:5]))
    bad = []
    for t, tb in boxes:
        for art in data_artists:
            db = art.get_window_extent(r)
            if isinstance(art, Line2D):
                db = db.padded(art.get_linewidth() * fig.dpi / 72.0 / 2)
            if (min(tb.x1, db.x1) - max(tb.x0, db.x0) > tol) and (min(tb.y1, db.y1) - max(tb.y0, db.y0) > tol):
                bad.append(t.get_text())
    gate(not bad, 'G3 no text overlaps a data mark, interval, strip cell or reference line', str(bad[:5]))
    out = [t.get_text() for t, tb in boxes if tb.x0 < 0 or tb.y0 < 0 or tb.x1 > Wpx or tb.y1 > Hpx]
    gate(not out, 'G4 all text inside the canvas', str(out))
    small = sorted({t.get_fontsize() for t, _ in boxes})
    gate(min(small) >= 5.7, 'G5 text sizes (pt) >= 5.7', str(small))


def check_pdf():
    import fitz
    d = fitz.open(TMP['pdf'])
    p = d[0]
    box = d.xref_get_key(p.xref, 'MediaBox')[1].strip('[]').split()      # exact text written in the PDF
    gate(len(d) == 1 and box[:3] == ['0', '0', '518.4'] and float(box[3]) <= 432.0,
         'G6 page 518.4 pt (7.2 in) wide, height <= 432 pt (6 in)', 'MediaBox [%s]' % ' '.join(box))
    fonts = p.get_fonts()
    gate(fonts and all(f[1] == 'ttf' and 'Arial' in f[3] for f in fonts),
         'G7 fonts embedded as TrueType (fonttype 42), Arial only',
         ', '.join(sorted({'%s (%s)' % (f[3], f[2]) for f in fonts})))
    sizes, sub = Counter(), []
    for blk in p.get_text('dict')['blocks']:
        for ln in blk.get('lines', []):
            for s in ln['spans']:
                sizes[round(s['size'], 2)] += 1
                if s['size'] < 5.7 - 1e-6:
                    sub.append(s['text'])
    gate(all(x.strip() == '2' for x in sub), 'G7 only the log2 subscripts are below 5.7 pt',
         'sizes %s; below 5.7: %s' % (dict(sorted(sizes.items())), sub))


def main():
    try:
        for k, p in IN.items():
            log('input %-20s sha256 %s  %s' % (k, sha(p), p))
        VC = build_source_data()
        fig, drawn, data_artists, num = draw(VC)
        check_drawing(fig, drawn, data_artists, num)
        fig.savefig(TMP['pdf'], format='pdf', metadata={'CreationDate': None, 'ModDate': None})
        fig.savefig(TMP['png'], format='png', dpi=200, metadata={'Software': None})
        plt.close(fig)
        check_pdf()
    except BaseException:
        for p in TMP.values():
            if os.path.exists(p):
                os.remove(p)
        raise
    for k in ('sd', 'pdf', 'png'):                  # every gate passed: publish
        os.replace(TMP[k], OUT[k])
        log('output sha256 %s  %s' % (sha(OUT[k]), os.path.basename(OUT[k])))
    log('script sha256 %s  %s' % (sha(os.path.abspath(__file__)), os.path.basename(__file__)))
    with open(OUT_LOG, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(LOG) + '\n')


if __name__ == '__main__':
    main()
