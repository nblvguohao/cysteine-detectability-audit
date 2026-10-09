"""Independent re-read: every numeric Source Data value of the new figures against a second copy of the number.

The second copy is read with different code (csv/json/markdown parsing written here, not common.py) from the table
that the manuscript cites for it (Supplemental Notes 4, 7, 12, 14, 20, 22; the positional-profile CSV; the matched JSONs;
the pLMSNOSite result CSVs), compared at the precision that table gives. Output: qa_logs/verify_values.csv and a
per-figure summary. Also lists every number printed on each new figure (PDF text layer) with the Source Data row(s)
it matches (printed_numbers.csv). The tables committed under qa/ were written by
    python3 posthoc_2026-10-09_figures/verify_values.py figures_final source_data_final posthoc_2026-10-09_figures/qa
run from the repository root.
Usage: python3 verify_values.py [figures_dir source_data_dir out_dir]   (defaults: figures_final, source_data_final,
       posthoc_2026-10-09_figures/build_output/qa_logs)
"""
import csv
import json
import pathlib
import re
import sys

import fitz

REPO = pathlib.Path(__file__).resolve().parent.parent      # repository root
SUP = REPO / 'supplemental'
SD32 = REPO / 'source_data_submitted'
FIG = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / 'figures_final'
SDD = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else REPO / 'source_data_final'
OUT = pathlib.Path(sys.argv[3]) if len(sys.argv) > 3 else REPO / 'posthoc_2026-10-09_figures' / 'build_output' / 'qa_logs'


def tables(path):
    """all markdown tables of a file as lists of dicts"""
    out, cur, hdr = [], None, None
    for ln in pathlib.Path(path).read_text(encoding='utf-8').splitlines():
        if ln.startswith('|'):
            cells = [c.strip().strip('`') for c in ln.strip().strip('|').split('|')]
            if cur is None:
                hdr, cur = cells, []
            elif set(''.join(cells)) <= set('-: '):
                continue
            else:
                cur.append(dict(zip(hdr, cells)))
        elif cur is not None:
            out.append(cur)
            cur = None
    if cur is not None:
        out.append(cur)
    return out


def f(s):
    return float(s.replace('−', '-').replace(',', '').replace('ᵈ', '').strip())


def ci(s):
    m = re.findall(r'[-−]?\d+\.?\d*', s.replace('ᵈ', ''))
    return [f(x) for x in m[:3]]


def sd(name):
    return list(csv.DictReader(open(SDD / name, encoding='utf-8')))


RES = []


def check(fig, panel, row, field, sdval, ref, ref_src, tol):
    ok = abs(float(sdval) - float(ref)) <= tol
    RES.append(dict(figure=fig, panel=panel, row=row, field=field, source_data_value=sdval, independent_value=ref,
                    independent_source=ref_src, tolerance=tol, match='yes' if ok else 'NO'))


# ------------------------------------------------------------------------------------------- Fig 2 vs Note 20 S20.1
names = {'qtrp_S1_ph5': 'QTRP S1 pH5', 'qtrp_S2_ph5': 'QTRP S2 pH5', 'qpers_sid_tierB': 'qPerS-SID tier B',
         'cysboost2019_human_sno_hela': 'Cys-BOOST HeLa', 'cysboost2019_human_sno_shsy5y': 'Cys-BOOST SH-SY5Y',
         'natcomm2023_ath_sno': 'FAT-switch', 'abiotech2025_ath_sno': 'PAT-switch', 'fps2020_ath_sulfenyl': 'YAP1C reporter'}
s201 = [t for t in tables(SUP / 'Supplemental_Note_20_fig5a_absolute_discrimination.md') if 'VIS10 global' in t[0]][0]
col = {'auc_global_VIS10': 'VIS10 global', 'auc_within_protein_VIS10': 'VIS10 within protein',
       'auc_within_protein_DIG25': 'DIG25 within protein'}
for r in sd('Source_Data_Fig2_negative_sets.csv'):
    m = re.match(r'(auc_\w+?)(_lo|_hi)?$', r['field'])
    if not m or '|' not in r['row']:
        continue
    cid, neg = r['row'].split('|')
    t = [x for x in s201 if x['Cohort'] == names[cid] and x['Neg. set'] == neg][0]
    k = {None: 0, '_lo': 1, '_hi': 2}[m.group(2)]
    check('Fig2', r['panel'], r['row'], r['field'], r['value'], ci(t[col[m.group(1)]])[k],
          'Note 20 Table S20.1', 0.00051)

# ------------------------------------------------------------------------------------------- Fig 3
pos = {}
with open(REPO / 'artifacts12/positional_profile/results/positional_kr_profile_public_2026-09-17.csv') as fh:
    for r in csv.DictReader(fh):
        if r['dataset_id'] == 'natcomm2023_ath_sno':
            pos[int(r['offset'])] = r
mj = {k: json.load(open(REPO / ('artifacts12/results/12_matched_%s.json' % k), encoding='utf-8')) for k in ('sno', 'so', 'ox')}
n4 = [t for t in tables(SUP / 'Supplemental_Note_4_four_protease_public_deposit.md') if 'distal 6-12' in t[0]][0]
n4n = [t for t in tables(SUP / 'Supplemental_Note_4_four_protease_public_deposit.md') if 'identified cysteines (+hydroxylamine)' in t[0]][0]
arm = {'AspN': 'AspN', 'CT': 'chymotrypsin', 'GluC': 'GluC', 'Trypsin': 'trypsin'}
for r in sd('Source_Data_Fig3_cleavage_geometry.csv'):
    if r['panel'] == 'a' and r['row'].startswith('offset') and r['field'].endswith('_printed'):
        o = int(r['row'].split()[1])
        check('Fig3', 'a', r['row'], r['field'], r['value'], round(float(pos[o]['ratio_vs_proteome']), 2), 'positional CSV, rounded to 2 dp', 1e-9)
    elif r['panel'] == 'a' and r['row'].startswith('offset'):
        o = int(r['row'].split()[1])
        check('Fig3', 'a', r['row'], r['field'], r['value'], pos[o][r['field']], 'positional CSV (csv module)', 1e-9)
    elif r['panel'] == 'b':
        k = {'S-nitrosylation': 'sno', 'S-sulfenylation': 'so', 'Reversible oxidation': 'ox'}[r['row'].split(',')[0]]
        cls = '碱性 KRH' if 'K/R/H' in r['row'] else '酸性 DE'
        st = r['field'].split('_')[1]
        check('Fig3', 'b', r['row'], r['field'], r['value'], mj[k]['classes'][cls][st]['z'], '12_matched JSON', 1e-9)
    elif r['panel'] in ('c', 'd') and '|' in r['row']:
        code, rule, bg, _ = r['row'].split('|')
        t = [x for x in n4 if x['arm'] == arm[code] and x['rule'] == rule and x['background'] == bg][0]
        k = {'estimate': 0, 'ci_low': 1, 'ci_high': 2}[r['field']]
        check('Fig3', r['panel'], r['row'], r['field'], r['value'], ci(t['distal 6-12'])[k], 'Note 4 cleavage table', 1e-6)
    elif r['field'] == 'identified_cysteines':
        nm = {'AspN': 'AspN', 'Chymotrypsin': 'chymotrypsin', 'GluC': 'GluC', 'Trypsin': 'trypsin'}[r['row']]
        t = [x for x in n4n if x['arm'] == nm][0]
        check('Fig3', 'c', r['row'], r['field'], r['value'], t['identified cysteines (+hydroxylamine)'], 'Note 4 input table', 0)

# ------------------------------------------------------------------------------------------- Fig 4
n7 = [t for t in tables(SUP / 'Supplemental_Note_7_artefact4_public_deposit.md') if 'wilson_low' in t[0]][0]
n7p = [x for x in n7 if x['definition'] == 'PRIMARY_leading_both_sides' and x['reading'] == 'CAM_only'][0]
n4c = [t for t in tables(SUP / 'Supplemental_Note_4_four_protease_public_deposit.md') if 'coincidence' in t[0]][0]
coh = {r['dataset_id']: r for r in csv.DictReader(open(SD32 / 'Source_Data_text_artefact4_public_cohorts.csv'))}
fig3_32 = list(csv.DictReader(open(SD32 / 'Source_Data_Fig3_search_space.csv', encoding='utf-8')))
for r in sd('Source_Data_Fig4_removed_comparisons.csv'):
    if r['panel'] == 'a' and r['row'] == 'PXD048216':
        key = {'share': 'p', 'wilson_low': 'wilson_low', 'wilson_high': 'wilson_high', 'identified': 'observed_cysteines',
               'with_site': 'matched'}[r['field']]
        check('Fig4', 'a', r['row'], r['field'], r['value'], n7p[key], 'Note 7 coincidence table', 1e-9)
    elif r['panel'] == 'a' and r['row'].startswith('PXD063463'):
        nm = r['row'].split()[1]
        t = [x for x in n4c if x['arm'] == nm][0]
        v = ci(t['coincidence'] + ' ' + t['95 per cent interval'])
        k = {'estimate': 0, 'ci_low': 1, 'ci_high': 2}[r['field']]
        check('Fig4', 'a', r['row'], r['field'], r['value'], v[k], 'Note 4 coincidence table', 0.00051)
    elif r['panel'] == 'a' and r['row'] in coh:
        check('Fig4', 'a', r['row'], r['field'], r['value'], coh[r['row']][r['field']], 'artefact4 cohorts CSV', 1e-9)
    elif r['panel'] in ('b', 'c') and r['row'] != 'printed constants' and r['field'] not in (
            'mass_limit_da_printed', 'percent_printed'):
        p = 'a' if r['panel'] == 'b' else 'd'
        h = [x for x in fig3_32 if x['panel'] == p and x['row'] == r['row'] and x['field'] == r['field']]
        if len(h) == 1:
            check('Fig4', r['panel'], r['row'], r['field'], r['value'], h[0]['value'], 'v3.2 Source Data Fig 3 (verbatim)', 1e-9)

# ------------------------------------------------------------------------------------------- Fig 5 / S4 vs Note 14
n14 = tables(SUP / 'Supplemental_Note_14_claim_retest_tally.md')
s141 = {x['Claim']: x for x in [t for t in n14 if 'Random control [95% CI]' in t[0]][0]}
s142 = {x['Claim']: x for x in [t for t in n14 if 'Matched / baseline' in t[0]][0]}
tally = [t for t in n14 if 'null broken by control' in t[0]][0]
cols = {'baseline': 'Baseline [95% CI]', 'matched': 'Matched [95% CI]', 'random_control': 'Random control [95% CI]'}
for fname, fid in (('Source_Data_Fig5_claim_retests.csv', 'Fig5'), ('Source_Data_FigS4_protein_level_claims.csv', 'FigS4')):
    for r in sd(fname):
        m = re.match(r'(baseline|matched|random_control)_(log2_or|ci_low|ci_high)$', r['field'])
        if m and r['row'] in s141 and not (fid == 'Fig5' and r['panel'] == 'b'):
            v = ci(s141[r['row']][cols[m.group(1)]])[{'log2_or': 0, 'ci_low': 1, 'ci_high': 2}[m.group(2)]]
            check(fid, r['panel'], r['row'], r['field'], r['value'], v, 'Note 14 Table S14.1 (4 dp)', 0.00006)
        elif fid == 'Fig5' and r['panel'] == 'b' and r['row'] in s142 and r['field'] != 'attribute_class':
            k = {'baseline_log2_or': 'Baseline', 'matched_log2_or': 'Matched', 'matched_over_baseline': 'Matched / baseline',
                 'random_control_log2_or': 'Random control', 'random_over_baseline': 'Random / baseline'}[r['field']]
            check(fid, 'b', r['row'], r['field'], r['value'], f(s142[r['row']][k]), 'Note 14 Table S14.2', 1e-9)
        elif fid == 'Fig5' and r['panel'] == 'c' and r['field'] == 'site_level_claims' and r['row'] != 'All re-tested claims':
            site = [x for x in tally if x['Claims'] == 'Site level'][0]
            check(fid, 'c', r['row'], r['field'], r['value'], site[r['row'].lower().replace('null broken by control',
                  'null broken by control')], 'Note 14 Tally table (site level)', 0)
        elif fid == 'Fig5' and r['panel'] == 'c' and r['row'] == 'All re-tested claims':
            lvl = 'Site level' if r['field'] == 'site_level_claims' else 'Protein level'
            check(fid, 'c', r['row'], r['field'], r['value'], [x for x in tally if x['Claims'] == lvl][0]['total'],
                  'Note 14 Tally table (total)', 0)

# ------------------------------------------------------------------------------------------- Fig 6
n22 = {x['Model or comparison']: x for x in [t for t in tables(SUP / 'Supplemental_Note_22_plmsnosite_rescored.md')
                                            if 'Model or comparison' in t[0]][0]}
row22 = {'pLMSNOSite': 'pLMSNOSite, released models (probabilities)', 'VIS10': 'VIS10', 'DIG25': 'DIG25',
         'pLMSNOSite_call_0.5': 'pLMSNOSite, calls thresholded at 0.5ᶜ'}
n12 = [t for t in tables(SUP / 'Supplemental_Note_12_public_detectability_only.md') if 'recovery_ratio' in t[0]][0]
f7 = list(csv.DictReader(open(SD32 / 'Source_Data_Fig7_self_audit_public.csv', encoding='utf-8')))
for r in sd('Source_Data_Fig6_predictors.csv'):
    if r['panel'] == 'a' and r['row'] in row22:
        if r['field'].startswith('all_sites_auroc'):
            k = {'all_sites_auroc': 0, 'all_sites_auroc_ci_low': 1, 'all_sites_auroc_ci_high': 2}[r['field']]
            check('Fig6', 'a', r['row'], r['field'], r['value'], ci(n22[row22[r['row']]]['AUROC, all test sites'])[k],
                  'Note 22 Table S22.1 (3 dp)', 0.00051)
        else:
            k = {'within_protein_auroc_value': 0, 'within_protein_auroc_ci_low': 1, 'within_protein_auroc_ci_high': 2}[r['field']]
            check('Fig6', 'a', r['row'], r['field'], r['value'],
                  ci(n22[row22[r['row']]]['Within-protein AUROC, mean of per-protein AUROCsᵃ'])[k], 'Note 22 Table S22.1', 0.00051)
    elif r['panel'] == 'a' and ' | ' in r['row'] and r['field'] in ('difference', 'ci_low', 'ci_high'):
        scope, comp = r['row'].split(' | ')
        rr = n22['pLMSNOSite − %s, paired difference' % comp.split(' - ')[1]]
        c = 'AUROC, all test sites' if scope == 'All sites' else 'Within-protein AUROC, mean of per-protein AUROCsᵃ'
        check('Fig6', 'a', r['row'], r['field'], r['value'], ci(rr[c])[{'difference': 0, 'ci_low': 1, 'ci_high': 2}[r['field']]],
              'Note 22 Table S22.1', 0.00051)
    elif r['panel'] == 'b' and r['row'].endswith(('PXD044043', 'PXD072089')):
        t = [x for x in n12 if x['cohort'] == r['row']][0]
        k = {'within_protein_auc_detect_only': 'wp_auc_detect_only', 'within_protein_auc_full': 'wp_auc_full',
             'recovery_share': 'recovery_ratio', 'recovery_share_ci_low': 'ci_low', 'recovery_share_ci_high': 'ci_high'}[r['field']]
        check('Fig6', 'b', r['row'], r['field'], r['value'], t[k], 'Note 12 table (LightGBM ranking member)', 1e-9)

# ------------------------------------------------------------------------------------------- Fig 1c, S7 vs v3.2 Source Data
f1 = list(csv.DictReader(open(SD32 / 'Source_Data_Fig1_overview.csv', encoding='utf-8')))
for fname, fid, p in (('Source_Data_Fig1_overview.csv', 'Fig1', 'c'), ('Source_Data_FigS7_literature_survey.csv', 'FigS7', 'b')):
    for r in sd(fname):
        h = [x for x in f1 if x['panel'] == p and x['row'] == r['row']]
        if h and r['field'] in ('classified_analyses', 'records'):
            check(fid, r['panel'], r['row'], r['field'], r['value'], h[0]['value'], 'v3.2 Source Data Fig 1%s' % p, 0)

OUT.mkdir(parents=True, exist_ok=True)
with open(OUT / 'verify_values.csv', 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=list(RES[0].keys()), lineterminator='\n')
    w.writeheader()
    w.writerows(RES)
summ = {}
for r in RES:
    s = summ.setdefault(r['figure'], [0, 0])
    s[0] += 1
    s[1] += r['match'] == 'NO'
for k, (n, bad) in summ.items():
    print('%-6s values checked %4d  mismatches %d' % (k, n, bad))
for r in RES:
    if r['match'] == 'NO':
        print('MISMATCH', r)

# ------------------------------------------------------------------------------------------- printed numbers
rows_out = []
for stem in ('Fig1_overview', 'Fig2_negative_sets', 'Fig3_cleavage_geometry', 'Fig4_removed_comparisons',
             'Fig5_claim_retests', 'Fig6_predictors', 'FigS4_protein_level_claims', 'FigS7_literature_survey'):
    page = fitz.open(FIG / (stem + '.pdf'))[0]
    sdr = sd('Source_Data_%s.csv' % stem)
    for b in page.get_text('dict')['blocks']:
        for ln in b.get('lines', []):
            txt = ''.join(s['text'] for s in ln['spans'])
            for tok in re.findall(r'(?<![A-Za-z\d.])[-−]?\d[\d,]*\.?\d*', txt):
                t = tok.replace('−', '-').replace(',', '').rstrip('.')
                dec = len(t.split('.')[1]) if '.' in t else 0
                v = float(t)
                hits = []
                for r in sdr:
                    try:
                        x = float(r['value'].replace(',', ''))
                    except ValueError:
                        continue
                    if abs(round(x, dec) - v) < 10 ** (-dec) * 0.51:
                        hits.append('%s/%s/%s' % (r['panel'], r['row'], r['field']))
                rows_out.append(dict(figure=stem, printed_text=txt.strip(), token=tok, source_data_rows='; '.join(hits[:3]) + (
                    ' (+%d more)' % (len(hits) - 3) if len(hits) > 3 else ''), matched='yes' if hits else 'no'))
with open(OUT / 'printed_numbers.csv', 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()), lineterminator='\n')
    w.writeheader()
    w.writerows(rows_out)
print('printed numeric tokens:', len(rows_out), 'unmatched:', sum(r['matched'] == 'no' for r in rows_out))
