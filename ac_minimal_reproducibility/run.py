#!/usr/bin/env python3
"""Fit the frozen core estimator on separately supplied sites, FASTA and folds."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sites', type=Path, required=True)
    parser.add_argument('--fasta', type=Path, required=True)
    parser.add_argument('--folds', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='New output directory')
    parser.add_argument('--definitions', nargs='+', choices=['strict_mc_v1', 'probe_compatible_v1'],
                        default=['strict_mc_v1', 'probe_compatible_v1'])
    parser.add_argument('--feature-sets', nargs='+', choices=['VIS10', 'DIG25', 'CPL15'], default=['VIS10'])
    parser.add_argument('--bootstrap-reps', type=int, default=2000)
    args = parser.parse_args(argv)
    if args.bootstrap_reps < 20:
        parser.error('--bootstrap-reps must be at least 20; use 2000 for the scientific analysis')
    if args.output.exists():
        parser.error('Output already exists; use a fresh directory to preserve previous evidence')
    for path in (args.sites, args.fasta, args.folds):
        if not path.is_file():
            parser.error('Missing separately supplied input: ' + str(path))
    for key in ['OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        os.environ[key] = '1'
    os.environ['OMP_NUM_THREADS'] = '2'
    os.environ['CUDA_VISIBLE_DEVICES'] = ''
    sys.path.insert(0, str(HERE/'code'))
    import numpy as np
    import pandas as pd
    from threadpoolctl import threadpool_limits
    import run_closure as closure
    import check_boundaries
    from prepare_inputs import fasta

    columns = ['acc', 'site', 'unit', 'is_pos', 'in_A', 'in_B']
    sites = pd.read_csv(args.sites, sep='\t', float_precision='round_trip')
    assert list(sites.columns) == columns, 'Sites must contain exactly the documented six columns'
    assert len(sites) and not sites.isna().any().any(), 'Empty or missing site identities/memberships'
    assert not sites.duplicated(['acc', 'site']).any(), 'Duplicate accession-position identities'
    assert all(sites[c].isin([0, 1]).all() for c in ['is_pos', 'in_A', 'in_B']), 'Non-binary membership'
    assert ((sites.is_pos + sites.in_A + sites.in_B) > 0).all(), 'Row outside all comparison sets'
    assert not ((sites.is_pos == 1) & ((sites.in_A == 1) | (sites.in_B == 1))).any(), 'Positive/reference overlap'
    assert np.isfinite(sites.site).all() and (sites.site == sites.site.astype(int)).all(), 'Invalid positions'
    seqs = fasta(args.fasta)
    for row in sites.itertuples():
        assert row.acc in seqs and 1 <= row.site <= len(seqs[row.acc]), 'Unmapped position'
        assert seqs[row.acc][int(row.site)-1] == 'C', 'Position is not cysteine'
    # Read identities and fixed folds only. Prediction/feature columns are never consumed.
    table = pd.read_csv(args.folds, sep='\t', usecols=columns+['new_fold'], float_precision='round_trip')
    assert sites.equals(table[columns]), 'Fold-table row identities/memberships/order differ from sites'
    assert table.new_fold.isin(range(5)).all() and set(table.new_fold) == set(range(5)), 'Expected five fixed folds'
    assert table.groupby('unit').new_fold.nunique().max() == 1, 'A unit crosses held-out folds'
    folds = table.new_fold.to_numpy()
    y = sites.is_pos.to_numpy()
    for k in range(5):
        train, test = np.flatnonzero(folds != k), np.flatnonzero(folds == k)
        assert not set(sites.unit.iloc[train]) & set(sites.unit.iloc[test]), 'Train/test group overlap'
        for b in ['A', 'B']:
            selected = train[(y[train] == 1) | sites['in_'+b].to_numpy()[train].astype(bool)]
            assert set(y[selected]) == {0, 1}, 'Training fold lacks a class'
    for b in ['A', 'B']:
        assert set(y[(y == 1) | sites['in_'+b].to_numpy().astype(bool)]) == {0, 1}, 'Evaluation lacks a class'

    out = args.output.resolve()
    inp = out/'inputs/prepared/provided'
    inp.mkdir(parents=True)
    for name in ['logs', 'results']:
        (out/name).mkdir()
    shutil.copy2(args.sites, inp/'sites.tsv')
    shutil.copy2(args.fasta, inp/'reference.fasta')
    # Scores, if present in the supplied folds file, are omitted from this saved input.
    table[columns+['new_fold']].to_csv(out/'inputs/fixed_folds.tsv', sep='\t', index=False)
    closure.ROOT = out
    check_boundaries.ROOT = out

    class FixedFolds:
        def __init__(self, n_splits=5):
            assert n_splits == 5

        def split(self, X, labels, groups):
            assert len(X) == len(sites) and np.array_equal(labels, y)
            assert np.array_equal(groups, sites.unit.to_numpy())
            for k in range(5):
                yield np.flatnonzero(folds != k), np.flatnonzero(folds == k)

    closure.GroupKFold = FixedFolds
    sources = [HERE/'run.py', HERE/'requirements.txt', *sorted((HERE/'code').glob('*.py'))]
    record = dict(adapter='external_inputs_fixed_folds_v1',
                  code={str(p.relative_to(HERE)): sha(p) for p in sources},
                  inputs={name: sha(out/name) for name in ['inputs/prepared/provided/sites.tsv',
                          'inputs/prepared/provided/reference.fasta', 'inputs/fixed_folds.tsv']},
                  environment=closure.environment(), definitions=args.definitions,
                  feature_sets=args.feature_sets, bootstrap_reps=args.bootstrap_reps)
    identity = hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()
    closure.write_json(out/'logs/execution_identity.json', dict(identity=identity, **record))
    with threadpool_limits(limits=1, user_api='blas'):
        check_boundaries.check()
        for definition in args.definitions:
            for feature_set in args.feature_sets:
                closure.run_task('provided', definition, feature_set, identity, args.bootstrap_reps)
    closure.write_json(out/'logs/run_complete.json', dict(identity=identity,
        tasks=len(args.definitions)*len(args.feature_sets),
        models=10*len(args.definitions)*len(args.feature_sets)))
    print('Outputs:', out)


if __name__ == '__main__':
    main()
