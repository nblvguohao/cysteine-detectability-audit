import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
import joblib

HERE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('minimal_adapter', HERE/'run.py')
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        rows = []
        for i in range(10):
            for position, positive, ref_b in [(4, 1, 0), (8, 0, 1), (12, 0, 0)]:
                rows.append(dict(acc='synthetic_'+str(i), site=position, unit='group_'+str(i),
                                 is_pos=positive, in_A=1-positive, in_B=ref_b))
        self.sites = pd.DataFrame(rows)
        self.sites_path = self.root/'sites.tsv'
        self.sites.to_csv(self.sites_path, sep='\t', index=False)
        self.folds = self.sites.copy()
        self.folds['new_fold'] = np.repeat(np.arange(10) % 5, 3)
        # These unused columns must not enter feature construction or fitting.
        self.folds['new_score_A'] = np.nan
        self.folds['new_score_B'] = 'ignored'
        self.fold_path = self.root/'folds.tsv'
        self.folds.to_csv(self.fold_path, sep='\t', index=False)
        self.fasta_path = self.root/'reference.fasta'
        self.fasta_path.write_text(''.join('>synthetic_'+str(i)+'\nAAACAAACAAACAAA\n' for i in range(10)))
        self.output = self.root/'new_output'

    def run_adapter(self):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            adapter.main(['--sites', str(self.sites_path), '--fasta', str(self.fasta_path),
                          '--folds', str(self.fold_path), '--output', str(self.output),
                          '--definitions', 'strict_mc_v1', '--bootstrap-reps', '20'])

    def test_fit_predictions_and_fixed_folds(self):
        self.run_adapter()
        folder = self.output/'results/strict_mc_v1/provided/VIS10'
        features = pd.read_csv(folder/'features.tsv.gz', sep='\t', float_precision='round_trip').to_numpy()
        oof = pd.read_csv(folder/'new_oof.tsv.gz', sep='\t', float_precision='round_trip')
        np.testing.assert_array_equal(oof.new_fold, self.folds.new_fold)
        for fold in range(5):
            for background in ['A', 'B']:
                model = joblib.load(folder/f'fold_{fold}_{background}.joblib')
                saved = np.load(folder/f'fold_{fold}_{background}.npz')
                scores = model.predict_proba(features[saved['test_rows']])[:, 1]
                np.testing.assert_array_equal(scores, saved['scores'])
                np.testing.assert_array_equal(scores, oof['new_score_'+background].to_numpy()[saved['test_rows']])
        consumed = pd.read_csv(self.output/'inputs/fixed_folds.tsv', sep='\t')
        self.assertEqual(list(consumed.columns), list(self.sites.columns)+['new_fold'])
        self.assertEqual(json.loads((self.output/'logs/run_complete.json').read_text())['models'], 10)

    def test_reject_group_crossing_folds(self):
        self.folds.loc[0, 'new_fold'] = 1
        self.folds.to_csv(self.fold_path, sep='\t', index=False)
        with self.assertRaisesRegex(AssertionError, 'unit crosses'):
            self.run_adapter()
        self.assertFalse(self.output.exists())

    def test_reject_duplicate_sites(self):
        pd.concat([self.sites, self.sites.iloc[:1]]).to_csv(self.sites_path, sep='\t', index=False)
        with self.assertRaisesRegex(AssertionError, 'Duplicate'):
            self.run_adapter()
        self.assertFalse(self.output.exists())

    def test_reject_non_cysteine(self):
        self.sites.loc[0, 'site'] = 1
        self.sites.to_csv(self.sites_path, sep='\t', index=False)
        with self.assertRaisesRegex(AssertionError, 'not cysteine'):
            self.run_adapter()
        self.assertFalse(self.output.exists())

    def test_preserve_existing_output(self):
        self.output.mkdir()
        marker = self.output/'marker'
        marker.write_text('preserve')
        with self.assertRaises(SystemExit):
            self.run_adapter()
        self.assertEqual(marker.read_text(), 'preserve')


if __name__ == '__main__':
    unittest.main()
