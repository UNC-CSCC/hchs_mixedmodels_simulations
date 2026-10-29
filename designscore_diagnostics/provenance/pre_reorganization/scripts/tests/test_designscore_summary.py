"""Check that missing fits and missing intervals remain visible in summaries."""
import csv
import importlib.util
import io
from contextlib import redirect_stdout
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


class SummaryDenominators(unittest.TestCase):
    def test_boundary_and_missing_fit_use_distinct_denominators(self):
        script = Path(__file__).parents[1] / 'evaluation/summarize_designscore.py'
        spec = importlib.util.spec_from_file_location('designscore_summary', script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            truth = root / 'data/population_coefficients/population_coef_oct2026_hhsub.csv'
            truth.parent.mkdir(parents=True)
            truth.write_text('term,estimate\nvar_subid,175\n')
            for name in ('scen3_designscore.do', 'variance_designscore.mata'):
                path = root / 'scripts/stata' / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('test fixture')
            output = root / 'data/estimates/oct2026/scen3_designscore/exact'
            output.mkdir(parents=True)
            for sim, point, usable in ((1, 170, 1), (2, 180, 0)):
                common = dict(sim_id=sim, population='oct2026', probability_mode='exact',
                    fit_status='OK' if usable else 'BOUNDARY_POINT_ONLY',
                    se_status=0 if usable else 6)
                coefficient = dict(term='var_subid', estimate=point, std_error=4 if usable else '',
                    conf_low=160 if usable else '', conf_high=180 if usable else '',
                    inference_usable=usable, **common)
                audit = dict(**common, boundary_mode=0 if usable else 1, kkt_error=1e-9,
                    step1_score_error=1e-7, jacobian_condition=100)
                module.write(output / f'scen3_designscore_coef_{sim}.csv', [coefficient])
                module.write(output / f'scen3_designscore_audit_{sim}.csv', [audit])
            with patch('sys.argv', ['summary', '--repo', str(root), '--mode', 'exact', '--last', '3']):
                with redirect_stdout(io.StringIO()):
                    module.main()
            folder = root / 'data/summaries/oct2026'
            row = module.read(folder / 'coverage_scen3_designscore_exact.csv')[0]
            self.assertEqual(int(row['n_expected']), 3)
            self.assertEqual(int(row['n_reps']), 2)
            self.assertEqual(int(row['n_valid_intervals']), 1)
            self.assertEqual(int(row['n_missing_points']), 1)
            self.assertEqual(int(row['n_missing_intervals']), 2)
            self.assertEqual(float(row['empirical_bias']), 0)
            self.assertEqual(float(row['coverage_rate']), 1)
            self.assertAlmostEqual(float(row['coverage_over_all_requested']), 1/3)
            audit = module.read(folder / 'audit_scen3_designscore_exact.csv')
            self.assertEqual(audit[-1]['fit_status'], 'MISSING_OUTPUT')


if __name__ == '__main__':
    unittest.main()
