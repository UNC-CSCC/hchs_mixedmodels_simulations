"""Non-model tests: diagnostics, run guards and shared pooling formulas.

No Stata fit, sampler, imputation generation, or real data is used. Temporary
trace fixtures are numerical arrays, not draws from the implemented model.
"""
from __future__ import annotations
import contextlib
import csv
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
TEST_TMP = ROOT/'.statatmp/mlmi_code_tests'


@contextlib.contextmanager
def temporary_directory():
    TEST_TMP.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=TEST_TMP) as path:
        if not Path(path).resolve().is_relative_to(ROOT.resolve()):
            raise ValueError('Test temporary directory escaped diagnostic root')
        yield path


def load(name, path):
    spec = importlib.util.spec_from_file_location(name,path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


check = load('mlmi_check_test',ROOT/'code/evaluation/check_multilevel_mi.py')
run = load('mlmi_run_test',ROOT/'code/execution/run_designscore_mlmi.py')
summary = load('mlmi_summary_test',ROOT/'code/evaluation/summarize_designscore_mlmi.py')


class Diagnostics(unittest.TestCase):
    def fixture(self):
        # Fixed deterministic permutations with equal empirical distributions.
        base = np.sin(np.arange(512)*1.61803398875)
        return np.array([base[(np.arange(512)*k+13*k)%512] for k in (3,7,11,17)])

    def test_average_tied_ranks(self):
        np.testing.assert_array_equal(check.ranks([[4,1],[1,9]]),[[3,1.5],[1.5,4]])

    def test_rhat_known_value(self):
        # W=5/3; B=8; var_plus=13/4; Rhat=sqrt(39/20).
        self.assertAlmostEqual(check.rhat(np.array([[0,1,2,3],[2,3,4,5.]])),np.sqrt(39/20))

    def test_location_and_scale_nonconvergence(self):
        base = self.fixture()
        good = check.diagnose(base)
        self.assertLess(good['rhat'],1.02)
        self.assertGreater(good['bulk_ess'],400)
        shifted = base.copy(); shifted[0] += 8
        scaled = base.copy(); scaled[0] *= 100
        self.assertGreater(check.diagnose(shifted)['rhat'],1.1)
        self.assertGreater(check.diagnose(scaled)['rhat'],1.1)

    def test_rank_invariance(self):
        a = check.diagnose(self.fixture())
        b = check.diagnose(self.fixture()*7+42)
        for key in a:
            self.assertAlmostEqual(a[key],b[key])

    def test_stuck_or_invalid_chains_fail(self):
        result = check.diagnose(np.ones((4,100)))
        self.assertFalse(np.isfinite(result['rhat']))
        self.assertEqual(result['bulk_ess'],0)
        with self.assertRaises(ValueError):
            check.diagnose(np.full((4,100),np.nan))

    def test_trace_validation_requires_complete_iterations(self):
        config = dict(chains=4,imputations=4,burn=20,between=100,every=1,max_rhat=1.05,min_ess=20)
        with temporary_directory() as td:
            for chain in range(1,5):
                with (Path(td)/f'trace_1_chain{chain}.csv').open('w',newline='') as stream:
                    writer = csv.writer(stream)
                    writer.writerow(['chain','iteration','warmup','log_joint'])
                    for i in range(1,121):
                        writer.writerow([chain,i,int(i<=20),np.sin(i*chain*1.618)])
            result = check.check_traces(td,1,config)
            self.assertEqual(result['postburn_trace_draws_per_chain'],100)
            path = Path(td)/'trace_1_chain2.csv'
            path.write_text('\n'.join(path.read_text().splitlines()[:-1])+'\n')
            with self.assertRaisesRegex(ValueError,'Incomplete'):
                check.check_traces(td,1,config)


class Pipeline(unittest.TestCase):
    def test_plan_launches_nothing_and_creates_nothing(self):
        with patch.object(run.subprocess,'run',side_effect=AssertionError('Process launched')), \
             patch.object(run,'freeze_inputs',side_effect=AssertionError('Manifest written')), \
             contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(run.main([]),0)
        self.assertIn('"models_launched": 0',out.getvalue())

    def test_old_labels_and_unequal_chains_rejected(self):
        for args in (['--run-label','mi_s1_v2_m5'],['--imputations','5'],['--every','0']):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                run.options(args)

    def test_fitting_requires_saved_validated_imputations(self):
        a = run.options(['--stage','fit'])
        with temporary_directory() as td:
            with self.assertRaisesRegex(ValueError,'No saved multilevel'):
                run.validate_and_promote(ROOT,Path(td),1,a,check)

    def test_existing_pool_function_is_reused_exactly(self):
        self.assertIs(summary.pool,summary._legacy.pool)

    def test_saved_cells_reject_leakage_or_changed_observations(self):
        import pandas as pd
        raw_rows, saved_rows = [], []
        for person in (1,2):
            row = dict(bgid=1,hhid=1,subid=person)
            for visit in (1,2,3):
                missing = int((person,visit) in ((1,2),(2,3)))
                raw_rows.append(dict(bgid=1,hhid=1,subid=person,v_num=visit,
                                     y_gfr=100+person+visit,y_bmi=20+person+visit,
                                     miss_ind_mar=missing,x6=visit-1))
                row[f'mar_{visit}']=missing
                row[f'x6_{visit}']=visit-1
                for outcome,base in (('y_gfr',100),('y_bmi',20)):
                    value=base+person+visit
                    row[f'{outcome}_{visit}']=np.nan if missing else value
                    if visit>1:
                        for m in (1,2):
                            row[f'_{m}_{outcome}_{visit}']=value+m if missing else value
            saved_rows.append(row)
        raw, saved = pd.DataFrame(raw_rows),pd.DataFrame(saved_rows)
        with temporary_directory() as td:
            root = Path(td)/'fake_parent/fake_submodule/designscore_diagnostics'
            with patch.object(pd,'read_csv',return_value=raw),patch.object(pd,'read_stata',return_value=saved):
                self.assertTrue(check.verify_imputations(root,'fixture.dta',1,2)['passed'])
                changed=saved.copy();changed.loc[1,'_1_y_gfr_2']+=1
                with patch.object(pd,'read_stata',return_value=changed), self.assertRaisesRegex(ValueError,'modified observed'):
                    check.verify_imputations(root,'fixture.dta',1,2)
                leaked=saved.copy();leaked.loc[0,'y_gfr_2']=103
                with patch.object(pd,'read_stata',return_value=leaked), self.assertRaisesRegex(ValueError,'leaked'):
                    check.verify_imputations(root,'fixture.dta',1,2)


if __name__=='__main__':
    unittest.main(verbosity=2)
