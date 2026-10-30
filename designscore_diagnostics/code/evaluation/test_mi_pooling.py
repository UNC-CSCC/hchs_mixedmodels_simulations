"""Small numerical and failure-handling checks; does not run analysis models."""
import math
import unittest
from summarize_designscore_mi import pool, t_quantiles


class PoolingTests(unittest.TestCase):
    def test_rubin_hand_calculation(self):
        result = pool([1, 2, 3], [.5, .5, .5])
        self.assertEqual(result['estimate'], 2)
        self.assertEqual(result['within_variance'], .25)
        self.assertEqual(result['between_variance'], 1)
        self.assertAlmostEqual(result['total_variance'], 19/12)
        self.assertAlmostEqual(result['mi_df'], 2.8203125)
        self.assertAlmostEqual(result['critical_value'], 3.3002841783, places=9)

    def test_zero_between_imputation_variance(self):
        result = pool([4]*5, [2]*5)
        self.assertEqual(result['std_error'], 2)
        self.assertEqual(result['mi_df'], math.inf)
        self.assertAlmostEqual(result['conf_low'], 4-1.959963984540054*2)

    def test_log_pooling_changes_the_estimand_scale_consistently(self):
        result = pool([1, 4, 16], [.2, .8, 3.2], log_scale=True)
        self.assertAlmostEqual(result['estimate'], 4)
        self.assertEqual(result['arithmetic_mean_estimate'], 7)
        self.assertAlmostEqual(result['within_variance'], .04)
        self.assertAlmostEqual(math.log(result['conf_low']) + math.log(result['conf_high']), 2*math.log(4))
        self.assertGreater(result['conf_low'], 0)

    def test_rejects_incomplete_or_boundary_inputs(self):
        cases = [([1], [1], False), ([1, 2], [1], False),
                 ([1, math.nan], [1, 1], False), ([1, 2], [0, 1], False),
                 ([0, 2], [1, 1], True)]
        for values, errors, log_scale in cases:
            with self.subTest(values=values, errors=errors, log_scale=log_scale):
                with self.assertRaises(ValueError):
                    pool(values, errors, log_scale=log_scale)

    def test_tail_quantiles(self):
        for actual, expected in zip(t_quantiles([1, 10, math.inf]),
                                    [12.7062047361747, 2.22813885198627, 1.95996398454005]):
            self.assertAlmostEqual(actual, expected, places=11)


if __name__ == '__main__':
    unittest.main()
