"""Independent scalar algebra checks. No data models, fitting or Stata calls."""
from __future__ import annotations

import importlib.util
import itertools
import math
from pathlib import Path
import tempfile
import unittest


def summaries(residual, selected, p, pair):
    nhat = sum(1 / p[i] for i in selected)
    bhat = sum(residual[i] ** 2 / p[i] for i in selected)
    qhat = bhat + 2 * sum(residual[i] * residual[j] / pair[i, j]
                          for i, j in itertools.combinations(selected, 2))
    return nhat, bhat, qhat


def scores(a, b, e, size, stats):
    n, B, Q = stats
    d = b + e / 3
    D = d + a * size
    return (.5 * (Q / D**2 - n / D),
            .5 * (B / d**2 - a * (2*d + a*size)*Q / (d*d*D*D)
                  - n/d + a*n/(d*D)))


def objective(a, b, e, size, stats):
    n, B, Q = stats
    d = b + e / 3
    D = d + a * size
    return -.5 * (n*math.log(d) + n/size*math.log(D/d) + B/d - a*Q/(d*D))


def inverse(matrix):
    """Small dense reference inverse, independent of random-intercept identities."""
    n = len(matrix)
    rows = [list(row) + [float(i == j) for j in range(n)]
            for i, row in enumerate(matrix)]
    for i in range(n):
        pivot = rows[i][i]
        rows[i] = [x/pivot for x in rows[i]]
        for j in range(n):
            if j != i:
                factor = rows[j][i]
                rows[j] = [x-factor*y for x, y in zip(rows[j], rows[i])]
    return [row[n:] for row in rows]


class DesignScoreIdentities(unittest.TestCase):
    def test_census_matches_dense_gaussian_scores(self):
        r = [-7., 12., 4., -2.]
        a, b, e = 22., 175., 95.
        d = b + e/3
        inv = inverse([[a + (d if i == j else 0) for j in range(4)] for i in range(4)])
        vr = [sum(x*y for x, y in zip(row, r)) for row in inv]
        dense_h = .5 * (sum(vr)**2 - sum(map(sum, inv)))
        dense_s = .5 * (sum(x*x for x in vr) - sum(inv[i][i] for i in range(4)))
        p = [1.] * 4
        pair = {ij: 1. for ij in itertools.combinations(range(4), 2)}
        actual = scores(a, b, e, 4, summaries(r, range(4), p, pair))
        self.assertAlmostEqual(actual[0], dense_h, places=12)
        self.assertAlmostEqual(actual[1], dense_s, places=12)

    def test_objective_derivatives_equal_scores(self):
        args = (22., 175., 95., 4, (3.2, 350., -80.))
        analytic = scores(*args)
        for j in range(2):
            plus, minus = list(args), list(args)
            plus[j] += 1e-3
            minus[j] -= 1e-3
            numerical = (objective(*plus)-objective(*minus))/.002
            self.assertAlmostEqual(analytic[j], numerical, places=10)

    def test_bernoulli_expectation_including_empty_household(self):
        r = [-20., 8., 17., 5.]
        p = [.3575, .55, .3575, .55]
        pair = {(i, j): p[i]*p[j] for i, j in itertools.combinations(range(4), 2)}
        expected = [0., 0.]
        for mask in itertools.product((0, 1), repeat=4):
            probability = math.prod(p[i] if mask[i] else 1-p[i] for i in range(4))
            selected = [i for i in range(4) if mask[i]]
            val = scores(22., 175., 95., 4, summaries(r, selected, p, pair))
            expected = [x+probability*y for x, y in zip(expected, val)]
        census = scores(22., 175., 95., 4, (4, sum(x*x for x in r), sum(r)**2))
        for x, y in zip(expected, census):
            self.assertAlmostEqual(x, y, places=12)

    def test_exact_fixed_size_age_sampling_expectation(self):
        # Two younger and two older HH members, drawn from larger stage-3 pools.
        # Enumerate the full pools, so same-age pair dependence is genuine.
        r = [-20., 8., 17., 5.]
        p = [3/6, 3/6, 3/5, 3/5]
        pair = {(i, j): (3*2/(6*5) if j < 2 else 3*2/(5*4) if i >= 2
                         else p[i]*p[j]) for i, j in itertools.combinations(range(4), 2)}
        expected = [0., 0.]
        choices = list(itertools.product(itertools.combinations(range(6), 3),
                                         itertools.combinations(range(5), 3)))
        for younger, older in choices:
            selected = [i for i in (0, 1) if i in younger]
            selected += [i+2 for i in (0, 1) if i in older]
            val = scores(22., 175., 95., 4, summaries(r, selected, p, pair))
            expected = [x+y/len(choices) for x, y in zip(expected, val)]
        census = scores(22., 175., 95., 4, (4, sum(x*x for x in r), sum(r)**2))
        for x, y in zip(expected, census):
            self.assertAlmostEqual(x, y, places=12)

    def test_singleton_keeps_inverse_selection_weight(self):
        n, B, Q = summaries([10.], [0], [.3575], {})
        self.assertAlmostEqual(n, 1/.3575)
        self.assertAlmostEqual(B, 100/.3575)
        self.assertEqual(B, Q)


class FrameReader(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).parents[1] / 'preparation/prepare_designscore_frame.py'
        spec = importlib.util.spec_from_file_location('designscore_frame', path)
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def read(self, rows):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'population.csv'
            path.write_text('BGid,hhid,subid,strat,v_num,hisp_strat,age_base\n' + rows)
            return self.module.read_households(path)

    def test_counts_age_groups_and_does_not_need_outcomes(self):
        rows = ''.join(f'1,1,{subject},1,{visit},TRUE,{age}\n'
                       for subject, age in [(1, 30), (2, 60)] for visit in (1, 2, 3))
        result = self.read(rows)[1, 1]
        self.assertEqual((result['m_h'], result['n_younger'], result['n_older']), (2, 1, 1))

    def test_rejects_incomplete_subject(self):
        with self.assertRaises(ValueError):
            self.read('1,1,1,1,1,TRUE,30\n1,1,1,1,2,TRUE,30\n')

    def test_rejects_duplicate_subject(self):
        rows = ''.join(f'1,1,1,1,{v},TRUE,30\n' for v in (1, 2, 3, 1, 2, 3))
        with self.assertRaises(ValueError):
            self.read(rows)


if __name__ == '__main__':
    unittest.main()
