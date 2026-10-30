"""Read-only validation of saved multilevel MI draws; never fits a model.

Rank-normalized folded/split R-hat; bulk and 5%/95% tail ESS using the
initial-positive, initial-monotone autocorrelation sequence. These diagnostics
are necessary screening conditions, not a proof of convergence or valid MI.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np


def ranks(values):
    x = np.asarray(values, dtype=float).ravel()
    order = np.argsort(x, kind='stable')
    sorted_x = x[order]
    start = np.r_[0, np.flatnonzero(np.diff(sorted_x)) + 1]
    end = np.r_[start[1:], len(x)]
    result = np.empty(len(x))
    for lo, hi in zip(start, end):
        result[order[lo:hi]] = (lo + hi + 1) / 2
    return result.reshape(np.shape(values))


def rank_normalize(values):
    # Blom normal scores, with average ranks for ties.
    r = ranks(values)
    probabilities = (r - 3/8) / (r.size + 1/4)
    return np.array([NormalDist().inv_cdf(p) for p in probabilities.ravel()]).reshape(r.shape)


def split(values):
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or x.shape[0] < 2 or x.shape[1] < 40 or not np.isfinite(x).all():
        raise ValueError('Need at least two finite chains with 40 draws each')
    half = x.shape[1] // 2
    return np.concatenate((x[:, :half], x[:, -half:]), axis=0)


def variance_parts(x):
    n = x.shape[1]
    within = np.var(x, axis=1, ddof=1).mean()
    between = n * np.var(x.mean(axis=1), ddof=1)
    return within, (n - 1)/n * within + between/n


def rhat(x):
    within, total = variance_parts(x)
    return max(1., math.sqrt(total/within)) if within > 0 else math.inf


def ess(x):
    """Multi-chain ESS, capped at the number of stored draws (conservative)."""
    chains, n = x.shape
    within, total = variance_parts(x)
    if within <= 0 or total <= 0:
        return 0.
    z = x - x.mean(axis=1, keepdims=True)
    nfft = 1 << (2*n - 1).bit_length()
    f = np.fft.rfft(z, n=nfft, axis=1)
    acov = np.fft.irfft(f * f.conjugate(), n=nfft, axis=1)[:, :n] / n
    rho = 1 - (within - acov.mean(axis=0)) / total
    rho[0] = 1
    pairs = []
    for k in range(0, n-1, 2):
        pair = rho[k] + rho[k+1]
        if pair <= 0:
            break
        pairs.append(min(pair, pairs[-1]) if pairs else pair)
    tau = max(1., -1 + 2*sum(pairs))
    return min(float(chains*n), chains*n/tau)


def diagnose(values):
    x = split(values)
    z = rank_normalize(x)
    folded = rank_normalize(np.abs(x - np.median(x)))
    tail = [ess((x <= np.quantile(x, q)).astype(float)) for q in (.05, .95)]
    return dict(rhat=max(rhat(z), rhat(folded)), bulk_ess=ess(z), tail_ess=min(tail))


def check_traces(directory, sim, config):
    arrays, columns = [], None
    steps = config['burn'] + config['imputations']//config['chains'] * config['between']
    expected_iterations = list(range(config['every'], steps+1, config['every']))
    for chain in range(1, config['chains']+1):
        path = Path(directory)/f'trace_{sim}_chain{chain}.csv'
        with path.open(newline='', encoding='utf-8-sig') as stream:
            reader = csv.DictReader(stream)
            header = reader.fieldnames
            if not header or len(set(header)) != len(header) or header[:4] != ['chain','iteration','warmup','log_joint']:
                raise ValueError('Invalid trace header')
            rows = [{k: float(v) for k, v in row.items()} for row in reader]
        if columns is None:
            columns = header[3:]
        if columns != header[3:] or [r['iteration'] for r in rows] != expected_iterations:
            raise ValueError('Incomplete or inconsistent trace')
        if any(r['chain'] != chain or r['warmup'] != int(r['iteration'] <= config['burn']) for r in rows):
            raise ValueError('Trace chain/warmup mismatch')
        arrays.append([[r[k] for k in columns] for r in rows if not r['warmup']])
    draws = np.asarray(arrays)
    if not np.isfinite(draws).all():
        raise ValueError('Nonfinite sampler trace')
    records = []
    for j, name in enumerate(columns):
        result = diagnose(draws[:, :, j])
        passed = bool(result['rhat'] <= config['max_rhat'] and
                      result['bulk_ess'] >= config['min_ess'] and result['tail_ess'] >= config['min_ess'])
        records.append(dict(parameter=name, **result, passed=passed))
    return dict(passed=all(r['passed'] for r in records), chains=config['chains'],
                postburn_trace_draws_per_chain=draws.shape[1], parameters=records)


def verify_imputations(root, path, sim, m):
    """Check all saved cells and paired covariance diagnostics against inputs.

    Complete hidden outcomes are used ONLY here for post-imputation diagnostics,
    never by the imputer. No acceptance gate is based on matching their moments.
    """
    import pandas as pd
    root = Path(root)
    d = pd.read_stata(path, convert_categoricals=False)
    source = root.parents[1]/'data/Oct2026/sample_miss00pct'/f'samplemiss_{sim}.csv'
    raw = pd.read_csv(source)
    raw.columns = raw.columns.str.lower()
    ids = ['bgid', 'hhid', 'subid']
    if d.duplicated(ids).any() or raw.duplicated(ids+['v_num']).any():
        raise ValueError('Duplicate subject/visit identifier')
    index = pd.MultiIndex.from_frame(d[ids])
    wide = raw.pivot(index=ids, columns='v_num', values=['y_gfr','y_bmi','miss_ind_mar','x6']).reindex(index)
    if wide.isna().any().any() or len(raw) != len(d)*3:
        raise ValueError('Incomplete/mismatched source subjects')
    for v in (1, 2, 3):
        flag = wide['miss_ind_mar'][v].to_numpy().astype(bool)
        if not np.array_equal(d[f'mar_{v}'].to_numpy(), flag):
            raise ValueError('MAR indicators changed')
        if not np.allclose(d[f'x6_{v}'], wide['x6'][v], rtol=0, atol=1e-10):
            raise ValueError('Visit times changed')
        for outcome in ('y_gfr','y_bmi'):
            actual = d[f'{outcome}_{v}'].to_numpy()
            truth = wide[outcome][v].to_numpy()
            if not np.array_equal(np.isnan(actual), flag):
                raise ValueError('Flagged outcomes leaked into m=0, or observed cells missing')
            if not np.allclose(actual[~flag], truth[~flag], rtol=0, atol=1e-10):
                raise ValueError('Observed source outcome changed')
            if v > 1:
                for imp in range(1, m+1):
                    completed = d[f'_{imp}_{outcome}_{v}'].to_numpy()
                    if not np.isfinite(completed).all() or not np.array_equal(completed[~flag], actual[~flag]):
                        raise ValueError('Incomplete imputation or modified observed values')
    expected = {f'_{imp}_{outcome}_{v}' for imp in range(1,m+1) for outcome in ('y_gfr','y_bmi') for v in (2,3)}
    actual = {c for c in d if c.startswith('_') and '_y_' in c}
    if actual != expected:
        raise ValueError('Saved imputation count/columns mismatch')
    report = dict(passed=True, subjects=len(d), imputations=m, observed_values_preserved=True,
                  missing_values_masked_in_m0=True, all_imputations_complete=True)
    coef = root/'data/estimates/oct2026/exact'/f'scen3_designscore_coef_{sim}.csv'
    if coef.exists():
        b = pd.read_csv(coef).set_index('term').estimate
        mu = np.full(len(d), b['(Intercept)'])
        for term in ('age_strat_new','x12','x17','x18'):
            mu += d[term].to_numpy()*b[term]
        mu = mu[:,None] + d[['x6_1','x6_2','x6_3']].to_numpy()*b['x6']
        hh = d.gl_hhid.to_numpy(dtype=int)
        counts = np.bincount(hh)
        npairs = np.sum(counts*(counts-1))
        def moments(y, bmi):
            r = y - mu - bmi*b['y_bmi']
            centered = r - r.mean(axis=0, keepdims=True)
            offdiag = (centered[:,0]*centered[:,1] + centered[:,0]*centered[:,2] + centered[:,1]*centered[:,2]).mean()/3
            means = r.mean(axis=1); means -= means.mean()
            sums = np.bincount(hh, weights=means)
            sums2 = np.bincount(hh, weights=means**2)
            return dict(household_crossproduct=float(np.sum(sums*sums-sums2)/npairs) if npairs else None,
                        same_subject_crossvisit_covariance=float(offdiag),
                        within_person_mean_square=float(np.mean(np.sum((r-r.mean(axis=1)[:,None])**2,axis=1)/2)),
                        contrast_v1_v2=float(np.mean((r[:,0]-r[:,1])**2)/2),
                        contrast_v1_v3=float(np.mean((r[:,0]-r[:,2])**2)/2),
                        contrast_v2_v3=float(np.mean((r[:,1]-r[:,2])**2)/2))
        original = moments(wide['y_gfr'].to_numpy(), wide['y_bmi'].to_numpy())
        imputed = []
        for imp in range(1,m+1):
            y = d[['y_gfr_1',f'_{imp}_y_gfr_2',f'_{imp}_y_gfr_3']].to_numpy()
            bmi = d[['y_bmi_1',f'_{imp}_y_bmi_2',f'_{imp}_y_bmi_3']].to_numpy()
            imputed.append(dict(imputation=imp, **moments(y,bmi)))
        report['covariance_diagnostics'] = dict(description='Unweighted descriptive moments using fixed saved complete-data coefficients; not refitted variance components or acceptance criteria.', complete=original, imputed=imputed)
    return report
