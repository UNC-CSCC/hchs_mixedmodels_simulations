"""Summarize an explicitly requested design-score replicate range, including failures."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics


def read(path):
    with path.open(newline='', encoding='utf-8-sig') as source:
        return list(csv.DictReader(source))


def number(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return math.nan


def average(values):
    return statistics.mean(values) if values else math.nan


def write(path, rows):
    with path.open('w', newline='', encoding='utf-8') as target:
        writer = csv.DictWriter(target, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows({key: '' if isinstance(value, float) and not math.isfinite(value)
                         else value for key, value in row.items()} for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', '--repo', dest='repo', type=Path,
                        default=Path(__file__).resolve().parents[2],
                        help='Path to designscore_diagnostics (--repo is an alias)')
    parser.add_argument('--output-dir', type=Path,
                        help='Optional destination for summary-only verification or a separate report')
    parser.add_argument('--population', choices=('oct2026', 'dec2024'), default='oct2026')
    parser.add_argument('--mode', choices=('exact', 'bernoulli'), required=True)
    parser.add_argument('--first', type=int, default=1)
    parser.add_argument('--last', type=int, default=200)
    args = parser.parse_args()
    if not 1 <= args.first <= args.last:
        parser.error('Invalid replicate range')
    repo = args.repo.resolve()
    truth_path = repo / f'data/population_coefficients/population_coef_{args.population}_hhsub.csv'
    truth_rows = read(truth_path)
    truth = {list(row.values())[0]: number(list(row.values())[1]) for row in truth_rows}
    estimates = repo / f'data/estimates/{args.population}/{args.mode}'
    destination = args.output_dir.resolve() if args.output_dir else repo / f'summaries/{args.population}'
    destination.mkdir(parents=True, exist_ok=True)
    by_term = {term: [] for term in truth}
    audits = []
    checksums = {}
    for sim_id in range(args.first, args.last + 1):
        coefficient = estimates / f'scen3_designscore_coef_{sim_id}.csv'
        audit_path = estimates / f'scen3_designscore_audit_{sim_id}.csv'
        if not coefficient.exists() or not audit_path.exists():
            audits.append(dict(sim_id=sim_id, fit_status='MISSING_OUTPUT', se_status='',
                boundary_mode='', kkt_error='', step1_score_error='', jacobian_condition=''))
            continue
        rows, audit = read(coefficient), read(audit_path)
        assert len(audit) == 1 and len(rows) == len(truth), sim_id
        assert {row['term'] for row in rows} == set(truth), sim_id
        assert all(int(row['sim_id']) == sim_id and row['population'] == args.population
                   and row['probability_mode'] == args.mode for row in rows + audit), sim_id
        assert all(row['fit_status'] == audit[0]['fit_status']
                   and row['se_status'] == audit[0]['se_status'] for row in rows), sim_id
        assert all(number(row['std_error']) > 0
                   and number(row['conf_low']) < number(row['conf_high'])
                   for row in rows if row['inference_usable'] == '1'), sim_id
        checksums[coefficient.name] = hashlib.sha256(coefficient.read_bytes()).hexdigest()
        audits.append({key: audit[0][key] for key in
            ('sim_id', 'fit_status', 'se_status', 'boundary_mode', 'kkt_error',
             'step1_score_error', 'jacobian_condition')})
        for row in rows:
            by_term[row['term']].append(row)
    expected = args.last - args.first + 1
    summaries = []
    for term in sorted(truth, key=lambda t: (t.startswith('var_') or t == 's2', t)):
        rows = by_term[term]
        points = [number(row['estimate']) for row in rows if math.isfinite(number(row['estimate']))]
        intervals = [row for row in rows if row['inference_usable'] == '1'
                     and all(math.isfinite(number(row[key])) for key in
                             ('estimate', 'std_error', 'conf_low', 'conf_high'))]
        true = truth[term]
        mean = average(points)
        bias = mean - true
        empirical_se = statistics.stdev(points) if len(points) > 1 else math.nan
        model_se = average([number(row['std_error']) for row in intervals])
        covered = [number(row['conf_low']) <= true <= number(row['conf_high']) for row in intervals]
        coverage = average(covered)
        summaries.append(dict(scenario=f'scen3_designscore_{args.mode}_{args.population}',
            term_type='variance_component' if term.startswith('var_') or term == 's2' else 'fixed_effect',
            term=term, n_reps=len(points), true_value=true, mean_estimate=mean,
            empirical_bias=bias, relative_bias=bias/true if true else math.nan,
            empirical_se=empirical_se, model_se=model_se,
            relative_se_diff=(model_se-empirical_se)/empirical_se if empirical_se > 0 else math.nan,
            coverage_rate=coverage,
            coverage_emp_se=average([abs(point-true) <= 1.959963984540054*empirical_se for point in points]),
            p_reject_h0=average([number(row['conf_low']) > 0 or number(row['conf_high']) < 0 for row in intervals]),
            n_expected=expected, n_missing_points=expected-len(points), n_valid_intervals=len(intervals),
            n_missing_intervals=expected-len(intervals),
            coverage_over_all_requested=sum(covered)/expected,
            bias_mcse=empirical_se/math.sqrt(len(points)) if points else math.nan,
            coverage_mcse=math.sqrt(coverage*(1-coverage)/len(intervals)) if intervals else math.nan,
            probability_mode=args.mode, se_method='stacked_stratified_psu_wr',
            first_replicate=args.first, last_replicate=args.last))
    stem = f'scen3_designscore_{args.mode}'
    write(destination / f'coverage_{stem}.csv', summaries)
    write(destination / f'audit_{stem}.csv', audits)
    provenance = dict(population=args.population, probability_mode=args.mode,
        first_replicate=args.first, last_replicate=args.last, expected_replicates=expected,
        coefficient_sha256=checksums,
        truth_sha256=hashlib.sha256(truth_path.read_bytes()).hexdigest(),
        source_sha256={str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in (repo/'code/stata/scen3_designscore.do',
                                   repo/'code/stata/variance_designscore.mata')})
    provenance['metadata_sha256'] = {
        str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (repo/f'data/design_frames/{args.population}_households.json',
                     repo/f'data/design_frames/{args.population}_subject_counts.csv')
        if path.exists()}
    (destination/f'provenance_{stem}.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    print(json.dumps(summaries[-3:], indent=2))
    print(f'Saved summaries and audit to {destination}')


if __name__ == '__main__':
    main()
