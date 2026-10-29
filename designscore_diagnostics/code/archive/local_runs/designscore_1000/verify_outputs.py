"""Final integrity checks for the completed local 1,000-sample run."""
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

repo = Path.cwd()
state = repo / '.statatmp/designscore_1000'
folder = repo / 'data/estimates/oct2026/scen3_designscore/exact'
summaries = repo / 'data/summaries/oct2026'

def rows(path):
    with path.open(newline='', encoding='utf-8-sig') as source:
        return list(csv.DictReader(source))

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

progress = json.loads((state / 'status.json').read_text())
assert progress['status'] == 'complete', progress['status']
assert len(progress['jobs']) == 16
assert all(j['stage'] == 'complete' and not j['failed_fits']
           for j in progress['jobs'].values())
prior = json.loads((summaries / 'archive_200/provenance_scen3_designscore_exact.json').read_text())
current = json.loads((summaries / 'provenance_scen3_designscore_exact.json').read_text())
assert current['expected_replicates'] == 1000
assert len(current['coefficient_sha256']) == 1000
for name, expected in prior['coefficient_sha256'].items():
    assert digest(folder / name) == expected, name
for name, expected in current['coefficient_sha256'].items():
    assert digest(folder / name) == expected, name
assert prior['source_sha256'] == current['source_sha256']
for name, expected in current['source_sha256'].items():
    assert digest(repo / name) == expected, name
assert digest(repo / 'scripts/stata/scen3.do') == '99d14072842388a820bd7dbc779e24e2b7d1703ad85f2719cf9b57ac45ebc293'
assert prior['truth_sha256'] == current['truth_sha256']
for name, expected in current['metadata_sha256'].items():
    assert digest(repo / name) == expected, name

counts = rows(repo / 'data/design_frames/oct2026_subject_counts.csv')
assert [int(r['sim_id']) for r in counts] == list(range(1, 1001))
assert all(r['ids_verified'] == 'TRUE' for r in counts)
assert counts[:200] == rows(repo / 'data/design_frames/oct2026_subject_counts_1_200.csv')
audits = []
for sim, count in enumerate(counts, 1):
    audit_rows = rows(folder / f'scen3_designscore_audit_{sim}.csv')
    assert len(audit_rows) == 1
    audit = audit_rows[0]
    audits.append(audit)
    assert int(audit['sim_id']) == sim
    assert audit['fit_status'] == 'OK' and audit['se_status'] == '0'
    assert audit['boundary_mode'] == '0'
    assert float(audit['kkt_error']) <= 1e-6
    assert float(audit['step1_score_error']) <= 1e-3
    assert float(audit['jacobian_condition']) <= 1e12
    assert int(audit['n_subjects']) == int(count['selected_younger']) + int(count['selected_older'])
    for age in ('younger', 'older'):
        n, total = int(count[f'selected_{age}']), int(count[f'eligible_{age}'])
        assert int(audit[f'eligible_{age}']) == total
        assert math.isclose(float(audit[f'p_{age}']), n / total, rel_tol=1e-12)
        assert math.isclose(float(audit[f'pair_{age}']), n*(n-1)/(total*(total-1)), rel_tol=1e-12)
    coefficients = rows(folder / f'scen3_designscore_coef_{sim}.csv')
    assert len(coefficients) == 10
    for row in coefficients:
        assert int(row['sim_id']) == sim and row['inference_usable'] == '1'
        assert all(math.isfinite(float(row[k])) for k in ('estimate', 'std_error', 'conf_low', 'conf_high'))
        assert float(row['std_error']) > 0
        assert float(row['conf_low']) < float(row['conf_high'])
for first in range(201, 1001, 50):
    manifest = rows(folder / f'local_run_{first}_{first+49}.csv')
    assert [int(r['sim_id']) for r in manifest] == list(range(first, first+50))
    assert all(r['return_code'] == '0' for r in manifest)
summary = rows(summaries / 'coverage_scen3_designscore_exact.csv')
assert len(summary) == 10
assert all(int(r['n_expected']) == 1000 and int(r['n_reps']) == 1000
           and int(r['n_valid_intervals']) == 1000
           and int(r['n_missing_intervals']) == 0 for r in summary)
result = dict(replicates=1000, coefficient_rows=10000, original_200_unchanged=True,
              original_scen3_unchanged=True, exact_count_records_verified=1000,
              fit_status=dict(Counter(a['fit_status'] for a in audits)),
              max_kkt_error=max(float(a['kkt_error']) for a in audits),
              max_step1_score_error=max(float(a['step1_score_error']) for a in audits),
              max_jacobian_condition=max(float(a['jacobian_condition']) for a in audits),
              n_psu=sorted({int(a['n_psu']) for a in audits}),
              min_subjects=min(int(a['n_subjects']) for a in audits),
              max_subjects=max(int(a['n_subjects']) for a in audits))
(state / 'verification.json').write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
