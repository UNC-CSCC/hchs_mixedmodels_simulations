"""Finalize a completed paired MI batch, optionally waiting for its summary marker.

No analysis models are run here. The Stata check only reads saved imputations.
"""
import argparse
from datetime import datetime, timezone
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first', type=int, default=1)
    parser.add_argument('--last', type=int, default=100)
    parser.add_argument('--run-label', default='mi_s1_v2_m5')
    parser.add_argument('--wait-for-summary', type=Path)
    parser.add_argument('--stata', default='C:/Program Files/StataNow19/StataBE-64.exe')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    os.chdir(root)
    if args.wait_for_summary:
        started = time.monotonic()
        previous = None
        while True:
            try:
                marker = json.loads(args.wait_for_summary.read_text())
                expected = marker['last_replicate'] - marker['first_replicate'] + 1
                if marker['completed_replicates'] != expected:
                    raise RuntimeError('Batch finished with incomplete replicates; resolve them before finalizing')
                break
            except (FileNotFoundError, json.JSONDecodeError):
                pass
            n = len(list((root/'data/estimates/oct2026'/args.run_label/'pooled').glob('mi_coef_*.csv')))
            if n != previous:
                print(f'Waiting for completed batch: {n} pooled replicates available.', flush=True)
                previous = n
            if time.monotonic()-started > 43200:
                raise TimeoutError('Batch did not finish within 12 hours; fits were not stopped')
            time.sleep(30)
    spec = importlib.util.spec_from_file_location('summary', root/'code/evaluation/summarize_designscore_mi.py')
    summary = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(summary)
    dest, coverage, audits = summary.summarize(root, args.run_label, args.first, args.last, 5)
    assert len(audits) == args.last-args.first+1 and all(a['status']=='OK' for a in audits)
    manifest = json.loads((root/'provenance'/f'{args.run_label}_inputs.json').read_text())
    for sim in range(args.first,args.last+1):
        source = Path(manifest['sample_directory'])/f'samplemiss_{sim}.csv'
        assert summary.digest(source)==manifest['samples'][str(sim)]['sha256'], f'Changed sample: {sim}'
    for source, expected_hash in manifest['analysis_source_sha256'].items():
        assert summary.digest(root/source)==expected_hash, f'Changed analysis source: {source}'
    assert summary.digest(root/'data/design_frames/oct2026_households.csv')==manifest['frame_sha256']
    assert summary.digest(root/'data/design_frames/oct2026_subject_counts.csv')==manifest['counts_sha256']
    assert summary.digest(root.parent/'scripts/stata/scen4_mi.do')==manifest['original_mi_script_sha256']
    print('Consolidated paired estimates; checking pooling and input preservation.', flush=True)
    flags = subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0
    def check(command, filename):
        result = subprocess.run(command, text=True, capture_output=True, creationflags=flags)
        (dest/filename).write_text(result.stdout+result.stderr, encoding='utf-8')
        result.check_returncode()
        return result.stdout+result.stderr
    check([sys.executable, 'code/evaluation/test_mi_pooling.py'], 'pooling_tests.txt')
    rscript = os.environ.get('R_SCRIPT','C:/Program Files/R/R-4.5.0/bin/Rscript.exe')
    verification = check([rscript,'--vanilla','code/evaluation/verify_designscore_mi.R',str(root),
                          args.run_label,str(args.first),str(args.last)], 'independent_r_verification.txt')
    assert 'PASS: independently verified' in verification
    state = root/'.statatmp'/args.run_label
    wrapper = state/'verify_all_inputs.do'
    wrapper.write_text(f'cd "{root.as_posix()}"\ndo code/stata/verify_mi_inputs.do {args.first} {args.last} {args.run_label} 5\nexit, clear\n')
    startup = subprocess.STARTUPINFO()
    startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    startup.wShowWindow = subprocess.SW_HIDE
    env = dict(os.environ,STATATMP=str(state/'runtime_tmp'),TMP=str(state/'runtime_tmp'),TEMP=str(state/'runtime_tmp'))
    subprocess.run([args.stata,'/e','do',str(wrapper)],cwd=root/'logs/oct2026'/args.run_label,
                   env=env,startupinfo=startup,creationflags=flags,check=True,timeout=3600)
    input_log = root/'logs/oct2026'/args.run_label/f'verify_mi_inputs_{args.first}_{args.last}.log'
    text = input_log.read_text(errors='replace')
    assert f'MI_INPUT_VERIFICATION_COMPLETE first={args.first} last={args.last} M=5' in text
    # Match actual output lines, not the command echoed by Stata's do-file log.
    verified_ids = [int(v) for v in re.findall(r'^MI_INPUT_PASS sim=(\d+):',text,re.MULTILINE)]
    assert verified_ids == list(range(args.first,args.last+1))
    validation = dict(completed_at=datetime.now(timezone.utc).isoformat(),
                      requested_replicates=len(audits),completed_replicates=len(audits),
                      shared_imputations_verified=len(audits)*5,
                      independent_R_recalculation='PASS',pooling_unit_tests='PASS',
                      input_and_analysis_source_checksums='PASS',
                      observed_values_and_MAR_masking='PASS',input_log=input_log.relative_to(root).as_posix(),
                      input_log_sha256=summary.digest(input_log),
                      verification_source_sha256={p:summary.digest(root/p) for p in [
                          'code/stata/verify_mi_inputs.do','code/evaluation/verify_designscore_mi.R',
                          'code/evaluation/test_mi_pooling.py','code/execution/finalize_designscore_mi.py']})
    (dest/'validation_mi.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
    rows = []
    for term in ['var_hhid','var_subid','s2']:
        a = next(r for r in coverage if r['term']==term and 'original' in r['scenario'])
        b = next(r for r in coverage if r['term']==term and 'designscore' in r['scenario'])
        rows.append(f"| {term} | {a['true_value']:.3f} | {a['empirical_bias']:+.3f} | {a['coverage_rate']:.0%} | {b['empirical_bias']:+.3f} | {b['coverage_rate']:.0%} |")
    miss2 = sum(a['mean_missing_v2'] for a in audits)/len(audits)
    miss3 = sum(a['mean_missing_v3'] for a in audits)/len(audits)
    report = f'''# Oct2026 paired MI comparison: replicates {args.first}–{args.last}

All **{len(audits)} replicates** completed with five shared imputations each.
Both methods used the same {len(audits)*5} completed datasets. The source folder
is `sample_miss00pct`; the MAR flag still masks visits (mean visit-2 missingness
{miss2:.1%}, visit-3 missingness {miss3:.1%}).

| Parameter | Population target | Original MI bias | Original coverage | Design-score MI bias | Design-score coverage |
|---|---:|---:|---:|---:|---:|
'''+ '\n'.join(rows)+'''

Coverage is for 95% Rubin intervals, using log-scale pooling for positive
variance components. Full fixed-effect and variance summaries include empirical
SE, mean model SE, and Monte Carlo SEs. The paired comparison isolates differences
between the analysis methods given these shared imputations.

The imputation model is the existing scenario 4 s1/v2 wide chained-regression
model. It does not explicitly include household or PSU random effects. Therefore
this run does not establish that the imputation model preserves all multilevel
covariance. With 100 replicates, coverage has substantial Monte Carlo uncertainty.

Independent R calculations reproduce the pooled estimates and summary statistics.
Checks of all saved MI datasets confirm that observed values are preserved,
MAR-flagged values are masked before imputation, and all completed imputations
contain the required visits. See `validation_mi.json` for verification records.

- [Original scenario 3 MI summary](coverage_scen3_original_mi.csv)
- [Design-score MI summary](coverage_scen3_designscore_mi.csv)
- [Paired differences](paired_comparison.csv)
- [Pooled estimates](pooled_estimates.csv)
- [Replicate audit](audit_mi.csv)
- [Methods and code guide](../../../MI_COMPARISON_README.md)
'''
    (dest/'RESULTS.md').write_text(report,encoding='utf-8')
    print(f'FINALIZED: {dest}',flush=True)


if __name__=='__main__':
    main()
