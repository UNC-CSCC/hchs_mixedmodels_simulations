"""Verify and compare a completed multilevel MI batch; never runs models."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-label', required=True)
    ap.add_argument('--first', type=int, default=1)
    ap.add_argument('--last', type=int, default=100)
    ap.add_argument('--source-runs', nargs='+', help='Collect validated replicates from these runs in preference order')
    ap.add_argument('--rscript', type=Path, default=Path('C:/Program Files/R/R-4.5.0/bin/Rscript.exe'))
    a = ap.parse_args()
    baseline_count = 0
    for baseline_name in ('multilevel_legacy_source_baseline.json',
                          'multilevel_comparison_current_mi_baseline.json'):
        baseline = json.loads((ROOT/'provenance'/baseline_name).read_text())
        for relative, expected in baseline.items():
            if digest(ROOT/relative) != expected:
                raise ValueError(f'Existing MI source/result changed: {relative}')
        baseline_count += len(baseline)
    manifest_path = ROOT/'provenance'/f'{a.run_label}_inputs.json'
    runner = module(ROOT/'code/execution/run_designscore_mlmi.py', 'mlmi_runner_final')
    summary = module(ROOT/'code/evaluation/summarize_designscore_mlmi.py', 'mlmi_summary_final')
    if a.source_runs:
        if not a.run_label.startswith('mlmi_') or a.run_label in a.source_runs:
            raise ValueError('Use a distinct mlmi_ collection label')
        collector = module(ROOT/'code/evaluation/collect_designscore_mlmi.py','mlmi_collector_final')
        destination, coverage, audits, manifest = collector.collect(
            ROOT,a.run_label,a.source_runs,a.first,a.last)
        mcount = manifest['configuration']['imputations']
    else:
        manifest = json.loads(manifest_path.read_text())
        args = ['--stage', 'summarize', '--run-label', a.run_label,
                '--first', str(a.first), '--last', str(a.last)]
        for key, value in manifest['configuration'].items():
            args.extend(['--'+key.replace('_','-'), str(value)])
        config = runner.options(args)
        # Recheck sample hashes, exact design inputs, source files, and settings.
        runner.freeze_inputs(ROOT, config, summary)
        destination, coverage, audits = summary.summarize(
            ROOT, a.run_label, a.first, a.last, config.imputations)
        mcount = config.imputations
    if len(audits) != a.last-a.first+1 or any(row['status'] != 'OK' for row in audits):
        raise ValueError('Final comparison requires every requested replicate to finish successfully')
    si = subprocess.STARTUPINFO()
    si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    si.wShowWindow = subprocess.SW_HIDE
    rcheck = ROOT/'code/evaluation/verify_designscore_mlmi.R'
    with (destination/'independent_R_verification.txt').open('w',encoding='utf-8') as log:
        subprocess.run([str(a.rscript), '--vanilla', str(rcheck), str(ROOT), a.run_label,
                        str(a.first), str(a.last), str(mcount)],
                       cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True,
                       startupinfo=si, creationflags=subprocess.CREATE_NO_WINDOW)
    compare = ROOT/'code/evaluation/compare_multilevel_with_current_mi.py'
    subprocess.run([sys.executable, str(compare), '--new-label', a.run_label,
                    '--first', str(a.first), '--last', str(a.last)], cwd=ROOT, check=True,
                   startupinfo=si, creationflags=subprocess.CREATE_NO_WINDOW)
    report = dict(completed_utc=datetime.now(timezone.utc).isoformat(),
                  run_label=a.run_label, first=a.first, last=a.last,
                  imputations_per_replicate=mcount,
                  completed_replicates=len(audits),
                  identical_imputations_verified=True,
                  independent_R_recalculation='PASS',
                  unchanged_legacy_files_verified=baseline_count,
                  source_sha256={str(p.relative_to(ROOT)):digest(p)
                                 for p in [Path(__file__), rcheck, compare, manifest_path,
                                           ROOT/'code/evaluation/collect_designscore_mlmi.py',
                                           ROOT/'code/execution/complete_designscore_mlmi.py']},
                  summary_sha256={p.name:digest(p) for p in destination.iterdir()
                                  if p.is_file() and p.name!='final_verification.json'})
    runner.json_write(destination/'final_verification.json', report)
    print(json.dumps(dict(status='COMPLETE', summary=str(destination),
                          replicates=len(audits), M=mcount),indent=2))


if __name__ == '__main__':
    main()
