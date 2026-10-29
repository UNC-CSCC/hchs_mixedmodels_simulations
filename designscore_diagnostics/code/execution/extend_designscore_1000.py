"""Resume the local October design-score study from 200 to 1,000 replicates.

Each worker verifies exact sampling counts for its range, then runs Stata.
Existing fits are retained. A complete run updates the combined summary.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time


def rows(path):
    with path.open(newline='', encoding='utf-8-sig') as source:
        return list(csv.DictReader(source))


def digest(path, algorithm='sha256'):
    with path.open('rb') as source:
        return hashlib.file_digest(source, algorithm).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--stata', type=Path, default=Path('C:/Program Files/StataNow19/StataBE-64.exe'))
    parser.add_argument('--rscript', type=Path, default=Path('C:/Program Files/R/R-4.5.0/bin/Rscript.exe'))
    args = parser.parse_args()
    if not 1 <= args.workers <= 24:
        parser.error('workers must be between 1 and 24')
    repo = Path(__file__).resolve().parents[2]
    state = repo / '.statatmp/designscore_1000'
    state.mkdir(parents=True, exist_ok=True)
    temporary = state / 'runtime_tmp'
    temporary.mkdir(exist_ok=True)
    log_dir = repo / 'logs/local_runs/designscore_1000'
    log_dir.mkdir(parents=True, exist_ok=True)
    batch_logs = repo / 'logs/local_batches'
    batch_logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TMP=str(temporary), TEMP=str(temporary), TMPDIR=str(temporary),
               LANG='C', LC_ALL='C', PYTHONDONTWRITEBYTECODE='1')
    startup = subprocess.STARTUPINFO()
    startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    startup.wShowWindow = subprocess.SW_HIDE
    process_options = dict(cwd=repo, env=env, startupinfo=startup,
                           creationflags=subprocess.CREATE_NO_WINDOW)
    frame_dir = repo / 'data/design_frames'
    estimate_dir = repo / 'data/estimates/oct2026/exact'
    hchs_root = repo.parents[1]
    samples = hchs_root / 'data/Oct2026/sample_miss00pct'
    original_counts = rows(frame_dir / 'oct2026_subject_counts_1_200.csv')
    assert [int(row['sim_id']) for row in original_counts] == list(range(1, 201))
    inventory = sorted(int(p.stem.split('_')[-1]) for p in samples.glob('samplemiss_*.csv'))
    assert inventory == list(range(1, 1001)), 'Require exactly samples 1 through 1000'
    baseline = json.loads((repo/'summaries/oct2026/provenance_scen3_designscore_exact.json').read_text())
    relocation = json.loads((repo/'provenance/relocation_manifest.json').read_text(encoding='utf-8-sig'))
    moved = {entry['old_path']: entry for entry in relocation['moves']}
    for relative, expected in baseline['source_sha256'].items():
        relative = relative.replace('\\', '/')
        assert digest(repo/'provenance/pre_reorganization'/relative) == expected, relative
        entry = moved[relative]
        assert digest(repo/entry['new_path']) == entry['sha256_after'], f'Estimator changed since relocation: {relative}'
    for sim in range(1, 201):
        name = f'scen3_designscore_coef_{sim}.csv'
        expected = baseline['coefficient_sha256'][name]
        assert digest(estimate_dir/name) == expected, f'Prior fit changed: {name}'

    lock = threading.Lock()
    progress = dict(started_at=datetime.now(timezone.utc).isoformat(), status='verifying_previous_samples',
                    workers=args.workers, expected_replicates=1000, jobs={})

    def publish():
        temporary_file = state/'status.json.tmp'
        # Windows sync/indexing readers can briefly hold the destination open.
        # A transient progress-file lock must not abort a scientific fit.
        for attempt in range(40):
            try:
                temporary_file.write_text(json.dumps(progress, indent=2), encoding='utf-8')
                temporary_file.replace(state/'status.json')
                break
            except PermissionError:
                if attempt == 39:
                    print('Progress snapshot temporarily locked; continuing batch', flush=True)
                    break
                time.sleep(.25)

    def update(label, stage, **details):
        with lock:
            progress['jobs'][label] = dict(stage=stage, updated_at=datetime.now(timezone.utc).isoformat(), **details)
            publish()
        print(f'{label}: {stage}', flush=True)

    publish()
    # Reuse prior estimates only if their source samples are still identical.
    def check_previous(row):
        sim = int(row['sim_id'])
        assert digest(samples/f'samplemiss_{sim}.csv', 'md5').upper() == row['sample_md5'], sim
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(check_previous, original_counts))
    progress['status'] = 'running'
    progress['previous_200_samples_verified'] = True
    publish()

    def job(first, last):
        label = f'{first}_{last}'
        counts = frame_dir/f'oct2026_subject_counts_{label}.csv'
        update(label, 'recovering_counts')
        if not counts.exists():
            with (log_dir/f'counts_{label}.log').open('w') as log:
                result = subprocess.run([str(args.rscript), '--vanilla',
                    'code/preparation/recover_designscore_counts.R', str(hchs_root),
                    str(counts), str(last), str(first)], stdout=log, stderr=subprocess.STDOUT,
                    timeout=3600, **process_options)
                if result.returncode:
                    raise RuntimeError(f'Count recovery failed for {label}: return code {result.returncode}')
        records = rows(counts)
        assert [int(r['sim_id']) for r in records] == list(range(first, last+1)), label
        assert all(r['ids_verified'] == 'TRUE' for r in records), label
        wrapper = state/f'batch_{label}.do'
        counts_relative = counts.relative_to(repo).as_posix()
        wrapper.write_text('version 19\nset more off\n'
            f'cd "{repo.as_posix()}"\n'
            f'do code/stata/run_designscore_local.do {first} {last} exact "{counts_relative}"\n'
            'exit, clear\n', encoding='utf-8')
        update(label, 'fitting', verified_counts=len(records))
        stata_options = dict(process_options, cwd=batch_logs)
        process = subprocess.Popen([str(args.stata), '/e', 'do', str(wrapper)], **stata_options)
        update(label, 'fitting', verified_counts=len(records), pid=process.pid)
        result = process.wait(timeout=6*3600)
        manifest = rows(batch_logs/f'local_run_{label}.csv')
        assert [int(r['sim_id']) for r in manifest] == list(range(first, last+1)), label
        assert result == 0 and all(r['return_code'] in ('0', '430') for r in manifest), label
        failures = [int(r['sim_id']) for r in manifest if r['return_code'] != '0']
        update(label, 'complete', verified_counts=len(records), attempted=len(manifest), failed_fits=failures)
        return records

    combined = list(original_counts)
    failures = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(job, first, first+49): f'{first}_{first+49}'
                   for first in range(201, 1001, 50)}
        for future in as_completed(futures):
            label = futures[future]
            try:
                combined.extend(future.result())
            except Exception as error:
                update(label, 'error', error=str(error))
                failures.append(label)
    if failures:
        progress['status'] = 'failed'
        progress['failed_jobs'] = failures
        publish()
        raise RuntimeError(f'Incomplete ranges: {failures}')
    combined.sort(key=lambda row: int(row['sim_id']))
    assert [int(row['sim_id']) for row in combined] == list(range(1, 1001))
    destination = frame_dir/'oct2026_subject_counts.csv'
    staged = frame_dir/'oct2026_subject_counts.csv.tmp'
    with staged.open('w', newline='', encoding='utf-8') as target:
        writer = csv.DictWriter(target, fieldnames=list(original_counts[0]))
        writer.writeheader()
        writer.writerows(combined)
    staged.replace(destination)
    progress['status'] = 'summarizing'
    publish()
    with (log_dir/'summary.log').open('w') as log:
        subprocess.run([sys.executable, '-B', 'code/evaluation/summarize_designscore.py',
            '--root', str(repo), '--mode', 'exact', '--first', '1', '--last', '1000'],
            stdout=log, stderr=subprocess.STDOUT, check=True, **process_options)
    progress['status'] = 'complete'
    progress['completed_at'] = datetime.now(timezone.utc).isoformat()
    publish()
    print('DESIGNSCORE_1000_COMPLETE', flush=True)


if __name__ == '__main__':
    main()
