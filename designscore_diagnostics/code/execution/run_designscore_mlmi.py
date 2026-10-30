"""Separate multilevel MI pipeline. Default stage=plan launches no processes.

Explicit stages: impute (sampler + validation), fit (validated saved MI only),
all (both), summarize (saved paired estimates only). Never edits legacy code.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import subprocess


MODEL = 'joint6_psu_household_v1'
SOURCES = (
    'code/stata/multilevel_mi.mata', 'code/stata/impute_multilevel_mi.do',
    'code/stata/multilevel_mi_collapsed.mata',
    'code/stata/scen3_designscore_mi.do', 'code/stata/variance_designscore.mata',
    'code/stata/mi_rename_compat.do', 'code/evaluation/check_multilevel_mi.py',
    'code/evaluation/summarize_designscore_mlmi.py',
    'code/evaluation/summarize_designscore_mi.py',
    'code/execution/run_designscore_mlmi.py',
)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def digest(path, algorithm='sha256'):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def json_write(path, value):
    def clean(x):
        if isinstance(x, float) and not math.isfinite(x):
            return None
        if isinstance(x, dict):
            return {k: clean(v) for k,v in x.items()}
        if isinstance(x, (list, tuple)):
            return [clean(v) for v in x]
        return x
    path = Path(path)
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(clean(value), indent=2, allow_nan=False), encoding='utf-8')
    temp.replace(path)


def options(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage', choices=['plan','impute','fit','all','summarize'], default='plan')
    ap.add_argument('--first', type=int, default=1)
    ap.add_argument('--last', type=int, default=1)
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--run-label', default='mlmi_joint6_m20_v1')
    ap.add_argument('--imputations', type=int, default=20)
    ap.add_argument('--chains', type=int, default=4)
    ap.add_argument('--burn', type=int, default=1000)
    ap.add_argument('--between', type=int, default=200)
    ap.add_argument('--every', type=int, default=5, help='Trace interval, including warmup')
    ap.add_argument('--ridge', type=float, default=1e-4)
    ap.add_argument('--hh-prior', type=float, default=.05)
    ap.add_argument('--psu-prior', type=float, default=.01)
    ap.add_argument('--person-prior', type=float, default=.5)
    ap.add_argument('--max-rhat', type=float, default=1.01)
    ap.add_argument('--min-ess', type=float, default=400)
    ap.add_argument('--timeout', type=int, default=86400, help='Per Stata process, seconds')
    ap.add_argument('--stata', type=Path, default=Path('C:/Program Files/StataNow19/StataBE-64.exe'))
    a = ap.parse_args(argv)
    if not re.fullmatch(r'mlmi_[A-Za-z0-9_]+', a.run_label):
        ap.error('Use a distinct mlmi_... run label, with letters/digits/underscores only')
    if not (1 <= a.first <= a.last <= 1000 and 1 <= a.workers <= 16 and 2 <= a.chains <= 16):
        ap.error('Invalid replicate, worker or chain count')
    if a.imputations < a.chains or a.imputations % a.chains:
        ap.error('Imputations must be a positive multiple of chains')
    if min(a.burn,a.between,a.every,a.timeout) <= 0 or a.burn % a.every or a.between % a.every:
        ap.error('Positive iteration counts, with burn/between divisible by every, required')
    if a.imputations//a.chains*a.between//a.every < 100:
        ap.error('Need at least 100 post-burn trace observations per chain')
    numeric = [a.ridge,a.hh_prior,a.psu_prior,a.person_prior,a.max_rhat,a.min_ess]
    if not all(math.isfinite(v) and v > 0 for v in numeric) or a.max_rhat <= 1:
        ap.error('Invalid priors or convergence thresholds')
    return a


def settings(a):
    keys = ('imputations','chains','burn','between','every','ridge','hh_prior','psu_prior','person_prior','max_rhat','min_ess')
    return {k: getattr(a,k) for k in keys}


def freeze_inputs(root, a, summary):
    config = settings(a)
    counts = root/'data/design_frames/oct2026_subject_counts.csv'
    frame = root/'data/design_frames/oct2026_households.csv'
    samples = root.parents[1]/'data/Oct2026/sample_miss00pct'
    rows = {int(r['sim_id']): r for r in summary.read(counts)}
    spec = dict(model=MODEL, sampler='gibbs_interweaving_collapsed_hr_v3', configuration=config, sample_directory=str(samples),
                source_sha256={s: digest(root/s) for s in SOURCES},
                frame_sha256=digest(frame), counts_sha256=digest(counts),
                truth_sha256=digest(root/'data/population_coefficients/population_coef_oct2026_hhsub.csv'))
    path = root/'provenance'/f'{a.run_label}_inputs.json'
    previous = json.loads(path.read_text()) if path.exists() else None
    if previous and any(previous.get(k) != v for k,v in spec.items()):
        raise ValueError('Run inputs/configuration changed; use a new mlmi_ run label')
    sample_info = dict(previous.get('samples',{})) if previous else {}
    for sim in range(a.first,a.last+1):
        source = samples/f'samplemiss_{sim}.csv'
        md5 = digest(source,'md5').upper()
        if rows[sim]['ids_verified'] != 'TRUE' or md5 != rows[sim]['sample_md5']:
            raise ValueError(f'Sample {sim} differs from the exact design frame/count provenance')
        info = dict(md5=md5, sha256=digest(source), mi_seed=500000+sim,
                    chain_seeds=[500000+sim+1000000*c for c in range(1,a.chains+1)])
        if str(sim) in sample_info and sample_info[str(sim)] != info:
            raise ValueError(f'Sample {sim} changed')
        sample_info[str(sim)] = info
    if a.stage in ('fit','summarize') and previous is None:
        raise ValueError('No frozen multilevel MI input manifest; generate imputations first')
    spec['samples'] = sample_info
    path.parent.mkdir(parents=True,exist_ok=True)
    json_write(path,spec)
    return spec


def stata_run(a, root, state, logs, sim, stage, command, marker):
    wrapper = state/f'{stage}_{sim}.do'
    wrapper.write_text('version 19\nset more off\n'+
                      f'cd "{root.as_posix()}"\n'+command+'\nexit, clear\n',encoding='utf-8')
    si = subprocess.STARTUPINFO()
    si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    si.wShowWindow = subprocess.SW_HIDE
    tmp = state/f'tmp_{sim}'
    tmp.mkdir(exist_ok=True)
    env = dict(os.environ,STATATMP=str(tmp),TMP=str(tmp),TEMP=str(tmp))
    result = subprocess.run([str(a.stata),'/e','do',str(wrapper)],cwd=logs,env=env,
                            startupinfo=si,creationflags=subprocess.CREATE_NO_WINDOW,timeout=a.timeout)
    # Stata Windows /e can exit 0 even on a do-file error. Require completion
    # marker as an executed display line, not its echo inside a program body.
    outer = logs/f'{stage}_{sim}.log'
    text = outer.read_text(errors='replace') if outer.exists() else ''
    if result.returncode or not any(line.strip() == marker for line in text.splitlines()):
        raise RuntimeError(f'Stata {stage} failed for replicate {sim}; inspect {outer}')


def validate_and_promote(root, data, sim, a, checker):
    final = data/'imputations'/f'mi_{sim}.dta'
    candidate = data/'imputations'/f'candidate_{sim}.dta'
    path = final if final.exists() else candidate
    if not path.exists():
        raise ValueError(f'No saved multilevel imputations for sample {sim}')
    traces = checker.check_traces(data/'sampler',sim,settings(a))
    integrity = checker.verify_imputations(root,path,sim,a.imputations)
    report = dict(model=MODEL, sim_id=sim, configuration=settings(a),
                  passed=traces['passed'] and integrity['passed'],
                  imputation_sha256=digest(path),
                  trace_sha256={p.name:digest(p) for p in (data/'sampler').glob(f'trace_{sim}_chain*.csv')},
                  mcmc=traces, integrity=integrity,
                  checked_utc=datetime.now(timezone.utc).isoformat())
    gate = data/'sampler'/f'validation_{sim}.json'
    json_write(gate,report)
    if not report['passed']:
        raise ValueError(f'Chain checks failed for sample {sim}; candidate retained; no analysis launched. See {gate}')
    if not final.exists():
        candidate.replace(final)
    return final


def main(argv=None):
    a = options(argv)
    root = Path(__file__).resolve().parents[2]
    if a.stage == 'plan':
        print(json.dumps(dict(stage='plan',models_launched=0,run_label=a.run_label,
                              replicates=[a.first,a.last],configuration=settings(a)),indent=2))
        return 0
    summary = module(root/'code/evaluation/summarize_designscore_mlmi.py','mlmi_summary')
    checker = module(root/'code/evaluation/check_multilevel_mi.py','mlmi_check')
    freeze_inputs(root,a,summary)
    data = root/'data/estimates/oct2026'/a.run_label
    logs = root/'logs/oct2026'/a.run_label
    state = root/'.statatmp'/a.run_label
    for folder in [*(data/x for x in ('imputations','per_imputation','pooled','audits','sampler')),logs,state]:
        folder.mkdir(parents=True,exist_ok=True)
    if a.stage == 'summarize':
        destination,_,_ = summary.summarize(root,a.run_label,a.first,a.last,a.imputations)
        print(destination)
        return 0

    def job(sim):
        final = data/'imputations'/f'mi_{sim}.dta'
        candidate = data/'imputations'/f'candidate_{sim}.dta'
        if not final.exists() and not candidate.exists():
            if a.stage == 'fit':
                raise ValueError('fit stage requires saved imputations; it never generates them')
            # A partial trace means interrupted generation. Preserve it and
            # require a fresh label rather than silently replacing evidence.
            if list((data/'sampler').glob(f'*_{sim}_chain*.csv')):
                raise ValueError('Interrupted sampler files present; use a new run label')
            numbers = ' '.join(str(v) for v in (a.imputations,a.chains,a.burn,a.between,a.every,a.ridge,a.hh_prior,a.psu_prior,a.person_prior))
            command = f'do code/stata/impute_multilevel_mi.do {sim} "{root.as_posix()}" {a.run_label} {numbers}'
            print(f'Replicate {sim}: generating multilevel imputations',flush=True)
            stata_run(a,root,state,logs,sim,'impute',command,
                      f'MLMI_GENERATION_COMPLETE sim={sim} M={a.imputations}; candidate only, no analysis fitted.')
        final = validate_and_promote(root,data,sim,a,checker)
        print(f'Replicate {sim}: saved imputation and chain checks passed',flush=True)
        if a.stage in ('fit','all'):
            estimates = data/'per_imputation'/f'mi_estimates_{sim}.csv'
            audit = data/'audits'/f'mi_audit_{sim}.csv'
            seal = data/'sampler'/f'paired_fit_{sim}.json'
            sealed = json.loads(seal.read_text()) if seal.exists() else None
            expected = dict(imputation_sha256=digest(final))
            if estimates.exists() or audit.exists() or seal.exists():
                expected.update(estimates_sha256=digest(estimates),audit_sha256=digest(audit))
                if sealed != expected:
                    raise ValueError('Incomplete or changed paired output; use a new run label')
            else:
                command = f'do code/stata/scen3_designscore_mi.do {sim} "{root.as_posix()}" {a.run_label} {a.imputations}'
                stata_run(a,root,state,logs,sim,'fit',command,
                          f'PAIRED_MI_COMPLETE sim={sim} M={a.imputations}')
                if digest(final) != expected['imputation_sha256']:
                    raise ValueError('Shared imputation dataset changed during fitting')
                expected.update(estimates_sha256=digest(estimates),audit_sha256=digest(audit))
                json_write(seal,expected)
            pooled = summary.pool_replicate(root,a.run_label,sim,a.imputations)
            if any(row['fit_status'] != 'OK' for row in pooled):
                raise ValueError('Paired fits saved, but at least one pooled term lacks usable inference')
        return sim

    failures = []
    with ThreadPoolExecutor(max_workers=a.workers) as executor:
        futures = {executor.submit(job,sim):sim for sim in range(a.first,a.last+1)}
        for future in as_completed(futures):
            sim = futures[future]
            try:
                future.result()
                print(f'Replicate {sim}: {a.stage} complete',flush=True)
            except Exception as exc:
                failures.append(sim)
                json_write(data/'sampler'/f'failure_{sim}.json',dict(stage=a.stage,error=str(exc)))
                print(f'Replicate {sim}: FAILED: {exc}',flush=True)
    if a.stage in ('fit','all'):
        summary.summarize(root,a.run_label,a.first,a.last,a.imputations)
    return int(bool(failures))


if __name__ == '__main__':
    raise SystemExit(main())
