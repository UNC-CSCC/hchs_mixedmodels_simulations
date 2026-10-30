"""Extend failed MI chains while the primary batch runs, then verify all results.

No accepted analysis is rerun. Retries retain the same sampler/model and M;
only iteration spacing increases. Definitions/importing launch nothing.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time

ROOT=Path(__file__).resolve().parents[2]


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj)
    return obj


def read_json(path):
    return json.loads(path.read_text()) if path.exists() else None


def state(label,sim):
    data=ROOT/'data/estimates/oct2026'/label
    gate=read_json(data/'sampler'/f'validation_{sim}.json')
    failure=read_json(data/'sampler'/f'failure_{sim}.json')
    seal=data/'sampler'/f'paired_fit_{sim}.json'
    pooled=data/'pooled'/f'mi_coef_{sim}.csv'
    if gate and gate['passed'] and seal.exists() and pooled.exists():
        with pooled.open(newline='') as f:rows=list(csv.DictReader(f))
        if len(rows)==20 and all(r['fit_status']=='OK' for r in rows):return 'complete'
    if failure:
        if gate and not gate['passed'] and failure['error'].startswith('Chain checks failed'):
            return 'failed_chains'
        return 'needs_review'
    started=(data/'sampler'/f'scaling_{sim}.csv').exists()
    return 'active' if started else 'not_started'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--primary-label',required=True)
    ap.add_argument('--retry-labels',nargs='+',required=True)
    ap.add_argument('--between',nargs='+',type=int,required=True)
    ap.add_argument('--collection-label',required=True)
    ap.add_argument('--adopt-labels',nargs='*',default=[],help='Also accept completed results from these earlier retry runs')
    ap.add_argument('--adopt-active',nargs='*',type=int,default=[],help='Earlier retry samples that are still running externally')
    ap.add_argument('--workers',type=int,default=3)
    ap.add_argument('--first',type=int,default=1);ap.add_argument('--last',type=int,default=100)
    a=ap.parse_args()
    if len(a.retry_labels)!=len(a.between) or not 1<=a.workers<=8:
        ap.error('Provide one spacing per retry label and 1-8 workers')
    runner=module(ROOT/'code/execution/run_designscore_mlmi.py','retry_runner')
    summary=module(ROOT/'code/evaluation/summarize_designscore_mlmi.py','retry_summary')
    primary=read_json(ROOT/'provenance'/f'{a.primary_label}_inputs.json')
    config=primary['configuration']
    adopted={label:read_json(ROOT/'provenance'/f'{label}_inputs.json') for label in a.adopt_labels}
    if any(v<=config['between'] for v in a.between) or a.between!=sorted(set(a.between)):
        ap.error('Retry spacing must strictly increase beyond the primary run')
    lock=threading.Lock();original_freeze=runner.freeze_inputs
    def locked_freeze(*args,**kwargs):
        with lock:return original_freeze(*args,**kwargs)
    runner.freeze_inputs=locked_freeze
    def arguments(label,spacing,first,last):
        values=dict(config,between=spacing)
        result=['--stage','all','--run-label',label,'--first',str(first),'--last',str(last),'--workers','1']
        for k,v in values.items():result.extend(['--'+k.replace('_','-'),str(v)])
        return result
    # Register the complete sample set before concurrent per-replicate jobs.
    for label,spacing in zip(a.retry_labels,a.between):
        args=arguments(label,spacing,a.first,a.last)
        runner.freeze_inputs(ROOT,runner.options(args),summary)
    labels=[a.primary_label,*a.retry_labels]
    source_labels=[a.primary_label,*a.adopt_labels,*a.retry_labels]
    process=dict(process_id=os.getpid(),started_local=datetime.now().isoformat(),workers=a.workers,
                 primary=a.primary_label,retries=a.retry_labels,adopted=a.adopt_labels,adopt_active=a.adopt_active)
    (ROOT/'logs/oct2026'/a.primary_label/'completion_process.json').write_text(json.dumps(process,indent=2))
    print(json.dumps(process),flush=True)
    futures={};last_status=None;scheduler_errors={}
    with ThreadPoolExecutor(max_workers=a.workers) as executor:
        while True:
            for key,future in list(futures.items()):
                if future.done():
                    try:future.result()
                    except Exception as exc:
                        scheduler_errors[key[1]]=str(exc)
                        print(f'Retry scheduler error {key}: {exc}',flush=True)
                    del futures[key]
            complete={};review=[];exhausted=[];pending=0
            for sim in range(a.first,a.last+1):
                if sim in scheduler_errors:
                    review.append(sim);continue
                observed={label:state(label,sim) for label in source_labels}
                completed=[label for label in source_labels if observed[label]=='complete']
                if completed:
                    complete[sim]=completed[0];continue
                states=[observed[label] for label in labels]
                for label,manifest in adopted.items():
                    spacing=manifest['configuration']['between']
                    if spacing in a.between:
                        index=1+a.between.index(spacing)
                        value=observed[label]
                        if states[index]=='not_started' and (value in ('failed_chains','needs_review') or (value=='active' and sim in a.adopt_active)):
                            states[index]=value
                if 'needs_review' in states:
                    review.append(sim);continue
                if all(s=='failed_chains' for s in states):
                    exhausted.append(sim);continue
                # Only a recorded chain failure makes a longer run eligible.
                for j in range(1,len(labels)):
                    if states[j-1]=='failed_chains' and states[j]=='not_started':
                        key=(labels[j],sim)
                        if key not in futures:
                            print(f'Longer-chain retry: sample {sim}, {labels[j]}, spacing {a.between[j-1]}',flush=True)
                            futures[key]=executor.submit(runner.main,arguments(labels[j],a.between[j-1],sim,sim))
                        break
                pending+=1
            status=dict(completed=len(complete),requested=a.last-a.first+1,
                        completed_by_run={label:sum(s==label for s in complete.values()) for label in source_labels},
                        scheduled_retries=len(futures),needs_review=review,exhausted=exhausted)
            if status!=last_status:
                print(json.dumps(dict(local_time=datetime.now().isoformat(timespec='seconds'),**status)),flush=True)
                (ROOT/'logs/oct2026'/a.primary_label/'completion_progress.json').write_text(json.dumps(status,indent=2))
                last_status=status
            if len(complete)==a.last-a.first+1:
                break
            if not pending and not futures:
                raise RuntimeError('All runs ended but some replicates still require review')
            time.sleep(30)
    si=subprocess.STARTUPINFO();si.dwFlags|=subprocess.STARTF_USESHOWWINDOW;si.wShowWindow=subprocess.SW_HIDE
    command=[sys.executable,'code/execution/finalize_designscore_mlmi.py',
             '--run-label',a.collection_label,'--source-runs',*source_labels,
             '--first',str(a.first),'--last',str(a.last)]
    result=subprocess.run(command,cwd=ROOT,startupinfo=si,creationflags=subprocess.CREATE_NO_WINDOW)
    return result.returncode


if __name__=='__main__':
    raise SystemExit(main())
