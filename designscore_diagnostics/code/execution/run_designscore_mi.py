"""Run paired original/design-score MI analyses locally; preserve completed fits."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import threading
import time


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first',type=int,default=1);parser.add_argument('--last',type=int,default=100)
    parser.add_argument('--workers',type=int,default=8);parser.add_argument('--imputations',type=int,default=5)
    parser.add_argument('--run-label',default='mi_s1_v2_m5')
    parser.add_argument('--stata',type=Path,default=Path('C:/Program Files/StataNow19/StataBE-64.exe'))
    args=parser.parse_args()
    assert 1<=args.first<=args.last<=1000 and 1<=args.workers<=16 and args.imputations>=2
    assert args.run_label.replace('_','').isalnum()
    root=Path(__file__).resolve().parents[2]
    spec=importlib.util.spec_from_file_location('summary',root/'code/evaluation/summarize_designscore_mi.py')
    summary=importlib.util.module_from_spec(spec);spec.loader.exec_module(summary)
    data=root/'data/estimates/oct2026'/args.run_label
    logs=root/'logs/oct2026'/args.run_label
    state=root/'.statatmp'/args.run_label
    for folder in [*(data/p for p in ['imputations','per_imputation','pooled','audits']),logs,state/'runtime_tmp']:
        folder.mkdir(parents=True,exist_ok=True)
    samples=root.parents[1]/'data/Oct2026/sample_miss00pct'
    countfile=root/'data/design_frames/oct2026_subject_counts.csv'
    countrows={int(r['sim_id']):r for r in summary.read(countfile)}
    manifest={'sample_directory':str(samples),'imputation_model':'scen4_s1_v2','M':args.imputations,'samples':{}}
    for sim in range(args.first,args.last+1):
        source=samples/f'samplemiss_{sim}.csv'
        with source.open('rb') as f:
            md5=hashlib.file_digest(f,'md5').hexdigest().upper()
        assert md5==countrows[sim]['sample_md5'] and countrows[sim]['ids_verified']=='TRUE',sim
        manifest['samples'][str(sim)]={'md5':md5,'sha256':summary.digest(source),'mi_seed':100000+sim}
    manifest['frame_sha256']=summary.digest(root/'data/design_frames/oct2026_households.csv')
    manifest['counts_sha256']=summary.digest(countfile)
    manifest['original_mi_script_sha256']=summary.digest(root.parent/'scripts/stata/scen4_mi.do')
    manifest['analysis_source_sha256']={p:summary.digest(root/p) for p in [
        'code/stata/scen3_designscore_mi.do','code/stata/mi_rename_compat.do',
        'code/stata/variance_designscore.mata']}
    manifest_path=root/'provenance'/f'{args.run_label}_inputs.json'
    if manifest_path.exists():
        prior=json.loads(manifest_path.read_text())
        assert prior['M']==manifest['M'] and prior['imputation_model']==manifest['imputation_model']
        for key in ['frame_sha256','counts_sha256','original_mi_script_sha256','analysis_source_sha256']:
            if key in prior: assert prior[key]==manifest[key],f'Changed run input: {key}'
        for sim,item in manifest['samples'].items():
            if sim in prior['samples']: assert prior['samples'][sim]==item
        samples_merged={**prior['samples'],**manifest['samples']}
        manifest={**prior,**manifest,'samples':samples_merged}
    manifest_path.write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(f'Verified source samples {args.first}-{args.last}; launching up to {args.workers} Stata processes.',flush=True)
    startup=subprocess.STARTUPINFO();startup.dwFlags|=subprocess.STARTF_USESHOWWINDOW;startup.wShowWindow=subprocess.SW_HIDE
    env=dict(os.environ,STATATMP=str(state/'runtime_tmp'),TMP=str(state/'runtime_tmp'),TEMP=str(state/'runtime_tmp'),PYTHONDONTWRITEBYTECODE='1')
    progress={'started_at':datetime.now(timezone.utc).isoformat(),'run_label':args.run_label,'first':args.first,'last':args.last,'imputations':args.imputations,'jobs':{}}
    lock=threading.Lock()
    def update(sim,stage,**kwargs):
        with lock:
            progress['jobs'][str(sim)]={'stage':stage,**kwargs}
            tmp=state/'status.tmp';tmp.write_text(json.dumps(progress,indent=2),encoding='utf-8')
            for attempt in range(20):
                try:
                    tmp.replace(state/'status.json');break
                except PermissionError:
                    time.sleep(.1)
        print(f'Replicate {sim}: {stage} {kwargs}',flush=True)
    def job(sim):
        output=data/'per_imputation'/f'mi_estimates_{sim}.csv'
        audit=data/'audits'/f'mi_audit_{sim}.csv'
        if output.exists() and audit.exists():
            rows=summary.pool_replicate(root,args.run_label,sim,args.imputations)
            ok=all(r['fit_status']=='OK' for r in rows)
            update(sim,'existing' if ok else 'incomplete_inference')
            return ok
        wrapper=state/f'replicate_{sim}.do'
        wrapper.write_text('version 19\nset more off\n'+f'cd "{root.as_posix()}"\n'+f'do code/stata/scen3_designscore_mi.do {sim} "{root.as_posix()}" {args.run_label} {args.imputations}\nexit, clear\n',encoding='utf-8')
        started=time.monotonic()
        proc=subprocess.Popen([str(args.stata),'/e','do',str(wrapper)],cwd=logs,env=env,startupinfo=startup,creationflags=subprocess.CREATE_NO_WINDOW)
        update(sim,'running',pid=proc.pid)
        try:
            rc=proc.wait(timeout=7200)
        except subprocess.TimeoutExpired:
            proc.kill();proc.wait();update(sim,'timeout');return False
        if not output.exists() or not audit.exists():
            update(sim,'failed',return_code=rc,seconds=round(time.monotonic()-started,1));return False
        rows=summary.pool_replicate(root,args.run_label,sim,args.imputations)
        ok=rc==0 and all(r['fit_status']=='OK' for r in rows)
        update(sim,'complete' if ok else 'incomplete_inference',return_code=rc,seconds=round(time.monotonic()-started,1))
        return ok
    failed=[]
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures={executor.submit(job,sim):sim for sim in range(args.first,args.last+1)}
        for future in as_completed(futures):
            sim=futures[future]
            try:
                if not future.result():failed.append(sim)
            except Exception as exc:
                update(sim,'error',message=str(exc));failed.append(sim)
    dest,_,_=summary.summarize(root,args.run_label,args.first,args.last,args.imputations)
    print(json.dumps({'summaries':str(dest),'failed':sorted(failed)}),flush=True)
    return bool(failed)


if __name__=='__main__':
    raise SystemExit(main())
