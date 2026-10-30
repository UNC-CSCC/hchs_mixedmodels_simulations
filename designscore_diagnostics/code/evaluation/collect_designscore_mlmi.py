"""Collect complete paired replicates when chain lengths differ across runs.

The model, priors, imputation count, chain count and acceptance thresholds must
agree. No failed sampler is used, and no model is fitted by this module.
"""
from __future__ import annotations
from datetime import datetime, timezone
import importlib.util
import json
import math
from pathlib import Path
import statistics as st


def module(path, name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj)
    return obj


def collect(root, label, source_labels, first, last):
    root=Path(root)
    pooler=module(root/'code/evaluation/summarize_designscore_mlmi.py','collection_pooler')
    comparison=module(root/'code/evaluation/compare_multilevel_with_current_mi.py','collection_comparison')
    base=comparison.base
    manifests={s:json.loads((root/'provenance'/f'{s}_inputs.json').read_text()) for s in source_labels}
    common=None
    for s,m in manifests.items():
        config={k:v for k,v in m['configuration'].items() if k not in ('burn','between','every')}
        target=(m['model'],m['sampler'],config,m['source_sha256'],m['frame_sha256'],m['counts_sha256'],m['truth_sha256'])
        if common is None:common=target
        elif target!=common:raise ValueError('Collection runs differ in model, priors, code or inference settings')
        for rel,h in m['source_sha256'].items():
            if base.digest(root/rel)!=h:raise ValueError(f'Changed source {rel}')
        for rel,key in [('data/design_frames/oct2026_households.csv','frame_sha256'),
                        ('data/design_frames/oct2026_subject_counts.csv','counts_sha256'),
                        ('data/population_coefficients/population_coef_oct2026_hhsub.csv','truth_sha256')]:
            if base.digest(root/rel)!=m[key]:raise ValueError(f'Changed design input {rel}')
    requested=last-first+1
    mcount=common[2]['imputations']
    selections={};rows=[];audits=[];hashes={}
    for sim in range(first,last+1):
        chosen=None
        # Prefer the first completed, valid run. Later runs are retries of
        # failed checks, never selected by estimates or proximity to truth.
        for source in source_labels:
            data=root/'data/estimates/oct2026'/source
            gate=data/'sampler'/f'validation_{sim}.json'
            if not gate.exists() or not (data/'sampler'/f'paired_fit_{sim}.json').exists():continue
            if not json.loads(gate.read_text())['passed']:continue
            candidate=pooler.pool_replicate(root,source,sim,mcount)
            if any(r['fit_status']!='OK' for r in candidate):continue
            chosen=source;break
        if chosen is None:raise ValueError(f'Replicate {sim} has no complete validated paired result')
        manifest=manifests[chosen]
        sample=Path(manifest['sample_directory'])/f'samplemiss_{sim}.csv'
        if base.digest(sample)!=manifest['samples'][str(sim)]['sha256']:
            raise ValueError(f'Sample {sim} changed after imputation')
        selections[str(sim)]=dict(run_label=chosen,configuration=manifest['configuration'])
        for r in candidate:rows.append(dict(r,source_run=chosen))
        data=root/'data/estimates/oct2026'/chosen
        a=base.read(data/'audits'/f'mi_audit_{sim}.csv')
        audits.append(dict(sim_id=sim,status='OK',n_imputations=mcount,paired_inputs_verified=1,
            max_kkt_error=max(float(x['kkt_error']) for x in a),
            max_step1_score_error=max(float(x['step1_score_error']) for x in a),
            max_jacobian_condition=max(float(x['jacobian_condition']) for x in a),
            mean_missing_v2=st.mean(float(x['missing_v2']) for x in a),
            mean_missing_v3=st.mean(float(x['missing_v3']) for x in a),source_run=chosen))
        for p in [data/'per_imputation'/f'mi_estimates_{sim}.csv',data/'pooled'/f'mi_coef_{sim}.csv',
                  data/'audits'/f'mi_audit_{sim}.csv',data/'sampler'/f'validation_{sim}.json',
                  data/'sampler'/f'paired_fit_{sim}.json',data/'imputations'/f'mi_{sim}.dta']:
            hashes[p.relative_to(root).as_posix()]=base.digest(p)
    truth={r['term']:float(r['estimate']) for r in base.read(root/'data/population_coefficients/population_coef_oct2026_hhsub.csv')}
    comparison.validate(rows,first,last,truth)
    coverage=comparison.summarize(rows,truth,f'multilevel_MI_M{mcount}')
    for r in coverage:
        r.update(scenario=f"scen3_{r['method']}_mlmi_joint6_oct2026",n_expected=requested,
                 n_missing_points=0,n_valid_intervals=requested,n_missing_intervals=0,
                 coverage_over_all_requested=r['coverage_rate'],n_imputations=mcount,
                 first_replicate=first,last_replicate=last,
                 term_type='variance_component' if r['term'].startswith('var_') or r['term']=='s2' else 'fixed_effect')
        values=[x for x in rows if x['method']==r['method'] and x['term']==r['term']]
        r.update(relative_se_diff=(r['model_se']-r['empirical_se'])/r['empirical_se'],
                 coverage_emp_se=st.mean(abs(x['estimate']-r['true_value'])<=base.NORMAL_CRITICAL*r['empirical_se'] for x in values),
                 p_reject_h0=st.mean(not(x['conf_low']<=0<=x['conf_high']) for x in values),
                 pooling_scale=values[0]['pooling_scale'],mean_arithmetic_estimate=st.mean(x['arithmetic_mean_estimate'] for x in values))
    paired=[]
    for term,t in truth.items():
        original={r['sim_id']:r for r in rows if r['term']==term and r['method']=='original'}
        design={r['sim_id']:r for r in rows if r['term']==term and r['method']=='designscore'}
        d=[design[i]['estimate']-original[i]['estimate'] for i in range(first,last+1)]
        c=[int(design[i]['conf_low']<=t<=design[i]['conf_high'])-int(original[i]['conf_low']<=t<=original[i]['conf_high']) for i in range(first,last+1)]
        paired.append(dict(term=term,n_pairs=requested,mean_designscore_minus_original=st.mean(d),
                           paired_difference_mcse=st.stdev(d)/math.sqrt(requested),
                           coverage_difference=st.mean(c),coverage_difference_mcse=st.stdev(c)/math.sqrt(requested)))
    dest=root/'summaries/oct2026'/f'{label}_{first}_{last}'
    for method in base.METHODS:
        base.write(dest/f'coverage_scen3_{method}_mi.csv',[r for r in coverage if r['method']==method])
    base.write(dest/'pooled_estimates.csv',rows);base.write(dest/'audit_mi.csv',audits)
    base.write(dest/'paired_comparison.csv',paired)
    base.write(dest/'replicate_source_runs.csv',[dict(sim_id=int(i),source_run=s['run_label'],**s['configuration']) for i,s in selections.items()])
    collection=dict(model=common[0],sampler=common[1],configuration=common[2],
                    replicate_sources=selections,source_sha256=common[3],
                    source_input_manifests={s:dict(path=f'provenance/{s}_inputs.json',sha256=base.digest(root/'provenance'/f'{s}_inputs.json')) for s in source_labels},
                    selection_rule='First completed validated run in supplied order; longer chains only for failed validation, never chosen using estimates')
    (root/'provenance'/f'{label}_inputs.json').write_text(json.dumps(collection,indent=2),encoding='utf-8')
    provenance=dict(run_label=label,first_replicate=first,last_replicate=last,n_imputations=mcount,
                    completed_replicates=requested,created_utc=datetime.now(timezone.utc).isoformat(),
                    replicate_sources=selections,output_sha256=hashes,
                    collection_code_sha256=base.digest(Path(__file__)),
                    note='Same model, priors and M across replicates; some chains extended to meet unchanged convergence criteria. Source data are referenced in place, not copied.')
    (dest/'provenance_mi.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    return dest,coverage,audits,collection
