"""Pool multilevel MI outputs with the unchanged legacy Rubin formulas.

No sampler or analysis model is run. Labels, seeds and provenance identify the
new imputer. Estimates are accepted only with matching saved-data/fit seals.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import math
from pathlib import Path
import statistics as st

_spec = importlib.util.spec_from_file_location(
    'legacy_mi_pooling', Path(__file__).with_name('summarize_designscore_mi.py'))
_legacy = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_legacy)
read, write, number, digest = _legacy.read, _legacy.write, _legacy.number, _legacy.digest
pool, t_quantiles = _legacy.pool, _legacy.t_quantiles
METHODS, NORMAL_CRITICAL = _legacy.METHODS, _legacy.NORMAL_CRITICAL


def validate_multilevel_output(root, label, sim, m):
    root = Path(root)
    if not label.startswith('mlmi_') or '/' in label or '\\' in label:
        raise ValueError('A distinct multilevel run label is required')
    data = root/'data/estimates/oct2026'/label
    manifest = json.loads((root/'provenance'/f'{label}_inputs.json').read_text())
    if manifest['model'] != 'joint6_psu_household_v1' or manifest['configuration']['imputations'] != m:
        raise ValueError('Multilevel input manifest mismatch')
    for path, expected in manifest['source_sha256'].items():
        if digest(root/path) != expected:
            raise ValueError(f'Changed source: {path}')
    gate = json.loads((data/'sampler'/f'validation_{sim}.json').read_text())
    if not gate['passed'] or gate['configuration'] != manifest['configuration']:
        raise ValueError('No passing multilevel MI validation for these settings')
    current = digest(data/'imputations'/f'mi_{sim}.dta')
    if current != gate['imputation_sha256']:
        raise ValueError('Saved imputations changed after validation')
    for name, expected in gate['trace_sha256'].items():
        if digest(data/'sampler'/name) != expected:
            raise ValueError('Sampler trace changed after validation')
    seal = json.loads((data/'sampler'/f'paired_fit_{sim}.json').read_text())
    if seal != dict(imputation_sha256=current,
                    estimates_sha256=digest(data/'per_imputation'/f'mi_estimates_{sim}.csv'),
                    audit_sha256=digest(data/'audits'/f'mi_audit_{sim}.csv')):
        raise ValueError('Paired fits do not match the validated shared imputations')


def pool_replicate(root, label, sim, m):
    validate_multilevel_output(root, label, sim, m)
    directory = Path(root)/'data/estimates/oct2026'/label
    records = read(directory/'per_imputation'/f'mi_estimates_{sim}.csv')
    audits = read(directory/'audits'/f'mi_audit_{sim}.csv')
    truth = read(Path(root)/'data/population_coefficients/population_coef_oct2026_hhsub.csv')
    terms = [r['term'] for r in truth]
    expected = {(method, imp, term) for method in METHODS for imp in range(1, m+1) for term in terms}
    keys = [(r['method'], int(r['imputation']), r['term']) for r in records]
    if set(keys) != expected or len(keys) != len(expected):
        raise ValueError(f'Replicate {sim}: incomplete or duplicate imputation estimates')
    if not all(int(r['sim_id']) == sim for r in records + audits):
        raise ValueError('Replicate identifier mismatch')
    if sorted(int(r['imputation']) for r in audits) != list(range(1,m+1)):
        raise ValueError('Incomplete audit')
    for imp in range(1,m+1):
        sigs = {r['input_signature'] for r in records if int(r['imputation'])==imp}
        sigs.update(r['input_signature'] for r in audits if int(r['imputation'])==imp)
        if len(sigs)!=1 or not next(iter(sigs)):
            raise ValueError('The two analyses did not use identical completed input')
        for term in terms[:7] + ['s2']:
            pair = [r for r in records if int(r['imputation'])==imp and r['term']==term]
            if not math.isclose(number(pair[0]['estimate']), number(pair[1]['estimate']), rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError('Shared first-stage point estimates disagree')
    output=[]
    for method in METHODS:
        for term in terms:
            selected = sorted((r for r in records if r['method']==method and r['term']==term), key=lambda r:int(r['imputation']))
            base=dict(sim_id=sim, population='oct2026', method=method, term=term,
                      n_imputations_requested=m, n_imputations_available=len(selected))
            vals=[number(r['estimate']) for r in selected]
            ses=[number(r['std_error']) for r in selected]
            valid=all(r['inference_usable']=='1' for r in selected)
            try:
                if not valid:
                    raise ValueError('A within-imputation fit lacks usable inference')
                # Replace provisional critical values together after all terms.
                result=pool(vals,ses,log_scale=term.startswith('var_') or term=='s2',critical_value=NORMAL_CRITICAL)
                status='OK'
            except ValueError:
                result={k:math.nan for k in ['estimate','std_error','conf_low','conf_high','arithmetic_mean_estimate','within_variance','between_variance','total_variance','mi_df','fraction_variance_between','critical_value']}
                result['pooling_scale']='log' if term.startswith('var_') or term=='s2' else 'identity'
                status='INCOMPLETE_INFERENCE'
            output.append({**base,**result,'fit_status':status,'inference_usable':int(status=='OK'),
                           'imputation_model':'joint6_psu_household_v1', 'mi_seed':500000+sim,
                           'within_se_method':'household_cluster_original' if method=='original' else 'stacked_stratified_psu_wr'})
    usable = [row for row in output if row['fit_status']=='OK']
    for row, critical in zip(usable,t_quantiles([row['mi_df'] for row in usable]) if usable else []):
        qbar = math.log(row['estimate']) if row['pooling_scale']=='log' else row['estimate']
        sd = math.sqrt(row['total_variance'])
        lower, upper = qbar-critical*sd, qbar+critical*sd
        row['conf_low'] = math.exp(lower) if row['pooling_scale']=='log' else lower
        row['conf_high'] = math.exp(upper) if row['pooling_scale']=='log' else upper
        row['critical_value'] = critical
    write(directory/'pooled'/f'mi_coef_{sim}.csv', output)
    return output


def summarize(root, label, first, last, m):
    root=Path(root)
    data=root/'data/estimates/oct2026'/label
    destination=root/'summaries/oct2026'/f'{label}_{first}_{last}'
    truth={r['term']:number(r['estimate']) for r in read(root/'data/population_coefficients/population_coef_oct2026_hhsub.csv')}
    allrows=[]; auditrows=[]; hashes={}; repmap={}
    for sim in range(first,last+1):
        path=data/'per_imputation'/f'mi_estimates_{sim}.csv'
        if not path.exists() or not (data/'audits'/f'mi_audit_{sim}.csv').exists():
            auditrows.append(dict(sim_id=sim,status='MISSING_OUTPUT',n_imputations=0,paired_inputs_verified=0,max_kkt_error='',max_step1_score_error='',max_jacobian_condition='',mean_missing_v2='',mean_missing_v3=''))
            continue
        rows=pool_replicate(root,label,sim,m)
        allrows.extend(rows); repmap[sim]={(r['method'],r['term']):r for r in rows}
        a=read(data/'audits'/f'mi_audit_{sim}.csv')
        auditrows.append(dict(sim_id=sim,status='OK' if all(r['fit_status']=='OK' for r in rows) else 'INCOMPLETE_INFERENCE',n_imputations=m,paired_inputs_verified=1,max_kkt_error=max(number(r['kkt_error']) for r in a),max_step1_score_error=max(number(r['step1_score_error']) for r in a),max_jacobian_condition=max(number(r['jacobian_condition']) for r in a),mean_missing_v2=st.mean(number(r['missing_v2']) for r in a),mean_missing_v3=st.mean(number(r['missing_v3']) for r in a)))
        for item in [path,data/'audits'/f'mi_audit_{sim}.csv',data/'pooled'/f'mi_coef_{sim}.csv',data/'imputations'/f'mi_{sim}.dta']:
            hashes[item.relative_to(root).as_posix()]=digest(item)
    expected=last-first+1
    coverage=[]; paired=[]
    for method in METHODS:
        for term,true in truth.items():
            rows=[r for r in allrows if r['method']==method and r['term']==term and r['fit_status']=='OK']
            n=len(rows)
            values=[r['estimate'] for r in rows]
            mean=st.mean(values) if n else math.nan
            empse=st.stdev(values) if n>1 else math.nan
            modelse=st.mean(r['std_error'] for r in rows) if n else math.nan
            covered=[r['conf_low']<=true<=r['conf_high'] for r in rows]
            cov=st.mean(covered) if n else math.nan
            coverage.append(dict(scenario=f'scen3_{method}_mlmi_joint6_oct2026',term_type='variance_component' if term.startswith('var_') or term=='s2' else 'fixed_effect',term=term,n_reps=n,true_value=true,mean_estimate=mean,empirical_bias=mean-true,relative_bias=(mean-true)/true if true else math.nan,empirical_se=empse,model_se=modelse,relative_se_diff=(modelse-empse)/empse if empse>0 else math.nan,coverage_rate=cov,coverage_emp_se=st.mean(abs(v-true)<=1.95996398454*empse for v in values) if n>1 else math.nan,p_reject_h0=st.mean(not(r['conf_low']<=0<=r['conf_high']) for r in rows) if n else math.nan,n_expected=expected,n_missing_points=expected-n,n_valid_intervals=n,n_missing_intervals=expected-n,coverage_over_all_requested=sum(covered)/expected,bias_mcse=empse/math.sqrt(n) if n>1 else math.nan,coverage_mcse=math.sqrt(cov*(1-cov)/n) if n else math.nan,n_imputations=m,pooling_scale='log' if term.startswith('var_') or term=='s2' else 'identity',mean_arithmetic_estimate=st.mean(r['arithmetic_mean_estimate'] for r in rows) if n else math.nan,first_replicate=first,last_replicate=last))
    for term,true in truth.items():
        pairs=[(d[('original',term)],d[('designscore',term)]) for d in repmap.values() if d[('original',term)]['fit_status']=='OK' and d[('designscore',term)]['fit_status']=='OK']
        n=len(pairs)
        delta=[b['estimate']-a['estimate'] for a,b in pairs]
        cd=[int(b['conf_low']<=true<=b['conf_high'])-int(a['conf_low']<=true<=a['conf_high']) for a,b in pairs]
        paired.append(dict(term=term,n_pairs=n,mean_designscore_minus_original=st.mean(delta) if n else math.nan,paired_difference_mcse=st.stdev(delta)/math.sqrt(n) if n>1 else math.nan,coverage_difference=st.mean(cd) if n else math.nan,coverage_difference_mcse=st.stdev(cd)/math.sqrt(n) if n>1 else math.nan))
    for method in METHODS:
        write(destination/f'coverage_scen3_{method}_mi.csv',[r for r in coverage if r['scenario'].startswith(f'scen3_{method}_')])
    write(destination/'paired_comparison.csv',paired)
    write(destination/'audit_mi.csv',auditrows)
    if allrows:
        write(destination/'pooled_estimates.csv',allrows)
    provenance=dict(run_label=label,first_replicate=first,last_replicate=last,n_imputations=m,
                    imputation_model='Joint six-outcome Gaussian MI with PSU and household random effects; configuration and chain seeds in input manifest',
                    missingness='y_gfr and y_bmi masked using miss_ind_mar; baseline 00pct',
                    pooling='Rubin rules; identity for fixed effects, log for positive variances; finite-M t degrees of freedom, asymptotic complete-data variance; both methods identical rules',
                    completed_replicates=sum(r['status']=='OK' for r in auditrows),
                    output_sha256=hashes,
                    source_sha256={p.relative_to(root).as_posix():digest(p) for p in [root/'code/stata/scen3_designscore_mi.do',root/'code/stata/mi_rename_compat.do',root/'code/stata/variance_designscore.mata',Path(__file__),root/'code/execution/run_designscore_mlmi.py',root/'code/stata/multilevel_mi.mata',root/'code/stata/impute_multilevel_mi.do',root/'code/evaluation/check_multilevel_mi.py',root/'code/evaluation/summarize_designscore_mi.py']},
                    truth_sha256=digest(root/'data/population_coefficients/population_coef_oct2026_hhsub.csv'))
    manifest=root/'provenance'/f'{label}_inputs.json'
    if manifest.exists():
        provenance['multilevel_configuration']=json.loads(manifest.read_text())['configuration']
        provenance['input_manifest']=manifest.relative_to(root).as_posix()
        provenance['input_manifest_sha256']=digest(manifest)
    (destination/'provenance_mi.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    return destination,coverage,auditrows


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    ap.add_argument('--run-label',default='mlmi_joint6_m20_v1')
    ap.add_argument('--first',type=int,default=1);ap.add_argument('--last',type=int,default=100)
    ap.add_argument('--imputations',type=int,default=20)
    a=ap.parse_args()
    dest,_,audits=summarize(a.root,a.run_label,a.first,a.last,a.imputations)
    print(json.dumps({'summaries':str(dest),'complete':sum(r['status']=='OK' for r in audits),'requested':len(audits)}))
