"""Pool a paired MI run and summarize explicitly requested simulation replicates.

Both methods use Rubin pooling with finite-imputation t degrees of freedom.
Positive variances are pooled consistently on the log scale. No failed
imputation is silently dropped. No model is fitted by this module.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics as st
import subprocess

NORMAL_CRITICAL = st.NormalDist().inv_cdf(.975)


def t_quantiles(dfs):
    """Use installed R's qt for numerical tail inversion, once per replicate."""
    rscript = os.environ.get('R_SCRIPT', 'C:/Program Files/R/R-4.5.0/bin/Rscript.exe')
    options = {'creationflags': subprocess.CREATE_NO_WINDOW} if os.name == 'nt' else {}
    command = [rscript, '--vanilla', '-e',
               'options(digits=17); d<-as.numeric(commandArgs(TRUE)); cat(qt(.975,d),sep="\\n")',
               *('Inf' if not math.isfinite(df) else str(df) for df in dfs)]
    result = subprocess.run(command, capture_output=True, text=True, check=True, **options)
    values = [float(v) for v in result.stdout.split()]
    if len(values) != len(dfs):
        raise ValueError('R returned an unexpected number of t quantiles')
    return values

METHODS = ('original', 'designscore')


def read(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def write(path, records):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    with tmp.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        for record in records:
            writer.writerow({k: '' if isinstance(v, float) and not math.isfinite(v) else v
                             for k, v in record.items()})
    tmp.replace(path)


def pool(values, errors, *, log_scale=False, critical_value=None):
    if len(values) < 2 or len(errors) != len(values):
        raise ValueError('At least two paired estimates and SEs are required')
    if not all(math.isfinite(v) for v in values) or not all(math.isfinite(s) and s > 0 for s in errors):
        raise ValueError('Incomplete estimates or invalid standard errors')
    if log_scale and min(values) <= 0:
        raise ValueError('Log pooling requires positive variance estimates')
    q = [math.log(v) for v in values] if log_scale else values
    u = [(s / v)**2 for s, v in zip(errors, values)] if log_scale else [s*s for s in errors]
    m = len(values)
    qbar, ubar, between = st.mean(q), st.mean(u), st.variance(q)
    inflation = (1 + 1/m) * between
    total = ubar + inflation
    df = (m-1) * (1 + ubar/inflation)**2 if inflation > 0 else math.inf
    critical = critical_value if critical_value is not None else t_quantiles([df])[0]
    sd = math.sqrt(total)
    estimate = math.exp(qbar) if log_scale else qbar
    lower, upper = qbar-critical*sd, qbar+critical*sd
    return dict(estimate=estimate, std_error=estimate*sd if log_scale else sd,
                conf_low=math.exp(lower) if log_scale else lower,
                conf_high=math.exp(upper) if log_scale else upper,
                pooling_scale='log' if log_scale else 'identity',
                arithmetic_mean_estimate=st.mean(values), within_variance=ubar,
                between_variance=between, total_variance=total, mi_df=df,
                fraction_variance_between=inflation/total, critical_value=critical)


def pool_replicate(root, label, sim, m):
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
                           'imputation_model':'scen4_s1_v2', 'mi_seed':100000+sim,
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
            coverage.append(dict(scenario=f'scen3_{method}_mi_s1_v2_oct2026',term_type='variance_component' if term.startswith('var_') or term=='s2' else 'fixed_effect',term=term,n_reps=n,true_value=true,mean_estimate=mean,empirical_bias=mean-true,relative_bias=(mean-true)/true if true else math.nan,empirical_se=empse,model_se=modelse,relative_se_diff=(modelse-empse)/empse if empse>0 else math.nan,coverage_rate=cov,coverage_emp_se=st.mean(abs(v-true)<=1.95996398454*empse for v in values) if n>1 else math.nan,p_reject_h0=st.mean(not(r['conf_low']<=0<=r['conf_high']) for r in rows) if n else math.nan,n_expected=expected,n_missing_points=expected-n,n_valid_intervals=n,n_missing_intervals=expected-n,coverage_over_all_requested=sum(covered)/expected,bias_mcse=empse/math.sqrt(n) if n>1 else math.nan,coverage_mcse=math.sqrt(cov*(1-cov)/n) if n else math.nan,n_imputations=m,pooling_scale='log' if term.startswith('var_') or term=='s2' else 'identity',mean_arithmetic_estimate=st.mean(r['arithmetic_mean_estimate'] for r in rows) if n else math.nan,first_replicate=first,last_replicate=last))
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
                    imputation_model='Existing scen4 s1/v2 wide chained linear regression; default burn-in 10; seed 100000+replicate',
                    missingness='y_gfr and y_bmi masked using miss_ind_mar; baseline 00pct',
                    pooling='Rubin rules; identity for fixed effects, log for positive variances; finite-M t degrees of freedom, asymptotic complete-data variance; both methods identical rules',
                    completed_replicates=sum(r['status']=='OK' for r in auditrows),
                    output_sha256=hashes,
                    source_sha256={p.relative_to(root).as_posix():digest(p) for p in [root/'code/stata/scen3_designscore_mi.do',root/'code/stata/mi_rename_compat.do',root/'code/stata/variance_designscore.mata',Path(__file__),root/'code/execution/run_designscore_mi.py']},
                    truth_sha256=digest(root/'data/population_coefficients/population_coef_oct2026_hhsub.csv'))
    manifest=root/'provenance'/f'{label}_inputs.json'
    if manifest.exists():
        provenance['input_manifest']=manifest.relative_to(root).as_posix()
        provenance['input_manifest_sha256']=digest(manifest)
    (destination/'provenance_mi.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    return destination,coverage,auditrows


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    ap.add_argument('--run-label',default='mi_s1_v2_m5')
    ap.add_argument('--first',type=int,default=1);ap.add_argument('--last',type=int,default=100)
    ap.add_argument('--imputations',type=int,default=5)
    a=ap.parse_args()
    dest,_,audits=summarize(a.root,a.run_label,a.first,a.last,a.imputations)
    print(json.dumps({'summaries':str(dest),'complete':sum(r['status']=='OK' for r in audits),'requested':len(audits)}))
