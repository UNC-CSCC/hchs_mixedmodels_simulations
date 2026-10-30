"""Compare matched saved MI results; never imputes or fits an analysis model.

The five-draw sensitivity subset is fixed before inspecting results: the first
retained draw of each of four chains, then the second draw of chain 1. Both
analysis methods always use the same selected completed datasets.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import importlib.util
import json
import math
from pathlib import Path
import statistics as st


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('current_mi_summary', ROOT/'code/evaluation/summarize_designscore_mi.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def five_draw_sensitivity(root, label, first, last, m, chains=4):
    if chains != 4 or m % chains or m//chains < 2:
        raise ValueError('The prespecified five-draw sensitivity requires four chains with >=2 draws each')
    per_chain = m//chains
    selected = [1 + c*per_chain for c in range(chains)] + [2]
    terms = [r['term'] for r in base.read(root/'data/population_coefficients/population_coef_oct2026_hhsub.csv')]
    result = []
    input_manifest=root/'provenance'/f'{label}_inputs.json'
    sources=json.loads(input_manifest.read_text()).get('replicate_sources',{}) if input_manifest.exists() else {}
    for sim in range(first,last+1):
        actual_label=sources.get(str(sim),{}).get('run_label',label)
        source = base.read(root/'data/estimates/oct2026'/actual_label/'per_imputation'/f'mi_estimates_{sim}.csv')
        for imp in selected:
            sigs = {r['input_signature'] for r in source if int(r['imputation'])==imp}
            if len(sigs)!=1 or not next(iter(sigs)):
                raise ValueError('Shared imputation signatures disagree')
        for method in base.METHODS:
            for term in terms:
                rows = [r for r in source if r['method']==method and r['term']==term and int(r['imputation']) in selected]
                if len(rows)!=5 or len({r['imputation'] for r in rows})!=5 or any(r['inference_usable']!='1' for r in rows):
                    raise ValueError('Incomplete five-draw subset')
                pooled = base.pool([float(r['estimate']) for r in rows], [float(r['std_error']) for r in rows],
                                   log_scale=term.startswith('var_') or term=='s2', critical_value=base.NORMAL_CRITICAL)
                result.append(dict(sim_id=sim,method=method,term=term,**pooled,fit_status='OK',
                                   n_imputations_requested=5,selected_imputations=';'.join(map(str,selected))))
    for start in range(0,len(result),100):
        chunk = result[start:start+100]
        for row, critical in zip(chunk,base.t_quantiles([r['mi_df'] for r in chunk])):
            q = math.log(row['estimate']) if row['pooling_scale']=='log' else row['estimate']
            delta = critical*math.sqrt(row['total_variance'])
            row['critical_value']=critical
            row['conf_low']=math.exp(q-delta) if row['pooling_scale']=='log' else q-delta
            row['conf_high']=math.exp(q+delta) if row['pooling_scale']=='log' else q+delta
    return result,selected


def normalize(rows):
    for row in rows:
        row = dict(row)
        row['sim_id']=int(row['sim_id'])
        for key in ('estimate','std_error','conf_low','conf_high'):
            row[key]=float(row[key])
        yield row


def validate(rows, first, last, terms):
    expected={(sim,method,term) for sim in range(first,last+1) for method in base.METHODS for term in terms}
    keys=[(r['sim_id'],r['method'],r['term']) for r in rows]
    if set(keys)!=expected or len(keys)!=len(expected):
        raise ValueError('Comparison requires every requested replicate and parameter, exactly once')
    if any(r['fit_status']!='OK' or not all(math.isfinite(r[k]) for k in ('estimate','std_error','conf_low','conf_high')) for r in rows):
        raise ValueError('Comparison contains incomplete inference')


def summarize(rows, truth, name):
    out=[]
    for method in base.METHODS:
        for term,t in truth.items():
            r=[x for x in rows if x['method']==method and x['term']==term]
            values=[x['estimate'] for x in r]
            mean=st.mean(values); n=len(r); empirical=st.stdev(values)
            coverage=st.mean(x['conf_low']<=t<=x['conf_high'] for x in r)
            out.append(dict(imputation=name,method=method,term=term,n_reps=n,true_value=t,
                            mean_estimate=mean,empirical_bias=mean-t,relative_bias=(mean-t)/t if t else '',
                            empirical_se=empirical,model_se=st.mean(x['std_error'] for x in r),
                            coverage_rate=coverage,bias_mcse=empirical/math.sqrt(n),
                            coverage_mcse=math.sqrt(coverage*(1-coverage)/n)))
    return out


def paired_changes(old, new, truth, name):
    old={(r['sim_id'],r['method'],r['term']):r for r in old}
    out=[]
    for method in base.METHODS:
        for term,t in truth.items():
            pairs=[(old[(r['sim_id'],method,term)],r) for r in new if r['method']==method and r['term']==term]
            d=[b['estimate']-a['estimate'] for a,b in pairs]
            c=[int(b['conf_low']<=t<=b['conf_high'])-int(a['conf_low']<=t<=a['conf_high']) for a,b in pairs]
            out.append(dict(comparison=name,method=method,term=term,n_pairs=len(pairs),
                            mean_new_minus_current=st.mean(d),difference_mcse=st.stdev(d)/math.sqrt(len(d)),
                            coverage_difference=st.mean(c),coverage_difference_mcse=st.stdev(c)/math.sqrt(len(c))))
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--new-label',required=True)
    ap.add_argument('--current-label',default='mi_s1_v2_m5')
    ap.add_argument('--first',type=int,default=1); ap.add_argument('--last',type=int,default=100)
    a=ap.parse_args()
    root=ROOT
    manifest=json.loads((root/'provenance'/f'{a.new_label}_inputs.json').read_text())
    m=manifest['configuration']['imputations']; chains=manifest['configuration']['chains']
    dest=root/'summaries/oct2026'/f'{a.new_label}_{a.first}_{a.last}'
    oldpath=root/'summaries/oct2026'/f'{a.current_label}_{a.first}_{a.last}'/'pooled_estimates.csv'
    newpath=dest/'pooled_estimates.csv'
    old=list(normalize(base.read(oldpath))); new=list(normalize(base.read(newpath)))
    truth={r['term']:float(r['estimate']) for r in base.read(root/'data/population_coefficients/population_coef_oct2026_hhsub.csv')}
    validate(old,a.first,a.last,truth); validate(new,a.first,a.last,truth)
    small,selected=five_draw_sensitivity(root,a.new_label,a.first,a.last,m,chains)
    validate(small,a.first,a.last,truth)
    summary=(summarize(old,truth,'current_wide_MI_M5')+summarize(new,truth,f'multilevel_MI_M{m}')+
             summarize(small,truth,'multilevel_MI_M5_subset'))
    changes=paired_changes(old,new,truth,f'multilevel_M{m}_minus_current_M5')+paired_changes(old,small,truth,'multilevel_M5_minus_current_M5')
    base.write(dest/'comparison_with_current_mi.csv',summary)
    base.write(dest/'paired_changes_from_current_mi.csv',changes)
    base.write(dest/'pooled_multilevel_m5_sensitivity.csv',small)
    # Aggregate read-only residual diagnostics; each replicate gets equal weight.
    covariance=[]; chain_summary=[]
    reference_path=root/'data/run_artifacts'/f'{a.current_label}_covariance_reference_{a.first}_{a.last}.json'
    reference=None;reference_rows={}
    if reference_path.exists():
        reference=json.loads(reference_path.read_text())
        if (reference['label'],reference['first'],reference['last'])!=(a.current_label,a.first,a.last):
            raise ValueError('Current-MI covariance reference mismatch')
        for rel,expected in reference['source_sha256'].items():
            if base.digest(root/rel)!=expected:raise ValueError(f'Changed covariance reference input: {rel}')
        reference_rows={r['sim_id']:r for r in reference['records']}
        if set(reference_rows)!=set(range(a.first,a.last+1)):raise ValueError('Incomplete covariance reference')
    for sim in range(a.first,a.last+1):
        actual_label=manifest.get('replicate_sources',{}).get(str(sim),{}).get('run_label',a.new_label)
        gate=json.loads((root/'data/estimates/oct2026'/actual_label/'sampler'/f'validation_{sim}.json').read_text())
        if not gate['passed']: raise ValueError('A replicate failed multilevel validation')
        pars=gate['mcmc']['parameters']
        chain_summary.append(dict(sim_id=sim,max_rhat=max(p['rhat'] for p in pars),
                                  min_bulk_ess=min(p['bulk_ess'] for p in pars),min_tail_ess=min(p['tail_ess'] for p in pars),
                                  source_run=actual_label,burn=gate['configuration']['burn'],between=gate['configuration']['between']))
        cov=gate['integrity'].get('covariance_diagnostics')
        if cov:
            for term,full in cov['complete'].items():
                if full is not None:
                    completed=st.mean(r[term] for r in cov['imputed'])
                    row=dict(sim_id=sim,moment=term,complete=full,imputed=completed,change=completed-full)
                    if reference:
                        old_cov=reference_rows[sim]
                        if not math.isclose(full,old_cov['complete'][term],rel_tol=1e-9,abs_tol=1e-9):
                            raise ValueError('Covariance checks used different complete-data references')
                        current=old_cov['imputed_mean'][term]
                        row.update(current_mi=current,current_mi_change=current-full,new_minus_current=completed-current)
                    covariance.append(row)
    base.write(dest/'chain_diagnostics_summary.csv',chain_summary)
    if covariance: base.write(dest/'covariance_preservation.csv',covariance)
    rows={ (r['imputation'],r['method'],r['term']):r for r in summary }
    old_sub=rows[('current_wide_MI_M5','designscore','var_subid')]
    new_sub=rows[(f'multilevel_MI_M{m}','designscore','var_subid')]
    old_s2=rows[('current_wide_MI_M5','designscore','s2')]
    new_s2=rows[(f'multilevel_MI_M{m}','designscore','s2')]
    old_original=rows[('current_wide_MI_M5','original','var_subid')]
    new_original=rows[(f'multilevel_MI_M{m}','original','var_subid')]
    lines=['# Multilevel MI versus current MI','',
           f"With multilevel MI (M={m}), design-score subject-variance bias changes from {old_sub['empirical_bias']:+.2f} to {new_sub['empirical_bias']:+.2f}, and residual-variance bias from {old_s2['empirical_bias']:+.2f} to {new_s2['empirical_bias']:+.2f}. Original scenario 3 subject-variance bias changes from {old_original['empirical_bias']:+.2f} to {new_original['empirical_bias']:+.2f}.",'',
           f'Matched Oct2026 00pct samples {a.first}-{a.last}. Both GFR and BMI were masked using miss_ind_mar. Each pair of analysis methods used identical imputations within its imputation method.','',
           f'The multilevel run uses M={m}; existing wide chained MI uses M=5 and was not rerun. The prespecified M=5 sensitivity uses multilevel imputation indices {selected}. This subset changes pooling only; it does not refit models.','',
           'The new imputer jointly models the six GFR/BMI measurements with shared PSU and household effects and an unrestricted subject-trajectory covariance. It also adds continuous baseline age, full stratum indicators, and the specified age/stratum interactions. This comparison evaluates that complete imputation specification; it does not isolate the contribution of random effects alone. The analysis estimators and population targets are unchanged.','',
           '| Analysis | Component | Target | Current MI mean (bias) | Multilevel MI mean (bias) | Coverage: current / multilevel |',
           '|---|---|---:|---:|---:|---:|']
    for method in base.METHODS:
        for term in ('var_hhid','var_subid','s2'):
            o=rows[('current_wide_MI_M5',method,term)]; n=rows[(f'multilevel_MI_M{m}',method,term)]
            lines.append(f"| {method} | {term} | {n['true_value']:.2f} | {o['mean_estimate']:.2f} ({o['empirical_bias']:+.2f}) | {n['mean_estimate']:.2f} ({n['empirical_bias']:+.2f}) | {o['coverage_rate']:.0%} / {n['coverage_rate']:.0%} |")
    sensitivity=[rows[('multilevel_MI_M5_subset','designscore',term)] for term in ('var_hhid','var_subid','s2')]
    lines += ['', 'For the prespecified five-imputation subset, design-score biases are '+
              ', '.join(f"{r['term']}: {r['empirical_bias']:+.2f} (coverage {r['coverage_rate']:.0%})" for r in sensitivity)+
              '. The CSV comparison includes both methods and all parameters.']
    if reference:
        lines += ['', '## Residual covariance check', '',
                  'These unweighted residual moments use the same saved complete-data coefficients and source covariates. They help check covariance preservation but are not fitted variance components. Each sample receives equal weight.', '',
                  '| Moment | Complete data | Current MI | Multilevel MI |', '|---|---:|---:|---:|']
        for term,title in [('household_crossproduct','Household cross-person product'),
                           ('same_subject_crossvisit_covariance','Same-subject cross-visit covariance'),
                           ('within_person_mean_square','Within-person mean square')]:
            selected_rows=[r for r in covariance if r['moment']==term]
            lines.append(f"| {title} | {st.mean(r['complete'] for r in selected_rows):.2f} | {st.mean(r['current_mi'] for r in selected_rows):.2f} | {st.mean(r['imputed'] for r in selected_rows):.2f} |")
    lines += ['', '## Checks and interpretation','',
              f"All requested samples have complete paired results and passed the recorded input/chain checks. Largest R-hat: {max(r['max_rhat'] for r in chain_summary):.4f}; smallest bulk ESS: {min(r['min_bulk_ess'] for r in chain_summary):.1f}; smallest tail ESS: {min(r['min_tail_ess'] for r in chain_summary):.1f}.",'',
              'Chain lengths and source runs are recorded per replicate in `chain_diagnostics_summary.csv`. Any longer run retains the same model, priors, imputation count, and convergence thresholds; failed candidate imputations are excluded.','',
              'Coverage uses the same Rubin formulas for both analyses, with log-scale variance pooling and finite-imputation t intervals. The five-imputation sensitivity and paired Monte Carlo standard errors are in the accompanying CSV files. With 100 replicates, coverage differences retain appreciable Monte Carlo uncertainty.','',
              'Covariance-preservation diagnostics are unweighted descriptive calculations using fixed, previously saved complete-data coefficients. They are not additional fitted variance components and were not used to tune imputations to hidden simulation outcomes. Passing chain diagnostics does not establish that the imputation model is correctly specified.','',
              '## Files','',
              '- `comparison_with_current_mi.csv`: all parameters, both methods, M=20 and M=5 comparisons.',
              '- `paired_changes_from_current_mi.csv`: paired changes and their Monte Carlo standard errors.',
              '- `pooled_multilevel_m5_sensitivity.csv`: replicate-level five-imputation sensitivity.',
              '- `chain_diagnostics_summary.csv`: sampler checks for all requested replicates.',
              '- `covariance_preservation.csv`: complete versus imputed residual moments.','']
    (dest/'COMPARISON_WITH_CURRENT_MI.md').write_text('\n'.join(lines),encoding='utf-8')
    provenance=dict(created_utc=datetime.now(timezone.utc).isoformat(),new_label=a.new_label,current_label=a.current_label,
                    first=a.first,last=a.last,selected_five_draws=selected,selection_rule='first draw from each of four chains, then second draw of chain 1; fixed before result inspection',
                    source_hashes={str(p.relative_to(root)):base.digest(p) for p in [oldpath,newpath,Path(__file__),root/'provenance'/f'{a.new_label}_inputs.json',*([reference_path] if reference else [])]})
    (dest/'comparison_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    print(dest/'COMPARISON_WITH_CURRENT_MI.md')


if __name__=='__main__':
    main()
