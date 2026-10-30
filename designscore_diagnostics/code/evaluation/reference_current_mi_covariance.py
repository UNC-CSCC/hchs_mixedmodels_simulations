"""Read-only covariance reference for existing MI; never imputes or fits."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import statistics as st
import numpy as np
import pandas as pd


ROOT=Path(__file__).resolve().parents[2]


def digest(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def reference_moments(data,raw,b,m):
    ids=['bgid','hhid','subid']
    if data.duplicated(ids).any() or raw.duplicated(ids+['v_num']).any():
        raise ValueError('Duplicate input identifiers')
    index=pd.MultiIndex.from_frame(data[ids])
    wide=raw.pivot(index=ids,columns='v_num',values=['y_gfr','y_bmi','miss_ind_mar','x6']).reindex(index)
    if wide.isna().any().any() or len(raw)!=3*len(data):raise ValueError('Source subject mismatch')
    rounding={}
    for v in (1,2,3):
        flag=wide['miss_ind_mar'][v].to_numpy().astype(bool)
        if not np.array_equal(data[f'mar_{v}'].to_numpy(),flag):raise ValueError('Changed MAR flag')
        for outcome in ('y_gfr','y_bmi'):
            actual=data[f'{outcome}_{v}'].to_numpy();original=wide[outcome][v].to_numpy()
            if not np.array_equal(np.isnan(actual),flag):raise ValueError('Wrong missingness mask')
            # Legacy Stata imports stored outcomes as float32. Verify the exact
            # expected rounding, rather than relaxing the multilevel gate.
            if actual.dtype==np.float32:
                matches=np.array_equal(actual[~flag],original[~flag].astype(np.float32))
            else:
                matches=np.allclose(actual[~flag],original[~flag],rtol=0,atol=1e-10)
            if not matches:raise ValueError('Observed value differs beyond its original storage precision')
            rounding[f'{outcome}_{v}']=float(np.max(np.abs(actual[~flag]-original[~flag])))
            if v>1:
                for imp in range(1,m+1):
                    y=data[f'_{imp}_{outcome}_{v}'].to_numpy()
                    if not np.isfinite(y).all() or not np.array_equal(y[~flag],actual[~flag]):
                        raise ValueError('Legacy observed values changed across imputations')
    # Use the same unrounded source covariates for both imputation methods.
    base=raw[raw.v_num==1].set_index(ids).reindex(index)
    mu=np.full(len(data),b['(Intercept)'])+(base.age_base.to_numpy()>=45)*b['age_strat_new']
    for term in ('x12','x17','x18'):mu+=base[term].to_numpy()*b[term]
    mu=mu[:,None]+wide['x6'].to_numpy()*b['x6']
    hh=data.gl_hhid.to_numpy(dtype=int);counts=np.bincount(hh)
    pairs=np.sum(counts*(counts-1))
    def moments(y,bmi):
        r=y-mu-bmi*b['y_bmi'];centered=r-r.mean(axis=0,keepdims=True)
        cross=(centered[:,0]*centered[:,1]+centered[:,0]*centered[:,2]+centered[:,1]*centered[:,2]).mean()/3
        means=r.mean(axis=1);means-=means.mean()
        sums=np.bincount(hh,weights=means);sums2=np.bincount(hh,weights=means**2)
        return dict(household_crossproduct=float(np.sum(sums*sums-sums2)/pairs) if pairs else None,
                    same_subject_crossvisit_covariance=float(cross),
                    within_person_mean_square=float(np.mean(np.sum((r-r.mean(axis=1)[:,None])**2,axis=1)/2)),
                    contrast_v1_v2=float(np.mean((r[:,0]-r[:,1])**2)/2),
                    contrast_v1_v3=float(np.mean((r[:,0]-r[:,2])**2)/2),
                    contrast_v2_v3=float(np.mean((r[:,1]-r[:,2])**2)/2))
    complete=moments(wide['y_gfr'].to_numpy(),wide['y_bmi'].to_numpy())
    imputed=[]
    for imp in range(1,m+1):
        y=data[['y_gfr_1',f'_{imp}_y_gfr_2',f'_{imp}_y_gfr_3']].to_numpy()
        bmi=data[['y_bmi_1',f'_{imp}_y_bmi_2',f'_{imp}_y_bmi_3']].to_numpy()
        imputed.append(moments(y,bmi))
    return dict(complete=complete,imputed=imputed,observed_rounding_max=rounding)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--label',default='mi_s1_v2_m5')
    ap.add_argument('--first',type=int,default=1);ap.add_argument('--last',type=int,default=100)
    ap.add_argument('--imputations',type=int,default=5)
    a=ap.parse_args()
    checker_path=ROOT/'code/evaluation/check_multilevel_mi.py'
    records=[];hashes={}
    for sim in range(a.first,a.last+1):
        source=ROOT/'data/estimates/oct2026'/a.label/'imputations'/f'mi_{sim}.dta'
        before=digest(source)
        coef=ROOT/'data/estimates/oct2026/exact'/f'scen3_designscore_coef_{sim}.csv'
        data=pd.read_stata(source,convert_categoricals=False)
        raw=pd.read_csv(ROOT.parents[1]/'data/Oct2026/sample_miss00pct'/f'samplemiss_{sim}.csv')
        raw.columns=raw.columns.str.lower()
        c=reference_moments(data,raw,pd.read_csv(coef).set_index('term').estimate,a.imputations)
        if digest(source)!=before:raise ValueError('Existing MI integrity changed')
        records.append(dict(sim_id=sim,complete=c['complete'],
                            imputed_mean={k:st.mean(row[k] for row in c['imputed']) for k,v in c['complete'].items() if v is not None},
                            observed_rounding_max=c['observed_rounding_max']))
        hashes[source.relative_to(ROOT).as_posix()]=before
        hashes[coef.relative_to(ROOT).as_posix()]=digest(coef)
        if sim%10==0:print(f'Read-only covariance reference: {sim}/{a.last}',flush=True)
    output=ROOT/'data/run_artifacts'/f'{a.label}_covariance_reference_{a.first}_{a.last}.json'
    result=dict(label=a.label,first=a.first,last=a.last,imputations=a.imputations,records=records,
                source_sha256=hashes,code_sha256={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),checker_path]})
    output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(output)


if __name__=='__main__':main()
