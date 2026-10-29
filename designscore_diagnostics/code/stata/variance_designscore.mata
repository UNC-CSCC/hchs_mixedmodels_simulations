/*
   Scenario 3 design-weighted census variance scores, Stata/Mata 19.
   Load with: do code/stata/variance_designscore.mata
   Loading defines functions only; it never fits a model.

   beta order: intercept age_strat_new x12 x17 x18 x6 y_bmi.
   z order: y followed by those seven design columns.
   The caller supplies complete, sorted three-visit subjects.
*/
version 19
mata:

struct ds_data {
    // meta: PSU, stratum, full HH size, nominal W_H, exact W_H,
    //       sampled HH size, sum(raw W_sub), HT estimate of HH size
    real matrix meta, z1, B1, C1, B2, Q2
    real scalar scale
}

struct ds_solution {
    real rowvector variance, score
    real scalar objective, mode, kkt, attempts, accepted
}

struct ds_data scalar ds_build(real matrix z, real matrix design,
                               real scalar py, real scalar po,
                               real scalar pyy, real scalar poo,
                               real scalar pyo)
{
    struct ds_data scalar D
    real matrix panels, zs, zm, zz, b1, c1, b2, q2, dev
    real rowvector sw, pr, age, rr, ix
    real scalar h, i, j, first, last, nh, wh, pij, p, block

    // design: gl_hhid gl_psu strat m_h w_bghh p_h w_sub age_strat_new
    panels = panelsetup(design[.,1], 1)
    nh = rows(panels)
    D.meta = J(nh,8,0); D.z1 = J(nh,8,0)
    D.B1 = J(8*nh,8,0); D.C1 = J(8*nh,8,0)
    D.B2 = J(8*nh,8,0); D.Q2 = J(8*nh,8,0)
    for (h=1; h<=nh; h++) {
        first=panels[h,1]; last=panels[h,2]
        wh=(last-first+1)/3
        if (wh!=floor(wh)) _error(3498, "Incomplete subject in household")
        zm=J(wh,8,0); sw=pr=age=J(1,wh,0)
        b1=c1=b2=q2=J(8,8,0)
        for (i=1; i<=wh; i++) {
            j=first+3*(i-1)
            zs=z[|j,1 \ j+2,8|]
            rr=colsum(zs)/3
            zm[i,.]=rr
            sw[i]=design[j,7]; age[i]=design[j,8]
            p=(age[i]==0 ? py : po); pr[i]=p
            zz=rr'*rr
            dev=zs-J(3,1,1)*rr
            b1=b1+sw[i]*zz
            c1=c1+sw[i]*(dev'*dev)
            b2=b2+zz/p
            D.z1[h,.]=D.z1[h,.]+sw[i]*rr
        }
        q2=b2
        for (i=1; i<=wh; i++) for (j=i+1; j<=wh; j++) {
            if (age[i]!=age[j]) pij=pyo
            else pij=(age[i]==0 ? pyy : poo)
            if (pij<=0 | pij>=.) _error(3498, "Invalid observed-pair probability")
            q2=q2+(zm[i,.]'*zm[j,.]+zm[j,.]'*zm[i,.])/pij
        }
        D.meta[h,.]=(design[first,2], design[first,3], design[first,4],
                     design[first,5], 1/design[first,6], wh, sum(sw), sum(1:/pr))
        block=8*(h-1)+1; ix=block..(block+7)
        D.B1[ix,.]=b1; D.C1[ix,.]=c1
        D.B2[ix,.]=b2; D.Q2[ix,.]=q2
    }
    D.scale=sum(D.meta[.,5])
    return(D)
}

real matrix ds_quadratics(struct ds_data scalar D, real rowvector beta)
{
    real colvector q
    real matrix out
    real scalar h, k
    q=(1,-beta)'
    out=J(rows(D.meta),2,0)
    for (h=1; h<=rows(D.meta); h++) {
        k=8*(h-1)+1
        out[h,1]=q'*D.B2[|k,1 \ k+7,8|]*q
        out[h,2]=q'*D.Q2[|k,1 \ k+7,8|]*q
    }
    return(out)
}

real matrix ds_variance_scores(real rowvector v, real scalar e,
                              real matrix meta, real matrix moments)
{
    real scalar a, b, d
    real colvector M, N, W, B, Q, den, ua, ub
    a=v[1]; b=v[2]; d=b+e/3
    M=meta[.,3]; W=meta[.,5]; N=meta[.,8]
    B=moments[.,1]; Q=moments[.,2]; den=d:+a:*M
    ua=.5:*W:*(Q:/den:^2-N:/den)
    ub=.5:*W:*(B:/d^2-a:*(2*d:+a:*M):*Q:/(d^2:*den:^2)
                   -N:/d+a:*N:/(d:*den))
    return((ua,ub))
}

real scalar ds_objective(real rowvector v, real scalar e,
                        real matrix meta, real matrix moments)
{
    real scalar a, b, d
    real colvector M, N, W, den
    a=v[1]; b=v[2]; d=b+e/3
    M=meta[.,3]; W=meta[.,5]; N=meta[.,8]; den=d:+a:*M
    // HT estimate of the census conditional log likelihood, up to constants.
    return(-.5*sum(W:*(N:*ln(d)+(N:/M):*ln(den:/d)
              +moments[.,1]:/d-a:*moments[.,2]:/(d:*den))))
}

real rowvector ds_decode(real rowvector p, real scalar mode)
{
    if (mode==0) return(exp(p))
    if (mode==1) return((0,exp(p[1])))
    return((exp(p[1]),0))
}

void ds_evaluator(real scalar todo, real rowvector p, real scalar e,
                  real matrix meta, real matrix moments, real scalar mode,
                  real scalar scale, real scalar value, real rowvector g,
                  real matrix H)
{
    real rowvector v
    // Numerical guards; zero-variance boundaries are evaluated separately.
    if (max(p)>25 | min(p)<-25 | hasmissing(p)) {
        value=.
        return
    }
    v=ds_decode(p,mode)
    value=ds_objective(v,e,meta,moments)/scale
    if (todo>=1) {
        g=colsum(ds_variance_scores(v,e,meta,moments)):/scale
        if (mode==0) g=g:*v
        else if (mode==1) g=g[2]*v[2]
        else g=g[1]*v[1]
    }
}

real scalar ds_kkt(real rowvector v, real rowvector score,
                   real scalar e, real scalar scale)
{
    real rowvector s
    real scalar j
    // Use raw variance scores, so a tiny log variance cannot fake convergence.
    s=score:*(v:+e/3):/scale
    for (j=1; j<=2; j++) {
        if (v[j]==0) s[j]=max((0,s[j]))
        else s[j]=abs(s[j])
    }
    return(max(s))
}

struct ds_solution scalar ds_solve(struct ds_data scalar D, real rowvector beta,
                                  real scalar e, real rowvector initial)
{
    struct ds_solution scalar ans
    transmorphic S
    real matrix moments, factors
    real rowvector p, v, sc, start
    real scalar mode, j, ec, val, check

    moments=ds_quadratics(D,beta)
    ans.variance=ans.score=J(1,2,.)
    ans.objective=.; ans.mode=.; ans.kkt=.
    ans.attempts=0; ans.accepted=0
    factors=(1,1 \ .25,2 \ 4,.5)
    for (mode=0; mode<=2; mode++) for (j=1; j<=rows(factors); j++) {
        start=initial:*factors[j,.]
        start=(max((start[1],1e-4)),max((start[2],1e-4)))
        if (mode==0) p=ln(start)
        else if (mode==1) p=ln(start[2])
        else p=ln(start[1])
        S=optimize_init()
        optimize_init_evaluator(S,&ds_evaluator())
        optimize_init_evaluatortype(S,"d1")
        optimize_init_which(S,"max")
        optimize_init_technique(S,"bfgs")
        optimize_init_params(S,p)
        optimize_init_argument(S,1,e)
        optimize_init_argument(S,2,D.meta)
        optimize_init_argument(S,3,moments)
        optimize_init_argument(S,4,mode)
        optimize_init_argument(S,5,D.scale)
        optimize_init_tracelevel(S,"none")
        optimize_init_conv_maxiter(S,300)
        optimize_init_conv_ptol(S,1e-10)
        optimize_init_conv_vtol(S,1e-12)
        ec=_optimize(S)
        ans.attempts=ans.attempts+1
        if (ec!=0 | !optimize_result_converged(S)) continue
        p=optimize_result_params(S)
        if (hasmissing(p) | min(p)<=-24.9 | max(p)>=24.9) continue
        v=ds_decode(p,mode)
        sc=colsum(ds_variance_scores(v,e,D.meta,moments))
        val=ds_objective(v,e,D.meta,moments)/D.scale
        check=ds_kkt(v,sc,e,D.scale)
        if (check>=. | check>1e-6 | val>=.) continue
        ans.accepted=ans.accepted+1
        if (ans.objective>=. | val>ans.objective+1e-10 |
            (abs(val-ans.objective)<=1e-10 & mode>ans.mode)) {
            ans.variance=v; ans.score=sc; ans.objective=val
            ans.mode=mode; ans.kkt=check
        }
    }
    // The corner is also a valid constrained candidate.
    v=(0,0); sc=colsum(ds_variance_scores(v,e,D.meta,moments))
    val=ds_objective(v,e,D.meta,moments)/D.scale
    check=ds_kkt(v,sc,e,D.scale)
    if (check<. & check<=1e-6 & val<.) {
        ans.accepted=ans.accepted+1
        if (ans.objective>=. | val>=ans.objective-1e-10) {
            ans.variance=v; ans.score=sc; ans.objective=val
            ans.mode=3; ans.kkt=check
        }
    }
    return(ans)
}

real matrix ds_stacked_scores(real rowvector p, struct ds_data scalar D)
{
    // p: beta(7), log(a1), log(b1), log(e), log(global scale), log(a2), log(b2)
    real scalar h, k, a, b, e, c, d, A, B, C, R, den, ua, ub, ue
    real colvector q, bq, cq
    real rowvector beta, meanbeta, withinbeta, sumbeta
    real matrix out, blockb, blockc, moments
    beta=p[1..7]; q=(1,-beta)'
    a=exp(p[8]); b=exp(p[9]); e=exp(p[10]); c=exp(p[11]); d=b+e/3
    out=J(rows(D.meta),13,0)
    for (h=1; h<=rows(D.meta); h++) {
        k=8*(h-1)+1
        blockb=D.B1[|k,1 \ k+7,8|]; blockc=D.C1[|k,1 \ k+7,8|]
        bq=blockb*q; cq=blockc*q
        A=c*D.meta[h,7]; B=c*(q'*bq)
        C=c*D.z1[h,.]*q; R=c*(q'*cq); den=d+a*A
        ua=.5*(C^2/den^2-A/den)
        ub=.5*(B/d^2-a*(2*d+a*A)*C^2/(d^2*den^2)-A/d+a*A/(d*den))
        ue=ub/3+.5*(R/e^2-2*A/e)
        withinbeta=c*cq[2..8]'/e
        meanbeta=c*bq[2..8]'/d
        sumbeta=a*C*c*D.z1[h,2..8]/(d*den)
        out[h,1..7]=D.meta[h,4]*(withinbeta+meanbeta-sumbeta)
        out[h,8..10]=D.meta[h,4]*(ua,ub,ue)
        // Include the sample-dependent global weight normalization equation.
        out[h,11]=D.meta[h,6]-A
    }
    moments=ds_quadratics(D,beta)
    out[.,12..13]=ds_variance_scores(exp(p[12..13]),e,D.meta,moments)
    return(out)
}

real matrix ds_sandwich(real rowvector p, struct ds_data scalar D,
                       real scalar status, real scalar step1check,
                       real scalar jac_condition, real scalar npsu)
{
    real matrix raw, A, meat, V, invA, plus, minus, psus, psuscore, dev
    real rowvector scale, pp, pm, total, ref
    real colvector ix, strata, jj
    real scalar j, s, delta, n, K

    status=0; step1check=.; jac_condition=.; npsu=.
    raw=ds_stacked_scores(p,D)
    if (hasmissing(raw)) {
        status=1
        return(J(13,13,.))
    }
    // Compare mean score to mean absolute household score; catches wrong
    // extraction/parameterization of the fitted Step 1 before reporting SEs.
    total=colsum(raw); ref=colsum(abs(raw))
    for (j=1; j<=13; j++) if (ref[j]<1) ref[j]=1
    step1check=max(abs(total[1..11]):/ref[1..11])
    if (step1check>1e-3) {
        status=2
        return(J(13,13,.))
    }
    scale=sqrt(colsum(raw:^2)/rows(raw))
    for (j=1; j<=13; j++) if (scale[j]<1e-10) scale[j]=1e-10
    A=J(13,13,0)
    for (j=1; j<=13; j++) {
        delta=1e-5*max((1,abs(p[j])))
        pp=pm=p; pp[j]=p[j]+delta; pm[j]=p[j]-delta
        plus=ds_stacked_scores(pp,D); minus=ds_stacked_scores(pm,D)
        A[.,j]=-(colsum(plus-minus):/scale)'/(2*delta)
    }
    if (hasmissing(A) | rank(A)<13) {
        status=3
        return(J(13,13,.))
    }
    jac_condition=cond(A)
    if (jac_condition>=. | jac_condition>1e12) {
        status=3
        return(J(13,13,.))
    }
    psus=uniqrows(sort(D.meta[.,1..2],(2,1)))
    K=rows(psus); npsu=K
    psuscore=J(K,13,0)
    for (j=1; j<=K; j++) {
        ix=selectindex(D.meta[.,1]:==psus[j,1])
        psuscore[j,.]=colsum(raw[ix,.]):/scale
    }
    // Ultimate-cluster, with-replacement approximation: stratified BG scores.
    // This is not an exact multistage/FPC variance estimator.
    meat=J(13,13,0); strata=uniqrows(sort(psus[.,2],1))
    for (s=1; s<=rows(strata); s++) {
        jj=selectindex(psus[.,2]:==strata[s]); n=rows(jj)
        if (n<2) {
            status=4
            return(J(13,13,.))
        }
        dev=psuscore[jj,.]-J(n,1,1)*(colsum(psuscore[jj,.])/n)
        meat=meat+n/(n-1)*(dev'*dev)
    }
    invA=luinv(A); V=invA*meat*invA'; V=(V+V')/2
    if (hasmissing(V) | min(diagonal(V))<0) {
        status=5
        return(J(13,13,.))
    }
    return(V)
}

void ds_run(string scalar zvars, string scalar designvars,
            string scalar initialmatrix, string scalar probabilities,
            string scalar resultmatrix, string scalar diagnosticmatrix)
{
    struct ds_data scalar D
    struct ds_solution scalar fit
    real rowvector initial, probs, p
    real matrix result, diagnostics, V, z, design
    real scalar e, c, a1, b1, j, idx, v, se, crit, status, check, condition, npsu

    initial=st_matrix(initialmatrix); probs=st_matrix(probabilities)
    if (rows(initial)!=1 | cols(initial)!=11 | rows(probs)!=1 | cols(probs)!=5)
        _error(3498,"Invalid Step 1/probability matrix dimensions")
    if (hasmissing(initial) |
        min(initial[8..11])<=0) _error(3498,"Invalid Step 1 parameters")
    if (hasmissing(probs) | min(probs)<=0 | max(probs)>1)
        _error(3498,"Invalid inclusion probabilities")
    e=initial[10]; c=initial[11]; a1=initial[8]; b1=initial[9]
    z=st_data(.,tokens(zvars)); design=st_data(.,tokens(designvars))
    if (hasmissing(z) | hasmissing(design)) _error(3498,"Missing design-score inputs")
    D=ds_build(z,design,probs[1],probs[2],probs[3],probs[4],probs[5])
    fit=ds_solve(D,initial[1..7],e,(a1,b1))
    if (fit.mode>=.) _error(430,"No admissible design-score solution; no estimates exported")
    result=J(10,4,.)
    result[1..7,1]=initial[1..7]'
    result[8..10,1]=(fit.variance,e)'
    status=6; check=condition=npsu=.
    if (fit.mode==0) {
        p=(initial[1..7],ln(a1),ln(b1),ln(e),ln(c),ln(fit.variance))
        V=ds_sandwich(p,D,status,check,condition,npsu)
        if (status==0) {
            crit=invnormal(.975)
            for (j=1; j<=7; j++) {
                se=sqrt(V[j,j]); result[j,2..4]=(se,p[j]-crit*se,p[j]+crit*se)
            }
            for (j=8; j<=10; j++) {
                // Output order: household, subject, residual variance.
                if (j==8) idx=12
                else if (j==9) idx=13
                else idx=10
                se=sqrt(V[idx,idx]); v=result[j,1]
                result[j,2..4]=(v*se,exp(ln(v)-crit*se),exp(ln(v)+crit*se))
            }
            if (hasmissing(result)) {
                status=7
                result[.,2..4]=J(10,3,.)
            }
        }
    }
    diagnostics=(fit.mode,fit.kkt,fit.objective,fit.attempts,fit.accepted,
                 status,check,condition,npsu,rows(D.meta),sum(D.meta[.,6]),fit.score)
    st_matrix(resultmatrix,result)
    st_matrix(diagnosticmatrix,diagnostics)
}

end
