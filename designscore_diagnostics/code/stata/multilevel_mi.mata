/*
  Joint Gaussian MI for six-wide trajectories, Stata/Mata 19.
  Loading this file defines functions ONLY. No sampler or model is run.

  Y_i = X_i B + U_household[i] + V_PSU[i] + E_i.
  U ~ N_6(0,H), V ~ N_6(0,G), E ~ N_6(0,R), mutually independent.
  All three covariance matrices are unrestricted. R contains the remaining
  between-person and within-person trajectory covariance, NOT just s2.
  Standardized-data priors: C ~ IW(q+2, prior_C * I(q));
  B | R ~ MN(0, I(p)/ridge, R). See MULTILEVEL_MI_README.md.
*/
version 19
mata:

struct mlmi_state {
    real matrix Y, observed, raw, X, B, U, V, H, G, R, P, LP
    real matrix hp, gp
    real colvector hi, gi, pattern, rowid
    real rowvector center, scale, prior
    real scalar n, p, q, ridge, iteration
    real rowvector nc_accept
    real matrix hr_proposal, hr_initial, hr_history
    real rowvector hr_accept
    real scalar warmup, hr_count
    string rowvector ynames, xnames
}

real matrix mlmi_sym(real matrix A) return((A+A')/2)

real matrix mlmi_chol(real matrix A)
{
    real matrix L
    L=cholesky(mlmi_sym(A))
    if (hasmissing(L)) _error(3498,"MLMI covariance is not positive definite")
    return(L)
}

real matrix mlmi_inv(real matrix A)
{
    real matrix L
    L=mlmi_chol(A)
    return(mlmi_sym(invsym(A)))
}

// Deterministic IW transform, separated for algebra-only tests.
// T is a LOWER Bartlett factor: T*T' ~ Wishart(df,I).
real matrix mlmi_iw_transform(real matrix S, real matrix T)
{
    real matrix L, F
    L=mlmi_chol(S)
    F=L*luinv(T')
    return(mlmi_sym(F*F'))
}

real matrix mlmi_iw(real matrix S, real scalar df)
{
    real matrix T
    real scalar i, j, q
    q=rows(S)
    if (df<=q-1) _error(3498,"Invalid inverse-Wishart degrees of freedom")
    T=J(q,q,0)
    for (i=1; i<=q; i++) {
        T[i,i]=sqrt(rchi2(1,1,df-i+1))
        for (j=1; j<i; j++) T[i,j]=rnormal(1,1,0,1)
    }
    return(mlmi_iw_transform(S,T))
}

// Conditional missing-column regression and covariance. No random draws here.
void mlmi_conditional(real matrix R, real rowvector obs, real rowvector mis,
                      real matrix K, real matrix C)
{
    if (!cols(mis)) _error(3498,"Empty missing block")
    if (cols(obs)) {
        K=mlmi_inv(R[obs,obs])*R[obs,mis]
        C=mlmi_sym(R[mis,mis]-R[mis,obs]*K)
    }
    else {
        K=J(0,cols(mis),0)
        C=R[mis,mis]
    }
}

// Group-level posterior covariance shared by groups of the same size.
real matrix mlmi_group_cov(real scalar count, real matrix Ri, real matrix Ci)
    return(mlmi_inv(count*Ri+Ci))

real matrix mlmi_group_draw(real matrix residual, real matrix panels,
                            real matrix R, real matrix C)
{
    real matrix sums, out, Ri, Ci, Q, L
    real colvector counts, sizes, ix
    real scalar j, q
    q=cols(residual); sums=panelsum(residual,panels)
    counts=panels[.,2]-panels[.,1]:+1
    sizes=uniqrows(sort(counts,1)); out=J(rows(panels),q,.)
    Ri=mlmi_inv(R); Ci=mlmi_inv(C)
    for (j=1; j<=rows(sizes); j++) {
        ix=selectindex(counts:==sizes[j])
        Q=mlmi_group_cov(sizes[j],Ri,Ci); L=mlmi_chol(Q)
        out[ix,.]=sums[ix,.]*Ri*Q+rnormal(rows(ix),q,0,1)*L'
    }
    return(out)
}

// Noncentered interweaving leaves the joint model/priors unchanged. Given
// Z=U*inv(L'), propose raw lower-Cholesky entries with a symmetric Gaussian
// random walk. The covariance is based only on Z and R (fixed in this block).
// The target includes the C=L*L' Jacobian, and no U-density determinant:
// the U-to-Z Jacobian cancels it exactly.
real scalar mlmi_nc_target(real matrix L, real matrix T, real matrix D,
                           real matrix Ri, real scalar prior)
{
    real matrix C
    real scalar q, i, out
    q=rows(L)
    if (min(diagonal(L))<=0 | hasmissing(L)) return(-1e300)
    C=L*L'
    if (hasmissing(cholesky(C))) return(-1e300)
    out=-.5*trace(Ri*(L*T*L'-2*L*D))
    out=out-.5*((2*q+3)*2*sum(ln(diagonal(L)))+prior*trace(mlmi_inv(C)))
    for (i=1; i<=q; i++) out=out+(q+1-i)*ln(L[i,i])
    return(out)
}

void mlmi_interweave(real matrix U, real matrix C, real matrix residual,
                     real matrix panels, real matrix R, real scalar prior,
                     real scalar accepted)
{
    real matrix L, Z, T, D, Ri, Q, proposal, candidate
    real colvector counts, pick, delta
    real scalar q, i, j, k, step, oldtarget, newtarget
    q=cols(U); L=mlmi_chol(C); Z=U*luinv(L')
    counts=panels[.,2]-panels[.,1]:+1
    T=quadcross(Z,counts,Z); D=quadcross(Z,panelsum(residual,panels))
    Ri=mlmi_inv(R); Q=T#Ri
    pick=J(q*(q+1)/2,1,.); k=0
    for (i=1; i<=q; i++) for (j=1; j<=i; j++) pick[++k]=i+(j-1)*q
    proposal=(2.38/sqrt(rows(pick)))*mlmi_chol(mlmi_inv(Q[pick,pick]))
    oldtarget=mlmi_nc_target(L,T,D,Ri,prior); accepted=0
    for (step=1; step<=10; step++) {
        delta=proposal*rnormal(rows(pick),1,0,1)
        candidate=L; k=0
        for (i=1; i<=q; i++) for (j=1; j<=i; j++) candidate[i,j]=L[i,j]+delta[++k]
        newtarget=mlmi_nc_target(candidate,T,D,Ri,prior)
        if (newtarget<. & ln(runiform(1,1))<newtarget-oldtarget) {
            L=candidate; oldtarget=newtarget; accepted=accepted+1
        }
    }
    C=mlmi_sym(L*L'); U=Z*L'
}

end
do code/stata/multilevel_mi_collapsed.mata
mata:

struct mlmi_state scalar mlmi_initialize(string scalar yvars,
        string scalar xvars, string scalar hhvar, string scalar psuvar,
        real scalar ridge, real rowvector prior)
{
    struct mlmi_state scalar S
    real matrix raw, candidate, orth, z, ids
    real colvector ix
    real rowvector keep, proposal_sd
    real scalar j, i, k, d, sd, norm
    raw=st_data(.,tokens(yvars)); S.ynames=tokens(yvars)
    S.n=rows(raw); S.q=cols(raw); S.ridge=ridge; S.prior=prior
    if (S.q!=6 | S.n<10 | ridge<=0 | min(prior)<=0 | hasmissing(prior))
        _error(3498,"Invalid MLMI dimensions or priors")
    S.center=S.scale=J(1,S.q,.)
    S.raw=raw; S.observed=raw; S.Y=raw; S.pattern=J(S.n,1,0)
    S.rowid=st_data(.,"gl_subid")
    for (j=1; j<=S.q; j++) {
        ix=selectindex(raw[.,j]:<.)
        if (rows(ix)<2) _error(3498,"Outcome has fewer than two observations")
        S.center[j]=mean(raw[ix,j]); S.scale[j]=sqrt(variance(raw[ix,j]))
        if (S.scale[j]<=0 | S.scale[j]>=.) _error(3498,"Constant outcome")
        S.observed[.,j]=(raw[.,j]:-S.center[j])/S.scale[j]
        S.Y[.,j]=S.observed[.,j]
        ix=selectindex(raw[.,j]:>=.)
        if (rows(ix)) S.Y[ix,j]=rnormal(rows(ix),1,0,1)
        S.pattern=S.pattern+(raw[.,j]:>=.)*2^(j-1)
    }
    candidate=st_data(.,tokens(xvars))
    if (hasmissing(candidate)) _error(3498,"MLMI currently requires complete X")
    // Standardize X, then remove constants and exact linear dependencies.
    // B has a proper ridge prior, but redundant columns still obscure mixing.
    S.X=J(S.n,1,1); orth=S.X/sqrt(S.n); keep=J(1,0,.)
    for (j=1; j<=cols(candidate); j++) {
        sd=sqrt(variance(candidate[.,j]))
        if (sd<=1e-12) continue
        z=(candidate[.,j]:-mean(candidate[.,j]))/sd
        raw=z-orth*(orth'*z)
        raw=raw-orth*(orth'*raw)
        norm=sqrt(sum(raw:^2))
        if (norm<1e-8*sqrt(S.n)) continue
        orth=orth,raw/norm; S.X=S.X,z; keep=keep,j
    }
    S.xnames=("intercept",tokens(xvars)[keep]); S.p=cols(S.X)
    S.P=mlmi_inv(quadcross(S.X,S.X)+ridge*I(S.p)); S.LP=mlmi_chol(S.P)
    ids=st_data(.,(hhvar,psuvar))
    if (hasmissing(ids)) _error(3498,"Missing cluster identifier")
    S.hp=panelsetup(ids[.,1],1); S.gp=panelsetup(ids[.,2],1)
    if (rows(uniqrows(sort(ids[.,1],1)))!=rows(S.hp) |
        rows(uniqrows(sort(ids[.,2],1)))!=rows(S.gp))
        _error(3498,"Sort subjects by PSU, household, subject before MI")
    S.hi=S.gi=J(S.n,1,.)
    for (j=1; j<=rows(S.hp); j++) {
        ix=(S.hp[j,1]..S.hp[j,2])'
        if (min(ids[ix,2])!=max(ids[ix,2])) _error(3498,"Household crosses PSUs")
        S.hi[ix]=J(rows(ix),1,j)
    }
    for (j=1; j<=rows(S.gp); j++) {
        ix=(S.gp[j,1]..S.gp[j,2])'
        S.gi[ix]=J(rows(ix),1,j)
    }
    // Dispersed starts, without fitting complete or unmasked data.
    S.R=exp(rnormal(1,1,0,.4))*I(S.q)
    S.H=.1*exp(rnormal(1,1,0,.7))*I(S.q)
    S.G=.05*exp(rnormal(1,1,0,.7))*I(S.q)
    S.B=rnormal(S.p,S.q,0,.1)
    S.U=rnormal(rows(S.hp),S.q,0,1)*mlmi_chol(S.H)'
    S.V=rnormal(rows(S.gp),S.q,0,1)*mlmi_chol(S.G)'
    S.iteration=0
    S.nc_accept=J(1,4,0)
    d=S.q*(S.q+1)/2
    proposal_sd=J(1,2*d,.); k=0
    for(i=1;i<=S.q;i++) for(j=1;j<=i;j++) {
        k=k+1
        proposal_sd[k]=(i==j ? .10 : .010)
        proposal_sd[d+k]=(i==j ? .020 : .004)
    }
    S.hr_initial=diag(proposal_sd:^2)
    S.hr_proposal=(2.38/sqrt(2*d))*mlmi_chol(S.hr_initial)
    S.hr_history=J(500,2*d,.)
    S.hr_count=0; S.warmup=0; S.hr_accept=J(1,2,0)
    st_local("mlmi_kept_x",invtokens(S.xnames))
    st_matrix("MLMI_scaling",S.center\S.scale)
    return(S)
}

void mlmi_step(struct mlmi_state scalar S)
{
    struct mlmi_hr_stats scalar HR
    real matrix LR, mu, E, K, C, condmean
    real colvector ix
    real rowvector mis, obs, code
    real scalar pat, j, df, accepted
    K=C=J(0,0,.)
    LR=mlmi_chol(S.R)
    S.B=S.P*quadcross(S.X,S.Y-S.U[S.hi,.]-S.V[S.gi,.]) +
        S.LP*rnormal(S.p,S.q,0,1)*LR'
    mu=S.X*S.B
    S.U=mlmi_group_draw(S.Y-mu-S.V[S.gi,.],S.hp,S.R,S.H)
    S.V=mlmi_group_draw(S.Y-mu-S.U[S.hi,.],S.gp,S.R,S.G)
    df=S.q+2
    S.H=mlmi_iw(S.prior[1]*I(S.q)+quadcross(S.U,S.U),df+rows(S.U))
    S.G=mlmi_iw(S.prior[2]*I(S.q)+quadcross(S.V,S.V),df+rows(S.V))
    accepted=0
    mlmi_interweave(S.U,S.H,S.Y-mu-S.V[S.gi,.],S.hp,S.R,S.prior[1],accepted)
    S.nc_accept[1]=S.nc_accept[1]+accepted; S.nc_accept[2]=S.nc_accept[2]+10
    mlmi_interweave(S.V,S.G,S.Y-mu-S.U[S.hi,.],S.gp,S.R,S.prior[2],accepted)
    S.nc_accept[3]=S.nc_accept[3]+accepted; S.nc_accept[4]=S.nc_accept[4]+10
    mu=mu+S.U[S.hi,.]+S.V[S.gi,.]
    E=S.Y-mu
    // p and ridge*B'B come from the matrix-normal B|R prior.
    S.R=mlmi_iw(S.prior[3]*I(S.q)+quadcross(E,E)+
        S.ridge*quadcross(S.B,S.B),df+S.n+S.p)
    // Marginal H/R kernel integrates U out. Immediately redraw U before any
    // other conditional, restoring the correct joint augmented distribution.
    E=S.Y-S.X*S.B-S.V[S.gi,.]
    HR=mlmi_hr_build(E,S.hp,S.B,S.ridge,S.prior[1],S.prior[3])
    mlmi_hr_transition(S.H,S.R,HR,S.hr_proposal,accepted)
    S.hr_accept[1]=S.hr_accept[1]+accepted; S.hr_accept[2]=S.hr_accept[2]+20
    S.U=mlmi_group_draw(E,S.hp,S.R,S.H)
    mu=S.X*S.B+S.U[S.hi,.]+S.V[S.gi,.]
    // Finite warmup adaptation only. Use a recent window so dispersed starting
    // values do not determine the proposal retained for posterior sampling.
    if(S.iteration<S.warmup) {
        S.hr_count=S.hr_count+1
        S.hr_history[mod(S.hr_count-1,500)+1,.]=(mlmi_cov_pack(S.H),mlmi_cov_pack(S.R))
        if(S.hr_count>=100 & mod(S.hr_count,25)==0) {
            C=.9*variance(S.hr_history[1..min((S.hr_count,500)),.])+.1*S.hr_initial
            S.hr_proposal=(2.38/sqrt(cols(C)))*mlmi_chol(C)
        }
    }
    for (pat=1; pat<2^S.q; pat++) {
        ix=selectindex(S.pattern:==pat)
        if (!rows(ix)) continue
        code=J(1,S.q,0)
        for (j=1; j<=S.q; j++) code[j]=mod(floor(pat/2^(j-1)),2)
        mis=selectindex(code); obs=selectindex(!code)
        mlmi_conditional(S.R,obs,mis,K,C)
        condmean=mu[ix,mis]
        if (cols(obs)) condmean=condmean+(S.observed[ix,obs]-mu[ix,obs])*K
        S.Y[ix,mis]=condmean+rnormal(rows(ix),cols(mis),0,1)*mlmi_chol(C)'
    }
    if (hasmissing(S.Y)) _error(3498,"Incomplete MLMI draw")
    S.iteration=S.iteration+1
}

real scalar mlmi_logdet(real matrix A)
    return(2*sum(ln(diagonal(mlmi_chol(A)))))

real scalar mlmi_logjoint(struct mlmi_state scalar S)
{
    real matrix E, Ri, Hi, Gi
    real scalar out, df
    Ri=mlmi_inv(S.R); Hi=mlmi_inv(S.H); Gi=mlmi_inv(S.G); df=S.q+2
    E=S.Y-S.X*S.B-S.U[S.hi,.]-S.V[S.gi,.]
    out=-.5*((S.n+S.p)*mlmi_logdet(S.R)+
        trace((quadcross(E,E)+S.ridge*quadcross(S.B,S.B))*Ri))
    out=out-.5*(rows(S.U)*mlmi_logdet(S.H)+trace(quadcross(S.U,S.U)*Hi))
    out=out-.5*(rows(S.V)*mlmi_logdet(S.G)+trace(quadcross(S.V,S.V)*Gi))
    out=out-.5*((df+S.q+1)*(mlmi_logdet(S.R)+mlmi_logdet(S.H)+mlmi_logdet(S.G))+
        S.prior[1]*trace(Hi)+S.prior[2]*trace(Gi)+S.prior[3]*trace(Ri))
    return(out)
}

string scalar mlmi_trace_header(struct mlmi_state scalar S)
{
    string scalar out, name
    real scalar k, i, j
    out="chain,iteration,warmup,log_joint"
    for (k=1; k<=3; k++) {
        name=("hh","psu","person")[k]
        for (i=1; i<=S.q; i++) for (j=1; j<=i; j++)
            out=out+","+name+"_"+strofreal(i)+"_"+strofreal(j)
    }
    for (i=1; i<=S.p; i++) for (j=1; j<=S.q; j++)
        out=out+",beta_"+S.xnames[i]+"_"+strofreal(j)
    return(out)
}

real rowvector mlmi_trace_values(struct mlmi_state scalar S)
{
    real rowvector out
    real matrix A
    real scalar k, i, j
    out=mlmi_logjoint(S)
    for (k=1; k<=3; k++) {
        if (k==1) A=S.H
        else if (k==2) A=S.G
        else A=S.R
        for (i=1; i<=S.q; i++) for (j=1; j<=i; j++) out=out,A[i,j]
    }
    for (i=1; i<=S.p; i++) for (j=1; j<=S.q; j++) out=out,S.B[i,j]
    return(out)
}

void mlmi_advance(struct mlmi_state scalar S, real scalar steps,
        real scalar chain, real scalar burn, real scalar every, real scalar fh)
{
    real scalar t
    real rowvector values
    S.warmup=burn
    for (t=1; t<=steps; t++) {
        mlmi_step(S)
        if (mod(S.iteration,every)==0) {
            values=(chain,S.iteration,S.iteration<=burn,mlmi_trace_values(S))
            fput(fh,invtokens(strofreal(values,"%24.16g"),","))
        }
    }
}

// Writes ONLY originally missing cells. Observed values remain bit-for-bit fixed.
void mlmi_store(struct mlmi_state scalar S)
{
    real scalar j
    real colvector ix
    if (st_nobs()!=S.n) _error(3498,"Subject order changed during MI")
    if (any(st_data(.,"gl_subid"):!=S.rowid)) _error(3498,"Subject order changed during MI")
    if (any((st_data(.,S.ynames):!=S.raw):&(S.raw:<.)))
        _error(3498,"Observed input changed during MI")
    for (j=1; j<=S.q; j++) {
        ix=selectindex(S.observed[.,j]:>=.)
        if (rows(ix)) st_store(ix,S.ynames[j],S.Y[ix,j]*S.scale[j]:+S.center[j])
    }
}
end
