/* Additional exact-target transition for the same joint Gaussian MI model.
   Integrate household effects out, jointly update H/R, then the caller MUST
   immediately redraw household effects. Definitions only; no fitting on load.
   Load after multilevel_mi.mata. */
version 19
mata:
struct mlmi_hr_stats {
    real matrix within, quadratics, beta_penalty
    real colvector size, count
    real scalar n, nh, p, q, hprior, rprior
}

real rowvector mlmi_cov_pack(real matrix C)
{
    real matrix L
    real rowvector x
    real scalar i,j,k,q
    q=rows(C); L=mlmi_chol(C); x=J(1,q*(q+1)/2,.); k=0
    for(i=1;i<=q;i++) for(j=1;j<=i;j++) x[++k]=(i==j ? ln(L[i,j]) : L[i,j])
    return(x)
}

real matrix mlmi_cov_unpack(real rowvector x, real scalar q)
{
    real matrix L
    real scalar i,j,k
    L=J(q,q,0); k=0
    for(i=1;i<=q;i++) for(j=1;j<=i;j++) {
        k=k+1
        L[i,j]=(i==j ? exp(x[k]) : x[k])
    }
    return(mlmi_sym(L*L'))
}

real scalar mlmi_cov_logjac(real rowvector x, real scalar q)
{
    real scalar i,j,k,out
    k=0; out=0
    for(i=1;i<=q;i++) for(j=1;j<=i;j++) {
        k=k+1
        if(i==j) out=out+(q+2-i)*x[k]
    }
    return(out)
}

struct mlmi_hr_stats scalar mlmi_hr_build(real matrix E, real matrix panels,
              real matrix B, real scalar ridge, real scalar hprior, real scalar rprior)
{
    struct mlmi_hr_stats scalar T
    real matrix sums, Q
    real colvector ns, ix
    real scalar j, q
    q=cols(E); T.q=q; T.n=rows(E); T.nh=rows(panels); T.p=rows(B)
    T.hprior=hprior; T.rprior=rprior; T.beta_penalty=ridge*quadcross(B,B)
    sums=panelsum(E,panels); ns=panels[.,2]-panels[.,1]:+1
    Q=quadcross(sums,1:/ns,sums)
    T.within=mlmi_sym(quadcross(E,E)-Q)
    T.size=uniqrows(sort(ns,1)); T.count=J(rows(T.size),1,.)
    T.quadratics=J(rows(T.size)*q,q,.)
    for(j=1;j<=rows(T.size);j++) {
        ix=selectindex(ns:==T.size[j]); T.count[j]=rows(ix)
        T.quadratics[|(j-1)*q+1,1\j*q,q|]=quadcross(sums[ix,.],1:/ns[ix],sums[ix,.])
    }
    return(T)
}

real scalar mlmi_hr_target(real rowvector x, struct mlmi_hr_stats scalar T)
{
    real matrix H,R,A,HL,RL,AL,Q
    real scalar q,d,j,out
    q=T.q; d=q*(q+1)/2
    H=mlmi_cov_unpack(x[1..d],q); R=mlmi_cov_unpack(x[(d+1)..(2*d)],q)
    HL=cholesky(H); RL=cholesky(R)
    if(hasmissing(HL) | hasmissing(RL)) return(-1e300)
    out=-.5*((T.n-T.nh+T.p+2*q+3)*2*sum(ln(diagonal(RL)))+
        trace(invsym(R)*(T.within+T.beta_penalty+T.rprior*I(q)))+
        (2*q+3)*2*sum(ln(diagonal(HL)))+T.hprior*trace(invsym(H)))
    for(j=1;j<=rows(T.size);j++) {
        A=R+T.size[j]*H; AL=cholesky(A)
        if(hasmissing(AL)) return(-1e300)
        Q=T.quadratics[|(j-1)*q+1,1\j*q,q|]
        out=out-.5*(T.count[j]*2*sum(ln(diagonal(AL)))+trace(invsym(A)*Q))
    }
    out=out+mlmi_cov_logjac(x[1..d],q)+mlmi_cov_logjac(x[(d+1)..(2*d)],q)
    return(out)
}

// P is a fixed proposal factor within this block. Adaptation (if any) belongs
// to warmup only, outside this kernel. The raw theta proposal is symmetric.
void mlmi_hr_transition(real matrix H,real matrix R,struct mlmi_hr_stats scalar T,
                       real matrix P,real scalar accepted)
{
    real rowvector x,y
    real scalar a,b,step,d
    x=(mlmi_cov_pack(H),mlmi_cov_pack(R)); d=cols(x)/2
    a=mlmi_hr_target(x,T); accepted=0
    for(step=1;step<=20;step++) {
        y=x+rnormal(1,cols(x),0,1)*P'
        b=mlmi_hr_target(y,T)
        if(b<. & ln(runiform(1,1))<b-a) {
            x=y; a=b; accepted=accepted+1
        }
    }
    H=mlmi_cov_unpack(x[1..d],T.q); R=mlmi_cov_unpack(x[(d+1)..(2*d)],T.q)
}
end
