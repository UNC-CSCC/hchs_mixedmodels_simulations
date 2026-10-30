version 19
set more off
set linesize 255
* Definition loading and deterministic algebra ONLY. No data, model, sampler,
* imputations, random-number generation, or analysis command is executed here.
args diagnostic_root
cd `"`diagnostic_root'"'
log using "logs/multilevel_mi_code_check.log", text replace
do code/stata/mi_rename_compat.do
do code/stata/multilevel_mi.mata
mata:
MLMI_K=MLMI_C=J(0,0,.)
mlmi_conditional((4,2\2,9),(1),(2),MLMI_K,MLMI_C)
assert(abs(MLMI_K[1,1]-.5)<1e-12)
assert(abs(MLMI_C[1,1]-8)<1e-12)
assert(abs(mlmi_group_cov(2,.25,1/9)-18/11)<1e-12)
assert(abs(6*.25*mlmi_group_cov(2,.25,1/9)-27/11)<1e-12)
MLMI_IW=mlmi_iw_transform((4,2\2,5),(2,0\1,3))
assert(max(abs(MLMI_IW-(10/9,1/3\1/3,.5)))<1e-12)
assert(abs(mlmi_logdet((4,2\2,9))-ln(32))<1e-12)
real scalar mlmi_test_centered(real matrix L, real matrix Z, real matrix E,
                             real colvector group, real matrix R, real scalar prior)
{
    real matrix U, C, residual
    real scalar q, out, i
    q=cols(Z); U=Z*L'; C=L*L'; residual=E-U[group,.]
    out=-.5*trace(quadcross(residual,residual)*invsym(R))
    out=out-.5*(rows(Z)*mlmi_logdet(C)+trace(quadcross(U,U)*invsym(C)))
    out=out-.5*((2*q+3)*mlmi_logdet(C)+prior*trace(invsym(C)))
    out=out+rows(Z)*sum(ln(diagonal(L)))
    for(i=1;i<=q;i++) out=out+(q+1-i)*ln(L[i,i])
    return(out)
}
MLMI_Z=(1,.3\-.2,.7); MLMI_E=(1,2\.4,-.2\.3,.7\1,1\.6,.8)
MLMI_R=(1,.2\.2,2); MLMI_P=(1,2\3,5); MLMI_MAP=(1\1\2\2\2)
MLMI_T=quadcross(MLMI_Z,(2\3),MLMI_Z)
MLMI_D=quadcross(MLMI_Z,panelsum(MLMI_E,MLMI_P))
MLMI_L1=(.5,0\.2,.7); MLMI_L2=(.6,0\-.1,.8)
MLMI_NC_DIFF=mlmi_nc_target(MLMI_L2,MLMI_T,MLMI_D,invsym(MLMI_R),.05)-mlmi_nc_target(MLMI_L1,MLMI_T,MLMI_D,invsym(MLMI_R),.05)
MLMI_CP_DIFF=mlmi_test_centered(MLMI_L2,MLMI_Z,MLMI_E,MLMI_MAP,MLMI_R,.05)-mlmi_test_centered(MLMI_L1,MLMI_Z,MLMI_E,MLMI_MAP,MLMI_R,.05)
assert(abs(MLMI_NC_DIFF-MLMI_CP_DIFF)<1e-10)
real scalar mlmi_test_dense_hr(real matrix H,real matrix R,real matrix E,
                             real matrix panels,real matrix B,real scalar ridge)
{
    real scalar q,i,n,out
    real matrix V,e,eh,x
    q=cols(E); out=0
    for(i=1;i<=rows(panels);i++) {
        n=panels[i,2]-panels[i,1]+1
        eh=E[|panels[i,1],1\panels[i,2],q|]; e=vec(eh')
        V=I(n)#R+J(n,n,1)#H
        out=out-.5*(mlmi_logdet(V)+e'*invsym(V)*e)
    }
    out=out-.5*(rows(B)*mlmi_logdet(R)+ridge*trace(quadcross(B,B)*invsym(R)))
    out=out-.5*((2*q+3)*(mlmi_logdet(H)+mlmi_logdet(R))+.05*trace(invsym(H))+.5*trace(invsym(R)))
    out=out+mlmi_cov_logjac(mlmi_cov_pack(H),q)+mlmi_cov_logjac(mlmi_cov_pack(R),q)
    return(out)
}
MLMI_BT=(.1,.2\-.1,.4)
MLMI_HRT=mlmi_hr_build(MLMI_E,MLMI_P,MLMI_BT,.0001,.05,.5)
MLMI_H1=MLMI_L1*MLMI_L1'; MLMI_H2=MLMI_L2*MLMI_L2'
MLMI_R2=(1.3,.1\.1,1.8)
MLMI_HRX1=(mlmi_cov_pack(MLMI_H1),mlmi_cov_pack(MLMI_R))
MLMI_HRX2=(mlmi_cov_pack(MLMI_H2),mlmi_cov_pack(MLMI_R2))
assert(max(abs(mlmi_cov_unpack(mlmi_cov_pack(MLMI_R),2)-MLMI_R))<1e-12)
MLMI_COLLAPSED_DIFF=mlmi_hr_target(MLMI_HRX2,MLMI_HRT)-mlmi_hr_target(MLMI_HRX1,MLMI_HRT)
MLMI_DENSE_DIFF=mlmi_test_dense_hr(MLMI_H2,MLMI_R2,MLMI_E,MLMI_P,MLMI_BT,.0001)-mlmi_test_dense_hr(MLMI_H1,MLMI_R,MLMI_E,MLMI_P,MLMI_BT,.0001)
assert(abs(MLMI_COLLAPSED_DIFF-MLMI_DENSE_DIFF)<1e-10)
end
* Hand-written storage fixture only: checks the import format used by the
* wrapper. Values below are constants, not fitted/generated imputations.
clear
input byte pid byte fixture_m double y_gfr_1 double y_gfr_2 double y_gfr_3
1 0 100 . 80
2 0 110 95 .
1 1 100 90 80
2 1 110 95 85
1 2 100 91 80
2 2 110 95 86
end
mi import flong, m(fixture_m) id(pid) imputed(y_gfr_2 y_gfr_3) clear
mi convert wide, clear
quietly mi describe
assert r(M)==2
assert _1_y_gfr_2==90 & _2_y_gfr_2==91 if pid==1
assert _1_y_gfr_2==95 & _2_y_gfr_2==95 if pid==2
assert missing(y_gfr_2) if pid==1
assert y_gfr_1==100 if pid==1
display "MLMI_CODE_CHECK_PASS: definitions, deterministic algebra, and fixed storage fixture passed. No models run."
log close
