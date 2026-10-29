# Exact Gaussian ML: shared household intercept, age-specific subject intercepts,
# and optionally age-specific independent visit residuals. Fixed effects profiled.
# Three balanced visits permit exact aggregation by household age composition.
args <- commandArgs(TRUE)
out <- normalizePath(args[1], winslash="/", mustWork=TRUE)
root <- normalizePath(file.path(out,"../../../.."),winslash="/",mustWork=TRUE)
names_beta <- c("(Intercept)","age_strat_new","x12","x17","x18","x6","y_bmi")

decode <- function(theta, mode) {
  v <- exp(theta)
  switch(mode, common=c(v[1],v[2],v[2],v[3],v[3]),
    subject_age=c(v[1],v[2],v[3],v[4],v[4]),
    residual_age=c(v[1],v[2],v[2],v[3],v[4]),
    both_age=v)
}

objective <- function(theta, s, mode, details=FALSE, gradient=FALSE) {
  v <- decode(theta,mode); a<-v[1]; b<-v[2:3]; e<-v[4:5]; d<-b+e/3
  q <- matrix(0,8,8)
  for (g in 1:2) q <- q+s$mean[[g]]/d[g]+s$within[[g]]/e[g]
  ld <- sum(s$ng*(2*log(e)+log(3*d)))
  for (h in s$patterns) {
    factor <- 1+a*(h$n0/d[1]+h$n1/d[2])
    q <- q-a/factor*(h$c00/d[1]^2+h$c11/d[2]^2+(h$c01+t(h$c01))/(d[1]*d[2]))
    ld <- ld+h$households*log(factor)
  }
  # Profile GLS coefficients; retain Gaussian ML rather than REML to match truth.
  beta <- solve(q[1:7,1:7],q[1:7,8])
  rss <- q[8,8]-sum(q[1:7,8]*beta)
  value <- .5*(ld+rss+3*s$n*log(2*pi))
  if(gradient) {
    rr<-c(-beta,1)
    quadratic<-function(mat) sum(rr*(mat%*%rr))
    gm<-vapply(s$mean,quadratic,0.0)
    gw<-vapply(s$within,quadratic,0.0)
    ga<-0;gd<-s$ng/d-gm/d^2
    ge<-2*s$ng/e-gw/e^2
    for(h in s$patterns) {
      nh<-c(h$n0,h$n1);k<-sum(nh/d);f<-1+a*k
      c00<-quadratic(h$c00);c11<-quadratic(h$c11);c01<-quadratic(h$c01)
      bt<-c00/d[1]^2+c11/d[2]^2+2*c01/prod(d)
      ga<-ga+h$households*k/f-bt/f^2
      dk<--nh/d^2
      dbt<--2*c(c00/d[1]^3+c01/(d[1]^2*d[2]),c11/d[2]^3+c01/(d[2]^2*d[1]))
      gd<-gd+h$households*a*dk/f+a^2*dk*bt/f^2-a*dbt/f
    }
    ge<-ge+gd/3
    gv<-switch(mode,common=c(ga,sum(gd),sum(ge)),subject_age=c(ga,gd,sum(ge)),
      residual_age=c(ga,sum(gd),ge),both_age=c(ga,gd,ge))
    return(.5*gv*exp(theta)/s$n)
  }
  if (!details) return(value/s$n)
  beta[1] <- beta[1]+s$center[8]-sum(s$center[1:7]*beta)
  names(beta)<-names_beta
  list(variance=setNames(v,c("household","subject_younger","subject_older","residual_younger","residual_older")),
       beta=beta,negative_loglik=value,objective_per_subject=value/s$n)
}

small_stats <- function(data) {
  n<-nrow(data)/3; idx<-seq(1,nrow(data),3)
  older<-as.integer(data$age_base[idx]>=45)
  z<-cbind(1,as.integer(data$age_base>=45),as.matrix(data[,c("x12","x17","x18","x6","y_bmi","y_gfr")]))
  center<-colMeans(z);center[1]<-0;z<-sweep(z,2,center)
  z1<-z[idx,,drop=FALSE];z2<-z[idx+1,,drop=FALSE];z3<-z[idx+2,,drop=FALSE]
  m<-(z1+z2+z3)/3
  hh<-interaction(data$BGid[idx],data$hhid[idx],drop=TRUE)
  ng<-tabulate(older+1,2);cross<-within<-sums<-list()
  for(g in 0:1) {
    take<-older==g
    cross[[g+1]]<-crossprod(m[take,,drop=FALSE])
    within[[g+1]]<-crossprod((z1-m)[take,,drop=FALSE])+crossprod((z2-m)[take,,drop=FALSE])+crossprod((z3-m)[take,,drop=FALSE])
    sums[[g+1]]<-rowsum(m*take,hh,reorder=FALSE)
  }
  counts<-cbind(rowsum(as.integer(older==0),hh,reorder=FALSE),rowsum(as.integer(older==1),hh,reorder=FALSE))
  pats<-list()
  for(h in seq_len(nrow(counts))) pats[[h]]<-list(n0=counts[h,1],n1=counts[h,2],households=1,
    c00=crossprod(sums[[1]][h,,drop=FALSE]),c11=crossprod(sums[[2]][h,,drop=FALSE]),
    c01=crossprod(sums[[1]][h,,drop=FALSE],sums[[2]][h,,drop=FALSE]))
  list(n=n,ng=ng,center=center,mean=cross,within=within,patterns=pats)
}

validate_dense <- function(data) {
  st<-small_stats(data)
  x<-cbind(1,as.integer(data$age_base>=45),as.matrix(data[,c("x12","x17","x18","x6","y_bmi")]))
  hh<-interaction(data$BGid,data$hhid,drop=TRUE)
  sub<-interaction(hh,data$subid,drop=TRUE)
  age<-data$age_base>=45;y<-data$y_gfr
  for(v in list(c(25,100,180,60,110),c(40,120,120,57,57))) {
    vc<-v[1]*outer(hh,hh,"==")+outer(sub,sub,"==")*ifelse(age,v[3],v[2])+
      diag(ifelse(age,v[5],v[4]))
    ivx<-solve(vc,x);ivy<-solve(vc,y)
    beta<-solve(crossprod(x,ivx),crossprod(x,ivy))
    r<-y-x%*%beta
    direct<-.5*(as.numeric(determinant(vc,logarithm=TRUE)$modulus)+sum(r*solve(vc,r))+length(y)*log(2*pi))
    fast<-objective(log(v),st,"both_age",TRUE)
    stopifnot(abs(direct-fast$negative_loglik)<1e-6,max(abs(beta-fast$beta))<1e-6)
    numerical<-sapply(seq_along(v),function(j) {
      step<-rep(0,length(v));step[j]<-1e-5
      (objective(log(v)+step,st,"both_age")-objective(log(v)-step,st,"both_age"))/(2e-5)
    })
    stopifnot(max(abs(numerical-objective(log(v),st,"both_age",gradient=TRUE)))<1e-7)
  }
  cat("PASS: compressed likelihood and profiled fixed effects match direct dense covariance.\n")
}

fit_model <- function(s, mode, truth) {
  a<-truth[["var_hhid"]];b<-truth[["var_subid"]];e<-truth[["s2"]]
  start<-switch(mode,common=c(a,b,e),subject_age=c(a,b,b,e),residual_age=c(a,b,e,e),both_age=c(a,b,b,e,e))
  fn<-function(theta) objective(theta,s,mode)
  gr<-function(theta) objective(theta,s,mode,gradient=TRUE)
  fit<-nlminb(start=log(start),objective=fn,gradient=gr,lower=rep(log(1e-5),length(start)),
    upper=rep(log(1e4),length(start)),control=list(iter.max=500,eval.max=1200,rel.tol=1e-12,x.tol=1e-9))
  check<-optim(par=log(start)+seq(-.15,.15,length.out=length(start)),fn=fn,gr=gr,
    method="BFGS",control=list(reltol=1e-14,maxit=1000))
  # Relative likelihood stopping can be premature when census differences are tiny.
  # Polish each independently initialized fit using the analytic score and its Jacobian.
  polish<-function(theta) {
    for(it in 1:20) {
      g<-gr(theta)
      if(max(abs(g))<1e-11) break
      hess<-sapply(seq_along(theta),function(j) {
        delta<-rep(0,length(theta));delta[j]<-1e-4
        (gr(theta+delta)-gr(theta-delta))/2e-4
      })
      hess<-(hess+t(hess))/2
      stopifnot(min(eigen(hess,symmetric=TRUE,only.values=TRUE)$values)>0)
      step<-solve(hess,g);accepted<-FALSE
      for(fraction in 2^-(0:20)) {
        trial<-theta-fraction*step
        if(max(abs(gr(trial)))<max(abs(g)) && fn(trial)<=fn(theta)+1e-12) {
          theta<-trial;accepted<-TRUE;break
        }
      }
      if(!accepted) stop("Newton score polishing failed")
    }
    stopifnot(max(abs(gr(theta)))<1e-9)
    theta
  }
  final<-polish(check$par);independent<-polish(fit$par)
  agreement<-max(abs(exp(final)/exp(independent)-1))
  stopifnot(check$convergence==0,abs(fn(final)-fn(independent))<1e-9,agreement<1e-6)
  result<-objective(final,s,mode,TRUE)
  result$convergence<-list(nlminb=fit$message,BFGS=check$convergence,Newton="analytic score converged from both starts")
  result$max_absolute_log_variance_score_per_subject<-max(abs(gr(final)))
  result$second_optimizer_max_relative_difference<-agreement
  result$subject_older_younger_ratio<-result$variance[3]/result$variance[2]
  # A census contrast is descriptive. No superpopulation p-values are reported.
  result
}

all_results<-list()
for(release in c("Dec2024","Oct2026")) {
  cat("\nPOPULATION:",release,"\n")
  s<-dget(file.path(out,paste0("sufficient_",release,".R")))
  validation<-read.delim(file.path(out,paste0("validation_subset_",release,".txt")))
  validate_dense(validation)
  tab<-read.csv(file.path(root,"hchs_mixedmodels_simulations/data/population_coefficients",
                         paste0("population_coef_",tolower(release),"_hhsub.csv")))
  truth<-setNames(tab$estimate,tab$term)
  results<-list()
  for(mode in c("common","subject_age","residual_age","both_age")) {
    ans<-fit_model(s,mode,truth)
    results[[mode]]<-ans
    cat(mode,"\n");print(ans$variance,digits=10)
    cat("Subject older/younger ratio:",ans$subject_older_younger_ratio,"\n")
  }
  baseline<-results$common$variance[c(1,2,4)]
  target<-truth[c("var_hhid","var_subid","s2")]
  error<-max(abs(baseline/target-1))
  beta_error<-max(abs(results$common$beta-truth[names_beta]))
  stopifnot(error<5e-5,beta_error<1e-3)
  cat("PASS: original census fit reproduced; maximum variance relative error =",error,
      "; maximum fixed-effect absolute error =",beta_error,"\n")
  results$validation<-list(subjects=s$n,households=s$nh,age_group_counts=s$ng,
    baseline_max_relative_error=error,baseline_beta_max_absolute_error=beta_error)
  all_results[[release]]<-results
  saveRDS(all_results,file.path(out,"age_variance_results.rds"))
  dput(all_results,file=file.path(out,"age_variance_results.R"))
}
cat("\nAll diagnostics completed. No Stata process was launched or modified.\n")
