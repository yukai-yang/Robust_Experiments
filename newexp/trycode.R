source("utilfuncs.R")

# example 1

iT = 100

ar = matrix( 0.5, 1, 1)
ma = matrix( 0.5, 1, 1)

# VARMA model
ar = matrix( c(0.5, 0.5, 0, 0.5), 2, 2)
ma = matrix( c(0.5, 0.5, 0, 0.5), 2, 2)

# VAR
ar = matrix( c(0.5, 0.5, 0, 0.5), 2, 2)
ma = matrix( 0, 2, 2)

# MA
ar = matrix( 0, 2, 2)
ma = matrix( c(0.5, 0, 0, 0.5), 2, 2)

# model -->
model = list(ar=ar,ma=ma)
n = nrow(model$ar)

# simulation
set.seed(1)

innov = matrix(rnorm(iT * n), iT, n)
zeta = matrix(0, iT, n)
zeta[iT%/%2,] = 10
zeta[1,] = 10;

tmp = sim_varma(model, iT, n, innov, zeta)
matplot(tmp[[3]], type='l')

#Log = capture.output({ ret = VARMA(tmp[[1]], p=1, q=1, include.mean=F, details=F) })
#ret$Phi; ret$Theta

# residuals
res0 = res_varma(model=model, mx=matrix(tmp$mx0,iT,n))
res1 = res_varma(model=model, mx=matrix(tmp$mx1,iT,n))
res2 = res_varma(model=model, mx=matrix(tmp$mx2,iT,n))
model
#estimation

subsamp = rep(1, iT)
subsamp[c(1,iT%/%2)] = 0

est0 = est_varma(tmp$mx0, model, subsamp, kappa=0);est0
est1 = est_varma(tmp$mx1, model, subsamp, kappa=0);est1
est2 = est_var(tmp$mx2, model, subsamp, kappa=0);est2

est1 = est_varma(tmp$mx1, model, subsamp, kappa=0);est1
est1 = est_varma(tmp$mx1, model, subsamp, kappa=2);est1
est1 = est_varma(tmp$mx1, model, subsamp, kappa=4);est1
