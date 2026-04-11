# functions used for experiments

# input:
#  model, a list of ar ma, list(ar=matrix(0,n,n), ma=matrix(0,n,n)) for white noise
#   ar and ma are matrics like (phi_1, phi_2, ...) n by np, (theta_1, theta_2, ...) n by nq
#  iT sample size
#  n #variables
#  innov, a iT by n matrix of errors
#  zeta, a iT by n matrix of outliers
# output:
#  mx, no outliers
#  my1, additive outliers
#  mv2, innovational outliers
sim_varma <- function(model, iT, n, innov, zeta){
  pp = ncol(model$ar)/n # AR lag length
  qq = ncol(model$ma)/n # MA lag length
  
  mx = matrix(0, pp+iT, n); my2 = mx
  innov = rbind(matrix(0,qq,n), innov)
  
  for(iter in 1:iT){
    mx[iter+pp,] = model$ar %*% c(t(mx[(iter+pp-1):(iter),])) + innov[iter+qq,] +
      model$ma %*% c(t(innov[(iter+qq-1):(iter),]))
    
    my2[iter+pp,] = model$ar %*% c(t(my2[(iter+pp-1):(iter),])) + innov[iter+qq,] +
      model$ma %*% c(t(innov[(iter+qq-1):(iter),])) + zeta[iter,]
  }
  
  mx = mx[(pp+1):(pp+iT),,drop=FALSE]
  my1 = mx + zeta
  my2 = my2[(pp+1):(pp+iT),,drop=FALSE]
  
  return(list(mx0 = mx, mx1 = my1, mx2 = my2))
}

# get the residuals based on the parameters
# model, a list of ar ma, list(ar=0, ma=0) for white noise
# mx the data iT by n
res_varma <- function(model, mx){
  iT = nrow(mx)
  n = ncol(mx)
  pp = ncol(model$ar)/n # AR lag length
  qq = ncol(model$ma)/n # MA lag length
  
  mx = rbind(matrix(0,pp,n), mx)
  me = matrix(0, qq+iT, n)
  
  for(iter in 1:iT){
    me[iter+qq,] = mx[iter+pp,] -  model$ar %*% c(t(mx[(iter+pp-1):(iter),])) -
      model$ma %*% c(t(me[(iter+qq-1):(iter),]))
  }
  
  me = me[(qq+1):(qq+iT),,drop=FALSE]
  return(me)
}

# iT, 500, 1000
# cont, contamination, no, AO, IO
# outlier, 0, 5, 10, 50, 100
# alpha, 0.01, 0.05, 0.1
# kappa, 0, 1, 2, 4
make_table <- function(atables, iT=1, cont=1:3, outlier=2:5, alpha=1:3,kappa=1:4){
  tables = NULL
  for(cter in cont){
    subtab = NULL
    for(oter in outlier){
      subtab = rbind(subtab,atables[cter,iT,alpha,oter,kappa])
    }
    tables = cbind(tables, subtab)
  }
  return(tables)
}



# estimate the varma model
# mx the data iT by n
# model, model specification and initial value for estimation
# subsamp, the vector of subsample, 0(trimmed) or 1(chosen)
# kappa
est_varma <- function(mx, model,subsamp, kappa=0){
  iT = length(subsamp)
  nsubsamp = subsamp
  tmp = which(subsamp==0)
  if(length(tmp)>0){
    tmp = c(sapply(0:kappa, function(vx) tmp+vx))
    nsubsamp[tmp[tmp<=iT]] = 0
  }
  iTT = sum(nsubsamp)
  nsubsamp = nsubsamp == 1
  
  n = ncol(mx)
  pp = ncol(model$ar)/n # AR lag length
  qq = ncol(model$ma)/n # MA lag length
  
  ftmp <- function(par){
    ar = matrix( par[1:(n*n*pp)], n, n*pp)
    ma = matrix( par[(n*n*pp+1):length(par)], n, n*qq)
    model = list(ar=ar,ma=ma)
    
    me = res_varma(model=model, mx=mx)
    return(det(crossprod(me[nsubsamp,])/iTT))
  }
  
  ret = c(c(model$ar),c(model$ma))
  ret = optim(par=ret, fn=ftmp)
  
  ar = matrix( ret$par[1:(n*n*pp)], n, n*pp)
  ma = matrix( ret$par[(n*n*pp+1):length(ret$par)], n, n*qq)
  
  return(list(ar=ar, ma=ma))
}

# estimate the var model
# mx the data iT by n
# model, model specification and initial value for estimation
# subsamp, the vector of subsample, 0(trimmed) or 1(chosen)
# kappa
est_var <- function(mx, model,subsamp, kappa=0){
  iT = length(subsamp)
  nsubsamp = subsamp
  tmp = which(subsamp==0)
  if(length(tmp)>0){
    tmp = c(sapply(0:kappa, function(vx) tmp+vx))
    nsubsamp[tmp[tmp<=iT]] = 0
  }
  iTT = sum(nsubsamp)
  nsubsamp = nsubsamp == 1
  
  n = ncol(mx)
  pp = ncol(model$ar)/n # AR lag length
  
  ftmp <- function(par){
    ar = matrix(par, n, n*pp)
    ma = matrix(0, n, n)
    model = list(ar=ar,ma=ma)
    
    me = res_varma(model=model, mx=mx)
    return(det(crossprod(me[nsubsamp,])/iTT))
  }
  
  ret = c(model$ar)
  ret = optim(par=ret, fn=ftmp)
  
  ar = matrix(ret$par, n, n*pp)
  
  return(list(ar=ar))
}


# estimate the vma model
# mx the data iT by n
# model, model specification and initial value for estimation
# subsamp, the vector of subsample, 0(trimmed) or 1(chosen)
# kappa
est_vma <- function(mx, model, subsamp, kappa=0){
  iT = length(subsamp)
  nsubsamp = subsamp
  tmp = which(subsamp==0)
  if(length(tmp)>0){
    tmp = c(sapply(0:kappa, function(vx) tmp+vx))
    nsubsamp[tmp[tmp<=iT]] = 0
  }
  iTT = sum(nsubsamp)
  nsubsamp = nsubsamp == 1
  
  n = ncol(mx)
  qq = ncol(model$ma)/n # MA lag length
  
  ftmp <- function(par){
    ma = matrix(par, n, n*qq)
    model = list(ar=matrix(0,n,n), ma=ma)
    
    me = res_varma(model=model, mx=mx)
    return(det(crossprod(me[nsubsamp,])/iTT))
  }
  
  ret = c(model$ma)
  ret = optim(par=ret, fn=ftmp)
  
  ma = matrix(ret$par, n, n*qq)
  
  return(list(ma=ma))
}


# ============================================================
# Helper: stack lagged rows in the same order as the original code
# ============================================================
# need lag_stack 


# ============================================================
# Get steady-state residuals based on the parameters
# model: list(ar = ..., ma = ...), where ar or ma may be NULL
# mx: iT by n data matrix
# mu: n-dimensional mean vector
# ============================================================
res_varmam <- function(model, mx, mu){
  iT = nrow(mx)
  n  = ncol(mx)
  
  has_ar = !is.null(model$ar)
  has_ma = !is.null(model$ma)
  
  pp = if(has_ar) ncol(model$ar) / n else 0
  qq = if(has_ma) ncol(model$ma) / n else 0
  
  mx_c = sweep(mx, 2, mu, FUN = "-")
  mx_c = rbind(matrix(0, pp, n), mx_c)
  me   = matrix(0, qq + iT, n)
  
  for(iter in 1:iT){
    ar_part = if(pp > 0){
      model$ar %*% lag_stack(mx_c, iter + pp - 1, pp)
    } else {
      rep(0, n)
    }
    
    ma_part = if(qq > 0){
      model$ma %*% lag_stack(me, iter + qq - 1, qq)
    } else {
      rep(0, n)
    }
    
    me[iter + qq, ] = mx_c[iter + pp, ] - ar_part - ma_part
  }
  
  me = me[(qq + 1):(qq + iT), , drop = FALSE]
  return(me)
}


# ============================================================
# Estimate steady-state VAR / VMA / VARMA model
# mx: iT by n data matrix
# model: list(ar = ..., ma = ...), where ar or ma may be NULL
# subsamp: vector of subsample indicators, 0(trimmed) or 1(chosen)
# kappa: patch length
#
# The function jointly estimates mu and the dynamic parameters.
# It automatically handles VAR, VMA, and VARMA according to whether
# model$ar and model$ma are NULL.
#
# Returns:
#   list(mu = ..., ar = ..., ma = ...)
# where ar is returned only if model$ar is not NULL,
# and ma is returned only if model$ma is not NULL.
# ============================================================
est_varmam <- function(mx, model, subsamp, kappa = 0){
  iT = length(subsamp)
  nsubsamp = subsamp
  
  tmp = which(subsamp == 0)
  if(length(tmp) > 0){
    tmp = c(sapply(0:kappa, function(vx) tmp + vx))
    nsubsamp[tmp[tmp <= iT]] = 0
  }
  
  iTT = sum(nsubsamp)
  nsubsamp = nsubsamp == 1
  
  n = ncol(mx)
  
  has_ar = !is.null(model$ar)
  has_ma = !is.null(model$ma)
  
  pp = if(has_ar) ncol(model$ar) / n else 0
  qq = if(has_ma) ncol(model$ma) / n else 0
  
  nar = n * n * pp
  nma = n * n * qq
  
  ftmp <- function(par){
    mu = par[1:n]
    
    pos = n
    
    ar = NULL
    if(has_ar){
      ar = matrix(par[(pos + 1):(pos + nar)], n, n * pp)
      pos = pos + nar
    }
    
    ma = NULL
    if(has_ma){
      ma = matrix(par[(pos + 1):(pos + nma)], n, n * qq)
      pos = pos + nma
    }
    
    model_now = list(ar = ar, ma = ma)
    me = res_varmam(model = model_now, mx = mx, mu = mu)
    
    return(det(crossprod(me[nsubsamp, , drop = FALSE]) / iTT))
  }
  
  init_mu = apply(mx, 2, median)
  init_par = c(init_mu)
  
  if(has_ar) init_par = c(init_par, c(model$ar))
  if(has_ma) init_par = c(init_par, c(model$ma))
  
  ret = optim(par = init_par, fn = ftmp)
  
  mu = ret$par[1:n]
  
  pos = n
  out = list(mu = mu)
  
  if(has_ar){
    out$ar = matrix(ret$par[(pos + 1):(pos + nar)], n, n * pp)
    pos = pos + nar
  }
  
  if(has_ma){
    out$ma = matrix(ret$par[(pos + 1):(pos + nma)], n, n * qq)
    pos = pos + nma
  }
  
  return(out)
}
