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
  
  mx = mx[(pp+1):(pp+iT),]
  my1 = mx + zeta
  my2 = my2[(pp+1):(pp+iT),]
  
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
  
  me = me[(qq+1):(qq+iT),]
  return(me)
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