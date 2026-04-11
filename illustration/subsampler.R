# ============================================================
# Helper: stack lagged rows in the same order as the original code
# ============================================================
lag_stack <- function(mat, end_row, nlags){
  if(nlags <= 0) return(numeric(0))
  idx = seq(end_row, end_row - nlags + 1)
  return(c(t(mat[idx, , drop = FALSE])))
}


# ============================================================
# Steady-state residuals for VARMA(p,q)
# x_t - mu = Phi(L)(x_{t-1} - mu) + e_t + Theta(L)e_{t-1}
# model$ar or model$ma may be NULL
# ============================================================
res_varma_mean <- function(model, mx, mu){
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
# Steady-state Huber estimation for VARMA / VAR / VMA
# model$ar or model$ma may be NULL
# ============================================================
est_varma_mean_huber <- function(mx, model, c = 1.345){
  iT = nrow(mx)
  n  = ncol(mx)
  
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
    
    model1 = list(ar = ar, ma = ma)
    me = res_varma_mean(model = model1, mx = mx, mu = mu)
    
    d = sqrt(rowSums(me^2))
    s = median(d) / 0.6745
    if(!is.finite(s) || s <= 1e-8) s = 1e-8
    
    u = d / s
    w = ifelse(u <= c, 1, c / u)
    
    me_w = me * sqrt(w)
    return(det(crossprod(me_w) / sum(w)))
  }
  
  init_mu = apply(mx, 2, median)
  init_par = c(init_mu)
  if(has_ar) init_par = c(init_par, c(model$ar))
  if(has_ma) init_par = c(init_par, c(model$ma))
  
  ret = optim(par = init_par, fn = ftmp)
  
  mu = ret$par[1:n]
  pos = n
  
  ar = NULL
  if(has_ar){
    ar = matrix(ret$par[(pos + 1):(pos + nar)], n, n * pp)
    pos = pos + nar
  }
  
  ma = NULL
  if(has_ma){
    ma = matrix(ret$par[(pos + 1):(pos + nma)], n, n * qq)
    pos = pos + nma
  }
  
  return(list(mu = mu, ar = ar, ma = ma))
}


# ============================================================
# Huber-skip outlier detection for VARMA / VAR / VMA
# model_init$ar or model_init$ma may be NULL
# ============================================================
huber_skip_varma_idx <- function(mx, p, q, alpha, c = 1.345, model_init = NULL, max_iter = 20){
  iT = nrow(mx)
  n  = ncol(mx)
  
  ntrim = floor(alpha * iT)
  
  if(is.null(model_init)){
    model_init = list(
      ar = if(p > 0) matrix(0, n, n * p) else NULL,
      ma = if(q > 0) matrix(0, n, n * q) else NULL
    )
  }
  
  subsamp = rep(1, iT)
  
  for(iter in 1:max_iter){
    fit = est_varma_mean_huber(mx = mx, model = model_init, c = c)
    me  = res_varma_mean(model = list(ar = fit$ar, ma = fit$ma), mx = mx, mu = fit$mu)
    
    d = sqrt(rowSums(me^2))
    s = median(d) / 0.6745
    if(!is.finite(s) || s <= 1e-8) s = 1e-8
    
    u = d / s
    ord = order(u, decreasing = TRUE)
    
    new_subsamp = rep(1, iT)
    if(ntrim > 0){
      new_subsamp[ord[1:ntrim]] = 0
    }
    
    if(all(new_subsamp == subsamp)){
      subsamp = new_subsamp
      break
    }
    
    subsamp = new_subsamp
    model_init = list(ar = fit$ar, ma = fit$ma)
  }
  
  return(subsamp)
}


# ============================================================
# Iterative feasible patch-removal procedure
# Works for VARMA / VAR / VMA
# Assumes est_varmam() and res_varmam() support NULL
# ============================================================
iter_huber_patch_varma <- function(mx, H_init, p, q, alpha, kappa = 0,
                                   c = 1.345, max_iter = 20){
  iT = nrow(mx)
  n  = ncol(mx)
  
  H_old = H_init
  
  model_start = list(
    ar = if(p > 0) matrix(0, n, n * p) else NULL,
    ma = if(q > 0) matrix(0, n, n * q) else NULL
  )
  
  fit_old = est_varmam(mx = mx, model = model_start, subsamp = H_old, kappa = kappa)
  
  for(iter in 1:max_iter){
    me = res_varmam(model = fit_old, mx = mx, mu = fit_old$mu)
    
    d = sqrt(rowSums(me^2))
    s = median(d) / 0.6745
    if(!is.finite(s) || s <= 1e-8) s = 1e-8
    
    u = d / s
    ord = order(u, decreasing = TRUE)
    
    ntrim = floor(alpha * iT)
    H_new = rep(1, iT)
    if(ntrim > 0){
      H_new[ord[1:ntrim]] = 0
    }
    
    fit_new = est_varmam(mx = mx, model = fit_old, subsamp = H_new, kappa = kappa)
    
    if(all(H_new == H_old)){
      return(list(model = fit_new, subsamp = H_new, iter = iter))
    }
    
    H_old = H_new
    fit_old = fit_new
  }
  
  return(list(model = fit_old, subsamp = H_old, iter = max_iter))
}