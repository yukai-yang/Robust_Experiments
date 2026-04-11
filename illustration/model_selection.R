# cannot be run directly
# should be sourced by other script

# source dependencies
source("utilfuncs.R")
source("subsampler.R")

# ------------------------------------------------------------
# grid of (p, q)
# ------------------------------------------------------------
pq_grid = matrix(
  c(
    1, 0,
    2, 0,
    3, 0,
    1, 1,
    2, 1,
    3, 1,
    1, 2,
    2, 2,
    3, 2,
    1, 3,
    2, 3,
    3, 3
  ),
  ncol = 2,
  byrow = TRUE
)
colnames(pq_grid) = c("p", "q")

# trimming and patch parameter
alpha = 0.1
kappa = 5

# ------------------------------------------------------------
# load data
# ------------------------------------------------------------
source("getrivers.R")
# data: time-indexed series
# mx  : plain matrix used for estimation

plot(data)

# ------------------------------------------------------------
# helper: count retained observations after patch removal
# ------------------------------------------------------------
count_keep_after_kappa <- function(subsamp, kappa){
  iT = length(subsamp)
  nsubsamp = subsamp
  tmp = which(subsamp == 0)
  
  if(length(tmp) > 0){
    tmp = c(sapply(0:kappa, function(vx) tmp + vx))
    nsubsamp[tmp[tmp <= iT]] = 0
  }
  
  sum(nsubsamp == 1)
}

# ------------------------------------------------------------
# helper: Gaussian average AIC from fitted model and subset
# ------------------------------------------------------------
compute_aic_avg <- function(mx, fit, keep_idx = NULL){
  me = res_varmam(
    model = list(
      ar = fit$ar,
      ma = fit$ma
    ),
    mx = mx,
    mu = fit$mu
  )
  
  if(is.null(keep_idx)){
    e_use = me
  } else {
    e_use = me[keep_idx, , drop = FALSE]
  }
  
  T_eff = nrow(e_use)
  d     = ncol(e_use)
  
  Sigma_hat = crossprod(e_use) / T_eff
  det_Sigma = determinant(Sigma_hat, logarithm = TRUE)
  
  loglik = -0.5 * T_eff * (
    d * log(2 * pi) +
      as.numeric(det_Sigma$modulus) +
      d
  )
  
  k_par = length(fit$mu)
  if(!is.null(fit$ar)) k_par = k_par + length(fit$ar)
  if(!is.null(fit$ma)) k_par = k_par + length(fit$ma)
  
  AIC_avg = -2 * loglik / T_eff + 2 * k_par / T_eff
  
  list(
    AIC_avg = AIC_avg,
    loglik = loglik,
    T_eff = T_eff,
    k_par = k_par,
    Sigma_hat = Sigma_hat
  )
}

# ------------------------------------------------------------
# helper: build retained subset after patch removal
# ------------------------------------------------------------
build_keep_idx <- function(subsamp, kappa){
  keep_idx = subsamp
  tmp = which(subsamp == 0)
  
  if(length(tmp) > 0){
    tmp = c(sapply(0:kappa, function(vx) tmp + vx))
    keep_idx[tmp[tmp <= length(subsamp)]] = 0
  }
  
  keep_idx == 1
}

# ------------------------------------------------------------
# store results
# ------------------------------------------------------------
results = vector("list", nrow(pq_grid))

# ------------------------------------------------------------
# loop over (p, q)
# ------------------------------------------------------------
for(i in seq_len(nrow(pq_grid))){
  
  p = pq_grid[i, 1]
  q = pq_grid[i, 2]
  
  cat("Running (p, q) =", p, q, "\n")
  
  # ----------------------------------------------------------
  # initial model
  # ----------------------------------------------------------
  model_init = list(
    ar = if(p > 0) matrix(0, ncol(mx), ncol(mx) * p) else NULL,
    ma = if(q > 0) matrix(0, ncol(mx), ncol(mx) * q) else NULL
  )
  
  # ----------------------------------------------------------
  # Huber-skip screening
  # ----------------------------------------------------------
  idx_huber = huber_skip_varma_idx(
    mx = mx,
    p = p,
    q = q,
    alpha = alpha,
    c = 1.345,
    model_init = model_init,
    max_iter = 20
  )
  
  # ----------------------------------------------------------
  # one-shot robust estimate with patch removal
  # ----------------------------------------------------------
  fit_patch = est_varmam(
    mx = mx,
    model = model_init,
    subsamp = idx_huber,
    kappa = kappa
  )
  
  keep_idx = build_keep_idx(idx_huber, kappa)
  aic_patch = compute_aic_avg(mx = mx, fit = fit_patch, keep_idx = keep_idx)
  
  # ----------------------------------------------------------
  # full-sample estimate without trimming or patch removal
  # ----------------------------------------------------------
  fit_full = est_varmam(
    mx = mx,
    model = model_init,
    subsamp = rep(1, nrow(mx)),
    kappa = 0
  )
  
  aic_full = compute_aic_avg(mx = mx, fit = fit_full, keep_idx = NULL)
  
  # ----------------------------------------------------------
  # sample sizes
  # ----------------------------------------------------------
  n_total      = nrow(mx)
  n_keep_huber = sum(idx_huber == 1)
  n_keep_kapp  = count_keep_after_kappa(idx_huber, kappa)
  
  # ----------------------------------------------------------
  # store result
  # ----------------------------------------------------------
  results[[i]] = data.frame(
    p = p,
    q = q,
    AIC_avg_patch = aic_patch$AIC_avg,
    AIC_avg_full  = aic_full$AIC_avg,
    sample_size   = n_total,
    huber_keep    = n_keep_huber,
    patch_keep    = n_keep_kapp
  )
}

# ------------------------------------------------------------
# combine all results
# ------------------------------------------------------------
results_df = do.call(rbind, results)

print(results_df)