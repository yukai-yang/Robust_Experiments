#library(quantmod)
source("utilfuncs.R")
source("subsampler.R")

# ------------------------------------------------------------
# choose data window
# ------------------------------------------------------------
#from = "2015-01-01"; to = "2025-12-31"

# ------------------------------------------------------------
# choose VARMA orders and trimming level
# ------------------------------------------------------------
p = 2
q = 0
alpha = 0.1

source("getrivers.R")

varnames = colnames(data)

# quick plot
plot(data)

# optional initial model
model_init = list(
  ar = if(p > 0) matrix(0, ncol(mx), ncol(mx) * p) else NULL,
  ma = if(q > 0) matrix(0, ncol(mx), ncol(mx) * q) else NULL
)

# ------------------------------------------------------------
# run Huber-skip screening
# ------------------------------------------------------------
idx_huber = huber_skip_varma_idx(
  mx = mx,
  p = p,
  q = q,
  alpha = alpha,
  c = 1.345,
  model_init = model_init,
  max_iter = 20
)

# check number of discarded observations
sum(idx_huber == 0)

# dates flagged by Huber-skip
dates_huber = index(data)[idx_huber == 0]

# ------------------------------------------------------------
# inspect flagged dates
# ------------------------------------------------------------
head(dates_huber, 20)

# ------------------------------------------------------------
# build xts indicators for plotting
# ------------------------------------------------------------
flag_huber = xts(1 - idx_huber, order.by = index(data))
colnames(flag_huber) = "Huber_flag"

# ------------------------------------------------------------
# simple plots of returns with flagged dates
# ------------------------------------------------------------
tt    = index(data)
y_sp  = as.numeric(data[, varnames[1]])
y_nas = as.numeric(data[, varnames[2]])

flag_huber = which(idx_huber == 0)

par(mfrow = c(2, 1))

plot(tt, y_sp, type = "l", main = "S&P 500 log returns (Huber flags)",
     xlab = "", ylab = varnames[1])
abline(v = tt[flag_huber], lwd = 1)
points(tt[flag_huber], y_sp[flag_huber], pch = 19, cex = 0.5)

plot(tt, y_nas, type = "l", main = "NASDAQ log returns (Huber flags)",
     xlab = "", ylab = varnames[2])
abline(v = tt[flag_huber], lwd = 1)
points(tt[flag_huber], y_nas[flag_huber], pch = 19, cex = 0.5)

par(mfrow = c(1, 1))

# ------------------------------------------------------------
# optional: estimate patch-based VARMA
# ------------------------------------------------------------

fit_huber_k0 = est_varmam(mx = mx, model = model_init, subsamp = idx_huber, kappa = 0)
fit_huber_k1 = est_varmam(mx = mx, model = model_init, subsamp = idx_huber, kappa = 1)
fit_huber_k2 = est_varmam(mx = mx, model = model_init, subsamp = idx_huber, kappa = 2)
fit_huber_k5 = est_varmam(mx = mx, model = model_init, subsamp = idx_huber, kappa = 5)

fit_iter_k0 = iter_huber_patch_varma(mx = mx, H_init = idx_huber, p = p, q = q, alpha = alpha, kappa = 0)
fit_iter_k1 = iter_huber_patch_varma(mx = mx, H_init = idx_huber, p = p, q = q, alpha = alpha, kappa = 1)
fit_iter_k2 = iter_huber_patch_varma(mx = mx, H_init = idx_huber, p = p, q = q, alpha = alpha, kappa = 2)
fit_iter_k5 = iter_huber_patch_varma(mx = mx, H_init = idx_huber, p = p, q = q, alpha = alpha, kappa = 5)

# ------------------------------------------------------------
# optional: compare retained sample sizes after patch removal
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

sapply(c(0, 1, 2, 5), function(k) count_keep_after_kappa(idx_huber, k))

# ============================================================
# USER INPUT FOR FORECASTING EXPERIMENT
# ============================================================
i0 = 1
iN = 500
kappa = 2

n_workers = 6   # choose how many CPUs to use

source("forecast_experiment.R")

# ============================================================
# OUTPUT
# ============================================================
iter = 15
idx_huber[iter]
actual[iter,]
fc_full[iter,]
fc_robust_k0[iter,]
fc_robust_kapp[iter,]

sqrt(mse_full)
sqrt(mse_robust_k0)
sqrt(mse_robust_kapp)
