rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# =========================================================
# 0. source your own functions
# =========================================================
source("utilfuncs.R")

# =========================================================
# 1. basic setup
# =========================================================
iT    <- 16L
nvar  <- 1L
t0    <- 4L
shock <- 5

phi   <- 0.8
theta <- 0.8

set.seed(1)
innov <- matrix(rnorm(iT * nvar), nrow = iT, ncol = nvar)

zeta <- matrix(0, nrow = iT, ncol = nvar)
zeta[t0, 1] <- shock

# =========================================================
# 2. helper: construct models exactly in utilfuncs.R format
# =========================================================
make_model <- function(ar_coef = 0, ma_coef = 0, nvar = 1L) {
  ar_mat <- matrix(ar_coef, nrow = nvar, ncol = nvar)
  ma_mat <- matrix(ma_coef, nrow = nvar, ncol = nvar)
  
  list(ar = ar_mat, ma = ma_mat)
}

model_ar   <- make_model(ar_coef = phi,   ma_coef = 0,     nvar = nvar)
model_ma   <- make_model(ar_coef = 0,     ma_coef = theta, nvar = nvar)
model_arma <- make_model(ar_coef = phi,   ma_coef = theta, nvar = nvar)

# =========================================================
# 3. helper: validate model object for utilfuncs.R
# =========================================================
check_model <- function(model, nvar) {
  if (!is.matrix(model$ar)) {
    stop("model$ar must be a matrix.")
  }
  if (!is.matrix(model$ma)) {
    stop("model$ma must be a matrix.")
  }
  if (nrow(model$ar) != nvar) {
    stop("nrow(model$ar) must equal nvar.")
  }
  if (nrow(model$ma) != nvar) {
    stop("nrow(model$ma) must equal nvar.")
  }
  if ((ncol(model$ar) %% nvar) != 0) {
    stop("ncol(model$ar) must be a multiple of nvar.")
  }
  if ((ncol(model$ma) %% nvar) != 0) {
    stop("ncol(model$ma) must be a multiple of nvar.")
  }
  invisible(TRUE)
}

# =========================================================
# 4. simulate and compute residuals using your own code
# =========================================================
build_case <- function(model, title_text, innov, zeta, iT, nvar) {
  check_model(model, nvar)
  
  sim <- sim_varma(
    model = model,
    iT    = iT,
    n     = nvar,
    innov = innov,
    zeta  = zeta
  )
  
  mx0 <- sim$mx0
  mx1 <- sim$mx1
  mx2 <- sim$mx2
  
  me0 <- res_varma(model, mx0)
  me1 <- res_varma(model, mx1)
  me2 <- res_varma(model, mx2)
  
  list(
    title = title_text,
    level = tibble(
      t     = seq_len(iT),
      clean = as.numeric(mx0[, 1]),
      AO    = as.numeric(mx1[, 1]),
      IO    = as.numeric(mx2[, 1])
    ),
    resid = tibble(
      t     = seq_len(iT),
      clean = as.numeric(me0[, 1]),
      AO    = as.numeric(me1[, 1]),
      IO    = as.numeric(me2[, 1])
    )
  )
}

ar_dat   <- build_case(model_ar,   "AR level",   innov, zeta, iT, nvar)
ma_dat   <- build_case(model_ma,   "MA level",   innov, zeta, iT, nvar)
arma_dat <- build_case(model_arma, "ARMA level", innov, zeta, iT, nvar)

# =========================================================
# 5. plotting helpers
# =========================================================
expand_range <- function(x, pad = 0.08) {
  rx <- range(x, na.rm = TRUE)
  rr <- diff(rx)
  if (rr == 0) rr <- 1
  c(rx[1] - pad * rr, rx[2] + pad * rr)
}

make_two_series_plot <- function(dat_list, contam = c("AO", "IO"), t0 = 5L) {
  contam <- match.arg(contam)
  
  all_levels <- c(
    dat_list$ar$level$clean,   dat_list$ar$level[[contam]],
    dat_list$ma$level$clean,   dat_list$ma$level[[contam]],
    dat_list$arma$level$clean, dat_list$arma$level[[contam]]
  )
  
  all_resids <- c(
    dat_list$ar$resid$clean,   dat_list$ar$resid[[contam]],
    dat_list$ma$resid$clean,   dat_list$ma$resid[[contam]],
    dat_list$arma$resid$clean, dat_list$arma$resid[[contam]]
  )
  
  ylim_level <- expand_range(all_levels, 0.08)
  ylim_resid <- expand_range(all_resids, 0.08)
  
  lab2 <- paste(contam, "contaminated")
  
  col_map   <- c(clean = "black", contam = "grey45")
  lt_map    <- c(clean = "solid", contam = "dashed")
  shape_map <- c(clean = 16, contam = 1)
  
  base_guides <- guides(
    colour   = guide_legend(order = 1, title = NULL, override.aes = list(linewidth = 0.8)),
    linetype = "none",
    shape    = "none"
  )
  
  make_level_panel <- function(df, title_text) {
    plot_df <- tibble(
      t      = df$t,
      clean  = df$clean,
      contam = df[[contam]]
    ) |>
      pivot_longer(-t, names_to = "series", values_to = "value")
    
    ggplot(plot_df, aes(x = t, y = value,
                        colour = series, linetype = series, shape = series)) +
      annotate(
        "rect",
        xmin = t0 - 0.35, xmax = t0 + 0.35,
        ymin = -Inf, ymax = Inf,
        fill = "grey80", alpha = 0.18
      ) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 2.1, stroke = 0.8, fill = "white") +
      scale_colour_manual(
        values = col_map,
        breaks = c("clean", "contam"),
        labels = c("uncontaminated", lab2),
        name = NULL
      ) +
      scale_linetype_manual(
        values = lt_map,
        breaks = c("clean", "contam"),
        labels = c("uncontaminated", lab2),
        name = NULL
      ) +
      scale_shape_manual(
        values = shape_map,
        breaks = c("clean", "contam"),
        labels = c("uncontaminated", lab2),
        name = NULL
      ) +
      scale_x_continuous(breaks = c(4, 8, 12)) +
      coord_cartesian(ylim = ylim_level) +
      labs(title = title_text, x = NULL, y = NULL) +
      base_guides +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )
  }
  
  make_resid_panel <- function(df) {
    plot_df <- tibble(
      t      = df$t,
      clean  = df$clean,
      contam = df[[contam]]
    ) |>
      pivot_longer(-t, names_to = "series", values_to = "value")
    
    ggplot(plot_df, aes(x = t, y = value,
                        colour = series, linetype = series, shape = series)) +
      annotate(
        "rect",
        xmin = t0 - 0.35, xmax = t0 + 0.35,
        ymin = -Inf, ymax = Inf,
        fill = "grey80", alpha = 0.18
      ) +
      geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey45") +
      geom_segment(
        aes(xend = t, y = 0, yend = value),
        linewidth = 0.65,
        position = position_dodge(width = 0.18)
      ) +
      geom_point(
        size = 2.1, stroke = 0.8, fill = "white",
        position = position_dodge(width = 0.18)
      ) +
      scale_colour_manual(
        values = col_map,
        breaks = c("clean", "contam"),
        labels = c("uncontaminated", lab2),
        name = NULL
      ) +
      scale_linetype_manual(
        values = lt_map,
        breaks = c("clean", "contam"),
        labels = c("uncontaminated", lab2),
        name = NULL
      ) +
      scale_shape_manual(
        values = shape_map,
        breaks = c("clean", "contam"),
        labels = c("uncontaminated", lab2),
        name = NULL
      ) +
      scale_x_continuous(breaks = c(4, 8, 12)) +
      coord_cartesian(ylim = ylim_resid) +
      labs(title = "residuals", x = NULL, y = NULL) +
      base_guides +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )
  }
  
  p11 <- make_level_panel(dat_list$ar$level,   dat_list$ar$title)
  p12 <- make_resid_panel(dat_list$ar$resid)
  
  p21 <- make_level_panel(dat_list$ma$level,   dat_list$ma$title)
  p22 <- make_resid_panel(dat_list$ma$resid)
  
  p31 <- make_level_panel(dat_list$arma$level, dat_list$arma$title)
  p32 <- make_resid_panel(dat_list$arma$resid)
  
  ((p11 | p12) /
      (p21 | p22) /
      (p31 | p32)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# =========================================================
# 6. build AO and IO figures
# =========================================================
dat_list <- list(
  ar   = ar_dat,
  ma   = ma_dat,
  arma = arma_dat
)

fig_ao <- make_two_series_plot(dat_list, contam = "AO", t0 = t0)
fig_io <- make_two_series_plot(dat_list, contam = "IO", t0 = t0)

# =========================================================
# 7. print and save
# =========================================================
print(fig_ao)
print(fig_io)

ggsave(
  filename = "outlier_effect_AO.png",
  plot = fig_ao,
  width = 10,
  height = 10,
  dpi = 320
)

ggsave(
  filename = "outlier_effect_IO.png",
  plot = fig_io,
  width = 10,
  height = 10,
  dpi = 320
)