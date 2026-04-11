# run illustration.R first, then run this

# ============================================================
# Convert fitted VARMA-type model to LaTeX equations
# in steady-state form:
# z_{i,t} = x_{i,t} - mu_i
# with epsilon_{i,t} placed at the end of the RHS
# Supports VAR / VMA / VARMA
# fit must be a list with components mu, ar, ma
# ============================================================

num_to_tex <- function(x, digits = 4){
  s <- format(round(x, digits), nsmall = digits, trim = TRUE, scientific = FALSE)
  s <- sub("^-0\\.?0*$", "0", s)
  s
}

signed_term <- function(coef, tex_var, digits = 4){
  if(abs(coef) < 10^(-digits)) return("")
  
  coef_abs <- num_to_tex(abs(coef), digits)
  
  if(coef < 0){
    paste0(" - ", coef_abs, " ", tex_var)
  } else {
    paste0(" + ", coef_abs, " ", tex_var)
  }
}

steady_varma_to_latex <- function(
    fit,
    x_name = "x",
    z_name = "z",
    innov_name = "\\varepsilon",
    digits = 4,
    include_equation_env = TRUE
){
  if(is.null(fit$mu)){
    stop("fit must contain component 'mu'.")
  }
  
  n <- length(fit$mu)
  
  pp <- if(!is.null(fit$ar)) ncol(fit$ar) / n else 0
  qq <- if(!is.null(fit$ma)) ncol(fit$ma) / n else 0
  
  eqs <- character(n)
  
  for(i in seq_len(n)){
    mu_i <- num_to_tex(fit$mu[i], digits)
    
    lhs <- paste0(
      z_name, "_{", i, ", t} = ",
      x_name, "_{", i, ", t} - ", mu_i
    )
    
    rhs <- ""
    
    # --------------------------------------------------------
    # AR part in terms of z_{t-lag}
    # --------------------------------------------------------
    if(pp > 0){
      for(lag in seq_len(pp)){
        block <- fit$ar[, ((lag - 1) * n + 1):(lag * n), drop = FALSE]
        for(j in seq_len(n)){
          tex_var <- paste0(z_name, "_{", j, ", t-", lag, "}")
          rhs <- paste0(rhs, signed_term(block[i, j], tex_var, digits))
        }
      }
    }
    
    # --------------------------------------------------------
    # MA part in terms of epsilon_{t-lag}
    # --------------------------------------------------------
    if(qq > 0){
      for(lag in seq_len(qq)){
        block <- fit$ma[, ((lag - 1) * n + 1):(lag * n), drop = FALSE]
        for(j in seq_len(n)){
          tex_var <- paste0(innov_name, "_{", j, ", t-", lag, "}")
          rhs <- paste0(rhs, signed_term(block[i, j], tex_var, digits))
        }
      }
    }
    
    # --------------------------------------------------------
    # innovation at time t always at the end
    # --------------------------------------------------------
    rhs <- paste0(rhs, " + ", innov_name, "_{", i, ", t}")
    
    # remove leading " + " if RHS starts with a positive term
    rhs <- sub("^ \\+ ", "", rhs)
    
    eqs[i] <- paste0(lhs, " = ", rhs)
  }
  
  body <- paste(eqs, collapse = " \\\\\n")
  
  if(include_equation_env){
    out <- paste0(
      "\\begin{equation}\n",
      "\\begin{aligned}\n",
      body, "\n",
      "\\end{aligned}\n",
      "\\end{equation}\n"
    )
  } else {
    out <- body
  }
  
  out
}

cat("% Full-sample estimate\n")
cat(steady_varma_to_latex(fit_fullsample, digits = 3))
cat("\n\n")

cat("% Huber-skip estimate without patch removal\n")
cat(steady_varma_to_latex(fit_huber_k0, digits = 3))
cat("\n\n")

cat("% Huber-skip estimate with patch removal (iterative, kappa = 2)\n")
cat(steady_varma_to_latex(fit_iter_k2$model, digits = 3))
cat("\n\n")

cat("% Huber-skip estimate with patch removal (iterative, kappa = 5)\n")
cat(steady_varma_to_latex(fit_iter_k5$model, digits = 3))
cat("\n")