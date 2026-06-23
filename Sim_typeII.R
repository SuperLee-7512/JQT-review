# Type II right-censoring Weibull simulation

rm(list = ls())

# =============================================================================
# Function definitions
# =============================================================================

loglik_type2 <- function(beta, a, x, n) {
  if (beta <= 0) return(-Inf)
  r <- length(x)
  delta <- x - a
  z <- beta * delta
  r_i <- exp(z)
  ## r failures + n-r censored at x_r
  ell <- r * log(beta) + beta * sum(delta) - sum(r_i) - (n - r) * r_i[r]
  ell
}

## 2. Profile in beta: a_hat(beta)

a_hat_beta_type2 <- function(beta, x, n) {
  r <- length(x)
  bx <- beta * x
  m <- max(bx)
  w <- exp(bx - m)
  ## sum_{i=1}^{r-1} exp(beta x_i) + (n - r + 1) exp(beta x_r)
  if (r > 1) {
    sum_s <- sum(w[1:(r - 1)]) + (n - r + 1) * w[r]
  } else {
    sum_s <- (n) * w[1]
  }
  sum_s <- exp(m) * sum_s
  (log(sum_s) - log(r)) / beta
}

profile_ll_beta_type2 <- function(beta, x, n) {
  if (beta <= 0) return(-Inf)
  a_hat <- a_hat_beta_type2(beta, x, n)
  loglik_type2(beta, a_hat, x, n)
}

mle_beta_type2 <- function(x, n) {
  f <- function(logb) {
    beta <- exp(logb)
    -profile_ll_beta_type2(beta, x, n)
  }
  opt <- optimize(f, interval = c(log(1e-3), log(1e3)))
  beta_hat <- exp(opt$minimum)
  ll_max   <- -opt$objective
  a_hat    <- a_hat_beta_type2(beta_hat, x, n)
  list(beta_hat = beta_hat, a_hat = a_hat, ll_max = ll_max)
}

## 3. Canonical parameter phi(beta, a) under Type II

phi_type2 <- function(beta, a, x, n) {
  r <- length(x)
  delta <- x - a
  z <- beta * delta
  r_i <- exp(z)
  
  if (r > 1) {
    S0_fail <- sum(1 - r_i[1:(r - 1)])
    S1_fail <- sum((1 - r_i[1:(r - 1)]) * z[1:(r - 1)])
  } else {
    S0_fail <- 0
    S1_fail <- 0
  }
  
  d_last <- 1 - (n - r + 1) * r_i[r]
  
  phi1 <- beta * (S0_fail + d_last)
  phi2 <- beta * (S1_fail + d_last * z[r])
  
  c(phi1, phi2)
}

## numerical Jacobian of phi with respect to (beta, a)

jacobian_phi_type2 <- function(beta, a, x, n, h = 1e-5) {
  theta <- c(beta, a)
  fun <- function(th) {
    phi_type2(th[1], th[2], x, n)
  }
  k <- 2
  p <- 2
  J <- matrix(NA_real_, nrow = k, ncol = p)
  for (j in 1:p) {
    step <- h * max(1, abs(theta[j]))
    th_plus  <- theta
    th_minus <- theta
    th_plus[j]  <- theta[j] + step
    th_minus[j] <- theta[j] - step
    if (j == 1 && th_minus[1] <= 0) {
      f0 <- fun(theta)
      fp <- fun(th_plus)
      J[, j] <- (fp - f0) / step
    } else {
      fp <- fun(th_plus)
      fm <- fun(th_minus)
      J[, j] <- (fp - fm) / (2 * step)
    }
  }
  J
}

## 4. Observed information J(beta, a) for Type II

J_type2 <- function(beta, a, x, n) {
  if (beta <= 0) return(matrix(NA_real_, 2, 2))
  r <- length(x)
  delta <- x - a
  z <- beta * delta
  r_i <- exp(z)
  
  S_r    <- sum(r_i)
  S_d_r  <- sum(delta * r_i)
  S_d2_r <- sum(delta^2 * r_i)
  
  J_bb <- r / beta^2 + S_d2_r + (n - r) * delta[r]^2 * r_i[r]
  J_aa <- beta^2 * (S_r + (n - r) * r_i[r])
  J_ba <- r - (S_r + (n - r) * r_i[r]) -
    beta * (S_d_r + (n - r) * delta[r] * r_i[r])
  
  matrix(
    c(J_bb, J_ba,
      J_ba, J_aa),
    nrow = 2, byrow = TRUE
  )
}

## 5. r and r* for beta (Type II)

r_beta_type2 <- function(beta0, x, n, mle_obj = NULL) {
  if (beta0 <= 0) stop("beta0 must be positive")
  if (is.null(mle_obj)) mle_obj <- mle_beta_type2(x, n)
  
  beta_hat <- mle_obj$beta_hat
  ll_hat   <- mle_obj$ll_max
  ll0      <- profile_ll_beta_type2(beta0, x, n)
  
  sign(beta_hat - beta0) * sqrt(max(0, 2 * (ll_hat - ll0)))
}

rstar_beta_type2 <- function(beta0, x, n, mle_obj = NULL) {
  if (beta0 <= 0) stop("beta0 must be positive")
  if (is.null(mle_obj)) mle_obj <- mle_beta_type2(x, n)
  
  beta_hat <- mle_obj$beta_hat
  a_hat    <- mle_obj$a_hat
  ll_hat   <- mle_obj$ll_max
  ll0      <- profile_ll_beta_type2(beta0, x, n)
  
  r <- sign(beta_hat - beta0) * sqrt(max(0, 2 * (ll_hat - ll0)))
  if (!is.finite(r) || abs(r) < 1e-10) return(0)
  
  a0 <- a_hat_beta_type2(beta0, x, n)
  
  phi_hat <- phi_type2(beta_hat, a_hat, x, n)
  phi0    <- phi_type2(beta0,    a0,    x, n)
  phi_diff <- phi_hat - phi0
  
  Phi_hat <- jacobian_phi_type2(beta_hat, a_hat, x, n)
  Phi0    <- jacobian_phi_type2(beta0,    a0,    x, n)
  phi_lambda <- Phi0[, 2]
  
  J_hat <- J_type2(beta_hat, a_hat, x, n)
  J0    <- J_type2(beta0,    a0,    x, n)
  J_ll  <- J0[2, 2]
  
  if (!is.finite(J_ll) || J_ll <= 0) return(r)
  
  det_num  <- det(cbind(phi_diff, phi_lambda))
  det_den  <- det(Phi_hat)
  detJ_hat <- det(J_hat)
  
  if (!is.finite(det_num) || !is.finite(det_den) || !is.finite(detJ_hat)) return(r)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)
  
  q_val <- (det_num / det_den) * sqrt(detJ_hat / J_ll)
  if (!is.finite(q_val)) return(r)
  if (q_val * r <= 0) q_val <- abs(q_val) * sign(r)
  
  out <- r + (1 / r) * log(q_val / r)
  if (!is.finite(out)) out <- r
  out
}

## 6. r and r* for log scale a (Type II)

prof_a_type2 <- function(a, x, n) {
  f <- function(logb) {
    beta <- exp(logb)
    -loglik_type2(beta, a, x, n)
  }
  opt <- optimize(f, interval = c(log(1e-3), log(1e3)))
  beta_hat <- exp(opt$minimum)
  ll_max   <- -opt$objective
  list(beta_hat = beta_hat, ll_max = ll_max)
}

r_a_type2 <- function(a0, x, n, mle_obj) {
  beta_hat <- mle_obj$beta_hat
  a_hat    <- mle_obj$a_hat
  ll_hat   <- mle_obj$ll_max
  
  pa0 <- prof_a_type2(a0, x, n)
  ll0 <- pa0$ll_max
  
  sign(a_hat - a0) * sqrt(max(0, 2 * (ll_hat - ll0)))
}

rstar_a_type2 <- function(a0, x, n, mle_obj) {
  beta_hat <- mle_obj$beta_hat
  a_hat    <- mle_obj$a_hat
  ll_hat   <- mle_obj$ll_max
  
  pa0 <- prof_a_type2(a0, x, n)
  beta0 <- pa0$beta_hat
  ll0   <- pa0$ll_max
  
  r <- sign(a_hat - a0) * sqrt(max(0, 2 * (ll_hat - ll0)))
  if (!is.finite(r) || abs(r) < 1e-10) return(0)
  
  phi_hat <- phi_type2(beta_hat, a_hat, x, n)
  phi0    <- phi_type2(beta0,    a0,    x, n)
  phi_diff <- phi_hat - phi0
  
  Phi_hat <- jacobian_phi_type2(beta_hat, a_hat, x, n)
  Phi0    <- jacobian_phi_type2(beta0,    a0,    x, n)
  
  ## (psi, lambda) = (a, beta), T is permutation
  Tmat <- matrix(c(0, 1,
                   1, 0),
                 nrow = 2, byrow = TRUE)
  Phi_hat_pl <- Phi_hat %*% Tmat
  phi_lambda <- Phi0 %*% c(1, 0)
  
  J_hat <- J_type2(beta_hat, a_hat, x, n)
  J0    <- J_type2(beta0,    a0,    x, n)
  J_ll  <- J0[1, 1]
  
  if (!is.finite(J_ll) || J_ll <= 0) return(r)
  
  detJ_hat <- det(J_hat)
  if (!is.finite(detJ_hat)) return(r)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)
  
  det_num <- det(cbind(phi_diff, phi_lambda))
  det_den <- det(Phi_hat_pl)
  if (!is.finite(det_num) || !is.finite(det_den)) return(r)
  
  q_val <- (det_num / det_den) * sqrt(detJ_hat / J_ll)
  if (!is.finite(q_val)) return(r)
  if (q_val * r <= 0) q_val <- abs(q_val) * sign(r)
  
  out <- r + (1 / r) * log(q_val / r)
  if (!is.finite(out)) out <- r
  out
}

## 7. r and r* for log quantile zeta_p (Type II)

prof_zeta_type2 <- function(psi, x, n, p) {
  s_p <- -log(1 - p)
  w_p <- log(s_p)
  f <- function(logb) {
    beta <- exp(logb)
    a    <- psi - w_p / beta
    -loglik_type2(beta, a, x, n)
  }
  opt <- optimize(f, interval = c(log(1e-3), log(1e3)))
  beta_hat <- exp(opt$minimum)
  ll_max   <- -opt$objective
  list(beta_hat = beta_hat, ll_max = ll_max)
}

r_zeta_type2 <- function(psi0, x, n, p, mle_obj) {
  beta_hat <- mle_obj$beta_hat
  a_hat    <- mle_obj$a_hat
  ll_hat   <- mle_obj$ll_max
  
  s_p <- -log(1 - p)
  w_p <- log(s_p)
  psi_hat <- a_hat + w_p / beta_hat
  
  p0 <- prof_zeta_type2(psi0, x, n, p)
  ll0 <- p0$ll_max
  
  sign(psi_hat - psi0) * sqrt(max(0, 2 * (ll_hat - ll0)))
}

rstar_zeta_type2 <- function(psi0, x, n, p, mle_obj) {
  beta_hat <- mle_obj$beta_hat
  a_hat    <- mle_obj$a_hat
  ll_hat   <- mle_obj$ll_max
  
  s_p <- -log(1 - p)
  w_p <- log(s_p)
  psi_hat <- a_hat + w_p / beta_hat
  
  p0 <- prof_zeta_type2(psi0, x, n, p)
  beta0 <- p0$beta_hat
  ll0   <- p0$ll_max
  
  r <- sign(psi_hat - psi0) * sqrt(max(0, 2 * (ll_hat - ll0)))
  if (!is.finite(r) || abs(r) < 1e-10) return(0)
  
  a0 <- psi0 - w_p / beta0
  
  phi_hat <- phi_type2(beta_hat, a_hat, x, n)
  phi0    <- phi_type2(beta0,    a0,    x, n)
  phi_diff <- phi_hat - phi0
  
  Phi_hat <- jacobian_phi_type2(beta_hat, a_hat, x, n)
  Phi0    <- jacobian_phi_type2(beta0,    a0,    x, n)
  
  T_hat <- matrix(c(0,      1,
                    1,  w_p / beta_hat^2),
                  nrow = 2, byrow = TRUE)
  T0    <- matrix(c(0,      1,
                    1,  w_p / beta0^2),
                  nrow = 2, byrow = TRUE)
  
  Phi_hat_pl <- Phi_hat %*% T_hat
  t_lambda0 <- c(1, w_p / beta0^2)
  phi_lambda <- Phi0 %*% t_lambda0
  
  J_hat <- J_type2(beta_hat, a_hat, x, n)
  J0    <- J_type2(beta0,    a0,    x, n)
  J_hat_pl <- t(T_hat) %*% J_hat %*% T_hat
  J0_pl    <- t(T0)   %*% J0   %*% T0
  J_ll     <- J0_pl[2, 2]
  
  if (!is.finite(J_ll) || J_ll <= 0) return(r)
  
  detJ_hat <- det(J_hat_pl)
  if (!is.finite(detJ_hat)) return(r)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)
  
  det_num <- det(cbind(phi_diff, phi_lambda))
  det_den <- det(Phi_hat_pl)
  if (!is.finite(det_num) || !is.finite(det_den)) return(r)
  
  q_val <- (det_num / det_den) * sqrt(detJ_hat / J_ll)
  if (!is.finite(q_val)) return(r)
  if (q_val * r <= 0) q_val <- abs(q_val) * sign(r)
  
  out <- r + (1 / r) * log(q_val / r)
  if (!is.finite(out)) out <- r
  out
}

## 8. Simulation function for Type II censoring

sim_type2_grid <- function(
    beta_vec  = c(0.5, 1, 5),
    eta_vec   = c(2,   1, 0.5),
    n         = 10,
    r_vec     = c(3, 5, 8),
    p_quant   = 0.1,
    n_rep     = 100,
    conf_levels = c(0.9, 0.95),
    seed      = 123
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(beta_vec) == length(eta_vec))
  
  z_levels <- qnorm((1 + conf_levels) / 2)
  names(z_levels) <- as.character(conf_levels)
  
  res_list <- list()
  idx <- 1
  
  for (k in seq_along(beta_vec)) {
    beta_true <- beta_vec[k]
    eta_true  <- eta_vec[k]
    a_true    <- log(eta_true)
    
    s_p <- -log(1 - p_quant)
    w_p <- log(s_p)
    Q_true <- eta_true * s_p^(1 / beta_true)
    zeta_true <- log(Q_true)
    
    for (r_fail in r_vec) {
      ## storage: param x method x conf_level
      counts <- array(0L, dim = c(3, 2, length(conf_levels)),
                      dimnames = list(
                        param  = c("shape_beta", "log_scale_a", "log_quantile_zeta"),
                        method = c("signed_root", "modified_root"),
                        level  = as.character(conf_levels)
                      ))
      
      for (rep in seq_len(n_rep)) {
        y <- rweibull(n, shape = beta_true, scale = eta_true)
        y_ord <- sort(y)
        y_fail <- y_ord[1:r_fail]
        x_fail <- log(y_fail)
        
        mle_obj <- mle_beta_type2(x_fail, n)
        
        ## shape
        r_b     <- r_beta_type2(beta_true, x_fail, n, mle_obj)
        rstar_b <- rstar_beta_type2(beta_true, x_fail, n, mle_obj)
        
        ## log scale
        r_a     <- r_a_type2(a_true, x_fail, n, mle_obj)
        rstar_a <- rstar_a_type2(a_true, x_fail, n, mle_obj)
        
        ## log quantile
        r_z     <- r_zeta_type2(zeta_true, x_fail, n, p_quant, mle_obj)
        rstar_z <- rstar_zeta_type2(zeta_true, x_fail, n, p_quant, mle_obj)
        
        for (lev in seq_along(conf_levels)) {
          z <- z_levels[lev]
          lev_name <- names(z_levels)[lev]
          
          ## shape
          if (is.finite(r_b) && abs(r_b) <= z) {
            counts["shape_beta", "signed_root", lev_name] <-
              counts["shape_beta", "signed_root", lev_name] + 1L
          }
          if (is.finite(rstar_b) && abs(rstar_b) <= z) {
            counts["shape_beta", "modified_root", lev_name] <-
              counts["shape_beta", "modified_root", lev_name] + 1L
          }
          
          ## log scale
          if (is.finite(r_a) && abs(r_a) <= z) {
            counts["log_scale_a", "signed_root", lev_name] <-
              counts["log_scale_a", "signed_root", lev_name] + 1L
          }
          if (is.finite(rstar_a) && abs(rstar_a) <= z) {
            counts["log_scale_a", "modified_root", lev_name] <-
              counts["log_scale_a", "modified_root", lev_name] + 1L
          }
          
          ## log quantile
          if (is.finite(r_z) && abs(r_z) <= z) {
            counts["log_quantile_zeta", "signed_root", lev_name] <-
              counts["log_quantile_zeta", "signed_root", lev_name] + 1L
          }
          if (is.finite(rstar_z) && abs(rstar_z) <= z) {
            counts["log_quantile_zeta", "modified_root", lev_name] <-
              counts["log_quantile_zeta", "modified_root", lev_name] + 1L
          }
        }
        
        if (rep %% 100 == 0) {
          cat("Type II - beta =", beta_true,
              "eta =", eta_true, "n =", n,
              "r_fail =", r_fail,
              "rep =", rep, "of", n_rep, "\r")
          flush.console()
        }
      }
      cat("\n")
      
      ## turn counts into data.frame
      for (param in dimnames(counts)$param) {
        for (meth in dimnames(counts)$method) {
          for (lev_name in dimnames(counts)$level) {
            res_list[[idx]] <- data.frame(
              param   = param,
              beta    = beta_true,
              eta     = eta_true,
              n       = n,
              r_fail  = r_fail,
              conf_level = as.numeric(lev_name),
              method  = meth,
              coverage = counts[param, meth, lev_name] / n_rep
            )
            idx <- idx + 1
          }
        }
      }
    }
  }
  
  do.call(rbind, res_list)
}

# =============================================================================
# Plotting helper
# =============================================================================

plot_type2_coverage <- function(df, conf_level_target) {
  df_sub <- df %>% filter(conf_level == conf_level_target)
  
  ggplot(df_sub,
         aes(x = r, y = coverage,
             color = method_lab, group = method_lab)) +
    geom_hline(yintercept = conf_level_target, linetype = "dashed") +
    geom_point(size = 2) +
    geom_line() +
    facet_grid(
      param_lab ~ beta_eta,
      labeller = labeller(
        param_lab = label_parsed,
        beta_eta  = label_parsed
      )
    ) +
    scale_x_continuous(
      breaks = sort(unique(df_sub$r))
    ) +
    scale_y_continuous(
      limits = c(min(0.8, min(df_sub$coverage) - 0.02), 1.0),
      breaks = seq(0.8, 1.0, by = 0.05)
    ) +
    labs(
      x = "Number of observed failures r (n = 10)",
      y = "Coverage probability",
      color = "Method",
      title = paste0("Type II censoring, nominal ", 100 * conf_level_target, "% intervals")
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 10),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
}

# =============================================================================
# Interval length functions
# =============================================================================

step_root <- function(target_fun, init_val, step, direction = c("left", "right"),
                      max_steps = 500, tol = 1e-4) {
  direction <- match.arg(direction)
  sign_ref <- sign(target_fun(init_val))
  if (sign_ref == 0 || is.na(sign_ref)) return(init_val)
  
  step <- abs(step) * if (direction == "left") -1 else 1
  cur_val <- init_val + step
  f_prev <- target_fun(init_val)
  
  for (i in 1:max_steps) {
    f_cur <- tryCatch(target_fun(cur_val), error = function(e) NA_real_)
    if (is.na(f_cur)) {
      cur_val <- cur_val + step
      next
    }
    if (sign(f_cur) != sign_ref || abs(f_cur) < tol) {
      # Refine the root location by linear interpolation
      if (f_prev != f_cur) {
        root <- cur_val - step + (0 - f_prev) / (f_cur - f_prev) * step
        return(root)
      }
      return(cur_val)
    }
    f_prev <- f_cur
    cur_val <- cur_val + step
  }
  return(cur_val)
}

## ============================================================================
## 1. Interval lengths for the shape parameter beta
## ============================================================================
beta_length_signed <- function(x, n, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle_obj <- tryCatch(mle_beta_type2(x, n), error = function(e) NULL)
  if (is.null(mle_obj) || !is.finite(mle_obj$beta_hat) || mle_obj$beta_hat <= 0) return(NA_real_)
  
  target <- function(b) {
    if (b <= 0) return(NA_real_)
    r_val <- tryCatch(r_beta_type2(b, x, n, mle_obj), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  
  beta_hat <- mle_obj$beta_hat
  step <- max(0.01, beta_hat * 0.05)  # Step size is approximately 5% of the MLE
  lower <- step_root(target, beta_hat, step, "left", max_steps = 500)
  upper <- step_root(target, beta_hat, step, "right", max_steps = 500)
  upper - lower
}

beta_length_modified <- function(x, n, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle_obj <- tryCatch(mle_beta_type2(x, n), error = function(e) NULL)
  if (is.null(mle_obj) || !is.finite(mle_obj$beta_hat) || mle_obj$beta_hat <= 0) return(NA_real_)
  
  target <- function(b) {
    if (b <= 0) return(NA_real_)
    r_val <- tryCatch(rstar_beta_type2(b, x, n, mle_obj), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  
  beta_hat <- mle_obj$beta_hat
  step <- max(0.01, beta_hat * 0.05)
  lower <- step_root(target, beta_hat, step, "left", max_steps = 500)
  upper <- step_root(target, beta_hat, step, "right", max_steps = 500)
  upper - lower
}

## ============================================================================
## 2. Interval lengths for the log-scale parameter a
## ============================================================================
a_length_signed <- function(x, n, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle_obj <- tryCatch(mle_beta_type2(x, n), error = function(e) NULL)
  if (is.null(mle_obj) || !is.finite(mle_obj$a_hat)) return(NA_real_)
  
  target <- function(a) {
    r_val <- tryCatch(r_a_type2(a, x, n, mle_obj), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  
  a_hat <- mle_obj$a_hat
  step <- max(0.01, abs(a_hat) * 0.05)
  lower <- step_root(target, a_hat, step, "left", max_steps = 500)
  upper <- step_root(target, a_hat, step, "right", max_steps = 500)
  upper - lower
}

a_length_modified <- function(x, n, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle_obj <- tryCatch(mle_beta_type2(x, n), error = function(e) NULL)
  if (is.null(mle_obj) || !is.finite(mle_obj$a_hat)) return(NA_real_)
  
  target <- function(a) {
    r_val <- tryCatch(rstar_a_type2(a, x, n, mle_obj), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  
  a_hat <- mle_obj$a_hat
  step <- max(0.01, abs(a_hat) * 0.05)
  lower <- step_root(target, a_hat, step, "left", max_steps = 500)
  upper <- step_root(target, a_hat, step, "right", max_steps = 500)
  upper - lower
}

## ============================================================================
## 3. Interval lengths for the log-quantile parameter zeta_p
## ============================================================================
zeta_length_signed <- function(x, n, p_quant, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle_obj <- tryCatch(mle_beta_type2(x, n), error = function(e) NULL)
  if (is.null(mle_obj) || !is.finite(mle_obj$beta_hat) || !is.finite(mle_obj$a_hat)) return(NA_real_)
  
  s_p <- -log(1 - p_quant)
  w_p <- log(s_p)
  
  target <- function(zeta) {
    r_val <- tryCatch(r_zeta_type2(zeta, x, n, p_quant, mle_obj), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  
  zeta_hat <- mle_obj$a_hat + w_p / mle_obj$beta_hat
  step <- max(0.01, abs(zeta_hat) * 0.05)
  lower <- step_root(target, zeta_hat, step, "left", max_steps = 500)
  upper <- step_root(target, zeta_hat, step, "right", max_steps = 500)
  upper - lower
}

zeta_length_modified <- function(x, n, p_quant, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle_obj <- tryCatch(mle_beta_type2(x, n), error = function(e) NULL)
  if (is.null(mle_obj) || !is.finite(mle_obj$beta_hat) || !is.finite(mle_obj$a_hat)) return(NA_real_)
  
  s_p <- -log(1 - p_quant)
  w_p <- log(s_p)
  
  target <- function(zeta) {
    r_val <- tryCatch(rstar_zeta_type2(zeta, x, n, p_quant, mle_obj), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  
  zeta_hat <- mle_obj$a_hat + w_p / mle_obj$beta_hat
  step <- max(0.01, abs(zeta_hat) * 0.05)
  lower <- step_root(target, zeta_hat, step, "left", max_steps = 500)
  upper <- step_root(target, zeta_hat, step, "right", max_steps = 500)
  upper - lower
}

## ============================================================================
## 4. Length simulation for Type II right censoring
## ============================================================================
sim_length_type2 <- function(
    beta_vec  = c(0.5, 1, 5),
    eta_vec   = c(2, 1, 0.5),
    n         = 10,
    r_vec     = c(3, 5, 8),
    p_quant   = 0.1,
    n_rep     = 10000,
    conf_level = 0.95,
    seed      = 123
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(beta_vec) == length(eta_vec))
  
  res_list <- list()
  idx <- 1
  
  for (k in seq_along(beta_vec)) {
    beta_true <- beta_vec[k]
    eta_true  <- eta_vec[k]
    a_true    <- log(eta_true)
    
    s_p <- -log(1 - p_quant)
    w_p <- log(s_p)
    zeta_true <- a_true + w_p / beta_true
    
    for (r_fail in r_vec) {
      # Store interval lengths for each replication
      len_beta_s <- numeric(n_rep)
      len_beta_m <- numeric(n_rep)
      len_a_s    <- numeric(n_rep)
      len_a_m    <- numeric(n_rep)
      len_zeta_s <- numeric(n_rep)
      len_zeta_m <- numeric(n_rep)
      
      for (rep in seq_len(n_rep)) {
        # Generate data and keep the first r_fail failures
        y_full <- rweibull(n, shape = beta_true, scale = eta_true)
        y_ord <- sort(y_full)
        y_fail <- y_ord[1:r_fail]
        x_fail <- log(y_fail)
        
        len_beta_s[rep] <- beta_length_signed(x_fail, n, conf_level)
        len_beta_m[rep] <- beta_length_modified(x_fail, n, conf_level)
        len_a_s[rep]    <- a_length_signed(x_fail, n, conf_level)
        len_a_m[rep]    <- a_length_modified(x_fail, n, conf_level)
        len_zeta_s[rep] <- zeta_length_signed(x_fail, n, p_quant, conf_level)
        len_zeta_m[rep] <- zeta_length_modified(x_fail, n, p_quant, conf_level)
        
        if (rep %% 100 == 0) {
          cat(sprintf("beta=%.1f eta=%.1f n=%d r=%d rep=%d/%d\r",
                      beta_true, eta_true, n, r_fail, rep, n_rep))
          flush.console()
        }
      }
      cat("\n")
      
      # Average interval length
      avg_len <- function(v) mean(v[is.finite(v)], na.rm = TRUE)
      df_beta <- data.frame(
        param = "shape_beta", beta = beta_true, eta = eta_true,
        n = n, r_fail = r_fail,
        method = c("signed_root", "modified_root"),
        avg_length = c(avg_len(len_beta_s), avg_len(len_beta_m))
      )
      df_a <- data.frame(
        param = "log_scale_a", beta = beta_true, eta = eta_true,
        n = n, r_fail = r_fail,
        method = c("signed_root", "modified_root"),
        avg_length = c(avg_len(len_a_s), avg_len(len_a_m))
      )
      df_zeta <- data.frame(
        param = "log_quantile_zeta", beta = beta_true, eta = eta_true,
        n = n, r_fail = r_fail,
        method = c("signed_root", "modified_root"),
        avg_length = c(avg_len(len_zeta_s), avg_len(len_zeta_m))
      )
      
      res_list[[idx]] <- rbind(df_beta, df_a, df_zeta)
      idx <- idx + 1
    }
  }
  do.call(rbind, res_list)
}

# # =============================================================================
# # Simulation: coverage probabilities
# # =============================================================================
# 
# res_type2 <- sim_type2_grid(
#   beta_vec  = c(0.5, 1, 5),
#   eta_vec   = c(2,   1, 0.5),
#   n         = 10,
#   r_vec     = c(3, 5, 8),
#   p_quant   = 0.1,
#   n_rep     = 10000,
#   conf_levels = c(0.9, 0.95),
#   seed      = 123
# )
# print(res_type2)
# 
# # =============================================================================
# # Plot: coverage probabilities
# # =============================================================================
# 
# library(dplyr)
# library(ggplot2)
# 
# res_type2 <- res_type2 %>%
#   rename(r = r_fail) %>%   # adjust this name if your column is slightly different
#   mutate(
#     param_lab = factor(
#       param,
#       levels = c("shape_beta", "log_scale_a", "log_quantile_zeta"),
#       labels = c(
#         "Shape~beta",
#         "Log~scale~a",
#         "Log~quantile~zeta[0.1]"
#       )
#     ),
#     beta_eta = dplyr::case_when(
#       beta == 0.5 & eta == 2   ~ "beta==0.5*','~~eta==2",
#       beta == 1   & eta == 1   ~ "beta==1*','~~eta==1",
#       beta == 5   & eta == 0.5 ~ "beta==5*','~~eta==0.5",
#       TRUE ~ "beta==NA*','~~eta==NA"
#     ),
#     beta_eta = factor(
#       beta_eta,
#       levels = c(
#         "beta==0.5*','~~eta==2",
#         "beta==1*','~~eta==1",
#         "beta==5*','~~eta==0.5"
#       )
#     ),
#     method_lab = factor(
#       method,
#       levels = c("signed_root", "modified_root"),
#       labels = c("Signed root", "Modified root")
#     ),
#     conf_lab = factor(
#       conf_level,
#       levels = c(0.9, 0.95),
#       labels = c("90% intervals", "95% intervals")
#     )
#   )
# 
# # Build and save plots
# #p_type2_90 <- plot_type2_coverage(res_type2, conf_level_target = 0.90)
# p_type2_95 <- plot_type2_coverage(res_type2, conf_level_target = 0.95)
# 
# #print(p_type2_90)
# print(p_type2_95)
# 
# #ggsave("weibull_typeII_coverage_90.pdf", p_type2_90, width = 7, height = 6)
# # ggsave("weibull_typeII_coverage_95.pdf", p_type2_95,
# # width = 7, height = 6)

# =============================================================================
# Simulation: average interval lengths
# =============================================================================

res_length <- sim_length_type2(
  beta_vec  = c(0.5, 1, 5),
  eta_vec   = c(2,   1, 0.5),
  n         = 10,
  r_vec     = c(3, 5, 8),
  p_quant   = 0.1,
  n_rep     = 10000,        # Can be increased to 1000 or higher
  conf_level = 0.95,
  seed      = 123
)

library(dplyr)
library(ggplot2)
library(ggh4x)

plot_data <- res_length %>%
  mutate(
    param_lab = factor(param,
                       levels = c("shape_beta", "log_scale_a", "log_quantile_zeta"),
                       labels = c("Shape~beta", "Log~scale~a", "Log~quantile~zeta[0.1]")),
    beta_eta = case_when(
      beta == 0.5 & eta == 2   ~ "beta==0.5*','~~eta==2",
      beta == 1   & eta == 1   ~ "beta==1*','~~eta==1",
      beta == 5   & eta == 0.5 ~ "beta==5*','~~eta==0.5"
    ),
    beta_eta = factor(beta_eta,
                      levels = c("beta==0.5*','~~eta==2",
                                 "beta==1*','~~eta==1",
                                 "beta==5*','~~eta==0.5")),
    method_lab = factor(method,
                        levels = c("signed_root", "modified_root"),
                        labels = c("Signed root", "Modified root"))
  )

p_length <- ggplot(plot_data, aes(x = r_fail, y = avg_length,
                                  color = method_lab, group = method_lab)) +
  geom_point(size = 2) +
  geom_line() +
  # Independent y-axis settings
  facet_grid2(
    param_lab ~ beta_eta,
    scales = "free_y",        # Free y-axis scales
    independent = "y",        # Make the y-axis independent across panels
    labeller = labeller(
      param_lab = label_parsed,
      beta_eta = label_parsed
    )
  ) +
  # End independent y-axis settings
  scale_x_continuous(breaks = sort(unique(plot_data$r_fail))) +
  labs(x = "Number of observed failures r (n = 10)",
       y = "Average interval length",
       color = "Method") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(size = 10),
        legend.position = "bottom")

print(p_length)
