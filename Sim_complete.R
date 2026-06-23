# Length revision for complete Weibull samples

rm(list = ls())

###############################################################################
## Complete data - Weibull simulation
##
## Parameters:
##   - Shape beta
##   - Log-scale a = log(eta)
##   - Log-quantile zeta_p = log(Q_p)
##
## Methods:
##   - Signed root r
##   - Modified root r*
##   - Technometrics calibrated profile likelihood
###############################################################################


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

clikelevels <- function(N) {
  # Single-parameter constants from Diaz-Frances & Sapatinas (2010).
  cs90 <- 0.2585 - (1 / (0.775  + (1.641  * N)))
  cs95 <- 0.1465 - (1 / (2.757  + (2.082  * N)))
  cs99 <- 0.0362 - (1 / (18.294 + (5.229  * N)))
  csig <- c(cs90, cs95, cs99)
  names(csig) <- c("90%", "95%", "99%")
  csig
}

get_c_level <- function(N, conf_level) {
  lev_name <- paste0(round(100 * conf_level), "%")
  csig <- clikelevels(N)
  if (!lev_name %in% names(csig)) {
    stop("clikelevels only implemented for 90%, 95%, 99%")
  }
  unname(csig[lev_name])
}


# -----------------------------------------------------------------------------
# Shape parameter: profile in beta for complete data
# -----------------------------------------------------------------------------

profile_ll_beta <- function(beta, x) {
  if (beta <= 0) return(-Inf)
  n  <- length(x)
  bx <- beta * x
  m  <- max(bx)
  w  <- exp(bx - m)
  K1 <- exp(m) * sum(w)
  n * log(beta) + beta * sum(x) - n * log(K1)
}

a_hat_beta <- function(beta, x) {
  n  <- length(x)
  bx <- beta * x
  m  <- max(bx)
  w  <- exp(bx - m)
  K1 <- exp(m) * sum(w)
  (log(K1) - log(n)) / beta
}

delta_r_beta_a <- function(beta, a, x) {
  delta <- x - a
  r     <- exp(beta * delta)
  list(delta = delta, r = r)
}

phi_beta_a <- function(beta, a, x) {
  dr    <- delta_r_beta_a(beta, a, x)
  delta <- dr$delta
  r     <- dr$r
  S0    <- sum(1 - r)
  S1    <- sum((1 - r) * delta)
  phi1  <- beta * S0
  phi2  <- beta^2 * S1
  c(phi1, phi2)
}

phi_theta_beta_a <- function(beta, a, x) {
  dr    <- delta_r_beta_a(beta, a, x)
  delta <- dr$delta
  r     <- dr$r

  S0 <- sum(1 - r)
  S1 <- sum((1 - r) * delta)

  dS0_a <- sum(beta * r)
  dS0_b <- -sum(delta * r)

  dS1_a <- sum(beta * r * delta - (1 - r))
  dS1_b <- -sum(delta^2 * r)

  phi1_a <- beta^2 * sum(r)
  phi1_b <- S0 - beta * sum(delta * r)

  phi2_a <- beta^2 * dS1_a
  phi2_b <- 2 * beta * S1 + beta^2 * dS1_b

  matrix(
    c(phi1_b, phi1_a,
      phi2_b, phi2_a),
    nrow = 2, byrow = TRUE
  )
}

phi_lambda_beta_a <- function(beta, a, x) {
  dr    <- delta_r_beta_a(beta, a, x)
  delta <- dr$delta
  r     <- dr$r
  phi1_a <- beta^2 * sum(r)
  dS1_a  <- sum(beta * r * delta - (1 - r))
  phi2_a <- beta^2 * dS1_a
  c(phi1_a, phi2_a)
}

J_beta_a <- function(beta, a, x) {
  n     <- length(x)
  dr    <- delta_r_beta_a(beta, a, x)
  delta <- dr$delta
  r     <- dr$r

  J_bb <- n / beta^2 + sum(delta^2 * r)
  J_aa <- beta^2 * sum(r)
  J_ba <- n - sum(r) - beta * sum(delta * r)

  matrix(
    c(J_bb, J_ba,
      J_ba, J_aa),
    nrow = 2, byrow = TRUE
  )
}

mle_beta_profile <- function(x) {
  f <- function(logb) {
    beta <- exp(logb)
    -profile_ll_beta(beta, x)
  }
  opt <- optimize(f, interval = c(log(1e-3), log(1e3)))
  beta_hat <- exp(opt$minimum)
  ll_max   <- -opt$objective
  a_hat    <- a_hat_beta(beta_hat, x)
  list(beta_hat = beta_hat, a_hat = a_hat, ll_max = ll_max)
}

r_beta <- function(beta0, x, mle_obj = NULL) {
  if (beta0 <= 0) stop("beta0 must be positive")
  if (is.null(mle_obj)) mle_obj <- mle_beta_profile(x)
  beta_hat <- mle_obj$beta_hat
  ll_hat   <- mle_obj$ll_max
  ll0      <- profile_ll_beta(beta0, x)
  r <- sign(beta_hat - beta0) * sqrt(max(0, 2 * (ll_hat - ll0)))
  r
}

rstar_beta <- function(beta0, x, mle_obj = NULL) {
  if (beta0 <= 0) stop("beta0 must be positive")
  if (is.null(mle_obj)) mle_obj <- mle_beta_profile(x)

  beta_hat <- mle_obj$beta_hat
  a_hat    <- mle_obj$a_hat
  ll_hat   <- mle_obj$ll_max

  ll0 <- profile_ll_beta(beta0, x)
  r   <- sign(beta_hat - beta0) * sqrt(max(0, 2 * (ll_hat - ll0)))
  if (!is.finite(r) || abs(r) < 1e-10) return(0)

  a0 <- a_hat_beta(beta0, x)

  phi_hat    <- phi_beta_a(beta_hat, a_hat, x)
  phi_constr <- phi_beta_a(beta0, a0, x)
  phi_diff   <- phi_hat - phi_constr

  phi_lambda_constr <- phi_lambda_beta_a(beta0, a0, x)
  phi_theta_hat     <- phi_theta_beta_a(beta_hat, a_hat, x)

  J_hat    <- J_beta_a(beta_hat, a_hat, x)
  J_constr <- J_beta_a(beta0,    a0,    x)
  J_aa_constr <- J_constr[2, 2]

  if (!is.finite(J_aa_constr) || J_aa_constr <= 0) return(r)

  num_det  <- det(cbind(phi_diff, phi_lambda_constr))
  den_det  <- det(phi_theta_hat)
  detJ_hat <- det(J_hat)

  if (!is.finite(num_det) || !is.finite(den_det) || !is.finite(detJ_hat)) return(r)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)

  q_val <- (num_det / den_det) * sqrt(detJ_hat / J_aa_constr)
  if (!is.finite(q_val)) return(r)
  if (q_val * r <= 0) q_val <- abs(q_val) * sign(r)

  out <- r + (1 / r) * log(q_val / r)
  if (!is.finite(out)) out <- r
  out
}


# -----------------------------------------------------------------------------
# Log-Weibull parametrization: (mu, sigma) for scale and quantile
# -----------------------------------------------------------------------------

loglik_mu_sigma <- function(mu, sigma, x) {
  if (sigma <= 0) return(-Inf)
  n <- length(x)
  z <- (x - mu) / sigma
  -n * log(sigma) + sum(z) - sum(exp(z))
}

phi_mu_sigma <- function(mu, sigma, x) {
  z <- (x - mu) / sigma
  r <- exp(z)
  phi1 <- (1 / sigma) * sum(1 - r)
  phi2 <- (1 / sigma) * sum((1 - r) * z)
  c(phi1, phi2)
}

phi_mu_sigma_jac <- function(mu, sigma, x) {
  z <- (x - mu) / sigma
  r <- exp(z)
  n <- length(x)

  S_r   <- sum(r)
  S_z   <- sum(z)
  S_zr  <- sum(z * r)
  S_z2r <- sum(z^2 * r)

  # d phi1 / d mu, d phi1 / d sigma
  phi1_mu    <- S_r / sigma^2
  phi1_sigma <- (-n + S_r + S_zr) / sigma^2

  # d phi2 / d mu, d phi2 / d sigma
  phi2_mu    <- (-n + S_r + S_zr) / sigma^2
  phi2_sigma <- (-2 * S_z + 2 * S_zr + S_z2r) / sigma^2

  matrix(
    c(phi1_mu,    phi1_sigma,
      phi2_mu,    phi2_sigma),
    nrow = 2, byrow = TRUE
  )
}

J_mu_sigma <- function(mu, sigma, x) {
  if (sigma <= 0) return(matrix(NA_real_, 2, 2))
  n <- length(x)
  z <- (x - mu) / sigma
  r <- exp(z)

  J_mm <- (1 / sigma^2) * sum(r)
  J_ms <- (1 / sigma^2) * (n - sum(r) + sum(r * z))
  J_ss <- (1 / sigma^2) * (-n - 2 * sum(z) + sum(r * z * (2 + z)))

  matrix(
    c(J_mm, J_ms,
      J_ms, J_ss),
    nrow = 2, byrow = TRUE
  )
}


# -----------------------------------------------------------------------------
# Profile in mu (log scale) and r, r*
# -----------------------------------------------------------------------------

profile_ll_mu <- function(mu, x, sigma_ref) {
  f <- function(logsigma) {
    sigma <- exp(logsigma)
    -loglik_mu_sigma(mu, sigma, x)
  }
  logsig_ref <- log(sigma_ref)
  opt <- optimize(f, interval = c(logsig_ref - 5, logsig_ref + 5))
  -opt$objective
}

sigma_hat_mu <- function(mu, x, sigma_ref) {
  f <- function(logsigma) {
    sigma <- exp(logsigma)
    -loglik_mu_sigma(mu, sigma, x)
  }
  logsig_ref <- log(sigma_ref)
  opt <- optimize(f, interval = c(logsig_ref - 5, logsig_ref + 5))
  exp(opt$minimum)
}

r_mu <- function(mu0, x, mu_hat, sigma_hat) {
  lp_hat <- loglik_mu_sigma(mu_hat, sigma_hat, x)
  sigma0_hat <- sigma_hat_mu(mu0, x, sigma_hat)
  lp0 <- loglik_mu_sigma(mu0, sigma0_hat, x)
  sign(mu_hat - mu0) * sqrt(max(0, 2 * (lp_hat - lp0)))
}

rstar_mu <- function(mu0, x, mu_hat, sigma_hat) {
  lp_hat <- loglik_mu_sigma(mu_hat, sigma_hat, x)
  sigma0_hat <- sigma_hat_mu(mu0, x, sigma_hat)
  lp0 <- loglik_mu_sigma(mu0, sigma0_hat, x)
  r <- sign(mu_hat - mu0) * sqrt(max(0, 2 * (lp_hat - lp0)))
  if (!is.finite(r) || abs(r) < 1e-10) return(0)

  phi_hat    <- phi_mu_sigma(mu_hat, sigma_hat, x)
  phi_constr <- phi_mu_sigma(mu0, sigma0_hat, x)
  phi_diff   <- phi_hat - phi_constr

  Phi_hat    <- phi_mu_sigma_jac(mu_hat, sigma_hat, x)
  Phi_constr <- phi_mu_sigma_jac(mu0, sigma0_hat, x)
  phi_lambda <- Phi_constr[, 2]

  J_hat    <- J_mu_sigma(mu_hat, sigma_hat, x)
  J_constr <- J_mu_sigma(mu0,  sigma0_hat, x)
  J_ll     <- J_constr[2, 2]

  if (!is.finite(J_ll) || J_ll <= 0) return(r)

  num_det  <- det(cbind(phi_diff, phi_lambda))
  den_det  <- det(Phi_hat)
  detJ_hat <- det(J_hat)

  if (!is.finite(num_det) || !is.finite(den_det) || !is.finite(detJ_hat)) return(r)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)

  q_val <- (num_det / den_det) * sqrt(detJ_hat / J_ll)
  if (!is.finite(q_val)) return(r)
  if (q_val * r <= 0) q_val <- abs(q_val) * sign(r)

  out <- r + (1 / r) * log(q_val / r)
  if (!is.finite(out)) out <- r
  out
}


# -----------------------------------------------------------------------------
# Profile in zeta_p (log quantile) and r, r*
# -----------------------------------------------------------------------------

profile_ll_zeta <- function(zeta, x, w_p, sigma_ref) {
  f <- function(logsigma) {
    sigma <- exp(logsigma)
    mu    <- zeta - sigma * w_p
    -loglik_mu_sigma(mu, sigma, x)
  }
  logsig_ref <- log(sigma_ref)
  opt <- optimize(f, interval = c(logsig_ref - 5, logsig_ref + 5))
  -opt$objective
}

sigma_hat_zeta <- function(zeta, x, w_p, sigma_ref) {
  f <- function(logsigma) {
    sigma <- exp(logsigma)
    mu    <- zeta - sigma * w_p
    -loglik_mu_sigma(mu, sigma, x)
  }
  logsig_ref <- log(sigma_ref)
  opt <- optimize(f, interval = c(logsig_ref - 5, logsig_ref + 5))
  exp(opt$minimum)
}

r_zeta <- function(zeta0, x, mu_hat, sigma_hat, w_p) {
  lp_hat <- loglik_mu_sigma(mu_hat, sigma_hat, x)
  sigma0_hat <- sigma_hat_zeta(zeta0, x, w_p, sigma_hat)
  mu0 <- zeta0 - sigma0_hat * w_p
  lp0 <- loglik_mu_sigma(mu0, sigma0_hat, x)
  sign(mu_hat + sigma_hat * w_p - zeta0) *
    sqrt(max(0, 2 * (lp_hat - lp0)))
}

rstar_zeta <- function(zeta0, x, mu_hat, sigma_hat, w_p) {
  lp_hat <- loglik_mu_sigma(mu_hat, sigma_hat, x)
  sigma0_hat <- sigma_hat_zeta(zeta0, x, w_p, sigma_hat)
  mu0 <- zeta0 - sigma0_hat * w_p
  lp0 <- loglik_mu_sigma(mu0, sigma0_hat, x)
  r <- sign(mu_hat + sigma_hat * w_p - zeta0) *
    sqrt(max(0, 2 * (lp_hat - lp0)))
  if (!is.finite(r) || abs(r) < 1e-10) return(0)

  phi_hat    <- phi_mu_sigma(mu_hat, sigma_hat, x)
  phi_constr <- phi_mu_sigma(mu0, sigma0_hat, x)
  phi_diff   <- phi_hat - phi_constr

  Phi_hat    <- phi_mu_sigma_jac(mu_hat, sigma_hat, x)
  Phi_constr <- phi_mu_sigma_jac(mu0,  sigma0_hat, x)

  # Mapping (psi, lambda) = (zeta, sigma):
  # mu = psi - lambda * w_p, sigma = lambda.
  T_mat <- matrix(c(1, -w_p,
                    0,  1),
                  nrow = 2, byrow = TRUE)

  Phi_hat_pl <- Phi_hat %*% T_mat
  t_lambda <- c(-w_p, 1)
  phi_lambda <- Phi_constr %*% t_lambda

  J_hat_0    <- J_mu_sigma(mu_hat, sigma_hat, x)
  J_constr_0 <- J_mu_sigma(mu0,    sigma0_hat, x)
  J_hat <- t(T_mat) %*% J_hat_0 %*% T_mat
  J_constr <- t(T_mat) %*% J_constr_0 %*% T_mat
  J_ll <- J_constr[2, 2]

  if (!is.finite(J_ll) || J_ll <= 0) return(r)

  detJ_hat <- det(J_hat)
  if (!is.finite(detJ_hat)) return(r)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)

  num_det <- det(cbind(phi_diff, phi_lambda))
  den_det <- det(Phi_hat_pl)
  if (!is.finite(num_det) || !is.finite(den_det)) return(r)

  q_val <- (num_det / den_det) * sqrt(detJ_hat / J_ll)
  if (!is.finite(q_val)) return(r)
  if (q_val * r <= 0) q_val <- abs(q_val) * sign(r)

  out <- r + (1 / r) * log(q_val / r)
  if (!is.finite(out)) out <- r
  out
}


# -----------------------------------------------------------------------------
# Coverage simulation for complete data
# -----------------------------------------------------------------------------

sim_complete_grid <- function(
    beta_vec  = c(0.5, 1, 5),
    eta_vec   = c(2,   1, 0.5),
    n_vec     = c(5, 10, 20),
    p_quant   = 0.1,
    n_rep     = 1000,
    conf_level = 0.95,
    seed      = 123
) {
  if (!is.null(seed)) set.seed(seed)

  stopifnot(length(beta_vec) == length(eta_vec))

  z <- qnorm((1 + conf_level) / 2)
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

    for (n in n_vec) {
      # Technometrics c(n).
      c_level <- get_c_level(n, conf_level)
      log_c   <- log(c_level)

      # Counters.
      cover_beta_r     <- 0
      cover_beta_rstar <- 0
      cover_beta_tech  <- 0

      cover_a_r     <- 0
      cover_a_rstar <- 0
      cover_a_tech  <- 0

      cover_z_r     <- 0
      cover_z_rstar <- 0
      cover_z_tech  <- 0

      for (rep in seq_len(n_rep)) {
        y <- rweibull(n, shape = beta_true, scale = eta_true)
        x <- log(y)

        # MLE in (beta, a).
        mle_obj <- mle_beta_profile(x)
        beta_hat <- mle_obj$beta_hat
        a_hat    <- mle_obj$a_hat
        ll_beta_hat <- mle_obj$ll_max

        mu_hat    <- a_hat
        sigma_hat <- 1 / beta_hat

        # Shape.
        r_b     <- r_beta(beta_true, x, mle_obj)
        rstar_b <- rstar_beta(beta_true, x, mle_obj)
        ll_beta_true <- profile_ll_beta(beta_true, x)

        if (is.finite(r_b)     && abs(r_b)     <= z) cover_beta_r     <- cover_beta_r     + 1L
        if (is.finite(rstar_b) && abs(rstar_b) <= z) cover_beta_rstar <- cover_beta_rstar + 1L
        if (ll_beta_true - ll_beta_hat >= log_c)      cover_beta_tech  <- cover_beta_tech  + 1L

        # Log scale mu = log(eta).
        r_a     <- r_mu(a_true, x, mu_hat, sigma_hat)
        rstar_a <- rstar_mu(a_true, x, mu_hat, sigma_hat)

        lp_mu_hat  <- loglik_mu_sigma(mu_hat, sigma_hat, x)
        sigma_mu_true <- sigma_hat_mu(a_true, x, sigma_hat)
        lp_mu_true <- loglik_mu_sigma(a_true, sigma_mu_true, x)

        if (is.finite(r_a)     && abs(r_a)     <= z) cover_a_r     <- cover_a_r     + 1L
        if (is.finite(rstar_a) && abs(rstar_a) <= z) cover_a_rstar <- cover_a_rstar + 1L
        if (lp_mu_true - lp_mu_hat >= log_c)          cover_a_tech  <- cover_a_tech  + 1L

        # Log quantile zeta = log Q_p.
        r_z     <- r_zeta(zeta_true, x, mu_hat, sigma_hat, w_p)
        rstar_z <- rstar_zeta(zeta_true, x, mu_hat, sigma_hat, w_p)

        lp_z_hat <- lp_mu_hat
        sigma_z_true <- sigma_hat_zeta(zeta_true, x, w_p, sigma_hat)
        mu_z_true    <- zeta_true - sigma_z_true * w_p
        lp_z_true    <- loglik_mu_sigma(mu_z_true, sigma_z_true, x)

        if (is.finite(r_z)     && abs(r_z)     <= z) cover_z_r     <- cover_z_r     + 1L
        if (is.finite(rstar_z) && abs(rstar_z) <= z) cover_z_rstar <- cover_z_rstar + 1L
        if (lp_z_true - lp_z_hat >= log_c)           cover_z_tech  <- cover_z_tech  + 1L

        if (rep %% 100 == 0) {
          cat("Complete data - beta =", beta_true,
              "eta =", eta_true, "n =", n,
              "rep =", rep, "of", n_rep, "\r")
          flush.console()
        }
      }
      cat("\n")

      res_list[[idx]] <- data.frame(
        param   = "shape_beta",
        beta    = beta_true,
        eta     = eta_true,
        n       = n,
        method  = c("signed_root", "modified_root", "technometricsPL"),
        coverage = c(cover_beta_r, cover_beta_rstar, cover_beta_tech) / n_rep
      )
      idx <- idx + 1

      res_list[[idx]] <- data.frame(
        param   = "log_scale_a",
        beta    = beta_true,
        eta     = eta_true,
        n       = n,
        method  = c("signed_root", "modified_root", "technometricsPL"),
        coverage = c(cover_a_r, cover_a_rstar, cover_a_tech) / n_rep
      )
      idx <- idx + 1

      res_list[[idx]] <- data.frame(
        param   = "log_quantile_zeta",
        beta    = beta_true,
        eta     = eta_true,
        n       = n,
        method  = c("signed_root", "modified_root", "technometricsPL"),
        coverage = c(cover_z_r, cover_z_rstar, cover_z_tech) / n_rep
      )
      idx <- idx + 1
    }
  }

  do.call(rbind, res_list)
}


# -----------------------------------------------------------------------------
# Stepwise root-finding function
# -----------------------------------------------------------------------------

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
    if (sign(f_cur) != sign_ref || f_cur == 0) {
      # Linear interpolation.
      if (f_prev != f_cur) {
        root <- cur_val - step + (0 - f_prev) / (f_cur - f_prev) * step
        return(root)
      }
      return(cur_val)
    }
    f_prev <- f_cur
    cur_val <- cur_val + step
  }
  return(cur_val)  # Return the boundary value if no sign change is found.
}


# -----------------------------------------------------------------------------
# Interval lengths for the shape parameter beta
# -----------------------------------------------------------------------------

beta_ci_length_signed <- function(x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(NA_real_)
  beta_hat <- mle$beta_hat
  if (is.na(beta_hat) || beta_hat <= 0) return(NA_real_)

  target <- function(b) {
    if (b <= 0) return(NA_real_)
    r_val <- tryCatch(r_beta(b, x, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(0.01, beta_hat * 0.02)   # Step size is 2% of the MLE.
  beta_L <- step_root(target, beta_hat, step, "left", max_steps = 200)
  beta_U <- step_root(target, beta_hat, step, "right", max_steps = 200)
  beta_U - beta_L
}

beta_ci_length_modified <- function(x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(NA_real_)
  beta_hat <- mle$beta_hat
  if (is.na(beta_hat) || beta_hat <= 0) return(NA_real_)

  target <- function(b) {
    if (b <= 0) return(NA_real_)
    r_val <- tryCatch(rstar_beta(b, x, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(0.01, beta_hat * 0.02)
  beta_L <- step_root(target, beta_hat, step, "left", max_steps = 200)
  beta_U <- step_root(target, beta_hat, step, "right", max_steps = 200)
  beta_U - beta_L
}

beta_ci_length_tech <- function(x, conf_level = 0.95) {
  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(NA_real_)
  beta_hat <- mle$beta_hat
  a_hat    <- mle$a_hat
  if (is.na(beta_hat) || beta_hat <= 0) return(NA_real_)

  n <- length(x)
  c_level <- get_c_level(n, conf_level)
  if (is.na(c_level)) return(NA_real_)
  log_c <- log(c_level)

  J <- tryCatch(J_beta_a(beta_hat, a_hat, x), error = function(e) NULL)
  if (is.null(J)) return(NA_real_)
  J_inv <- solve(J)
  I_beta <- 1 / J_inv[1, 1]
  if (I_beta <= 0) return(NA_real_)

  term <- sqrt(-log_c / (2 * I_beta)) / beta_hat
  beta1 <- beta_hat * (1 - term)^2
  beta2 <- beta_hat * (1 + term)^2
  beta2 - beta1
}


# -----------------------------------------------------------------------------
# Interval lengths for the log-scale parameter a (mu)
# -----------------------------------------------------------------------------

mu_ci_length_signed <- function(x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(NA_real_)
  mu_hat    <- mle$a_hat
  sigma_hat <- 1 / mle$beta_hat
  if (is.na(sigma_hat) || sigma_hat <= 0) return(NA_real_)

  target <- function(mu) {
    r_val <- tryCatch(r_mu(mu, x, mu_hat, sigma_hat), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(0.01, sigma_hat * 0.02)   # Use the Gumbel scale as reference.
  mu_L <- step_root(target, mu_hat, step, "left", max_steps = 200)
  mu_U <- step_root(target, mu_hat, step, "right", max_steps = 200)
  mu_U - mu_L
}

mu_ci_length_modified <- function(x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(NA_real_)
  mu_hat    <- mle$a_hat
  sigma_hat <- 1 / mle$beta_hat
  if (is.na(sigma_hat) || sigma_hat <= 0) return(NA_real_)

  target <- function(mu) {
    r_val <- tryCatch(rstar_mu(mu, x, mu_hat, sigma_hat), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(0.01, sigma_hat * 0.02)
  mu_L <- step_root(target, mu_hat, step, "left", max_steps = 200)
  mu_U <- step_root(target, mu_hat, step, "right", max_steps = 200)
  mu_U - mu_L
}

mu_ci_length_tech <- function(x, conf_level = 0.95) {
  n <- length(x)
  c_level <- get_c_level(n, conf_level)
  if (is.na(c_level)) return(NA_real_)
  threshold <- log(c_level)

  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(NA_real_)
  mu_hat    <- mle$a_hat
  sigma_hat <- 1 / mle$beta_hat
  if (is.na(sigma_hat) || sigma_hat <= 0) return(NA_real_)

  lp_hat <- loglik_mu_sigma(mu_hat, sigma_hat, x)
  target <- function(mu) {
    lp_mu <- tryCatch(profile_ll_mu(mu, x, sigma_hat), error = function(e) NA_real_)
    if (is.na(lp_mu)) return(NA_real_)
    lp_mu - lp_hat - threshold
  }
  step <- max(0.01, sigma_hat * 0.02)
  mu_L <- step_root(target, mu_hat, step, "left", max_steps = 200)
  mu_U <- step_root(target, mu_hat, step, "right", max_steps = 200)
  mu_U - mu_L
}


# -----------------------------------------------------------------------------
# Interval lengths for the log-quantile parameter zeta_p
# -----------------------------------------------------------------------------

zeta_ci_length_signed <- function(x, w_p, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(NA_real_)
  mu_hat    <- mle$a_hat
  sigma_hat <- 1 / mle$beta_hat
  zeta_hat  <- mu_hat + sigma_hat * w_p
  if (is.na(sigma_hat) || sigma_hat <= 0) return(NA_real_)

  target <- function(zeta) {
    r_val <- tryCatch(r_zeta(zeta, x, mu_hat, sigma_hat, w_p), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(0.01, sigma_hat * 0.02)
  zeta_L <- step_root(target, zeta_hat, step, "left", max_steps = 200)
  zeta_U <- step_root(target, zeta_hat, step, "right", max_steps = 200)
  zeta_U - zeta_L
}

zeta_ci_length_modified <- function(x, w_p, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(NA_real_)
  mu_hat    <- mle$a_hat
  sigma_hat <- 1 / mle$beta_hat
  zeta_hat  <- mu_hat + sigma_hat * w_p
  if (is.na(sigma_hat) || sigma_hat <= 0) return(NA_real_)

  target <- function(zeta) {
    r_val <- tryCatch(rstar_zeta(zeta, x, mu_hat, sigma_hat, w_p), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(0.01, sigma_hat * 0.02)
  zeta_L <- step_root(target, zeta_hat, step, "left", max_steps = 200)
  zeta_U <- step_root(target, zeta_hat, step, "right", max_steps = 200)
  zeta_U - zeta_L
}

zeta_ci_length_tech <- function(x, w_p, conf_level = 0.95) {
  n <- length(x)
  c_level <- get_c_level(n, conf_level)
  if (is.na(c_level)) return(NA_real_)
  threshold <- log(c_level)

  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(NA_real_)
  mu_hat    <- mle$a_hat
  sigma_hat <- 1 / mle$beta_hat
  zeta_hat  <- mu_hat + sigma_hat * w_p
  if (is.na(sigma_hat) || sigma_hat <= 0) return(NA_real_)

  lp_hat <- loglik_mu_sigma(mu_hat, sigma_hat, x)
  target <- function(zeta) {
    lp_zeta <- tryCatch(profile_ll_zeta(zeta, x, w_p, sigma_hat), error = function(e) NA_real_)
    if (is.na(lp_zeta)) return(NA_real_)
    lp_zeta - lp_hat - threshold
  }
  step <- max(0.01, sigma_hat * 0.02)
  zeta_L <- step_root(target, zeta_hat, step, "left", max_steps = 200)
  zeta_U <- step_root(target, zeta_hat, step, "right", max_steps = 200)
  zeta_U - zeta_L
}


# -----------------------------------------------------------------------------
# Interval-length simulation for complete data
# -----------------------------------------------------------------------------

sim_length_grid <- function(
    beta_vec  = c(0.5, 1, 5),
    eta_vec   = c(2,   1, 0.5),
    n_vec     = c(5, 10, 20),
    p_quant   = 0.1,
    n_rep     = 200,
    conf_level = 0.95,
    seed      = 123
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(beta_vec) == length(eta_vec))
  w_p <- log(-log(1 - p_quant))

  res_list <- list()
  idx <- 1

  for (k in seq_along(beta_vec)) {
    beta_true <- beta_vec[k]
    eta_true  <- eta_vec[k]

    for (n in n_vec) {
      # Preallocate vectors.
      len_beta_r     <- numeric(n_rep)
      len_beta_rstar <- numeric(n_rep)
      len_beta_tech  <- numeric(n_rep)

      len_mu_r     <- numeric(n_rep)
      len_mu_rstar <- numeric(n_rep)
      len_mu_tech  <- numeric(n_rep)

      len_zeta_r     <- numeric(n_rep)
      len_zeta_rstar <- numeric(n_rep)
      len_zeta_tech  <- numeric(n_rep)

      for (rep in seq_len(n_rep)) {
        y <- rweibull(n, shape = beta_true, scale = eta_true)
        x <- log(y)

        len_beta_r[rep]     <- beta_ci_length_signed(x, conf_level)
        len_beta_rstar[rep] <- beta_ci_length_modified(x, conf_level)
        len_beta_tech[rep]  <- beta_ci_length_tech(x, conf_level)

        len_mu_r[rep]     <- mu_ci_length_signed(x, conf_level)
        len_mu_rstar[rep] <- mu_ci_length_modified(x, conf_level)
        len_mu_tech[rep]  <- mu_ci_length_tech(x, conf_level)

        len_zeta_r[rep]     <- zeta_ci_length_signed(x, w_p, conf_level)
        len_zeta_rstar[rep] <- zeta_ci_length_modified(x, w_p, conf_level)
        len_zeta_tech[rep]  <- zeta_ci_length_tech(x, w_p, conf_level)
      }

      avg_len <- function(v) mean(v, na.rm = TRUE)
      methods <- c("signed_root", "modified_root", "technometricsPL")

      df_beta <- data.frame(
        param = "shape_beta", beta = beta_true, eta = eta_true,
        n = n, method = methods,
        avg_length = c(avg_len(len_beta_r), avg_len(len_beta_rstar), avg_len(len_beta_tech))
      )
      df_mu <- data.frame(
        param = "log_scale_a", beta = beta_true, eta = eta_true,
        n = n, method = methods,
        avg_length = c(avg_len(len_mu_r), avg_len(len_mu_rstar), avg_len(len_mu_tech))
      )
      df_zeta <- data.frame(
        param = "log_quantile_zeta", beta = beta_true, eta = eta_true,
        n = n, method = methods,
        avg_length = c(avg_len(len_zeta_r), avg_len(len_zeta_rstar), avg_len(len_zeta_tech))
      )

      res_list[[idx]] <- rbind(df_beta, df_mu, df_zeta)
      idx <- idx + 1
    }
  }

  do.call(rbind, res_list)
}


# =============================================================================
# Main execution
# =============================================================================

# -----------------------------------------------------------------------------
# Coverage simulation
# -----------------------------------------------------------------------------

res_complete <- sim_complete_grid(
  beta_vec  = c(0.5, 1, 5),
  eta_vec   = c(2,   1, 0.5),
  n_vec     = c(5, 10, 20),
  p_quant   = 0.1,
  n_rep     = 10000,
  conf_level = 0.95,
  seed      = 123
)

print(res_complete)


# -----------------------------------------------------------------------------
# Coverage plot
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

res_complete <- res_complete %>%
  mutate(
    param_lab = factor(
      param,
      levels = c("shape_beta", "log_scale_a", "log_quantile_zeta"),
      labels = c(
        "Shape~beta",
        "Log~scale~a",
        "Log~quantile~zeta[0.1]"
      )
    ),
    beta_eta = dplyr::case_when(
      beta == 0.5 & eta == 2   ~ "beta==0.5*','~~eta==2",
      beta == 1   & eta == 1   ~ "beta==1*','~~eta==1",
      beta == 5   & eta == 0.5 ~ "beta==5*','~~eta==0.5",
      TRUE ~ "beta==NA*','~~eta==NA"
    ),
    beta_eta = factor(
      beta_eta,
      levels = c(
        "beta==0.5*','~~eta==2",
        "beta==1*','~~eta==1",
        "beta==5*','~~eta==0.5"
      )
    ),
    method_lab = factor(
      method,
      levels = c("signed_root", "modified_root", "technometricsPL"),
      labels = c("Signed root", "Modified root", "Profile")
    )
  )

p_cov <- ggplot(
  res_complete,
  aes(x = n, y = coverage, color = method_lab, group = method_lab)
) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(
    param_lab ~ beta_eta,
    labeller = labeller(
      param_lab = label_parsed,
      beta_eta  = label_parsed
    )
  ) +
  scale_x_continuous(breaks = c(5, 10, 20)) +
  scale_y_continuous(
    limits = c(0.85, 1.00),
    breaks = seq(0.85, 1.00, by = 0.05)
  ) +
  labs(
    x = "Sample size n",
    y = "Coverage probability",
    color = "Method"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  )

print(p_cov)

# ggsave("weibull_complete_coverage.pdf", p_cov, width = 7, height = 6)


# -----------------------------------------------------------------------------
# Interval-length simulation
# -----------------------------------------------------------------------------

res_length <- sim_length_grid(
  beta_vec  = c(0.5, 1, 5),
  eta_vec   = c(2,   1, 0.5),
  n_vec     = c(5, 10, 20),
  p_quant   = 0.1,
  n_rep     = 10000,        # Increase to 1000 if needed.
  conf_level = 0.95,
  seed      = 123
)

# save(res_length, file = "complete_buchang.Rdata")

library(ggplot2)
library(dplyr)

plot_data <- res_length %>%
  mutate(
    param_lab = factor(
      param,
      levels = c("shape_beta", "log_scale_a", "log_quantile_zeta"),
      labels = c("Shape~beta", "Log~scale~a", "Log~quantile~zeta[0.1]")
    ),
    beta_eta = case_when(
      beta == 0.5 & eta == 2   ~ "beta==0.5*','~~eta==2",
      beta == 1   & eta == 1   ~ "beta==1*','~~eta==1",
      beta == 5   & eta == 0.5 ~ "beta==5*','~~eta==0.5"
    ),
    beta_eta = factor(
      beta_eta,
      levels = c(
        "beta==0.5*','~~eta==2",
        "beta==1*','~~eta==1",
        "beta==5*','~~eta==0.5"
      )
    ),
    method_lab = factor(
      method,
      levels = c("signed_root", "modified_root", "technometricsPL"),
      labels = c("Signed root", "Modified root", "Profile")
    )
  )

p_len <- ggplot(plot_data, aes(x = n, y = avg_length, color = method_lab, group = method_lab)) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.5) +
  facet_grid2(
    param_lab ~ beta_eta,
    scales = "free",          # 行列均自由缩放
    independent = "all",      # 每个面板拥有独立坐标轴
    labeller = labeller(
      param_lab = label_parsed,
      beta_eta  = label_parsed
    )
  ) +
  scale_x_continuous(breaks = c(5, 10, 20)) +
  labs(x = "Sample size n", y = "Average interval length", color = "Method") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  )

print(p_len)
