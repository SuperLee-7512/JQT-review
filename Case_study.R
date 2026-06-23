# ============================================================
# Likelihood-based confidence intervals for Weibull ALT models
# ============================================================


# ------------------------------------------------------------
# Core likelihood, score-like quantities, and numerical derivatives
# ------------------------------------------------------------

alt_loglik <- function(theta, y, x) {
  mu0   <- theta[1]
  mu1   <- theta[2]
  sigma <- theta[3]

  if (!is.finite(mu0) || !is.finite(mu1) ||
      !is.finite(sigma) || sigma <= 0) {
    return(-1e20)
  }

  z <- (y - mu0 - mu1 * x) / sigma
  if (any(!is.finite(z))) return(-1e20)

  z <- pmin(z, 700)              # Avoid overflow in exp(z).
  r <- exp(z)
  if (any(!is.finite(r))) return(-1e20)

  ell <- sum(z - r) - length(y) * log(sigma)
  if (!is.finite(ell)) ell <- -1e20
  ell
}

alt_phi <- function(theta, y, x) {
  mu0   <- theta[1]
  mu1   <- theta[2]
  sigma <- theta[3]
  if (sigma <= 0) return(rep(NA_real_, 3))

  z <- (y - mu0 - mu1 * x) / sigma
  z <- pmin(z, 700)              # Avoid overflow in exp(z).
  r <- exp(z)

  S0 <- sum(1 - r)
  S1 <- sum(x * (1 - r))
  S2 <- sum(z * (1 - r))

  c(S0, S1, S2) / sigma
}

alt_phi_jac <- function(theta, y, x, h = 1e-5) {
  base <- alt_phi(theta, y, x)
  if (any(!is.finite(base))) return(matrix(NA_real_, 3, 3))

  J <- matrix(NA_real_, 3, 3)
  for (j in 1:3) {
    step <- h * max(1, abs(theta[j]))
    th_p <- theta; th_m <- theta
    th_p[j] <- theta[j] + step
    th_m[j] <- theta[j] - step

    if (j == 3 && th_m[3] <= 0) {
      fp <- alt_phi(th_p, y, x)
      J[, j] <- (fp - base) / step
    } else {
      fp <- alt_phi(th_p, y, x)
      fm <- alt_phi(th_m, y, x)
      J[, j] <- (fp - fm) / (2 * step)
    }
  }
  J
}

alt_J <- function(theta, y, x, eps = 1e-4) {
  p  <- length(theta)
  J  <- matrix(NA_real_, p, p)
  f0 <- alt_loglik(theta, y, x)

  for (k in 1:p) {
    hk   <- eps * max(1, abs(theta[k]))
    th_p <- theta; th_m <- theta
    th_p[k] <- theta[k] + hk
    th_m[k] <- theta[k] - hk

    f_p <- alt_loglik(th_p, y, x)
    f_m <- alt_loglik(th_m, y, x)

    H_kk <- (f_p - 2 * f0 + f_m) / hk^2
    J[k, k] <- -H_kk

    if (k < p) {
      for (l in (k + 1):p) {
        hl <- eps * max(1, abs(theta[l]))
        th_pp <- theta; th_pm <- theta
        th_mp <- theta; th_mm <- theta
        th_pp[k] <- theta[k] + hk; th_pp[l] <- theta[l] + hl
        th_pm[k] <- theta[k] + hk; th_pm[l] <- theta[l] - hl
        th_mp[k] <- theta[k] - hk; th_mp[l] <- theta[l] + hl
        th_mm[k] <- theta[k] - hk; th_mm[l] <- theta[l] - hl

        f_pp <- alt_loglik(th_pp, y, x)
        f_pm <- alt_loglik(th_pm, y, x)
        f_mp <- alt_loglik(th_mp, y, x)
        f_mm <- alt_loglik(th_mm, y, x)

        H_kl <- (f_pp - f_pm - f_mp + f_mm) / (4 * hk * hl)
        J[k, l] <- J[l, k] <- -H_kl
      }
    }
  }
  J
}


# ------------------------------------------------------------
# Maximum likelihood estimation and profile likelihood
# ------------------------------------------------------------

alt_mle <- function(y, x) {
  n <- length(y)
  fit_lm <- lm(y ~ x)
  mu0_start <- coef(fit_lm)[1]
  mu1_start <- coef(fit_lm)[2]
  res <- y - fitted(fit_lm)
  sigma_start <- sqrt(sum(res^2) / n)
  if (!is.finite(sigma_start) || sigma_start <= 0) sigma_start <- 1

  par_start <- c(mu0_start, mu1_start, log(sigma_start))

  obj <- function(par) {
    mu0   <- par[1]
    mu1   <- par[2]
    sigma <- exp(par[3])
    -alt_loglik(c(mu0, mu1, sigma), y, x)
  }

  opt <- optim(par_start, obj, method = "BFGS",
               control = list(maxit = 2000))

  mu0_hat   <- opt$par[1]
  mu1_hat   <- opt$par[2]
  sigma_hat <- exp(opt$par[3])
  ll_max    <- -opt$value

  list(mu0_hat = mu0_hat,
       mu1_hat = mu1_hat,
       sigma_hat = sigma_hat,
       ll_max = ll_max,
       convergence = opt$convergence)
}

alt_profile_psi <- function(psi0, y, x, x_use, p, mle_obj) {
  s_p <- -log(1 - p)
  w_p <- log(s_p)

  mu0_hat   <- mle_obj$mu0_hat
  mu1_hat   <- mle_obj$mu1_hat

  obj_lambda <- function(lambda) {
    mu0 <- lambda[1]
    mu1 <- lambda[2]
    sigma <- (psi0 - mu0 - mu1 * x_use) / w_p
    if (sigma <= 0) return(1e10)
    -alt_loglik(c(mu0, mu1, sigma), y, x)
  }

  opt <- optim(c(mu0_hat, mu1_hat), obj_lambda,
               method = "BFGS", control = list(maxit = 2000))

  mu0_0 <- opt$par[1]
  mu1_0 <- opt$par[2]
  sigma_0 <- (psi0 - mu0_0 - mu1_0 * x_use) / w_p

  theta_0 <- c(mu0_0, mu1_0, sigma_0)
  ll0     <- alt_loglik(theta_0, y, x)

  list(theta_0 = theta_0,
       lambda_0 = c(mu0_0, mu1_0),
       ll0 = ll0,
       convergence = opt$convergence)
}


# ------------------------------------------------------------
# Signed likelihood root and modified signed likelihood root
# ------------------------------------------------------------

alt_r_psi <- function(psi0, y, x, x_use, p, mle_obj = NULL) {
  if (is.null(mle_obj)) mle_obj <- alt_mle(y, x)
  mu0_hat   <- mle_obj$mu0_hat
  mu1_hat   <- mle_obj$mu1_hat
  sigma_hat <- mle_obj$sigma_hat
  ll_hat    <- mle_obj$ll_max

  s_p <- -log(1 - p)
  w_p <- log(s_p)
  psi_hat <- mu0_hat + mu1_hat * x_use + sigma_hat * w_p

  prof0 <- alt_profile_psi(psi0, y, x, x_use, p, mle_obj)
  ll0   <- prof0$ll0

  r <- sign(psi_hat - psi0) * sqrt(max(0, 2 * (ll_hat - ll0)))
  r
}

alt_rstar_psi <- function(psi0, y, x, x_use, p, mle_obj = NULL) {
  if (is.null(mle_obj)) mle_obj <- alt_mle(y, x)
  mu0_hat   <- mle_obj$mu0_hat
  mu1_hat   <- mle_obj$mu1_hat
  sigma_hat <- mle_obj$sigma_hat
  ll_hat    <- mle_obj$ll_max

  s_p <- -log(1 - p)
  w_p <- log(s_p)
  psi_hat <- mu0_hat + mu1_hat * x_use + sigma_hat * w_p

  prof0 <- alt_profile_psi(psi0, y, x, x_use, p, mle_obj)
  theta_0  <- prof0$theta_0
  lambda_0 <- prof0$lambda_0
  ll0      <- prof0$ll0

  r <- sign(psi_hat - psi0) * sqrt(max(0, 2 * (ll_hat - ll0)))
  if (!is.finite(r) || abs(r) < 1e-10) return(0)

  theta_hat <- c(mu0_hat, mu1_hat, sigma_hat)

  phi_hat <- alt_phi(theta_hat, y, x)
  phi_0   <- alt_phi(theta_0,   y, x)
  if (any(!is.finite(phi_hat)) || any(!is.finite(phi_0))) return(r)

  phi_diff <- phi_hat - phi_0

  Phi_hat <- alt_phi_jac(theta_hat, y, x)
  Phi_0   <- alt_phi_jac(theta_0,   y, x)
  if (any(!is.finite(Phi_hat)) || any(!is.finite(Phi_0))) return(r)

  # Transformation matrix T = d theta / d(psi, lambda1, lambda2)'.
  #   mu0   = lambda1
  #   mu1   = lambda2
  #   sigma = (psi - lambda1 - lambda2 * x_use) / w_p
  #
  # Columns: derivatives with respect to psi, lambda1, and lambda2.
  T_mat <- matrix(
    c(0,      1,            0,
      0,      0,            1,
      1/w_p, -1/w_p, -x_use / w_p),
    nrow = 3, byrow = TRUE
  )

  Phi_hat_pl <- Phi_hat %*% T_mat

  T_lambda <- T_mat[, 2:3, drop = FALSE]  # 3 x 2
  phi_lambda <- Phi_0 %*% T_lambda        # 3 x 2

  J_hat_theta <- alt_J(theta_hat, y, x)
  J_0_theta   <- alt_J(theta_0,   y, x)
  if (any(!is.finite(J_hat_theta)) || any(!is.finite(J_0_theta))) return(r)

  J_hat_pl <- t(T_mat) %*% J_hat_theta %*% T_mat
  J_0_pl   <- t(T_mat) %*% J_0_theta   %*% T_mat

  J_ll <- J_0_pl[2:3, 2:3, drop = FALSE]

  detJ_hat <- det(J_hat_pl)
  detJ_ll  <- det(J_ll)

  if (!is.finite(detJ_hat) || !is.finite(detJ_ll) || detJ_ll <= 0) return(r)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)

  num_mat <- cbind(phi_diff, phi_lambda)  # 3 x 3
  det_num <- det(num_mat)
  det_den <- det(Phi_hat_pl)

  if (!is.finite(det_num) || !is.finite(det_den)) return(r)

  q_val <- (det_num / det_den) * sqrt(detJ_hat / detJ_ll)
  if (!is.finite(q_val)) return(r)

  if (q_val * r <= 0) q_val <- abs(q_val) * sign(r)

  rstar <- r + (1 / r) * log(q_val / r)
  if (!is.finite(rstar)) rstar <- r
  rstar
}


# ------------------------------------------------------------
# Confidence interval helper functions
# ------------------------------------------------------------

find_ci_one <- function(psi_hat, r_fun, z_crit, sigma_hat) {
  step <- 4 * sigma_hat
  low  <- psi_hat - step
  up   <- psi_hat + step

  for (i in 1:20) {
    val <- r_fun(low)
    if (is.finite(val) && val > z_crit) break
    low <- low - step
  }
  for (i in 1:20) {
    val <- r_fun(up)
    if (is.finite(val) && val < -z_crit) break
    up <- up + step
  }

  val_low <- r_fun(low)
  val_up  <- r_fun(up)
  if (!is.finite(val_low) || !is.finite(val_up)) {
    return(c(NA_real_, NA_real_))
  }
  if (!(val_low > z_crit && val_up < -z_crit)) {
    return(c(NA_real_, NA_real_))
  }

  lower <- uniroot(function(psi) r_fun(psi) - z_crit, c(low, psi_hat))$root
  upper <- uniroot(function(psi) r_fun(psi) + z_crit, c(psi_hat, up))$root
  c(lower, upper)
}

find_ci_general <- function(psi_hat, r_fun, z_crit, step_scale = 2) {
  step <- step_scale * max(1, abs(psi_hat))
  low  <- psi_hat - step
  up   <- psi_hat + step

  for (i in 1:20) {
    val <- r_fun(low)
    if (is.finite(val) && val > z_crit) break
    low <- low - step
  }
  for (i in 1:20) {
    val <- r_fun(up)
    if (is.finite(val) && val < -z_crit) break
    up <- up + step
  }

  val_low <- r_fun(low)
  val_up  <- r_fun(up)
  if (!is.finite(val_low) || !is.finite(val_up)) {
    return(c(NA_real_, NA_real_))
  }
  if (!(val_low > z_crit && val_up < -z_crit)) {
    return(c(NA_real_, NA_real_))
  }

  lower <- uniroot(function(psi) r_fun(psi) - z_crit, c(low, psi_hat))$root
  upper <- uniroot(function(psi) r_fun(psi) + z_crit, c(psi_hat, up))$root
  c(lower, upper)
}


# ------------------------------------------------------------
# Parametric bootstrap confidence interval for psi
# ------------------------------------------------------------

boot_ci_psi <- function(y, x, x_use, p_quant, mle_obj, B = 199, conf_level = 0.95) {
  s_p <- -log(1 - p_quant)
  w_p <- log(s_p)
  unique_x <- sort(unique(x))
  n_per <- as.integer(table(x)[as.character(unique_x)])
  beta_hat <- 1 / mle_obj$sigma_hat
  psi_boot <- numeric(B)
  success <- 0

  for (b in 1:B) {
    y_b <- c()
    x_b <- c()
    for (j in seq_along(unique_x)) {
      xj <- unique_x[j]
      nj <- n_per[j]
      eta_j <- exp(mle_obj$mu0_hat + mle_obj$mu1_hat * xj)
      T_bj <- rweibull(nj, shape = beta_hat, scale = eta_j)
      y_b <- c(y_b, log(T_bj))
      x_b <- c(x_b, rep(xj, nj))
    }
    mle_b <- try(alt_mle(y_b, x_b), silent = TRUE)
    if (inherits(mle_b, "try-error") || mle_b$convergence != 0) next
    success <- success + 1
    psi_boot[success] <- mle_b$mu0_hat + mle_b$mu1_hat * x_use + mle_b$sigma_hat * w_p
  }
  if (success < 0.8 * B) return(c(NA_real_, NA_real_))
  psi_boot <- psi_boot[1:success]
  alpha <- 1 - conf_level
  quantile(psi_boot, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
}


# ------------------------------------------------------------
# Main simulation routine
# ------------------------------------------------------------

run_full_simulation <- function(
    mu0_true, mu1_true, sigma_true,
    x_levels      = c(28.84, 34.68, 38.02),
    n_per_level   = c(6, 6, 6),
    x_use         = 20,
    p_quant       = 0.1,
    conf_level    = 0.95,
    n_rep         = 1000,
    boot_B        = 199,
    seed          = 123
) {
  if (!is.null(seed)) set.seed(seed)

  s_p <- -log(1 - p_quant)
  w_p <- log(s_p)
  psi_true <- mu0_true + mu1_true * x_use + sigma_true * w_p
  beta_true <- 1 / sigma_true
  z_crit <- qnorm((1 + conf_level) / 2)

  psi_hat_vec <- numeric(n_rep)
  MSE_vec <- c()
  len_r <- len_rstar <- len_boot <- numeric(n_rep)
  cover_r <- cover_rstar <- cover_boot <- logical(n_rep)

  valid <- 0
  for (rep in 1:n_rep) {
    y_list <- list()
    x_list <- list()
    for (j in seq_along(x_levels)) {
      xj <- x_levels[j]
      nj <- n_per_level[j]
      eta_j <- exp(mu0_true + mu1_true * xj)
      Tj <- rweibull(nj, shape = beta_true, scale = eta_j)
      y_list[[j]] <- log(Tj)
      x_list[[j]] <- rep(xj, nj)
    }
    y <- unlist(y_list)
    x <- unlist(x_list)

    mle_obj <- try(alt_mle(y, x), silent = TRUE)
    if (inherits(mle_obj, "try-error") || mle_obj$convergence != 0) next
    valid <- valid + 1
    psi_hat <- mle_obj$mu0_hat + mle_obj$mu1_hat * x_use + mle_obj$sigma_hat * w_p
    psi_hat_vec[valid] <- psi_hat

    MSE_vec <- rbind(MSE_vec,
                     (c(mle_obj$mu0_hat, mle_obj$mu1_hat, mle_obj$sigma_hat) - c(mu0_true, mu1_true, sigma_true)) ** 2)

    r_fun <- function(psi) alt_r_psi(psi, y, x, x_use, p_quant, mle_obj)
    ci_r <- find_ci_one(psi_hat, r_fun, z_crit, mle_obj$sigma_hat)
    if (!any(is.na(ci_r))) {
      cover_r[valid] <- (psi_true >= ci_r[1] && psi_true <= ci_r[2])
      len_r[valid] <- ci_r[2] - ci_r[1]
    }

    rs_fun <- function(psi) alt_rstar_psi(psi, y, x, x_use, p_quant, mle_obj)
    ci_rs <- find_ci_one(psi_hat, rs_fun, z_crit, mle_obj$sigma_hat)
    if (!any(is.na(ci_rs))) {
      cover_rstar[valid] <- (psi_true >= ci_rs[1] && psi_true <= ci_rs[2])
      len_rstar[valid] <- ci_rs[2] - ci_rs[1]
    }

    ci_boot <- boot_ci_psi(y, x, x_use, p_quant, mle_obj, B = boot_B, conf_level)
    if (!any(is.na(ci_boot))) {
      cover_boot[valid] <- (psi_true >= ci_boot[1] && psi_true <= ci_boot[2])
      len_boot[valid] <- ci_boot[2] - ci_boot[1]
    }

    if (valid %% 100 == 0) cat("已完成", valid, "/", n_rep, "次有效模拟\r")
  }
  cat("\n")

  psi_hat_vec <- psi_hat_vec[1:valid]
  cover_r <- cover_r[1:valid]; cover_rstar <- cover_rstar[1:valid]; cover_boot <- cover_boot[1:valid]
  len_r <- len_r[1:valid]; len_rstar <- len_rstar[1:valid]; len_boot <- len_boot[1:valid]

  mse <- mean((psi_hat_vec - psi_true)^2, na.rm = TRUE)
  cov_r <- mean(cover_r, na.rm = TRUE) * 100
  cov_rs <- mean(cover_rstar, na.rm = TRUE) * 100
  cov_boot <- mean(cover_boot, na.rm = TRUE) * 100
  mean_len_r <- mean(len_r, na.rm = TRUE)
  mean_len_rs <- mean(len_rstar, na.rm = TRUE)
  mean_len_boot <- mean(len_boot, na.rm = TRUE)

  list(
    settings = list(mu0_true=mu0_true, mu1_true=mu1_true, sigma_true=sigma_true,
                    x_levels=x_levels, n_per_level=n_per_level, x_use=x_use,
                    p_quant=p_quant, conf_level=conf_level, n_rep=n_rep, boot_B=boot_B),
    results = data.frame(
      quantity = c("MSE", "Coverage_r(%)", "Coverage_rstar(%)", "Coverage_boot(%)",
                   "MeanLength_r", "MeanLength_rstar", "MeanLength_boot"),
      value = c(mse, cov_r, cov_rs, cov_boot, mean_len_r, mean_len_rs, mean_len_boot)
    ),
    raw = list(psi_hat = psi_hat_vec, len_r = len_r, len_rstar = len_rstar, len_boot = len_boot,
               cover_r = cover_r, cover_rstar = cover_rstar, cover_boot = cover_boot, MSE = MSE_vec)
  )
}


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

mu0_true   <- 19.72264
mu1_true   <- -0.4245567
sigma_true <- 0.22713
beta_true  <- 1 / sigma_true
x_levels    <- c(28.84, 34.68, 38.02)
n_per_level <- c(6, 6, 6)
x_use       <- 20
conf_level  <- 0.95
n_rep       <- 1000
seed        <- 123
p_grid      <- 0.1


# ------------------------------------------------------------
# Example run based on parameter estimates from the original data
# ------------------------------------------------------------

if (!exists("original_stress")) {
  original_stress <- c(28.84, 34.68, 38.02)
  original_lifetimes <- list(
    c(1637, 1658, 2437, 1709, 1267, 1785),
    c(132, 96, 76, 122, 115, 87),
    c(43, 22, 39, 42, 41, 37)
  )
  t_all <- unlist(original_lifetimes)
  x_all <- rep(original_stress, times = sapply(original_lifetimes, length))
  y_all <- log(t_all)
  mle_true <- alt_mle(y_all, x_all)
  mu0_true <- mle_true$mu0_hat
  mu1_true <- mle_true$mu1_hat
  sigma_true <- mle_true$sigma_hat
}

set.seed(123)
sim_result <- run_full_simulation(
  mu0_true = mu0_true,
  mu1_true = mu1_true,
  sigma_true = sigma_true,
  x_levels = original_stress,
  n_per_level = c(6, 6, 6),
  x_use = 20,           # Modify as needed, for example 3 for comparison with code1.
  p_quant = 0.1,
  conf_level = 0.95,
  n_rep = 10000,
  boot_B = 1000,
  seed = 123
)

print(sim_result$results)


# ------------------------------------------------------------
# Plot interval lengths
# ------------------------------------------------------------

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  df_len <- data.frame(
    Method = rep(c("r", "r*", "Bootstrap"), each = length(sim_result$raw$len_r)),
    Length = c(sim_result$raw$len_r, sim_result$raw$len_rstar, sim_result$raw$len_boot)
  )
  p <- ggplot(df_len, aes(x = Method, y = Length, fill = Method)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    labs(y = "Interval length", title = "") +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 10),
      legend.position = "bottom"
    )
} else {
  boxplot(sim_result$raw$len_r, sim_result$raw$len_rstar, sim_result$raw$len_boot,
          names = c("r", "r*", "Bootstrap"), ylab = "Interval length",
          main = "Interval lengths for p = 0.1")
}


ggsave("case_IT.pdf", p, width = 6, height = 4)
