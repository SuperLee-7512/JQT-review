
#########################################################
## 1. Likelihood for log Weibull ALT model
#########################################################

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
  
  z <- pmin(z, 700)              # avoid overflow
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
  z <- pmin(z, 700)
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

alt_mle <- function(y, x) {
  n   <- length(y)
  lm0 <- lm(y ~ x)
  mu0_start <- coef(lm0)[1]
  mu1_start <- coef(lm0)[2]
  res <- y - fitted(lm0)
  sigma_start <- sqrt(sum(res^2) / n)
  if (!is.finite(sigma_start) || sigma_start <= 0) sigma_start <- 1
  
  par_start <- c(mu0_start, mu1_start, log(sigma_start))
  
  obj <- function(par) {
    mu0   <- par[1]
    mu1   <- par[2]
    sigma <- exp(par[3])
    val   <- alt_loglik(c(mu0, mu1, sigma), y, x)
    if (!is.finite(val)) return(1e20)
    -val
  }
  
  opt <- optim(par_start, obj, method = "BFGS",
               control = list(maxit = 2000))
  
  list(
    mu0_hat   = opt$par[1],
    mu1_hat   = opt$par[2],
    sigma_hat = exp(opt$par[3]),
    ll_max    = -opt$value,
    convergence = opt$convergence
  )
}

#########################################################
## 2. Profile, r and r* for log quantile at use stress
#########################################################

alt_profile_psi <- function(psi0, y, x, x_use, p, mle_obj) {
  s_p <- -log(1 - p); w_p <- log(s_p)
  mu0_hat <- mle_obj$mu0_hat
  mu1_hat <- mle_obj$mu1_hat
  
  obj_lambda <- function(lambda) {
    mu0   <- lambda[1]
    mu1   <- lambda[2]
    sigma <- (psi0 - mu0 - mu1 * x_use) / w_p
    if (!is.finite(sigma) || sigma <= 0) return(1e20)
    val <- alt_loglik(c(mu0, mu1, sigma), y, x)
    if (!is.finite(val)) return(1e20)
    -val
  }
  
  opt <- try(
    optim(c(mu0_hat, mu1_hat), obj_lambda,
          method = "BFGS", control = list(maxit = 2000)),
    silent = TRUE
  )
  
  if (inherits(opt, "try-error")) {
    opt <- optim(c(mu0_hat, mu1_hat), obj_lambda,
                 method = "Nelder-Mead",
                 control = list(maxit = 4000))
  }
  
  mu0_0   <- opt$par[1]
  mu1_0   <- opt$par[2]
  sigma_0 <- (psi0 - mu0_0 - mu1_0 * x_use) / w_p
  
  if (!is.finite(sigma_0) || sigma_0 <= 0) {
    theta_0 <- rep(NA_real_, 3)
    ll0     <- -Inf
  } else {
    theta_0 <- c(mu0_0, mu1_0, sigma_0)
    ll0     <- alt_loglik(theta_0, y, x)
  }
  
  list(theta_0 = theta_0, ll0 = ll0)
}

alt_r_psi <- function(psi0, y, x, x_use, p, mle_obj) {
  mu0_hat   <- mle_obj$mu0_hat
  mu1_hat   <- mle_obj$mu1_hat
  sigma_hat <- mle_obj$sigma_hat
  ll_hat    <- mle_obj$ll_max
  
  s_p <- -log(1 - p); w_p <- log(s_p)
  psi_hat <- mu0_hat + mu1_hat * x_use + sigma_hat * w_p
  
  prof0 <- alt_profile_psi(psi0, y, x, x_use, p, mle_obj)
  ll0   <- prof0$ll0
  
  if (!is.finite(ll0)) {
    return(sign(psi_hat - psi0) * 1e6)
  }
  
  sign(psi_hat - psi0) * sqrt(pmax(0, 2 * (ll_hat - ll0)))
}

alt_rstar_psi <- function(psi0, y, x, x_use, p, mle_obj) {
  mu0_hat   <- mle_obj$mu0_hat
  mu1_hat   <- mle_obj$mu1_hat
  sigma_hat <- mle_obj$sigma_hat
  ll_hat    <- mle_obj$ll_max
  
  s_p <- -log(1 - p); w_p <- log(s_p)
  psi_hat <- mu0_hat + mu1_hat * x_use + sigma_hat * w_p
  
  prof0  <- alt_profile_psi(psi0, y, x, x_use, p, mle_obj)
  theta0 <- prof0$theta_0
  ll0    <- prof0$ll0
  
  r <- alt_r_psi(psi0, y, x, x_use, p, mle_obj)
  if (!all(is.finite(theta0)) || !is.finite(ll0) ||
      !is.finite(r) || abs(r) < 1e-10) {
    return(r)
  }
  
  theta_hat <- c(mu0_hat, mu1_hat, sigma_hat)
  
  phi_hat <- alt_phi(theta_hat, y, x)
  phi_0   <- alt_phi(theta0,   y, x)
  if (any(!is.finite(phi_hat)) || any(!is.finite(phi_0))) return(r)
  
  phi_diff <- phi_hat - phi_0
  
  Phi_hat <- alt_phi_jac(theta_hat, y, x)
  Phi_0   <- alt_phi_jac(theta0,   y, x)
  if (any(!is.finite(Phi_hat)) || any(!is.finite(Phi_0))) return(r)
  
  T_mat <- matrix(
    c(0,      1,            0,
      0,      0,            1,
      1/w_p, -1/w_p, -x_use / w_p),
    nrow = 3, byrow = TRUE
  )
  
  Phi_hat_pl <- Phi_hat %*% T_mat
  T_lambda   <- T_mat[, 2:3, drop = FALSE]
  phi_lambda <- Phi_0 %*% T_lambda
  
  J_hat_theta <- alt_J(theta_hat, y, x)
  J_0_theta   <- alt_J(theta0,   y, x)
  if (any(!is.finite(J_hat_theta)) || any(!is.finite(J_0_theta))) return(r)
  
  J_hat_pl <- t(T_mat) %*% J_hat_theta %*% T_mat
  J_0_pl   <- t(T_mat) %*% J_0_theta   %*% T_mat
  J_ll     <- J_0_pl[2:3, 2:3, drop = FALSE]
  
  detJ_hat <- det(J_hat_pl)
  detJ_ll  <- det(J_ll)
  if (!is.finite(detJ_hat) || !is.finite(detJ_ll) || detJ_ll <= 0) return(r)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)
  
  num_mat <- cbind(phi_diff, phi_lambda)
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

#########################################################
## 3. Helper to find one CI by r or r*
#########################################################

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
  
  lower <- uniroot(function(psi) r_fun(psi) - z_crit,
                   c(low, psi_hat))$root
  upper <- uniroot(function(psi) r_fun(psi) + z_crit,
                   c(psi_hat, up))$root
  
  c(lower, upper)
}

#########################################################
## 4. Case study with likelihood and bootstrap
#########################################################

run_alt_case_study <- function(stress_levels,
                               lifetimes_list,
                               use_stress,
                               p_grid     = seq(0.1, 0.3, length.out = 20),
                               conf_level = 0.95,
                               B_boot     = 1000,
                               log_scale  = TRUE,
                               seed       = 123) {
  
  t_vec <- unlist(lifetimes_list)
  x_vec <- rep(stress_levels, times = sapply(lifetimes_list, length))
  y_vec <- log(t_vec)
  
  mle <- alt_mle(y_vec, x_vec)
  if (mle$convergence != 0) stop("MLE did not converge")
  
  z_crit <- qnorm((1 + conf_level) / 2)
  s_p    <- -log(1 - p_grid)
  w_p    <- log(s_p)
  
  psi_hat <- mle$mu0_hat + mle$mu1_hat * use_stress + mle$sigma_hat * w_p
  
  K <- length(p_grid)
  psi_low_reg  <- psi_up_reg  <- numeric(K)
  psi_low_mod  <- psi_up_mod  <- numeric(K)
  
  for (k in seq_len(K)) {
    p  <- p_grid[k]
    ph <- psi_hat[k]
    
    r_fun  <- function(psi) alt_r_psi(psi,  y_vec, x_vec, use_stress, p, mle)
    rs_fun <- function(psi) alt_rstar_psi(psi, y_vec, x_vec, use_stress, p, mle)
    
    ci_reg <- find_ci_one(ph, r_fun,  z_crit, mle$sigma_hat)
    ci_mod <- find_ci_one(ph, rs_fun, z_crit, mle$sigma_hat)
    
    psi_low_reg[k] <- ci_reg[1]
    psi_up_reg[k]  <- ci_reg[2]
    psi_low_mod[k] <- ci_mod[1]
    psi_up_mod[k]  <- ci_mod[2]
  }
  
  ## parametric bootstrap
  if (!is.null(seed)) set.seed(seed)
  n_level  <- sapply(lifetimes_list, length)
  beta_hat <- 1 / mle$sigma_hat
  psi_boot <- matrix(NA_real_, nrow = B_boot, ncol = K)
  
  for (b in 1:B_boot) {
    y_b <- numeric(0)
    x_b <- numeric(0)
    for (j in seq_along(stress_levels)) {
      xj <- stress_levels[j]
      nj <- n_level[j]
      eta_j <- exp(mle$mu0_hat + mle$mu1_hat * xj)
      T_bj  <- rweibull(nj, shape = beta_hat, scale = eta_j)
      y_b   <- c(y_b, log(T_bj))
      x_b   <- c(x_b, rep(xj, nj))
    }
    mle_b <- alt_mle(y_b, x_b)
    if (mle_b$convergence != 0) next
    psi_boot[b, ] <- mle_b$mu0_hat + mle_b$mu1_hat * use_stress +
      mle_b$sigma_hat * w_p
  }
  
  psi_low_boot <- apply(psi_boot, 2, quantile,
                        probs = (1 - conf_level) / 2, na.rm = TRUE)
  psi_up_boot  <- apply(psi_boot, 2, quantile,
                        probs = 1 - (1 - conf_level) / 2, na.rm = TRUE)
  
  curve <- data.frame(
    p            = p_grid,
    psi_hat      = psi_hat,
    psi_low_reg  = psi_low_reg,
    psi_up_reg   = psi_up_reg,
    psi_low_mod  = psi_low_mod,
    psi_up_mod   = psi_up_mod,
    psi_low_boot = psi_low_boot,
    psi_up_boot  = psi_up_boot
  )
  
  if (!log_scale) {
    curve[, -1] <- exp(as.matrix(curve[, -1]))
  }
  
  list(curve = curve,
       mle   = mle,
       log_scale  = log_scale,
       use_stress = use_stress)
}

plot_alt_quantiles <- function(res) {
  curve <- res$curve
  p     <- curve$p
  ylab  <- if (res$log_scale) "log quantile at use stress"
  else "Lifetime at use stress"
  
  plot(p, curve$psi_hat, type = "l", lwd = 2,
       xlab = "probability p", ylab = ylab
  )
  
  lines(p, curve$psi_low_reg,  lty = 2, lwd = 2, col='blue')
  lines(p, curve$psi_up_reg,   lty = 2, lwd = 2, col='blue')
  
  lines(p, curve$psi_low_mod,  lty = 3, lwd = 2, col='red')
  lines(p, curve$psi_up_mod,   lty = 3, lwd = 2, col='red')
  
  lines(p, curve$psi_low_boot, lty = 4, lwd = 2, col='darkgrey')
  lines(p, curve$psi_up_boot,  lty = 4, lwd = 2, col='darkgrey')
  
  legend("bottomright",
         c("Point estimate",
           "Signed root",
           "Modified root ",
           "Bootstrap"),
         lwd = c(2, 2, 2, 2),
         lty = c(1, 2, 3, 4),
         col = c('black','blue','red','darkgrey'),
         bty = "n")
  
  #title(main = paste("Quantile function at use stress =", res$use_stress))
  title(main="")
}


weibull_prob_plot <- function(lifetimes_list, stress_levels,
                              main = "") {
  stopifnot(length(lifetimes_list) == length(stress_levels))
  groups <- lifetimes_list
  names(groups) <- stress_levels
  
  x_range <- range(log(unlist(groups)))
  y_range <- c(-4, 4)
  
  plot(NA, NA,
       xlim = c(2,8), 
       ylim = c(-4,4),
       xlab = "log lifetime",
       ylab = "Weibull probability scale",
       main = main)
  
  cols <- seq_along(groups)
  pchs <- seq_along(groups)
  i <- 1
  for (g in names(groups)) {
    t_g <- sort(groups[[g]])
    n_g <- length(t_g)
    p_g <- (1:n_g - 0.3) / (n_g + 0.4)
    x_g <- log(t_g)
    y_g <- log(-log(1 - p_g))
    points(x_g, y_g, col = cols[i], pch = pchs[i])
    fit_g <- lm(y_g ~ x_g)
    abline(fit_g, col = cols[i])
    i <- i + 1
  }
  legend("topleft",
         legend = paste("s =", names(groups)),
         col = cols, pch = pchs, bty = "n")
}




#########################################################
## Example ALT data 
#########################################################

stress_levels <- c(28.84, 
                   34.68, 
                   38.02) * 0.1

lifetimes_list <- list(
  c(1637, 1658, 2437, 1709, 1267, 1785),
  c(132, 96, 76, 122, 115, 87),
  c(43, 22, 39, 42, 41, 37)
)

use_stress <- 2


#########################################################
## Run on the example data
#########################################################

res_case <- run_alt_case_study(
  stress_levels   = stress_levels,
  lifetimes_list  = lifetimes_list,
  use_stress      = use_stress,
  p_grid          = seq(0.1, 0.4, length.out = 20),
  conf_level      = 0.95,
  B_boot          = 1000,
  log_scale       = TRUE,
  seed            = 123
)

par(mfrow = c(1, 2))

weibull_prob_plot(lifetimes_list, stress_levels,
                  main = "")
plot_alt_quantiles(res_case)


#########################################################
## 1. Likelihood for log Weibull ALT model
#########################################################

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
  
  z <- pmin(z, 700)              # avoid overflow
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
  z <- pmin(z, 700)
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

alt_mle <- function(y, x) {
  n   <- length(y)
  lm0 <- lm(y ~ x)
  mu0_start <- coef(lm0)[1]
  mu1_start <- coef(lm0)[2]
  res <- y - fitted(lm0)
  sigma_start <- sqrt(sum(res^2) / n)
  if (!is.finite(sigma_start) || sigma_start <= 0) sigma_start <- 1
  
  par_start <- c(mu0_start, mu1_start, log(sigma_start))
  
  obj <- function(par) {
    mu0   <- par[1]
    mu1   <- par[2]
    sigma <- exp(par[3])
    val   <- alt_loglik(c(mu0, mu1, sigma), y, x)
    if (!is.finite(val)) return(1e20)
    -val
  }
  
  opt <- optim(par_start, obj, method = "BFGS",
               control = list(maxit = 2000))
  
  list(
    mu0_hat   = opt$par[1],
    mu1_hat   = opt$par[2],
    sigma_hat = exp(opt$par[3]),
    ll_max    = -opt$value,
    convergence = opt$convergence
  )
}

#########################################################
## 2. Profile, r and r* for log quantile at use stress
#########################################################

alt_profile_psi <- function(psi0, y, x, x_use, p, mle_obj) {
  s_p <- -log(1 - p); w_p <- log(s_p)
  mu0_hat <- mle_obj$mu0_hat
  mu1_hat <- mle_obj$mu1_hat
  
  obj_lambda <- function(lambda) {
    mu0   <- lambda[1]
    mu1   <- lambda[2]
    sigma <- (psi0 - mu0 - mu1 * x_use) / w_p
    if (!is.finite(sigma) || sigma <= 0) return(1e20)
    val <- alt_loglik(c(mu0, mu1, sigma), y, x)
    if (!is.finite(val)) return(1e20)
    -val
  }
  
  opt <- try(
    optim(c(mu0_hat, mu1_hat), obj_lambda,
          method = "BFGS", control = list(maxit = 2000)),
    silent = TRUE
  )
  
  if (inherits(opt, "try-error")) {
    opt <- optim(c(mu0_hat, mu1_hat), obj_lambda,
                 method = "Nelder-Mead",
                 control = list(maxit = 4000))
  }
  
  mu0_0   <- opt$par[1]
  mu1_0   <- opt$par[2]
  sigma_0 <- (psi0 - mu0_0 - mu1_0 * x_use) / w_p
  
  if (!is.finite(sigma_0) || sigma_0 <= 0) {
    theta_0 <- rep(NA_real_, 3)
    ll0     <- -Inf
  } else {
    theta_0 <- c(mu0_0, mu1_0, sigma_0)
    ll0     <- alt_loglik(theta_0, y, x)
  }
  
  list(theta_0 = theta_0, ll0 = ll0)
}

alt_r_psi <- function(psi0, y, x, x_use, p, mle_obj) {
  mu0_hat   <- mle_obj$mu0_hat
  mu1_hat   <- mle_obj$mu1_hat
  sigma_hat <- mle_obj$sigma_hat
  ll_hat    <- mle_obj$ll_max
  
  s_p <- -log(1 - p); w_p <- log(s_p)
  psi_hat <- mu0_hat + mu1_hat * x_use + sigma_hat * w_p
  
  prof0 <- alt_profile_psi(psi0, y, x, x_use, p, mle_obj)
  ll0   <- prof0$ll0
  
  if (!is.finite(ll0)) {
    return(sign(psi_hat - psi0) * 1e6)
  }
  
  sign(psi_hat - psi0) * sqrt(pmax(0, 2 * (ll_hat - ll0)))
}

alt_rstar_psi <- function(psi0, y, x, x_use, p, mle_obj) {
  mu0_hat   <- mle_obj$mu0_hat
  mu1_hat   <- mle_obj$mu1_hat
  sigma_hat <- mle_obj$sigma_hat
  ll_hat    <- mle_obj$ll_max
  
  s_p <- -log(1 - p); w_p <- log(s_p)
  psi_hat <- mu0_hat + mu1_hat * x_use + sigma_hat * w_p
  
  prof0  <- alt_profile_psi(psi0, y, x, x_use, p, mle_obj)
  theta0 <- prof0$theta_0
  ll0    <- prof0$ll0
  
  r <- alt_r_psi(psi0, y, x, x_use, p, mle_obj)
  if (!all(is.finite(theta0)) || !is.finite(ll0) ||
      !is.finite(r) || abs(r) < 1e-10) {
    return(r)
  }
  
  theta_hat <- c(mu0_hat, mu1_hat, sigma_hat)
  
  phi_hat <- alt_phi(theta_hat, y, x)
  phi_0   <- alt_phi(theta0,   y, x)
  if (any(!is.finite(phi_hat)) || any(!is.finite(phi_0))) return(r)
  
  phi_diff <- phi_hat - phi_0
  
  Phi_hat <- alt_phi_jac(theta_hat, y, x)
  Phi_0   <- alt_phi_jac(theta0,   y, x)
  if (any(!is.finite(Phi_hat)) || any(!is.finite(Phi_0))) return(r)
  
  T_mat <- matrix(
    c(0,      1,            0,
      0,      0,            1,
      1/w_p, -1/w_p, -x_use / w_p),
    nrow = 3, byrow = TRUE
  )
  
  Phi_hat_pl <- Phi_hat %*% T_mat
  T_lambda   <- T_mat[, 2:3, drop = FALSE]
  phi_lambda <- Phi_0 %*% T_lambda
  
  J_hat_theta <- alt_J(theta_hat, y, x)
  J_0_theta   <- alt_J(theta0,   y, x)
  if (any(!is.finite(J_hat_theta)) || any(!is.finite(J_0_theta))) return(r)
  
  J_hat_pl <- t(T_mat) %*% J_hat_theta %*% T_mat
  J_0_pl   <- t(T_mat) %*% J_0_theta   %*% T_mat
  J_ll     <- J_0_pl[2:3, 2:3, drop = FALSE]
  
  detJ_hat <- det(J_hat_pl)
  detJ_ll  <- det(J_ll)
  if (!is.finite(detJ_hat) || !is.finite(detJ_ll) || detJ_ll <= 0) return(r)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)
  
  num_mat <- cbind(phi_diff, phi_lambda)
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

#########################################################
## 3. Helper to find one CI by r or r*
#########################################################

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
  
  lower <- uniroot(function(psi) r_fun(psi) - z_crit,
                   c(low, psi_hat))$root
  upper <- uniroot(function(psi) r_fun(psi) + z_crit,
                   c(psi_hat, up))$root
  
  c(lower, upper)
}

#########################################################
## 4. Run analysis for a single dataset (p single value possible)
#########################################################

run_alt_case_study <- function(stress_levels,
                               lifetimes_list,
                               use_stress,
                               p_grid,
                               conf_level = 0.95,
                               B_boot     = 1000,
                               log_scale  = TRUE,
                               seed       = NULL) {
  
  t_vec <- unlist(lifetimes_list)
  x_vec <- rep(stress_levels, times = sapply(lifetimes_list, length))
  y_vec <- log(t_vec)
  
  mle <- alt_mle(y_vec, x_vec)
  if (mle$convergence != 0) stop("MLE did not converge")
  
  z_crit <- qnorm((1 + conf_level) / 2)
  s_p    <- -log(1 - p_grid)
  w_p    <- log(s_p)
  
  psi_hat <- mle$mu0_hat + mle$mu1_hat * use_stress + mle$sigma_hat * w_p
  
  K <- length(p_grid)
  psi_low_reg  <- psi_up_reg  <- numeric(K)
  psi_low_mod  <- psi_up_mod  <- numeric(K)
  
  for (k in seq_len(K)) {
    p  <- p_grid[k]
    ph <- psi_hat[k]
    
    r_fun  <- function(psi) alt_r_psi(psi,  y_vec, x_vec, use_stress, p, mle)
    rs_fun <- function(psi) alt_rstar_psi(psi, y_vec, x_vec, use_stress, p, mle)
    
    ci_reg <- find_ci_one(ph, r_fun,  z_crit, mle$sigma_hat)
    ci_mod <- find_ci_one(ph, rs_fun, z_crit, mle$sigma_hat)
    
    psi_low_reg[k] <- ci_reg[1]
    psi_up_reg[k]  <- ci_reg[2]
    psi_low_mod[k] <- ci_mod[1]
    psi_up_mod[k]  <- ci_mod[2]
  }
  
  ## parametric bootstrap
  if (!is.null(seed)) set.seed(seed)
  n_level  <- sapply(lifetimes_list, length)
  beta_hat <- 1 / mle$sigma_hat
  psi_boot <- matrix(NA_real_, nrow = B_boot, ncol = K)
  
  for (b in 1:B_boot) {
    y_b <- numeric(0)
    x_b <- numeric(0)
    for (j in seq_along(stress_levels)) {
      xj <- stress_levels[j]
      nj <- n_level[j]
      eta_j <- exp(mle$mu0_hat + mle$mu1_hat * xj)
      T_bj  <- rweibull(nj, shape = beta_hat, scale = eta_j)
      y_b   <- c(y_b, log(T_bj))
      x_b   <- c(x_b, rep(xj, nj))
    }
    mle_b <- alt_mle(y_b, x_b)
    if (mle_b$convergence != 0) next
    psi_boot[b, ] <- mle_b$mu0_hat + mle_b$mu1_hat * use_stress +
      mle_b$sigma_hat * w_p
  }
  
  psi_low_boot <- apply(psi_boot, 2, quantile,
                        probs = (1 - conf_level) / 2, na.rm = TRUE)
  psi_up_boot  <- apply(psi_boot, 2, quantile,
                        probs = 1 - (1 - conf_level) / 2, na.rm = TRUE)
  
  curve <- data.frame(
    p            = p_grid,
    psi_hat      = psi_hat,
    psi_low_reg  = psi_low_reg,
    psi_up_reg   = psi_up_reg,
    psi_low_mod  = psi_low_mod,
    psi_up_mod   = psi_up_mod,
    psi_low_boot = psi_low_boot,
    psi_up_boot  = psi_up_boot
  )
  
  if (!log_scale) {
    curve[, -1] <- exp(as.matrix(curve[, -1]))
  }
  
  list(curve = curve,
       mle   = mle,
       log_scale  = log_scale,
       use_stress = use_stress)
}

simulate_grouped_original <- function(mu0, mu1, sigma, stress_levels, n_per) {
  lifetimes_list <- vector("list", length(stress_levels))
  for (i in seq_along(stress_levels)) {
    si <- stress_levels[i]
    ni <- n_per[i]
    mu_i <- mu0 + mu1 * si
    z <- -log(-log(runif(ni)))          # 标准 Gumbel
    lifetimes_list[[i]] <- exp(mu_i + sigma * z)
  }
  lifetimes_list
}

#########################################################
## 5. Simulation study: MSE, coverage, interval length
#########################################################

# ---- True parameters from original data ----
original_stress <- c(28.84, 34.68, 38.02) * 0.1
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

use_stress <- 2
p_target   <- 0.1
conf_level <- 0.95
z_crit     <- qnorm((1 + conf_level) / 2)

# True log quantile
s_p_target <- -log(1 - p_target)
w_p_target <- log(s_p_target)
psi_true   <- mu0_true + mu1_true * use_stress + sigma_true * w_p_target

# Simulation settings
n_sim <- 1000
B_boot_sim <- 1   # number of bootstrap replications inside each simulation

n_per_stress <- c(6, 6, 6)   # equal to original data

# Storage
psi_hat_vec <- numeric(n_sim)
cover_r     <- logical(n_sim)
cover_rstar <- logical(n_sim)
cover_boot  <- logical(n_sim)
len_r       <- numeric(n_sim)
len_rstar   <- numeric(n_sim)
len_boot    <- numeric(n_sim)

set.seed(123)  # for reproducibility

for (sim in 1:n_sim) {
  if (sim %% 100 == 0) cat("Simulation", sim, "/", n_sim, "\n")
  
  # Generate data under true model
  sim_lifetimes <- simulate_grouped_original(mu0_true, mu1_true, sigma_true,
                                             original_stress, n_per_stress)
  
  # Run analysis for p_target only
  res <- tryCatch(
    run_alt_case_study(
      stress_levels  = original_stress,
      lifetimes_list = sim_lifetimes,
      use_stress     = use_stress,
      p_grid         = p_target,   # single value
      conf_level     = conf_level,
      B_boot         = B_boot_sim,
      log_scale      = TRUE,
      seed           = NULL
    ),
    error = function(e) NULL
  )
  
  if (is.null(res)) {
    # if analysis fails, set NAs
    psi_hat_vec[sim] <- NA
    cover_r[sim]     <- NA
    cover_rstar[sim] <- NA
    cover_boot[sim]  <- NA
    len_r[sim]       <- NA
    len_rstar[sim]   <- NA
    len_boot[sim]    <- NA
    next
  }
  
  # Extract point estimate and CIs (log scale)
  curve <- res$curve
  psi_hat <- curve$psi_hat[1]
  
  low_r  <- curve$psi_low_reg[1]
  up_r   <- curve$psi_up_reg[1]
  low_rs <- curve$psi_low_mod[1]
  up_rs  <- curve$psi_up_mod[1]
  low_bt <- curve$psi_low_boot[1]
  up_bt  <- curve$psi_up_boot[1]
  
  psi_hat_vec[sim] <- psi_hat
  cover_r[sim]     <- (psi_true >= low_r  & psi_true <= up_r)
  cover_rstar[sim] <- (psi_true >= low_rs & psi_true <= up_rs)
  cover_boot[sim]  <- (psi_true >= low_bt & psi_true <= up_bt)
  len_r[sim]       <- up_r - low_r
  len_rstar[sim]   <- up_rs - low_rs
  len_boot[sim]    <- up_bt - low_bt
}

# ---- Summary results ----
# MSE of point estimator
mse_psi <- mean((psi_hat_vec - psi_true)^2, na.rm = TRUE)
cat("\nMSE of log-quantile estimator (p = 0.1):", mse_psi, "\n")

# Coverage rates
cov_r     <- mean(cover_r,     na.rm = TRUE) * 100
cov_rstar <- mean(cover_rstar, na.rm = TRUE) * 100
cov_boot  <- mean(cover_boot,  na.rm = TRUE) * 100
cat("\nCoverage rates (%):\n")
cat("  Signed root (r) :", cov_r, "\n")
cat("  Modified root (r*):", cov_rstar, "\n")
cat("  Bootstrap       :", cov_boot, "\n")

# Interval lengths boxplot
boxplot(len_r, len_rstar, len_boot,
        names = c("r", "r*", "Bootstrap"),
        ylab = "Interval length",
        main = paste0("Interval lengths for p = ", p_target))
abline(h = 0, lty = 2)

df_len <- data.frame(
  Method = rep(c("r", "r*", "Bootstrap"), each = n_sim),
  Length = c(len_r, len_rstar, len_boot)
)

ggplot(df_len, aes(x = Method, y = Length, fill = Method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("r" = "#66c2a5", "r*" = "#fc8d62", "Bootstrap" = "#8da0cb")) +
  labs(y = "Interval length") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

