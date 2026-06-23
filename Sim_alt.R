# ALT length-revise
rm(list = ls())

############################################################
## ALT simulation - Weibull regression (log-Weibull form)
## Model: Y_i = log T_i = mu0 + mu1 * x_i + sigma * W_i
##        W_i ~ EV(0, 1), independent
##
## Parameter: theta = (mu0, mu1, sigma)
## Quantile of interest at use stress x_use:
##   psi = zeta_p(x_use) = mu0 + mu1 * x_use + sigma * w_p
## where s_p = -log(1 - p), w_p = log(s_p)
##
## Methods: signed root r(psi), modified signed root r*(psi)
############################################################

########################################
## 1. Basic likelihood building blocks
########################################

## log-likelihood in theta = (mu0, mu1, sigma)
alt_loglik <- function(theta, y, x) {
  mu0   <- theta[1]
  mu1   <- theta[2]
  sigma <- theta[3]
  n     <- length(y)
  if (sigma <= 0) return(-1e30)  ## penalty for invalid sigma
  
  z <- (y - mu0 - mu1 * x) / sigma
  r <- exp(z)
  
  ell <- sum(z - r) - n * log(sigma)
  ell
}

## canonical parameter phi(theta) - 3 dimensional
alt_phi <- function(theta, y, x) {
  mu0   <- theta[1]
  mu1   <- theta[2]
  sigma <- theta[3]
  if (sigma <= 0) return(rep(NA_real_, 3))
  
  z <- (y - mu0 - mu1 * x) / sigma
  r <- exp(z)
  
  S0 <- sum(1 - r)
  S1 <- sum(x * (1 - r))
  S2 <- sum(z * (1 - r))
  
  phi1 <- S0 / sigma
  phi2 <- S1 / sigma
  phi3 <- S2 / sigma
  
  c(phi1, phi2, phi3)
}

## numerical Jacobian of phi wrt theta using central differences
alt_phi_jac <- function(theta, y, x, h = 1e-5) {
  k <- length(theta)  ## 3
  p <- length(theta)
  base <- alt_phi(theta, y, x)
  if (any(!is.finite(base))) {
    return(matrix(NA_real_, nrow = 3, ncol = 3))
  }
  J <- matrix(NA_real_, nrow = k, ncol = p)
  
  for (j in 1:p) {
    step <- h * max(1, abs(theta[j]))
    th_plus  <- theta
    th_minus <- theta
    
    th_plus[j]  <- theta[j] + step
    th_minus[j] <- theta[j] - step
    
    ## ensure sigma stays positive
    if (j == 3 && th_minus[3] <= 0) {
      ## forward difference
      f0 <- base
      fp <- alt_phi(th_plus, y, x)
      J[, j] <- (fp - f0) / step
    } else {
      fp <- alt_phi(th_plus, y, x)
      fm <- alt_phi(th_minus, y, x)
      J[, j] <- (fp - fm) / (2 * step)
    }
  }
  
  J
}

## numerical observed information J(theta) = -Hessian of log-likelihood
alt_J <- function(theta, y, x, eps = 1e-4) {
  p <- length(theta)  ## 3
  J <- matrix(NA_real_, nrow = p, ncol = p)
  f0 <- alt_loglik(theta, y, x)
  
  for (k in 1:p) {
    hk <- eps * max(1, abs(theta[k]))
    th_pk <- theta
    th_mk <- theta
    th_pk[k] <- theta[k] + hk
    th_mk[k] <- theta[k] - hk
    
    f_pk <- alt_loglik(th_pk, y, x)
    f_mk <- alt_loglik(th_mk, y, x)
    
    H_kk <- (f_pk - 2 * f0 + f_mk) / (hk^2)
    J[k, k] <- -H_kk
    
    if (k < p) {
      for (l in (k + 1):p) {
        hl <- eps * max(1, abs(theta[l]))
        th_pp <- theta
        th_pm <- theta
        th_mp <- theta
        th_mm <- theta
        
        th_pp[k] <- theta[k] + hk
        th_pp[l] <- theta[l] + hl
        
        th_pm[k] <- theta[k] + hk
        th_pm[l] <- theta[l] - hl
        
        th_mp[k] <- theta[k] - hk
        th_mp[l] <- theta[l] + hl
        
        th_mm[k] <- theta[k] - hk
        th_mm[l] <- theta[l] - hl
        
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

########################################
## 2. Unconstrained MLE
########################################

## We optimize over (mu0, mu1, log_sigma) for stability
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

########################################
## 3. Profile in psi (log quantile at use stress)
########################################

## psi = zeta_p(x_use) = mu0 + mu1 * x_use + sigma * w_p
## We treat (psi, lambda1, lambda2) = (psi, mu0, mu1)
## and eliminate sigma = (psi - mu0 - mu1 * x_use) / w_p

## For a given psi0, profile log-likelihood:
##   maximize ell(mu0, mu1, sigma) over (mu0, mu1) with sigma implied.

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

## signed root r(psi0)
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

########################################
## 4. Modified signed root r*(psi)
########################################

alt_rstar_psi <- function(psi0, y, x, x_use, p, mle_obj = NULL) {
  if (is.null(mle_obj)) mle_obj <- alt_mle(y, x)
  mu0_hat   <- mle_obj$mu0_hat
  mu1_hat   <- mle_obj$mu1_hat
  sigma_hat <- mle_obj$sigma_hat
  ll_hat    <- mle_obj$ll_max
  
  s_p <- -log(1 - p)
  w_p <- log(s_p)
  psi_hat <- mu0_hat + mu1_hat * x_use + sigma_hat * w_p
  
  ## signed root r
  prof0 <- alt_profile_psi(psi0, y, x, x_use, p, mle_obj)
  theta_0  <- prof0$theta_0
  lambda_0 <- prof0$lambda_0
  ll0      <- prof0$ll0
  
  r <- sign(psi_hat - psi0) * sqrt(max(0, 2 * (ll_hat - ll0)))
  if (!is.finite(r) || abs(r) < 1e-10) return(0)
  
  theta_hat <- c(mu0_hat, mu1_hat, sigma_hat)
  
  ## phi and phi_theta at hat and constrained
  phi_hat <- alt_phi(theta_hat, y, x)
  phi_0   <- alt_phi(theta_0,   y, x)
  if (any(!is.finite(phi_hat)) || any(!is.finite(phi_0))) return(r)
  
  phi_diff <- phi_hat - phi_0
  
  Phi_hat <- alt_phi_jac(theta_hat, y, x)
  Phi_0   <- alt_phi_jac(theta_0,   y, x)
  if (any(!is.finite(Phi_hat)) || any(!is.finite(Phi_0))) return(r)
  
  ## transformation matrix T = d theta / d(psi, lambda1, lambda2)'
  ##   mu0   = lambda1
  ##   mu1   = lambda2
  ##   sigma = (psi - lambda1 - lambda2 * x_use) / w_p
  ##
  ## columns: derivatives wrt psi, lambda1, lambda2
  T_mat <- matrix(
    c(0,      1,            0,
      0,      0,            1,
      1/w_p, -1/w_p, -x_use / w_p),
    nrow = 3, byrow = TRUE
  )
  
  ## phi_theta in (psi, lambda) coordinates at hat
  Phi_hat_pl <- Phi_hat %*% T_mat
  
  ## derivative wrt lambda = (lambda1, lambda2) at psi0
  T_lambda <- T_mat[, 2:3, drop = FALSE]  ## 3 x 2
  phi_lambda <- Phi_0 %*% T_lambda        ## 3 x 2
  
  ## observed information in theta and transformed to (psi, lambda)
  J_hat_theta <- alt_J(theta_hat, y, x)
  J_0_theta   <- alt_J(theta_0,   y, x)
  if (any(!is.finite(J_hat_theta)) || any(!is.finite(J_0_theta))) return(r)
  
  J_hat_pl <- t(T_mat) %*% J_hat_theta %*% T_mat
  J_0_pl   <- t(T_mat) %*% J_0_theta   %*% T_mat
  
  ## nuisance block is last 2 x 2
  J_ll <- J_0_pl[2:3, 2:3, drop = FALSE]
  
  detJ_hat <- det(J_hat_pl)
  detJ_ll  <- det(J_ll)
  
  if (!is.finite(detJ_hat) || !is.finite(detJ_ll) || detJ_ll <= 0) return(r)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)
  
  ## determinant in numerator: [phi_hat - phi_0, phi_lambda]
  num_mat <- cbind(phi_diff, phi_lambda)  ## 3 x 3
  det_num <- det(num_mat)
  det_den <- det(Phi_hat_pl)
  
  if (!is.finite(det_num) || !is.finite(det_den)) return(r)
  
  q_val <- (det_num / det_den) * sqrt(detJ_hat / detJ_ll)
  if (!is.finite(q_val)) return(r)
  
  ## sign adjustment
  if (q_val * r <= 0) q_val <- abs(q_val) * sign(r)
  
  rstar <- r + (1 / r) * log(q_val / r)
  if (!is.finite(rstar)) rstar <- r
  rstar
}

########################################
## 5. Simulation for ALT quantile
########################################

## Simulate one scenario:
##  - true parameters (mu0_true, mu1_true, sigma_true)
##  - stress levels x_levels and sample sizes n_per_level
##  - use stress x_use
##  - quantile p_quant
##  - conf_level (e.g. 0.95)
##  - n_rep replications
## Output: coverage of r and r* for psi

sim_alt_quantile <- function(
    mu0_true,
    mu1_true,
    sigma_true,
    x_levels      = c(-1, 0, 1),
    n_per_level   = c(10, 10, 10),
    x_use         = 0,
    p_quant       = 0.1,
    conf_level    = 0.95,
    n_rep         = 1000,
    seed          = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  if (length(x_levels) != length(n_per_level)) {
    stop("x_levels and n_per_level must have same length")
  }
  
  ## true shape and scale at each stress if you care
  beta_true <- 1 / sigma_true
  s_p <- -log(1 - p_quant)
  w_p <- log(s_p)
  
  ## true psi at use stress
  psi_true <- mu0_true + mu1_true * x_use + sigma_true * w_p
  
  z_crit <- qnorm((1 + conf_level) / 2)
  
  count_r     <- 0L
  count_rstar <- 0L
  
  for (rep in seq_len(n_rep)) {
    ## simulate data
    y_list <- list()
    x_list <- list()
    for (j in seq_along(x_levels)) {
      xj <- x_levels[j]
      nj <- n_per_level[j]
      eta_j <- exp(mu0_true + mu1_true * xj)
      Tj <- rweibull(nj, shape = beta_true, scale = eta_j)
      yj <- log(Tj)
      y_list[[j]] <- yj
      x_list[[j]] <- rep(xj, nj)
    }
    y <- unlist(y_list)
    x <- unlist(x_list)
    
    ## unconstrained MLE
    mle_obj <- alt_mle(y, x)
    if (mle_obj$convergence != 0) {
      next
    }
    
    ## compute r and r* at true psi
    r_val     <- alt_r_psi(psi_true, y, x, x_use, p_quant, mle_obj)
    rstar_val <- alt_rstar_psi(psi_true, y, x, x_use, p_quant, mle_obj)
    
    if (is.finite(r_val) && abs(r_val) <= z_crit) {
      count_r <- count_r + 1L
    }
    if (is.finite(rstar_val) && abs(rstar_val) <= z_crit) {
      count_rstar <- count_rstar + 1L
    }
    
    if (rep %% 100 == 0) {
      cat("ALT sim - rep", rep, "of", n_rep, "\r")
      flush.console()
    }
  }
  cat("\n")
  
  list(
    settings = list(
      mu0_true    = mu0_true,
      mu1_true    = mu1_true,
      sigma_true  = sigma_true,
      beta_true   = beta_true,
      x_levels    = x_levels,
      n_per_level = n_per_level,
      x_use       = x_use,
      p_quant     = p_quant,
      conf_level  = conf_level,
      n_rep       = n_rep
    ),
    coverage = c(
      signed_root   = count_r     / n_rep,
      modified_root = count_rstar / n_rep
    )
  )
}


########################################
## A. Profile, r and r* for parameters
##    psi in {mu0, mu1, sigma}
########################################

## psi_index: 1 = mu0, 2 = mu1, 3 = sigma

alt_profile_param <- function(psi_index, psi0, y, x, mle_obj) {
  mu0_hat   <- mle_obj$mu0_hat
  mu1_hat   <- mle_obj$mu1_hat
  sigma_hat <- mle_obj$sigma_hat
  
  if (psi_index == 1) {
    ## psi = mu0, nuisance = (mu1, sigma)
    mu0_fixed <- psi0
    par_start <- c(mu1_hat, log(sigma_hat))
    obj <- function(par) {
      mu1   <- par[1]
      sigma <- exp(par[2])
      -alt_loglik(c(mu0_fixed, mu1, sigma), y, x)
    }
    opt <- optim(par_start, obj, method = "BFGS",
                 control = list(maxit = 2000))
    mu1_0   <- opt$par[1]
    sigma_0 <- exp(opt$par[2])
    theta_0 <- c(mu0_fixed, mu1_0, sigma_0)
    
  } else if (psi_index == 2) {
    ## psi = mu1, nuisance = (mu0, sigma)
    mu1_fixed <- psi0
    par_start <- c(mu0_hat, log(sigma_hat))
    obj <- function(par) {
      mu0   <- par[1]
      sigma <- exp(par[2])
      -alt_loglik(c(mu0, mu1_fixed, sigma), y, x)
    }
    opt <- optim(par_start, obj, method = "BFGS",
                 control = list(maxit = 2000))
    mu0_0   <- opt$par[1]
    sigma_0 <- exp(opt$par[2])
    theta_0 <- c(mu0_0, mu1_fixed, sigma_0)
    
  } else if (psi_index == 3) {
    ## psi = sigma, nuisance = (mu0, mu1)
    sigma_fixed <- psi0
    if (sigma_fixed <= 0) {
      return(list(theta_0 = rep(NA_real_, 3),
                  ll0 = -Inf, convergence = 1))
    }
    par_start <- c(mu0_hat, mu1_hat)
    obj <- function(par) {
      mu0 <- par[1]
      mu1 <- par[2]
      -alt_loglik(c(mu0, mu1, sigma_fixed), y, x)
    }
    opt <- optim(par_start, obj, method = "BFGS",
                 control = list(maxit = 2000))
    mu0_0   <- opt$par[1]
    mu1_0   <- opt$par[2]
    theta_0 <- c(mu0_0, mu1_0, sigma_fixed)
    
  } else {
    stop("psi_index must be 1 (mu0), 2 (mu1) or 3 (sigma)")
  }
  
  ll0 <- alt_loglik(theta_0, y, x)
  
  list(theta_0 = theta_0,
       ll0     = ll0,
       convergence = opt$convergence)
}

## signed root r(psi0) for psi in {mu0, mu1, sigma}
alt_r_param <- function(psi_index, psi0, y, x, mle_obj = NULL) {
  if (is.null(mle_obj)) mle_obj <- alt_mle(y, x)
  if (mle_obj$convergence != 0) return(NA_real_)
  
  theta_hat <- c(mle_obj$mu0_hat, mle_obj$mu1_hat, mle_obj$sigma_hat)
  ll_hat    <- mle_obj$ll_max
  psi_hat   <- theta_hat[psi_index]
  
  prof0 <- alt_profile_param(psi_index, psi0, y, x, mle_obj)
  if (prof0$convergence != 0 || any(is.na(prof0$theta_0))) {
    return(NA_real_)
  }
  ll0 <- prof0$ll0
  
  r <- sign(psi_hat - psi0) * sqrt(max(0, 2 * (ll_hat - ll0)))
  r
}

## permutation matrix T for each psi_index
## theta = (mu0, mu1, sigma)
## new params = (psi, lambda1, lambda2)

alt_T_param <- function(psi_index) {
  if (psi_index == 1) {
    ## psi = mu0, lambda = (mu1, sigma)
    ## (psi, lambda1, lambda2) = (mu0, mu1, sigma)
    diag(3)
  } else if (psi_index == 2) {
    ## psi = mu1, lambda = (mu0, sigma)
    ## mu0 = lambda1, mu1 = psi, sigma = lambda2
    matrix(
      c(0, 1, 0,
        1, 0, 0,
        0, 0, 1),
      nrow = 3, byrow = TRUE
    )
  } else if (psi_index == 3) {
    ## psi = sigma, lambda = (mu0, mu1)
    ## mu0 = lambda1, mu1 = lambda2, sigma = psi
    matrix(
      c(0, 1, 0,
        0, 0, 1,
        1, 0, 0),
      nrow = 3, byrow = TRUE
    )
  } else {
    stop("psi_index must be 1, 2 or 3")
  }
}

## modified signed root r*(psi0) for psi in {mu0, mu1, sigma}
alt_rstar_param <- function(psi_index, psi0, y, x, mle_obj = NULL) {
  if (is.null(mle_obj)) mle_obj <- alt_mle(y, x)
  if (mle_obj$convergence != 0) return(NA_real_)
  
  theta_hat <- c(mle_obj$mu0_hat, mle_obj$mu1_hat, mle_obj$sigma_hat)
  ll_hat    <- mle_obj$ll_max
  psi_hat   <- theta_hat[psi_index]
  
  prof0 <- alt_profile_param(psi_index, psi0, y, x, mle_obj)
  theta_0 <- prof0$theta_0
  if (prof0$convergence != 0 || any(is.na(theta_0))) return(NA_real_)
  ll0 <- prof0$ll0
  
  ## signed root
  r <- sign(psi_hat - psi0) * sqrt(max(0, 2 * (ll_hat - ll0)))
  if (!is.finite(r) || abs(r) < 1e-10) return(r)
  
  ## phi and Jacobians
  phi_hat <- alt_phi(theta_hat, y, x)
  phi_0   <- alt_phi(theta_0,   y, x)
  if (any(!is.finite(phi_hat)) || any(!is.finite(phi_0))) return(r)
  
  phi_diff <- phi_hat - phi_0
  
  Phi_hat <- alt_phi_jac(theta_hat, y, x)
  Phi_0   <- alt_phi_jac(theta_0,   y, x)
  if (any(!is.finite(Phi_hat)) || any(!is.finite(Phi_0))) return(r)
  
  ## transformation T for this psi_index
  T_mat <- alt_T_param(psi_index)
  
  ## phi_theta in (psi, lambda) coords at hat
  Phi_hat_pl <- Phi_hat %*% T_mat
  
  ## derivative wrt lambda at psi0
  T_lambda   <- T_mat[, 2:3, drop = FALSE]  ## 3 x 2
  phi_lambda <- Phi_0 %*% T_lambda          ## 3 x 2
  
  ## observed information
  J_hat_theta <- alt_J(theta_hat, y, x)
  J_0_theta   <- alt_J(theta_0,   y, x)
  if (any(!is.finite(J_hat_theta)) || any(!is.finite(J_0_theta))) return(r)
  
  J_hat_pl <- t(T_mat) %*% J_hat_theta %*% T_mat
  J_0_pl   <- t(T_mat) %*% J_0_theta   %*% T_mat
  
  ## nuisance block: last 2x2
  J_ll <- J_0_pl[2:3, 2:3, drop = FALSE]
  
  detJ_hat <- det(J_hat_pl)
  detJ_ll  <- det(J_ll)
  if (!is.finite(detJ_hat) || !is.finite(detJ_ll) || detJ_ll <= 0) return(r)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)
  
  ## determinant in numerator: [phi_hat - phi_0, phi_lambda]
  num_mat <- cbind(phi_diff, phi_lambda)  ## 3 x 3
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


########################################
## B. ALT simulation for parameters
##    and 0.1-quantile at use stress
########################################

sim_alt_params_quantile <- function(
    beta_vec,        ## vector of true beta = 1/sigma
    mu0_vec,         ## matching true mu0
    mu1_vec,         ## matching true mu1
    x_levels,
    n_per_level = c(3, 4, 3),   ## total n = 10 (similar to Type II)
    x_use       = 0,
    p_quant     = 0.1,
    conf_level  = 0.95,
    n_rep       = 1000,
    seed        = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  if (length(beta_vec) != length(mu0_vec) ||
      length(beta_vec) != length(mu1_vec)) {
    stop("beta_vec, mu0_vec and mu1_vec must have the same length")
  }
  if (length(x_levels) != length(n_per_level)) {
    stop("x_levels and n_per_level must have the same length")
  }
  
  s_p <- -log(1 - p_quant)
  w_p <- log(s_p)
  z_crit <- qnorm((1 + conf_level) / 2)
  
  scenario_results <- list()
  
  for (j in seq_along(beta_vec)) {
    beta_true  <- beta_vec[j]
    sigma_true <- 1 / beta_true
    mu0_true   <- mu0_vec[j]
    mu1_true   <- mu1_vec[j]
    
    ## true psi (log 0.1-quantile) at use stress
    psi_true <- mu0_true + mu1_true * x_use + sigma_true * w_p
    
    ## counters
    count_valid <- 0L
    
    count_mu0_r     <- 0L
    count_mu0_rstar <- 0L
    count_mu1_r     <- 0L
    count_mu1_rstar <- 0L
    count_sig_r     <- 0L
    count_sig_rstar <- 0L
    count_q_r       <- 0L
    count_q_rstar   <- 0L
    
    for (rep in seq_len(n_rep)) {
      ## simulate data
      y_list <- list()
      x_list <- list()
      for (k in seq_along(x_levels)) {
        s_k <- x_levels[k]
        n_k <- n_per_level[k]
        eta_k <- exp(mu0_true + mu1_true * s_k)
        T_k <- rweibull(n_k, shape = beta_true, scale = eta_k)
        y_k <- log(T_k)
        y_list[[k]] <- y_k
        x_list[[k]] <- rep(s_k, n_k)
      }
      y <- unlist(y_list)
      x <- unlist(x_list)
      
      ## unconstrained MLE
      mle_obj <- alt_mle(y, x)
      if (mle_obj$convergence != 0) next  ## skip failed fit
      
      count_valid <- count_valid + 1L
      
      ## true values of parameters
      psi0_mu0   <- mu0_true
      psi0_mu1   <- mu1_true
      psi0_sigma <- sigma_true
      
      ## r and r* for mu0, mu1, sigma
      r_mu0     <- alt_r_param(1, psi0_mu0,   y, x, mle_obj)
      rstar_mu0 <- alt_rstar_param(1, psi0_mu0, y, x, mle_obj)
      
      r_mu1     <- alt_r_param(2, psi0_mu1,   y, x, mle_obj)
      rstar_mu1 <- alt_rstar_param(2, psi0_mu1, y, x, mle_obj)
      
      r_sig     <- alt_r_param(3, psi0_sigma, y, x, mle_obj)
      rstar_sig <- alt_rstar_param(3, psi0_sigma, y, x, mle_obj)
      
      ## r and r* for log 0.1-quantile at use stress
      r_q       <- alt_r_psi(psi_true, y, x, x_use, p_quant, mle_obj)
      rstar_q   <- alt_rstar_psi(psi_true, y, x, x_use, p_quant, mle_obj)
      
      if (is.finite(r_mu0) && abs(r_mu0) <= z_crit) {
        count_mu0_r <- count_mu0_r + 1L
      }
      if (is.finite(rstar_mu0) && abs(rstar_mu0) <= z_crit) {
        count_mu0_rstar <- count_mu0_rstar + 1L
      }
      
      if (is.finite(r_mu1) && abs(r_mu1) <= z_crit) {
        count_mu1_r <- count_mu1_r + 1L
      }
      if (is.finite(rstar_mu1) && abs(rstar_mu1) <= z_crit) {
        count_mu1_rstar <- count_mu1_rstar + 1L
      }
      
      if (is.finite(r_sig) && abs(r_sig) <= z_crit) {
        count_sig_r <- count_sig_r + 1L
      }
      if (is.finite(rstar_sig) && abs(rstar_sig) <= z_crit) {
        count_sig_rstar <- count_sig_rstar + 1L
      }
      
      if (is.finite(r_q) && abs(r_q) <= z_crit) {
        count_q_r <- count_q_r + 1L
      }
      if (is.finite(rstar_q) && abs(rstar_q) <= z_crit) {
        count_q_rstar <- count_q_rstar + 1L
      }
      
      if (rep %% 100 == 0) {
        cat("ALT scenario", j, "- rep", rep, "of", n_rep, "\r")
        flush.console()
      }
    }
    cat("\n")
    
    if (count_valid == 0L) {
      warning("No valid MLEs in scenario ", j)
      next
    }
    
    ## assemble results for scenario j
    df_j <- data.frame(
      scenario = j,
      beta_true = beta_true,
      mu0_true  = mu0_true,
      mu1_true  = mu1_true,
      sigma_true = sigma_true,
      param = rep(c("mu0", "mu1", "sigma", "log_quantile_zeta"),
                  each = 2),
      method = rep(c("signed_root", "modified_root"), times = 4),
      coverage = c(
        count_mu0_r     / count_valid,
        count_mu0_rstar / count_valid,
        count_mu1_r     / count_valid,
        count_mu1_rstar / count_valid,
        count_sig_r     / count_valid,
        count_sig_rstar / count_valid,
        count_q_r       / count_valid,
        count_q_rstar   / count_valid
      )
    )
    scenario_results[[j]] <- df_j
  }
  
  do.call(rbind, scenario_results)
}


###############################################################################
## Additional interval-length calculations and simulations
###############################################################################

## ============================================================================
## 0. Stepwise root-finding function
## ============================================================================
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
## 1. Two methods for mu0
## ============================================================================

mu0_ci_length_signed <- function(y, x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(alt_mle(y, x), error = function(e) NULL)
  if (is.null(mle) || mle$convergence != 0) return(NA_real_)
  psi_hat <- mle$mu0_hat
  if (is.na(psi_hat)) return(NA_real_)
  
  target <- function(psi) {
    r_val <- tryCatch(alt_r_param(1, psi, y, x, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(1e-4, abs(psi_hat) * 0.02)
  L <- step_root(target, psi_hat, step, "left", max_steps = 200)
  U <- step_root(target, psi_hat, step, "right", max_steps = 200)
  U - L
}

mu0_ci_length_modified <- function(y, x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(alt_mle(y, x), error = function(e) NULL)
  if (is.null(mle) || mle$convergence != 0) return(NA_real_)
  psi_hat <- mle$mu0_hat
  if (is.na(psi_hat)) return(NA_real_)
  
  target <- function(psi) {
    r_val <- tryCatch(alt_rstar_param(1, psi, y, x, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(1e-4, abs(psi_hat) * 0.02)
  L <- step_root(target, psi_hat, step, "left", max_steps = 200)
  U <- step_root(target, psi_hat, step, "right", max_steps = 200)
  U - L
}

## ============================================================================
## 2. Two methods for mu1
## ============================================================================

mu1_ci_length_signed <- function(y, x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(alt_mle(y, x), error = function(e) NULL)
  if (is.null(mle) || mle$convergence != 0) return(NA_real_)
  psi_hat <- mle$mu1_hat
  if (is.na(psi_hat)) return(NA_real_)
  
  target <- function(psi) {
    r_val <- tryCatch(alt_r_param(2, psi, y, x, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(1e-4, abs(psi_hat) * 0.02)
  L <- step_root(target, psi_hat, step, "left", max_steps = 200)
  U <- step_root(target, psi_hat, step, "right", max_steps = 200)
  U - L
}

mu1_ci_length_modified <- function(y, x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(alt_mle(y, x), error = function(e) NULL)
  if (is.null(mle) || mle$convergence != 0) return(NA_real_)
  psi_hat <- mle$mu1_hat
  if (is.na(psi_hat)) return(NA_real_)
  
  target <- function(psi) {
    r_val <- tryCatch(alt_rstar_param(2, psi, y, x, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(1e-4, abs(psi_hat) * 0.02)
  L <- step_root(target, psi_hat, step, "left", max_steps = 200)
  U <- step_root(target, psi_hat, step, "right", max_steps = 200)
  U - L
}

## ============================================================================
## 3. Two methods for sigma
## ============================================================================

sigma_ci_length_signed <- function(y, x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(alt_mle(y, x), error = function(e) NULL)
  if (is.null(mle) || mle$convergence != 0) return(NA_real_)
  psi_hat <- mle$sigma_hat
  if (is.na(psi_hat) || psi_hat <= 0) return(NA_real_)
  
  target <- function(psi) {
    if (psi <= 0) return(NA_real_)
    r_val <- tryCatch(alt_r_param(3, psi, y, x, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(1e-4, psi_hat * 0.02)
  L <- step_root(target, psi_hat, step, "left", max_steps = 200)
  U <- step_root(target, psi_hat, step, "right", max_steps = 200)
  L <- max(L, 1e-8)
  U - L
}

sigma_ci_length_modified <- function(y, x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(alt_mle(y, x), error = function(e) NULL)
  if (is.null(mle) || mle$convergence != 0) return(NA_real_)
  psi_hat <- mle$sigma_hat
  if (is.na(psi_hat) || psi_hat <= 0) return(NA_real_)
  
  target <- function(psi) {
    if (psi <= 0) return(NA_real_)
    r_val <- tryCatch(alt_rstar_param(3, psi, y, x, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(1e-4, psi_hat * 0.02)
  L <- step_root(target, psi_hat, step, "left", max_steps = 200)
  U <- step_root(target, psi_hat, step, "right", max_steps = 200)
  L <- max(L, 1e-8)
  U - L
}

## ============================================================================
## 4. Two methods for the log-quantile zeta_p at use stress
## ============================================================================

zeta_ci_length_signed <- function(y, x, x_use, p, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(alt_mle(y, x), error = function(e) NULL)
  if (is.null(mle) || mle$convergence != 0) return(NA_real_)
  s_p <- -log(1 - p)
  w_p <- log(s_p)
  psi_hat <- with(mle, mu0_hat + mu1_hat * x_use + sigma_hat * w_p)
  if (is.na(psi_hat)) return(NA_real_)
  
  target <- function(psi) {
    r_val <- tryCatch(alt_r_psi(psi, y, x, x_use, p, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(1e-4, abs(psi_hat) * 0.02)
  L <- step_root(target, psi_hat, step, "left", max_steps = 200)
  U <- step_root(target, psi_hat, step, "right", max_steps = 200)
  U - L
}

zeta_ci_length_modified <- function(y, x, x_use, p, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(alt_mle(y, x), error = function(e) NULL)
  if (is.null(mle) || mle$convergence != 0) return(NA_real_)
  s_p <- -log(1 - p)
  w_p <- log(s_p)
  psi_hat <- with(mle, mu0_hat + mu1_hat * x_use + sigma_hat * w_p)
  if (is.na(psi_hat)) return(NA_real_)
  
  target <- function(psi) {
    r_val <- tryCatch(alt_rstar_psi(psi, y, x, x_use, p, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(1e-4, abs(psi_hat) * 0.02)
  L <- step_root(target, psi_hat, step, "left", max_steps = 200)
  U <- step_root(target, psi_hat, step, "right", max_steps = 200)
  U - L
}

## ============================================================================
## 5. Simulation function for signed and modified roots
## ============================================================================
sim_alt_length_styled <- function(
    beta_vec     = c(0.5, 1, 5),
    mu0_vec      = log(c(2, 1, 0.5)),
    mu1_vec      = c(-0.5, -0.5, -0.5),
    x_levels     = c(0.5, 1, 2),
    n_per_level  = c(3, 4, 3),
    x_use        = 0,
    p_quant      = 0.1,
    conf_level   = 0.95,
    n_rep        = 500,
    seed         = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(beta_vec) == length(mu0_vec),
            length(beta_vec) == length(mu1_vec),
            length(x_levels) == length(n_per_level))
  
  params <- c("mu0", "mu1", "sigma", "log_quantile_zeta")
  methods <- c("signed", "modified")
  
  n_scenarios <- length(beta_vec)
  res_array <- array(NA_real_,
                     dim = c(n_scenarios, n_rep, length(params), length(methods)),
                     dimnames = list(scenario = 1:n_scenarios,
                                     rep = 1:n_rep,
                                     param = params,
                                     method = methods)
  )
  
  for (j in 1:n_scenarios) {
    beta_true  <- beta_vec[j]
    sigma_true <- 1 / beta_true
    mu0_true   <- mu0_vec[j]
    mu1_true   <- mu1_vec[j]
    
    for (rep in 1:n_rep) {
      y_all <- numeric(0)
      x_all <- numeric(0)
      for (k in seq_along(x_levels)) {
        s_k <- x_levels[k]
        n_k <- n_per_level[k]
        eta_k <- exp(mu0_true + mu1_true * s_k)
        T_k <- rweibull(n_k, shape = beta_true, scale = eta_k)
        y_all <- c(y_all, log(T_k))
        x_all <- c(x_all, rep(s_k, n_k))
      }
      
      res_array[j, rep, "mu0", "signed"]   <- mu0_ci_length_signed(y_all, x_all, conf_level)
      res_array[j, rep, "mu0", "modified"] <- mu0_ci_length_modified(y_all, x_all, conf_level)
      
      res_array[j, rep, "mu1", "signed"]   <- mu1_ci_length_signed(y_all, x_all, conf_level)
      res_array[j, rep, "mu1", "modified"] <- mu1_ci_length_modified(y_all, x_all, conf_level)
      
      res_array[j, rep, "sigma", "signed"]   <- sigma_ci_length_signed(y_all, x_all, conf_level)
      res_array[j, rep, "sigma", "modified"] <- sigma_ci_length_modified(y_all, x_all, conf_level)
      
      res_array[j, rep, "log_quantile_zeta", "signed"]   <- zeta_ci_length_signed(y_all, x_all, x_use, p_quant, conf_level)
      res_array[j, rep, "log_quantile_zeta", "modified"] <- zeta_ci_length_modified(y_all, x_all, x_use, p_quant, conf_level)
      
      if (rep %% 100 == 0) cat("Scenario", j, "- rep", rep, "\n")
    }
  }
  
  mean_array <- apply(res_array, c(1,3,4), mean, na.rm = TRUE)
  df_list <- list()
  for (j in 1:n_scenarios) {
    for (par in params) {
      for (meth in methods) {
        df_list[[length(df_list) + 1]] <- data.frame(
          scenario = j,
          beta_true = beta_vec[j],
          mu0_true  = mu0_vec[j],
          mu1_true  = mu1_vec[j],
          sigma_true = 1 / beta_vec[j],
          param = par,
          method = meth,
          avg_length = mean_array[j, par, meth],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, df_list)
}


############################################################
## Example run: three ALT scenarios
############################################################

beta_vec <- c(0.5, 1, 5)                      # shapes
mu0_vec  <- log(c(2, 1, 0.5))                 # eta(0) = 2, 1, 0.5
mu1_vec  <- c(-0.5, -0.5, -0.5)               # same acceleration slope

# mu0_vec   <- 19.72264
# mu1_vec   <- -0.4245567
# sigma_true <- 0.22713
# beta_vec  <- 1 / sigma_true   # approximately 4.403
# x_levels    <- c(28.84, 34.68, 38.02) * 0.1
# n_per_level <- c(6, 6, 6)
# x_use       <- 20 * 0.1
# conf_level  <- 0.95
# n_rep       <- 100              # Number of replications; increase to 1000 if needed
# seed        <- 123

res_alt_all <- sim_alt_params_quantile(
  beta_vec     = beta_vec,
  mu0_vec      = mu0_vec,
  mu1_vec      = mu1_vec,
  x_levels     = c(0.5, 1, 2),
  n_per_level  = c(3, 4, 3),   # total n = 10
  x_use        = 0,
  p_quant      = 0.1,
  conf_level   = 0.95,
  n_rep        = 10000,
  seed         = 123
)

print(res_alt_all)

## =========================================================
##  ALT coverage plot
##  res_alt_all has columns:
##  scenario, beta_true, mu0_true, mu1_true, sigma_true,
##  param, method, coverage
## =========================================================

library(dplyr)
library(ggplot2)

res_alt_plot <- res_alt_all %>%
  mutate(
    ## Row facets: parameter of interest (plotmath syntax)
    param_lab = factor(
      param,
      levels = c("mu0", "mu1", "sigma", "log_quantile_zeta"),
      labels = c(
        "mu[0]",
        "mu[1]",
        "beta",
        "Log~quantile~zeta[0.1]"
      )
    ),
    ## Legend labels
    method_lab = factor(
      method,
      levels = c("signed_root", "modified_root"),
      labels = c("Signed root", "Modified root")
    ),
    ## x axis: true shape beta, to be parsed as math
    beta_lab = factor(
      beta_true,
      levels = c(0.5, 1, 5),
      labels = c("beta==0.5", "beta==1", "beta==5")
    )
  )

p_alt <- ggplot(
  res_alt_plot,
  aes(x = beta_lab, y = coverage,
      color = method_lab, group = method_lab)
) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(
    param_lab ~ .,
    labeller = labeller(param_lab = label_parsed)
  ) +
  scale_x_discrete(
    labels = function(x) parse(text = x)
  ) +
  scale_y_continuous(
    limits = c(0.86, 0.97),
    breaks = seq(0.86, 0.96, by = 0.02)
  ) +
  labs(
    x = expression(True~shape~beta),
    y = "Coverage probability",
    color = "Method"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  )

## Show on screen
print(p_alt)

## Save to file
# ggsave("weibull_alt_coverage.pdf", p_alt,
#        width = 6, height = 6)


## ============================================================================
## 6. Run simulation and plot
## ============================================================================
res_len <- sim_alt_length_styled(
  beta_vec     = c(0.5, 1, 5),
  mu0_vec      = log(c(2, 1, 0.5)),
  mu1_vec      = c(-0.5, -0.5, -0.5),
  x_levels     = c(0.5, 1, 2),
  n_per_level  = c(3, 4, 3),
  x_use        = 0,
  p_quant      = 0.1,
  conf_level   = 0.95,
  n_rep        = 10000,
  seed         = 456
)

library(ggplot2)
library(dplyr)
library(ggh4x)

plot_len <- res_len %>%
  mutate(
    param_lab = factor(param,
                       levels = c("mu0", "mu1", "sigma", "log_quantile_zeta"),
                       labels = c("mu[0]", "mu[1]", "beta", "Log~quantile~zeta[0.1]")
    ),
    method_lab = factor(method,
                        levels = c("signed", "modified"),
                        labels = c("Signed root", "Modified root")
    ),
    beta_lab = factor(beta_true,
                      levels = c(0.5, 1, 5),
                      labels = c("beta==0.5", "beta==1", "beta==5")
    )
  )

p_len <- ggplot(plot_len,
                aes(x = beta_lab, y = avg_length,
                    color = method_lab, group = method_lab)) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(param_lab ~ ., scales = "free_y",
             labeller = labeller(param_lab = label_parsed)) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  labs(x = expression(True~shape~beta),
       y = "Average interval length",
       color = "Method") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(size = 10),
        legend.position = "bottom")

print(p_len)
