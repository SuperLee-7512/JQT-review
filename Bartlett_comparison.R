#length-revise-buchang
rm(list = ls())

###############################################
## Complete data - Weibull simulation
## Parameters: shape beta
## Methods: signed root r, modified root r*, Bartlett (r²·c)
###############################################

## --------- 0. Helper: step_root ---------

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

###############################################
## 1. Shape parameter - profile in beta (complete data)
###############################################

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

## canonical parameter in (beta, a)
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

###############################################
## 2. CI functions for shape_beta (three methods)
###############################################

ci_beta_signed <- function(x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(c(NA_real_, NA_real_))
  beta_hat <- mle$beta_hat
  if (is.na(beta_hat) || beta_hat <= 0) return(c(NA_real_, NA_real_))
  
  target <- function(b) {
    if (b <= 0) return(NA_real_)
    r_val <- tryCatch(r_beta(b, x, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  step <- max(0.01, beta_hat * 0.02)
  L <- step_root(target, beta_hat, step, "left",  max_steps = 200)
  U <- step_root(target, beta_hat, step, "right", max_steps = 200)
  c(L, U)
}

ci_beta_modified <- function(x, conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(c(NA_real_, NA_real_))
  beta_hat <- mle$beta_hat
  if (is.na(beta_hat) || beta_hat <= 0) return(c(NA_real_, NA_real_))
  
  target <- function(b) {
    if (b <= 0) return(NA_real_)
    rstar_val <- tryCatch(rstar_beta(b, x, mle), error = function(e) NA_real_)
    if (is.na(rstar_val)) return(NA_real_)
    abs(rstar_val) - z
  }
  step <- max(0.01, beta_hat * 0.02)
  L <- step_root(target, beta_hat, step, "left",  max_steps = 200)
  U <- step_root(target, beta_hat, step, "right", max_steps = 200)
  c(L, U)
}

ci_beta_bartlett <- function(x, conf_level = 0.95) {
  n <- length(x)
  R <- 1.56715 / n
  c_val <- 1 / (1 + R)
  chi2_thresh <- qchisq(conf_level, 1)
  
  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(c(NA_real_, NA_real_))
  beta_hat <- mle$beta_hat
  if (is.na(beta_hat) || beta_hat <= 0) return(c(NA_real_, NA_real_))
  
  target <- function(b) {
    if (b <= 0) return(NA_real_)
    r_val <- tryCatch(r_beta(b, x, mle), error = function(e) NA_real_)
    if (is.na(r_val)) return(NA_real_)
    r_val^2 * c_val - chi2_thresh
  }
  step <- max(0.01, beta_hat * 0.02)
  L <- step_root(target, beta_hat, step, "left",  max_steps = 200)
  U <- step_root(target, beta_hat, step, "right", max_steps = 200)
  c(L, U)
}

###############################################
## 3. Simulation function - coverage and length for shape_beta
##    (signed root, modified root, Bartlett)
###############################################

sim_beta_cp_length <- function(
    beta_vec  = c(0.5, 1, 5),
    eta_vec   = c(2,   1, 0.5),
    n_vec     = c(5, 10, 20),
    n_rep     = 1000,
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
    
    for (n in n_vec) {
      cover_signed   <- 0
      cover_modified <- 0
      cover_bartlett <- 0
      len_signed     <- numeric(n_rep)
      len_modified   <- numeric(n_rep)
      len_bartlett   <- numeric(n_rep)
      
      for (rep in seq_len(n_rep)) {
        y <- rweibull(n, shape = beta_true, scale = eta_true)
        x <- log(y)
        
        # signed root
        ci_s <- ci_beta_signed(x, conf_level)
        if (!any(is.na(ci_s))) {
          if (ci_s[1] <= beta_true && beta_true <= ci_s[2]) cover_signed <- cover_signed + 1L
          len_signed[rep] <- ci_s[2] - ci_s[1]
        } else {
          len_signed[rep] <- NA_real_
        }
        
        # modified root
        ci_m <- ci_beta_modified(x, conf_level)
        if (!any(is.na(ci_m))) {
          if (ci_m[1] <= beta_true && beta_true <= ci_m[2]) cover_modified <- cover_modified + 1L
          len_modified[rep] <- ci_m[2] - ci_m[1]
        } else {
          len_modified[rep] <- NA_real_
        }
        
        # Bartlett (r²·c)
        ci_b <- ci_beta_bartlett(x, conf_level)
        if (!any(is.na(ci_b))) {
          if (ci_b[1] <= beta_true && beta_true <= ci_b[2]) cover_bartlett <- cover_bartlett + 1L
          len_bartlett[rep] <- ci_b[2] - ci_b[1]
        } else {
          len_bartlett[rep] <- NA_real_
        }
        
        if (rep %% 100 == 0) {
          cat("beta =", beta_true, "eta =", eta_true, "n =", n,
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
        method  = c("signed_root", "modified_root", "bartlett"),
        coverage = c(cover_signed, cover_modified, cover_bartlett) / n_rep,
        avg_length = c(mean(len_signed, na.rm = TRUE), 
                       mean(len_modified, na.rm = TRUE),
                       mean(len_bartlett, na.rm = TRUE))
      )
      idx <- idx + 1
    }
  }
  
  do.call(rbind, res_list)
}

# ============================================================================
# Run simulation
# ============================================================================
res_beta <- sim_beta_cp_length(
  beta_vec  = c(0.5, 1, 5),
  eta_vec   = c(2,   1, 0.5),
  n_vec     = c(5, 10, 20),
  n_rep     = 1000,
  conf_level = 0.95,
  seed      = 123
)
print(res_beta)

# ============================================================================
# Plots
# ============================================================================
library(ggplot2)
library(dplyr)

plot_data <- res_beta %>%
  mutate(
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
                        levels = c("signed_root", "modified_root", "bartlett"),
                        labels = c("Signed root", "Modified root", "Bartlett"))
  )

# Coverage plot
p_cov <- ggplot(plot_data, aes(x = n, y = coverage, color = method_lab, group = method_lab)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(size = 2) + geom_line() +
  facet_wrap(~ beta_eta, labeller = label_parsed, ncol = 3) +
  scale_x_continuous(breaks = c(5, 10, 20)) +
  scale_y_continuous(limits = c(0.80, 1.00), breaks = seq(0.80, 1.00, 0.05)) +
  labs(x = "Sample size n", y = "Coverage probability", color = "Method") +
  theme_bw() + theme(legend.position = "bottom")
print(p_cov)

# Length plot
p_len <- ggplot(plot_data, aes(x = n, y = avg_length, color = method_lab, group = method_lab)) +
  geom_point(size = 2) + geom_line() +
  facet_wrap(~ beta_eta, labeller = label_parsed, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(5, 10, 20)) +
  labs(x = "Sample size n", y = "Average interval length", color = "Method") +
  theme_bw() + theme(legend.position = "bottom")
print(p_len)

plot_grid(p_cov, p_len, ncol = 1)

ggarrange(p_cov, p_len, ncol = 1, nrow = 2, common.legend = T)
