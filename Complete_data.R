# 精简版：比较三种方法的覆盖率和区间长度
rm(list = ls())

# 辅助函数
clikelevels <- function(N) {
  cs90 <- 0.2585 - (1 / (0.775 + (1.641 * N)))
  cs95 <- 0.1465 - (1 / (2.757 + (2.082 * N)))
  cs99 <- 0.0362 - (1 / (18.294 + (5.229 * N)))
  c(cs90, cs95, cs99)
}

get_c_level <- function(N, conf_level) {
  csig <- clikelevels(N)
  if (conf_level == 0.9) return(csig[1])
  if (conf_level == 0.95) return(csig[2])
  if (conf_level == 0.99) return(csig[3])
}

# 步进法求根
step_root <- function(target_fun, init_val, step, direction, max_steps = 500) {
  direction <- match.arg(direction, c("left", "right"))
  sign_ref <- sign(target_fun(init_val))
  if (sign_ref == 0 || is.na(sign_ref)) return(init_val)
  
  step <- abs(step) * if (direction == "left") -1 else 1
  cur_val <- init_val + step
  f_prev <- target_fun(init_val)
  
  for (i in 1:max_steps) {
    f_cur <- tryCatch(target_fun(cur_val), error = function(e) NA_real_)
    if (is.na(f_cur)) { cur_val <- cur_val + step; next }
    if (sign(f_cur) != sign_ref || f_cur == 0) {
      if (f_prev != f_cur) {
        return(cur_val - step + (0 - f_prev) / (f_cur - f_prev) * step)
      }
      return(cur_val)
    }
    f_prev <- f_cur
    cur_val <- cur_val + step
  }
  return(cur_val)
}

# 形状参数β的推断函数
profile_ll_beta <- function(beta, x) {
  if (beta <= 0) return(-Inf)
  n <- length(x)
  bx <- beta * x
  m <- max(bx)
  w <- exp(bx - m)
  K1 <- exp(m) * sum(w)
  n * log(beta) + beta * sum(x) - n * log(K1)
}

a_hat_beta <- function(beta, x) {
  n <- length(x)
  bx <- beta * x
  m <- max(bx)
  w <- exp(bx - m)
  K1 <- exp(m) * sum(w)
  (log(K1) - log(n)) / beta
}

mle_beta_profile <- function(x) {
  f <- function(logb) -profile_ll_beta(exp(logb), x)
  opt <- optimize(f, interval = c(log(1e-3), log(1e3)))
  beta_hat <- exp(opt$minimum)
  list(beta_hat = beta_hat, a_hat = a_hat_beta(beta_hat, x), ll_max = -opt$objective)
}

r_beta <- function(beta0, x, mle_obj) {
  if (beta0 <= 0) stop("beta0 must be positive")
  beta_hat <- mle_obj$beta_hat
  ll_hat <- mle_obj$ll_max
  ll0 <- profile_ll_beta(beta0, x)
  sign(beta_hat - beta0) * sqrt(max(0, 2 * (ll_hat - ll0)))
}

rstar_beta <- function(beta0, x, mle_obj) {
  if (beta0 <= 0) stop("beta0 must be positive")
  beta_hat <- mle_obj$beta_hat
  a_hat <- mle_obj$a_hat
  ll_hat <- mle_obj$ll_max
  
  ll0 <- profile_ll_beta(beta0, x)
  r <- sign(beta_hat - beta0) * sqrt(max(0, 2 * (ll_hat - ll0)))
  if (!is.finite(r) || abs(r) < 1e-10) return(0)
  
  a0 <- a_hat_beta(beta0, x)
  n <- length(x)
  
  # 简化版的J矩阵和phi计算
  dr <- function(beta, a) {
    delta <- x - a
    list(delta = delta, r = exp(beta * delta))
  }
  
  J_beta_a <- function(beta, a) {
    drv <- dr(beta, a)
    delta <- drv$delta
    r <- drv$r
    J_bb <- n / beta^2 + sum(delta^2 * r)
    J_aa <- beta^2 * sum(r)
    J_ba <- n - sum(r) - beta * sum(delta * r)
    matrix(c(J_bb, J_ba, J_ba, J_aa), 2, 2)
  }
  
  phi_beta_a <- function(beta, a) {
    drv <- dr(beta, a)
    delta <- drv$delta
    r <- drv$r
    S0 <- sum(1 - r)
    S1 <- sum((1 - r) * delta)
    c(beta * S0, beta^2 * S1)
  }
  
  phi_lambda <- function(beta, a) {
    drv <- dr(beta, a)
    r <- drv$r
    c(beta^2 * sum(r), beta^2 * sum(beta * drv$delta * r - (1 - r)))
  }
  
  phi_theta <- function(beta, a) {
    drv <- dr(beta, a)
    delta <- drv$delta
    r <- drv$r
    S0 <- sum(1 - r)
    S1 <- sum((1 - r) * delta)
    dS1_a <- sum(beta * r * delta - (1 - r))
    matrix(c(S0 - beta * sum(delta * r), beta^2 * sum(r),
             2 * beta * S1 + beta^2 * (-sum(delta^2 * r)), beta^2 * dS1_a), 2, 2)
  }
  
  phi_hat <- phi_beta_a(beta_hat, a_hat)
  phi_constr <- phi_beta_a(beta0, a0)
  phi_diff <- phi_hat - phi_constr
  phi_lambda_constr <- phi_lambda(beta0, a0)
  phi_theta_hat <- phi_theta(beta_hat, a_hat)
  J_hat <- J_beta_a(beta_hat, a_hat)
  J_constr <- J_beta_a(beta0, a0)
  J_aa_constr <- J_constr[2, 2]
  
  num_det <- det(cbind(phi_diff, phi_lambda_constr))
  den_det <- det(phi_theta_hat)
  detJ_hat <- det(J_hat)
  if (detJ_hat <= 0) detJ_hat <- abs(detJ_hat)
  
  q_val <- (num_det / den_det) * sqrt(detJ_hat / J_aa_constr)
  if (q_val * r <= 0) q_val <- abs(q_val) * sign(r)
  
  out <- r + (1 / r) * log(q_val / r)
  if (!is.finite(out)) out <- r
  out
}

# 置信区间长度计算
beta_ci_length <- function(x, method = "signed", conf_level = 0.95) {
  z <- qnorm((1 + conf_level) / 2)
  mle <- tryCatch(mle_beta_profile(x), error = function(e) NULL)
  if (is.null(mle)) return(NA_real_)
  beta_hat <- mle$beta_hat
  if (is.na(beta_hat) || beta_hat <= 0) return(NA_real_)
  
  target <- function(b) {
    if (b <= 0) return(NA_real_)
    r_val <- switch(method,
                    signed = tryCatch(r_beta(b, x, mle), error = function(e) NA_real_),
                    modified = tryCatch(rstar_beta(b, x, mle), error = function(e) NA_real_))
    if (is.na(r_val)) return(NA_real_)
    abs(r_val) - z
  }
  
  # Technometrics方法单独处理
  if (method == "tech") {
    n <- length(x)
    c_level <- get_c_level(n, conf_level)
    if (is.na(c_level)) return(NA_real_)
    log_c <- log(c_level)
    target <- function(b) {
      if (b <= 0) return(NA_real_)
      ll <- tryCatch(profile_ll_beta(b, x), error = function(e) NA_real_)
      if (is.na(ll)) return(NA_real_)
      ll - mle$ll_max - log_c
    }
  }
  
  step <- max(0.01, beta_hat * 0.02)
  L <- step_root(target, beta_hat, step, "left", max_steps = 200)
  U <- step_root(target, beta_hat, step, "right", max_steps = 200)
  U - L
}

# 模拟函数
sim_length_grid <- function(beta_vec = c(0.5, 1, 5), eta_vec = c(2, 1, 0.5),
                            n_vec = c(5, 10, 20), n_rep = 1000, 
                            conf_level = 0.95, seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(beta_vec) == length(eta_vec))
  
  res_list <- list()
  idx <- 1
  
  for (k in seq_along(beta_vec)) {
    beta_true <- beta_vec[k]
    eta_true <- eta_vec[k]
    
    for (n in n_vec) {
      len_signed <- numeric(n_rep)
      len_modified <- numeric(n_rep)
      len_tech <- numeric(n_rep)
      
      for (rep in seq_len(n_rep)) {
        y <- rweibull(n, shape = beta_true, scale = eta_true)
        x <- log(y)
        
        len_signed[rep] <- beta_ci_length(x, "signed", conf_level)
        len_modified[rep] <- beta_ci_length(x, "modified", conf_level)
        len_tech[rep] <- beta_ci_length(x, "tech", conf_level)
      }
      
      res_list[[idx]] <- data.frame(
        beta = beta_true, eta = eta_true, n = n,
        signed_root = mean(len_signed, na.rm = TRUE),
        modified_root = mean(len_modified, na.rm = TRUE),
        technometricsPL = mean(len_tech, na.rm = TRUE)
      )
      idx <- idx + 1
    }
  }
  
  do.call(rbind, res_list)
}

# 覆盖率模拟
sim_coverage_grid <- function(beta_vec = c(0.5, 1, 5), eta_vec = c(2, 1, 0.5),
                              n_vec = c(5, 10, 20), n_rep = 1000,
                              conf_level = 0.95, seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(beta_vec) == length(eta_vec))
  
  z <- qnorm((1 + conf_level) / 2)
  res_list <- list()
  idx <- 1
  
  for (k in seq_along(beta_vec)) {
    beta_true <- beta_vec[k]
    eta_true <- eta_vec[k]
    
    for (n in n_vec) {
      c_level <- get_c_level(n, conf_level)
      log_c <- log(c_level)
      cover_signed <- cover_modified <- cover_tech <- 0
      
      for (rep in seq_len(n_rep)) {
        y <- rweibull(n, shape = beta_true, scale = eta_true)
        x <- log(y)
        
        mle <- mle_beta_profile(x)
        beta_hat <- mle$beta_hat
        ll_hat <- mle$ll_max
        
        r_val <- r_beta(beta_true, x, mle)
        rstar_val <- rstar_beta(beta_true, x, mle)
        ll_true <- profile_ll_beta(beta_true, x)
        
        if (is.finite(r_val) && abs(r_val) <= z) cover_signed <- cover_signed + 1
        if (is.finite(rstar_val) && abs(rstar_val) <= z) cover_modified <- cover_modified + 1
        if (ll_true - ll_hat >= log_c) cover_tech <- cover_tech + 1
      }
      
      res_list[[idx]] <- data.frame(
        beta = beta_true, eta = eta_true, n = n,
        signed_root = cover_signed / n_rep,
        modified_root = cover_modified / n_rep,
        technometricsPL = cover_tech / n_rep
      )
      idx <- idx + 1
    }
  }
  
  do.call(rbind, res_list)
}

# 绘图对比
library(ggplot2)
library(tidyr)


# 长度图
len_plot <- len_res %>%
  pivot_longer(cols = c(signed_root, modified_root, technometricsPL),
               names_to = "method", values_to = "avg_length") %>%
  mutate(method = factor(method, levels = c("signed_root", "modified_root", "technometricsPL"),
                         labels = c("Signed root", "Modified root", "Technometrics"))) %>%
  ggplot(aes(x = n, y = avg_length, color = method, group = method)) +
  geom_point(size = 2) + geom_line() +
  facet_wrap(~ paste0("β=", beta, ", η=", eta), ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(5, 10, 20)) +
  labs(x = "Sample size n", y = "Average interval length", color = "Method") +
  theme_bw() + theme(legend.position = "bottom")

print(len_plot)