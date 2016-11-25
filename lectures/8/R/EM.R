# Author: José Fernado Moreno Gutiérrez

# load libraries and functions needed
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)

m <- 30 # number of variables taken for the model

# load data - Find it on https://1drv.ms/t/s!Ai0XbELt7PXquFJtWenDbxvksoXM ####
data <- read.table(file = "~/Documents/OneDrive/Documents/BGSE/First_Term/SMI/Datasets/synthetic_regression/synthetic_regression.txt",
                   nrow = 300)[,1:(m + 1)]

# Initial parameters
w0 <- rep(1,31)
q0 <- 1
nu <- 10
t <- as.vector(data[, 1])
X <- as.matrix(cbind(rep(1,nrow(data)), as.matrix(data[,2:31])))
N <- length(t)

mle_estimator_lm <- function (w, t_vector, phi) {
  if(is.vector(t_vector) == FALSE)
    t_vector <- as.vector(t_vector)
  if(is.matrix(phi) == FALSE)
    phi <- as.matrix(phi)
  M <- ncol(phi) 
  phi_w <- phi %*% w[1:M]
  Sig <- w[(M + 1):(M + 1)]
  sum(-(1 / 2) * log(2 * pi) - (1 / 2) * log(Sig ^ 2) - (1 / (2 * Sig ^ 2)) * (t_vector - phi_w) ^ 2)
}

mle_result <- function (mle_estimation_result, t_vector, phi) {
  phi        <- as.matrix(phi) 
  w          <- as.matrix(mle_estimation_result$par[-length(mle_estimation_result$par)])
  sigma      <- last(mle_estimation_result$par)
  variance   <- -solve(mle_estimation_result$hessian)
  w_se       <- sqrt(diag(variance))[1:length(w)]
  sigma_se   <- last(sqrt(diag(variance)))
  z_value    <- w/w_se
  p_value    <- 2 * (1 - pnorm(abs(z_value)))
  t_hat      <- phi %*% w 
  e          <- t_vector - t_hat
  e_st       <- e/sigma_se
  leverage   <- NULL 
  for(i in 1:nrow(phi)){
    leverage <- c(leverage, t(phi[i,]) %*% solve(t(phi) %*% phi) %*% phi[i,])
  }
  dev <- -2 * ((1 / 2) * log(2 * pi) - (1 / 2) * log(sigma ^ 2) - (1 / (2 * sigma ^ 2)) * (e ^ 2))
  # leverage is equivalent to diag(phi %*% solve(t(phi) %*% phi) %*% t(phi))
  # diagonal of the hat matrix
  return(list(w = w, sigma = sigma, variance = variance, w_se = w_se, sigma_se = sigma_se,
              z_value = z_value, p_value = p_value, t_hat = t_hat, e = e, e_st = e_st,
              leverage = leverage, dev = dev, n = nrow(phi)))
}

mle_plot <- function(mle_results){
  ci_plot_data <- data.frame(coeff_number = seq(1:nrow(mle_results$w))) %>% 
    bind_cols(data.frame(w = mle_results$w)) %>% 
    mutate(lower_bound = w - 1.96 * mle_results$w_se) %>% 
    mutate(upper_bound = w + 1.96 * mle_results$w_se) %>% 
    mutate(color_ci = ifelse(mle_results$p_value < 0.05, "blue", "red"))
  
  ci_plot <- ggplot(data = ci_plot_data, aes(x = as.factor(coeff_number), y = w)) + 
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
                  color = ci_plot_data$color_ci, width = 0.3) + 
    geom_point() + ggtitle("Maximum Likelihood Estimator") +
    labs(x = "Coefficient", y = "w, w +/- 1.96 s.e.") + theme_excel()
  
  de_plot_data <- data.frame(observation_number = seq(1:length(mle_results$dev))) %>% 
    bind_cols(data.frame(dev = mle_results$dev)) %>% 
    mutate(color_de = ifelse(dev > as.numeric(quantile(mle_results$dev, c(0.005, 0.995))[2]) | 
                             dev < as.numeric(quantile(mle_results$dev, c(0.005, 0.995))[1]), "red","black"))
  
  de_plot <- ggplot(data = de_plot_data, aes(x = observation_number, y = dev)) + 
    geom_point(color = de_plot_data$color_de) +  geom_hline(yintercept =  quantile(mle_results$dev, c(0.005, 0.995))[1]) +
    geom_hline(yintercept =  quantile(mle_results$dev, c(0.005, 0.995))[2]) + 
    ggtitle("Maximum Likelihood Estimator") +
    labs(x = "Observation", y = "Deviance Error") + theme_excel()
  
  return(list(ci_plot = ci_plot, de_plot = de_plot))
}

# mle_estimation
mle_estimation_result <- optim(runif(m + 2, 0, 1), mle_estimator_lm, phi = X, 
                               t_vector = t, method = "BFGS", 
                               control = list(trace = 1, maxit = 10000, fnscale = -1),
                               hessian = TRUE)

# mle_results
mle_results <- mle_result(mle_estimation_result, t, X)

# mle graphics
mle_graphics <- mle_plot(mle_results)
mle_graphics$ci_plot
mle_graphics$de_plot

# EM algorithm functions
e_step <- function(t, X, q, w, nu) {
  eta <- (nu + 1) / (nu + q * (t - X %*% w) ^ 2 - 2)
  return(eta)
}

m_step <- function(t, X, q0, w0, eta) {
  w <- solve(t(X) %*% diag(as.vector(eta)) %*% X, t(X) %*% diag(as.vector(eta)) %*% t)
  q <- (N / ( t(t - X %*% w) %*% diag(as.vector(eta)) %*% (t - X %*% w)))
  return(list(w = as.vector(w), q = as.numeric(q)))
}

rob_reg <- function(t, X, q0, w0, nu, max_iter = 100) {
  q <- q0
  w <- w0
  q_conv      <- NULL
  e_conv      <- NULL
  loglik_conv <- NULL
  for(i in 1:max_iter) {
    eta    <- e_step(t, X, q, w, nu)
    wq     <- m_step(t, X, q, w, eta)
    q      <- wq$q
    w      <- wq$w
    loglik <- (length(t) / 2) * log(q) - q/2 * (t((t - X %*% w)) * diag(eta)) %*% (t - X %*% w)
    q_conv[i]      <- q
    e_conv[i]      <- sum((t - X %*% w) ^ 2)
    loglik_conv[i] <- loglik
  }
  e   <- as.vector(t - X %*% w)
  se  <- as.vector(sqrt(diag(solve(q * t(as.vector((nu + 1) * (nu - 2 - q * e ^ 2) / ((nu + q * e ^ 2 - 2) ^ 2)) * X) %*% X))))
  dev <- NULL 
  for(i in 1:length(t)){
    dev[i] <- -2 * ((1 / 2) * log(q) - q/2 * (t[i] - X[i,] %*% w) * diag(eta) * (t[i] - X[i,] %*% w))
  }
  return(list(w = w, q = q, eta = eta, e = e, se = se, e_conv = e_conv, q_conv = q_conv, loglik_conv = loglik_conv, dev = dev))
}

# Execution
rob_reg_results <- rob_reg(t, X, q0, w0, nu)

rob_reg_plot <- function(rob_reg_results){
  ci_plot_data <- data.frame(coeff_number = seq(1:length(rob_reg_results$w))) %>% 
    bind_cols(data.frame(w = rob_reg_results$w)) %>% 
    mutate(lower_bound = w - 1.96 * rob_reg_results$se) %>% 
    mutate(upper_bound = w + 1.96 * rob_reg_results$se) %>% 
    mutate(color_ci = ifelse(upper_bound > 0 & lower_bound < 0, "red", "blue"))
  
  ci_plot <- ggplot(data = ci_plot_data, aes(x = as.factor(coeff_number), y = w)) + 
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
                  color = ci_plot_data$color_ci, width = 0.3) + 
    geom_point() + ggtitle("Robust Bayesian Estimator") +
    labs(x = "Coefficient", y = "w, w +/- 1.96 s.e.") + theme_excel()
  
  de_plot_data <- data.frame(observation_number = seq(1:length(rob_reg_results$dev))) %>% 
    bind_cols(data.frame(dev = rob_reg_results$dev)) %>% 
    mutate(color_de = ifelse(dev > as.numeric(quantile(rob_reg_results$dev, c(0.005, 0.995))[2]) | 
                             dev < as.numeric(quantile(rob_reg_results$dev, c(0.005, 0.995))[1]), "red","black"))
  
  de_plot <- ggplot(data = de_plot_data, aes(x = observation_number, y = dev)) + 
    geom_point(color = de_plot_data$color_de) +  geom_hline(yintercept =  quantile(rob_reg_results$dev, c(0.005, 0.995))[1]) +
    geom_hline(yintercept =  quantile(rob_reg_results$dev, c(0.005, 0.995))[2]) + 
    ggtitle("Robust Bayesian Estimator") +
    labs(x = "Observation", y = "Deviance Error") + theme_excel()
  
  tolerance <- 10 ^ -10
  conv_iter <- which((abs(diff(rob_reg_results$e_conv)) < tolerance & 
                      abs(diff(rob_reg_results$q_conv)) < tolerance &
                      abs(diff(rob_reg_results$loglik_conv)) < tolerance) == TRUE)[1]
  
  conv_plot_data <- data.frame(iter_number = seq(1:conv_iter), e_conv = rob_reg_results$e_conv[1:conv_iter], 
                               q_conv = rob_reg_results$q_conv[1:conv_iter],
                               loglik_conv = rob_reg_results$loglik_conv[1:conv_iter])
  
  q_conv_plot <- ggplot(data = conv_plot_data, aes(x = iter_number, y = q_conv)) +
    geom_point() + ggtitle("q convergence") +
    labs(x = "iteration", y = "q") + theme_excel()
  
  loglike_conv_plot <- ggplot(data = conv_plot_data, aes(x = iter_number, y = loglik_conv)) +
    geom_point() + ggtitle("loglik convergence") +
    labs(x = "iteration", y = "loglik") + theme_excel()
  
  e_conv_plot <- ggplot(data = conv_plot_data, aes(x = iter_number, y = e_conv)) +
    geom_point() + ggtitle("RSS convergence") +
    labs(x = "iteration", y = "RSS") + theme_excel()
  
  return(list(ci_plot = ci_plot, de_plot = de_plot, q_conv_plot = q_conv_plot, loglike_conv_plot = loglike_conv_plot,
              e_conv_plot = e_conv_plot))
}

# robust regression graphics
rob_reg_graphics <- rob_reg_plot(rob_reg_results)

rob_reg_graphics$ci_plot
rob_reg_graphics$de_plot
rob_reg_graphics$q_conv_plot
rob_reg_graphics$loglike_conv_plot
rob_reg_graphics$e_conv_plot
