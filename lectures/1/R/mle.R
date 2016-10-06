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

mle_results <- function (mle_estimation_result, t_vector, phi) {
  w          <- as.matrix(mle_estimation_result$par[-length(mle_estimation_result$par)])
  sigma      <- last(mle_estimation_result$par)
  variance   <- -solve(mle_estimation_result$hessian)
  w_se       <- sqrt(diag(w_variance))[1:length(w)]
  sigma_se   <- last(sqrt(diag(w_variance)))
  t_value    <- w/w_se
  p_value    <- 2 * (1 - pnorm(abs(t_value))) 
  e          <- t_vector - as.matrix(phi) %*% w
  return(list(w = w, sigma = sigma, variance = variance, w_se = w_se, sigma_se = sigma_se,
              t_value = t_value, p_value = p_value, e = e))
}

#summary(lm(as.matrix(t_vector)~as.matrix(phi)-1))
