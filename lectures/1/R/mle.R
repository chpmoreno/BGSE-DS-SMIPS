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
  w          <- as.matrix(mle_estimation_result$par[-length(mle_estimation_result$par)])
  sigma      <- last(mle_estimation_result$par)
  variance   <- -solve(mle_estimation_result$hessian)
  w_se       <- sqrt(diag(variance))[1:length(w)]
  sigma_se   <- last(sqrt(diag(variance)))
  t_value    <- w/w_se
  p_value    <- 2 * (1 - pnorm(abs(t_value))) 
  e          <- t_vector - as.matrix(phi) %*% w
  return(list(w = w, sigma = sigma, variance = variance, w_se = w_se, sigma_se = sigma_se,
              t_value = t_value, p_value = p_value, e = e))
}
expression(paste("Value is ", sigma,",", R^{2},'=0.6'))

mle_plots <- function(mle_results){
  w_plot_data <- data.frame(coeff_number = seq(1:nrow(mle_results$w))) %>% 
                 bind_cols(data.frame(w = mle_results$w)) %>% 
                 mutate(lower_bound = w - 1.96 * mle_results$w_se) %>% 
                 mutate(upper_bound = w + 1.96 * mle_results$w_se) %>% 
                 mutate(color_ci = ifelse(mle_results$p_value < 0.05, "green", "red"))
  ci_plot <- ggplot(data = w_plot_data, aes(x = as.factor(coeff_number), y = w)) + 
             geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
                           color = w_plot_data$color_ci, width = 0.3) + 
             geom_point() + ggtitle("w estimator and confidence intervals") +
             labs(x = "Coefficient", y = "w, w +/- 1.96 s.e.") + theme_economist()
  
}
summary(lm(as.matrix(t_vector)~as.matrix(phi)-1))
