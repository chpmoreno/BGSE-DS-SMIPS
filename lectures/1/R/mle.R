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
  z_value    <- w/w_se
  p_value    <- 2 * (1 - pnorm(abs(z_value)))
  t_hat      <- as.matrix(phi) %*% w 
  e          <- t_vector - t_hat
  e_st       <- e/sigma_se
  for(i in 1){
    
  } 
  return(list(w = w, sigma = sigma, variance = variance, w_se = w_se, sigma_se = sigma_se,
              z_value = z_value, p_value = p_value, t_hat = t_hat, e = e, e_st = e_st))
}

mle_plots <- function(mle_results){
  ci_plot_data <- data.frame(coeff_number = seq(1:nrow(mle_results$w))) %>% 
                  bind_cols(data.frame(w = mle_results$w)) %>% 
                  mutate(lower_bound = w - 1.96 * mle_results$w_se) %>% 
                  mutate(upper_bound = w + 1.96 * mle_results$w_se) %>% 
                  mutate(color_ci = ifelse(mle_results$p_value < 0.05, "green", "red"))
  
  ci_plot <- ggplot(data = ci_plot_data, aes(x = as.factor(coeff_number), y = w)) + 
             geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
                           color = w_plot_data$color_ci, width = 0.3) + 
             geom_point() + ggtitle("w estimation and confidence intervals") +
             labs(x = "Coefficient", y = "w, w +/- 1.96 s.e.") + theme_economist()
  
  st_res_plot_data <- data.frame(e_st = mle_results$e_st) %>% 
                       bind_cols(data.frame(t_hat = mle_results$t_hat)) 
                       
  st_res_plot <- ggplot(data = st_res_plot_data, aes(x = t_hat, y = e_st)) + 
                 geom_point()
  st_res_plot
                       
  
  
}
summary(lm(as.matrix(t_vector)~as.matrix(phi)-1))
