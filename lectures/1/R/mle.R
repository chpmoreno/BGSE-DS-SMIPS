mle_estimator_lm <- function (w, phi, t_vector) {
    if(is.vector(t_vector) == FALSE)
      t_vector <- as.vector(t_vector)
    if(is.matrix(phi) == FALSE)
      phi <- as.matrix(phi)
    M <- ncol(phi) 
    phi_w <- phi %*% w[1:M]
    Sig <- w[(M + 1):(M + 1)]
    sum(-(1 / 2) * log(2 * pi) - (1 / 2) * log(Sig ^ 2) - (1 / (2 * Sig ^ 2)) * (t_vector - phi_w) ^ 2)
}