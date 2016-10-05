library(readr)
library(dplyr)
library(ggplot2)


mle_estimator_lm <- function (w, phi, t) {
    if(is.vector(t) == FALSE)
      t <- as.vector(t)
    if(is.matrix(phi) == FALSE)
      phi <- as.matrix(phi)
    M <- ncol(phi) 
    phi_w <- phi %*% w[1:M]
    Sig <- par[(M + 1):(M + 1)]
    sum(-(1 / 2) * log(2 * pi) - (1 / 2) * log(Sig ^ 2) - (1 / (2 * Sig ^ 2)) * (t - phi_w) ^ 2)
}