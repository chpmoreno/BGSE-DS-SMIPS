# Author: José Fernado Moreno Gutiérrez

# load libraries and functions needed
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(stringr)

#load # load data - Find it on  ####

data_exercise <- read.table(file = "../Datasets/curve_data.txt")

data_plot <- ggplot(data = data_exercise, aes(x = x, y = t)) + 
             geom_point() + ggtitle("t / x scatter plot") + 
             theme_economist()

phix <- function(x, M, basis) {
  if(basis == "poly"){
    X <- matrix(NA, ncol = M + 1, nrow = length(x))
    for(i in 0:M) {
      X[ , i + 1] <- x ^ i
    }
  }
  if(basis == "Gauss") {
    means <- seq(0, 1, by = 1/M)
    X <- matrix(NA, ncol = M, nrow = length(x))
    for(i in 1:M) {
      X[, i] <- exp((-(x - means[i]) ^ 2) / 0.2)
    }
  }
  return(X)
}

post.params <- function(training_data, M, basis_type, phix, lambda, q){
  t_vector <- training_data$t
  x        <- training_data$x
  X        <- phix(x, M, basis_type)
  W        <- solve(lambda * diag(1, ncol(X), ncol(X)) + (t(X) %*% X)) %*% t(X) %*% t_vector
  Q        <- q * (lambda * diag(1, ncol(X), ncol(X)) + (t(X) %*% X))
}

