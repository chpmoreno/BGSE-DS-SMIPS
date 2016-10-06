# load libraries and functions needed
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
source("lectures/1/R/mle.R")

# initial parameters ####
m <- 30 # number of variables taken for the model

# load data ####
data_exercise <- read.table(file = "../Datasets/synthetic_regression/synthetic_regression.txt",
                                  nrow = 300)[,1:(m+1)]

# vector t (t_vector) and matrix phi
t_vector <- as.vector(data_exercise[ , "t"])
phi      <- cbind(rep(1, length(t_vector)), 
                  data_exercise[,-which(names(data_exercise) %in% c("t"))])
colnames(phi)[1] <- "const"

mle_estimation_result <- optim(runif(m + 2, 0, 1), mle_estimator_lm, phi = phi, 
                               t_vector = t_vector, method = "BFGS", 
                               control = list(trace = 1, maxit = 10000, fnscale = -1),
                               hessian = TRUE)

mle_results <- mle_result(mle_estimation_result, t_vector, phi)





