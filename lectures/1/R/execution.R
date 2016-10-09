# load libraries and functions needed
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
source("lectures/1/R/mle.R")

# initial parameters ####
m <- 30 # number of variables taken for the model

# load data - Find it on https://1drv.ms/t/s!Ai0XbELt7PXquFJtWenDbxvksoXM ####
data_exercise <- read.table(file = "../Datasets/synthetic_regression/synthetic_regression.txt",
                                  nrow = 300)[,1:(m + 1)]

# vector t (t_vector) and matrix phi
t_vector <- as.vector(data_exercise[ , "t"])
phi      <- cbind(rep(1, length(t_vector)), 
                  data_exercise[,-which(names(data_exercise) %in% c("t"))])
colnames(phi)[1] <- "const"

# mle_estimation
mle_estimation_result <- optim(runif(m + 2, 0, 1), mle_estimator_lm, phi = phi, 
                               t_vector = t_vector, method = "BFGS", 
                               control = list(trace = 1, maxit = 10000, fnscale = -1),
                               hessian = TRUE)
# mle_results
mle_results <- mle_result(mle_estimation_result, t_vector, phi)

# mle graphics
mle_graphics <- mle_plots(mle_results)

pdf(file = "lectures/1/Graphics/ci_plot.pdf", height = 7, width = 9)
  print(mle_graphics$ci_plot)
dev.off()

pdf(file = "lectures/1/Graphics/st_res_plot.pdf", height = 7, width = 9)
  print(mle_graphics$st_res_plot)
dev.off()

pdf(file = "lectures/1/Graphics/qq_plot.pdf", height = 7, width = 9)
print(mle_graphics$qq_plot)
dev.off()