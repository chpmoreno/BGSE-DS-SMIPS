# Author: José Fernado Moreno Gutiérrez

# load libraries and functions needed
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(stringr)
source("lectures/2/R/functions.R")

#load # load data - Find it on  ####

data_exercise <- read.table(file = "../Datasets/curve_data.txt")

plot_data <- ggplot(data = data_exercise, aes(x = x, y = t)) + 
             geom_point() + ggtitle("t / x scatter plot") + 
             theme_economist()
plot_data

x_plot          <- seq(min(data_exercise$x), max(data_exercise$x), by = 1/1000)
w_poly_plot     <- post.params(data_exercise, 9, "poly", phix, 2, (1/0.1) ^ 2)$w
t_poly_hat_plot <- phix(x_plot, M = 9, "poly") %*% w_poly_plot
w_gauss_plot     <- post.params(data_exercise, 9, "gauss", phix, 2, (1/0.1) ^ 2)$w
t_gauss_hat_plot <- phix(x_plot, M = 9, "gauss") %*% w_gauss_plot

data_plot_model_poly  <- cbind(as.data.frame(x_plot), as.data.frame(t_poly_hat_plot), rep("poly", length(x_plot)))
data_plot_model_gauss <- cbind(as.data.frame(x_plot), as.data.frame(t_gauss_hat_plot), rep("gauss", length(x_plot)))
data_plot_real        <- cbind(data_exercise, rep("observed", nrow(data_exercise)))

colnames(data_plot_model_poly)  <- c("x", "t", "basis")
colnames(data_plot_model_gauss) <- c("x", "t", "basis")
colnames(data_plot_real)        <- c("x", "t", "basis")

data_plot_tot       <- bind_rows(data_plot_model_poly, data_plot_model_gauss, data_plot_real)

plot_tot <- ggplot(data = filter(data_plot_tot, basis == "observed"), aes(x = x, y = t, color = basis)) + 
            geom_point(size = 3) +
            geom_line(data = filter(data_plot_tot, basis == "poly"), aes(x = x, y = t, color = basis),
                      size = 0.8) +
            geom_line(data = filter(data_plot_tot, basis == "gauss"), aes(x = x, y = t, color = basis),
                      size = 0.8) +
            scale_colour_manual(values = c("limegreen", "blue", "red"),
                                guide = guide_legend(title = NULL,
                                                     override.aes = list(linetype = c("solid", "blank", "solid")))) +
            ggtitle(expression(paste("t / x observed and estimated (M = 9, ", delta, " = 2, q = 100)"))) + 
            theme_economist()
plot_tot

pdf(file = "lectures/2/Graphics/plot_data.pdf", height = 7, width = 9)
  print(plot_data)
dev.off()

pdf(file = "lectures/2/Graphics/plot_estimations.pdf", height = 7, width = 9)
  print(plot_tot)
dev.off()
