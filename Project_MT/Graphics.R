library(ggplot2)
library(reshape2)
library(dplyr)
library(ggthemes)
require(gridExtra)

hypothesis_test_graphs <- function(m){
  alpha    <- seq(0, 1, 0.01)
  f1_alpha <- 1 - (1 - alpha) ^ (1/m)
  f2_alpha <- alpha / m
  f3_alpha <- alpha
  
  data_functions <- data.frame(alpha    = alpha,
                               f1_alpha = f1_alpha,
                               f2_alpha = f2_alpha,
                               f3_alpha = f3_alpha)
  data_functions <- reshape2::melt(data_functions, id.vars = alpha)
  greeks <- list(bquote(f[1](alpha) == 1 - (1 - alpha) ^ (1/m)), 
                 bquote(f[2](alpha) == alpha/m), 
                 bquote(f[2](alpha) == alpha))
  
  plot_tot <- ggplot(data = dplyr::filter(data_functions, variable == "f1_alpha"), 
                     aes(x = alpha, y = value, color = variable)) + 
              geom_line(size = 0.8) +
              geom_line(data = dplyr::filter(data_functions, variable == "f2_alpha"), 
                        aes(x = alpha, y = value, color = variable), size = 0.8) +
              geom_line(data = dplyr::filter(data_functions, variable == "f3_alpha"), 
                        aes(x = alpha, y = value, color = variable), size = 0.8) +
              scale_colour_manual(values = c("limegreen", "blue", "red"),
                                  labels = greeks,
                                  guide = guide_legend(title = NULL,
                                                       override.aes = list(linetype = c("solid", "solid", "solid")))) +
              labs(x = bquote(alpha), y = bquote(f(alpha))) +
              ggtitle(bquote(m == .(m))) +
              theme_economist()
  return(plot_tot)
}

grid.arrange(hypothesis_test_graphs(1), hypothesis_test_graphs(5),
             hypothesis_test_graphs(50), hypothesis_test_graphs(100),
             hypothesis_test_graphs(500), hypothesis_test_graphs(1000),
             ncol=2)