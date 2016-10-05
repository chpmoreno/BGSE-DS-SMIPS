# load libraries and functions needed
library(readr)
library(dplyr)
library(ggplot2)
source("lectures/1/R/mle.R")

# load data
data_exercise <- read_delim(file = "../Datasets/synthetic_regression/synthetic_regression.txt", 
                            delim = " ")
