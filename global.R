##The global.R script
load("Data/LM22.RData")
load("Data/TIL10.RData")
library(parallel)
library(ggplot2)
library(openxlsx)
.num.cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
