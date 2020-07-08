##The global.R script
load("Data/LM22.Rdata")
load("Data/TIL10.Rdata")
library(parallel)
library(ggplot2)
.num.cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
