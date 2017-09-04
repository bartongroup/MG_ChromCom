setwd("/cluster/gjb_lab/mgierlinski/projects/chromcomR/doc")
source("../../mylib/R/lib.R")
source("../R/lib.R")
binDir <- "../RData"

library(ggplot2)
library(gridExtra)
library(reshape2)
library(parallel)
library(methods)

args <- commandArgs(TRUE)
stopifnot(length(args) == 2)
set <- args[1]
name <- args[2]
stopifnot(!is.null(dataFile[[set]]))

print(paste("Fitting", set))

echr <- experimentalData(dataFile[[set]])
str(echr)

pars <- c3pars(
  t0 = 0,
  tau1 = 15,
  tau2 = 10,
  k1 = 0.05,
  k2 = 0.06
)
free <- c("tau1", "k1", "k2", "tau2")
chr <- cacheData(name, fitChr, echr, pars, free, nsim=10000, ntry=100, ncores=4, binDir=binDir)
