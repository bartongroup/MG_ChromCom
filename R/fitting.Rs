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
stopifnot(length(args) < 3)
set <- args[1]
name <- args[2]
t2ref <- args[3]
stopifnot(!is.null(dataFile[[set]]))
stopifnot(t2ref %in% c(0, 1))

nsim <- 1000
if(!is.null(args[4])) nsim <- args[4]

ntry <- 100
if(!is.null(args[5])) ntry <- args[5]


print(paste("Fitting", set))

echr <- experimentalData(dataFile[[set]])
str(echr)

pars <- c3pars(
  t0 = 0,
  tau1 = 15,
  tau2 = 10,
  k1 = 0.05,
  k2 = 0.06,
  t2ref = t2ref
)
free <- c("tau1", "k1", "k2", "tau2")
chr <- cacheData(name, fitChr, echr, pars, free, nsim=nsim, ntry=ntry, ncores=4, binDir=binDir)
