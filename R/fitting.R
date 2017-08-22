setwd("/cluster/gjb_lab/mgierlinski/projects/chromcomR/doc")

source("../../mylib/R/lib.R")
source("../R/lib.R")

library(ggplot2)
library(gridExtra)
library(reshape2)
library(parallel)
library(methods)

binDir <- "../RData"


echr.scr <- experimentalData(dataFile$scramble)
echr.d2 <- experimentalData(dataFile$NCAPD2)
echr.d3 <- experimentalData(dataFile$NCAPD3)
echr.smc <- experimentalData(dataFile$SMC2)


pars <- c3pars()
free <- c("tau", "k1", "k2", "dt2")

chr.scr <- cacheData("fitv_scramble_n10000_tau_k1_k2_dt2", fitChr, echr.scr, pars, free, nsim=10000, ntry=100, ncores=4, binDir=binDir)
chr.d2 <- cacheData("fitv_NCAPD2_n10000_tau_k1_k2_dt2", fitChr, echr.d2, pars, free, nsim=10000, ntry=100, ncores=4, binDir=binDir)
chr.d3 <- cacheData("fitv_NCAPD3_n10000_tau_k1_k2_dt2", fitChr, echr.d3, pars, free, nsim=10000, ntry=100, ncores=4, binDir=binDir)
chr.scm <- cacheData("fitv_SMC2_n10000_tau_k1_k2_dt2", fitChr, echr.smc, pars, free, nsim=10000, ntry=100, ncores=4, binDir=binDir)
