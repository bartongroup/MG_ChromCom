# Rscript fitting.Rs set name switch ncells ntry t0

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
stopifnot(length(args) >= 3)
set <- args[1]
name <- args[2]
t2ref <- as.integer(args[3])
stopifnot(!is.null(dataFile[[set]]))
stopifnot(t2ref %in% c(0, 1))

ncells <- 1000
if(!is.null(args[4])) ncells <- args[4]

ntry <- 100
if(!is.null(args[5])) ntry <- args[5]

t0 <- 0
if(!is.null(args[6])) t0 <- as.numeric(args[6])

print(paste("Fitting", set))

echr <- experimentalData(dataFile[[set]])
str(echr)

pars <- c3pars(
  t0 = t0,
  tau1 = 15,
  tau2 = 10,
  k1 = 0.05,
  k2 = 0.06,
  t2ref = t2ref
)
str(pars)

free <- c("tau1", "k1", "k2", "tau2")
chr <- cacheData(name, fitChr, echr, pars, free, ncells=ncells, ntry=ntry, ncores=8, binDir=binDir)
