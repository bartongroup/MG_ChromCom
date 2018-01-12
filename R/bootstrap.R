# Fit one bootstrap
# Rscript bootstrap.R set batch switch ncells ntry t0

setwd("/cluster/gjb_lab/mgierlinski/projects/chromcomR/doc")
source("../../mylib/R/lib.R")
source("../R/lib.R")
bootDir <- "../bootstrap/"

library(ggplot2)
library(gridExtra)
library(reshape2)
library(parallel)
library(methods)

args <- commandArgs(TRUE)
stopifnot(length(args) != 6)
set <- args[1]
batch <- args[2]
t2ref <- as.integer(args[3])
stopifnot(!is.null(dataFile[[set]]))
stopifnot(t2ref %in% c(0, 1))

ncells <- args[4]
ntry <- args[5]
t0 <- as.numeric(args[6])


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
fit <- fitChr(echr, pars, free, ncells=ncells, ntry=ntry, ncores=8, bootstrap=TRUE, binDir=binDir)

p <- t(as.matrix(do.call(c, fit$pars)))
file <- paste0(bootDir, set, "_", batch, ".pars")
write.table(p, file=file, col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)
