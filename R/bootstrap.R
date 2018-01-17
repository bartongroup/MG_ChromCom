# Fit one bootstrap
# Rscript bootstrap.R set root batch switch ncells ntry t0

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
stopifnot(length(args) == 7)
set <- args[1]
root <- args[2]
batch <- args[3]
t2ref <- as.integer(args[4])
stopifnot(!is.null(dataFile[[set]]))
stopifnot(t2ref %in% c(0, 1))

ncells <- args[5]
ntry <- args[6]
t0 <- as.numeric(args[7])


print(paste("Fitting", set, "batch", batch))

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
fit <- fitChr(echr, pars, free, ncells=ncells, ntry=ntry, ncores=8, bootstrap=TRUE)

p <- t(as.matrix(do.call(c, fit$pars)))
file <- paste0(bootDir, root, "_", set, "_", batch, ".pars")
write.table(p, file=file, col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)
