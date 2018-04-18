# Fit one bootstrap
# Rscript bootstrap.R infile outfile switch ncells ntry t0 npars

topDir <- "/cluster/gjb_lab/mgierlinski/projects/chromcomR/"
source(paste0(topDir, "/R/setup.R"))
source(paste0(topDir, "/R/lib.R"))

args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]
t2ref <- as.integer(args[3])
stopifnot(file.exists(infile))
stopifnot(t2ref %in% c(0, 1))

ncells <- args[4]
ntry <- args[5]
t0 <- as.numeric(args[6])
npars <- as.integer(args[7])

print(paste("Fitting", infile))

echr <- experimentalData(infile)
str(echr)

pars <- c3pars(
  t0 = t0,
  tau1 = 15,
  tau2 = 8,
  k1 = 0.05,
  k2 = 0.06,
  t2ref = t2ref
)
str(pars)

if(npars == 4) {
  free <- c("tau1", "k1", "k2", "tau2")
} else if(npars == 3) {
  free <- c("tau1", "k1", "k2")
} else {
  stop("Wrong number of parameters")
}



fit <- fitChr(echr, pars, free, ncells=ncells, ntry=ntry, ncores=8, bootstrap=TRUE)

p <- t(as.matrix(do.call(c, fit$pars)))
write.table(p, file=outfile, col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)
