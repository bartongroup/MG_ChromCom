# Rscript fitting.Rs set name switch ncells ntry t0 npars [tau2]

topDir <- "/cluster/gjb_lab/mgierlinski/projects/chromcomR/"

source(paste0(topDir, "/R/setup.R"))
source(paste0(topDir, "/R/lib.R"))

dataDir <- paste0(topDir, "data/")
fitDir <- paste0(topDir, "fits/")


args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(7, 8))
set <- args[1]
name <- args[2]
t2ref <- as.integer(args[3])
stopifnot(!is.null(dataFile[[set]]))
stopifnot(t2ref %in% c(0, 1))

ncells <- args[4]
ntry <- args[5]
t0 <- as.numeric(args[6])
npars <- as.integer(args[7])

if(length(args) == 7) {
  tau2 <- 8
} else {
  tau2 <- as.numeric(args[8])
}

print(paste("Fitting", set))

echr <- experimentalData(paste0(dataDir, dataFile[[set]]))
str(echr)

pars <- c3pars(
  t0 = t0,
  tau1 = 15,
  tau2 = tau2,
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

chr <- fitChr(echr, pars, free, ncells=ncells, ntry=ntry, ncores=8)
saveRDS(chr, file=paste0(fitDir, name, ".rds"))
