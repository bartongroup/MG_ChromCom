#' @import ggplot2
#' @import caTools
#' @import reshape2
#' @import gridExtra


data.colours <- c("", "_blue_", "_blueDark_g", "_blueDark_r",
                  "_brown_g", "_brown_r", "_brownDark_g", "_brownDark_r",
                  "_pink_", "_pinkDark_", "anaphase")
model.colours <- c(NA, "B", "B", "B", "B", "B", "B", "B", "P", "R", NA)
translateVector <- function(x) {
  model.colours[match(x, data.colours)]
}


#' Data files
#' @export
dataFile <- list(
  scramble = "../data/scramble_noncoloured.csv",
  NCAPD2 = "../data/NCAPD2_uncoloured.csv",
  NCAPD3 = "../data/NCAPD3_uncoloured.csv"
)

#' Simple theme for plotting
#' @export
simple_theme_grid <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(colour = "grey95"),
    axis.line = ggplot2::element_line(colour = "black")
  )



#' \code{ChrCom3} object constructor
#'
#'
#' \code{ChromCom3} will create an empty object with model parameters. Optionally, data can be added.
#'
#' @param pars A list of model parameters
#' @param time An optional time vector
#' @param cells An optional matrix with data
#' @param timepars Parameters for the time sequence
#' @param colours Colour names
#'
#' @return A \code{ChrCom3} object.
#' @export
ChromCom3 <- function(pars, time=NULL, cells=NULL, timepars=list(start=-140, stop=90, step=1), colours = c("B", "P", "R")) {
  stopifnot(is(pars, "c3pars"))
  if(is.null(time)) {
    start <- timepars$start
    stop <- timepars$stop
    step <- timepars$step
    time = seq(from=start, to=stop, by=step)
  } else {
    start <- time[1]
    stop <- time[length(time)]
    step <- time[2] - time[1]
  }
  n <- length(time)

  if(is.null(cells)) {
    cnt <- NULL
  } else {
    if(is.null(time)) stop("Need time")
    cnt <- cellCount(cells, time, colours)
  }

  obj <- list(
    colours = colours,
    timepars = list(
      start = start,
      stop = stop,
      step = step,
      n = n
    ),
    time = time,
    pars = pars,
    cells = cells,
    cnt = cnt
  )
  class(obj) <- append(class(obj), "ChromCom3")
  return(obj)
}


#' \code{c3pars} object constructor
#'
#' @param t1 Start time
#' @param dt2 Delay until P->R is possible
#' @param dt3 Delay until P>B is possible
#' @param k1 B->P rate
#' @param k2 P->R rate
#' @param k3 R->B rate
#' @param dummy Logical, if set empty list is returned
#'
#' @return Model parameters object
#' @export
c3pars <- function(
  t1 = -30,
  dt2 = 0,
  dt3 = 0,
  k1 = 0.04,
  k2 = 0.04,
  k3 = 0,
  dummy = FALSE
) {
  if(dummy) {
    pars = list()
  } else {
    pars <- list(t1=t1, dt2=dt2, dt3=dt3, k1=k1, k2=k2, k3=k3)
  }
  class(pars) <- append(class(pars), "c3pars")
  return(pars)
}


#' Generate transition times
#'
#' @param pars Parameters of the simulation (t1, dt2, k1, k2)
#'
#' @return A list with two transition times
transitionTimes <- function(pars) {
  BP <- ifelse(pars$k1 > 0, pars$t1 + rexp(1, pars$k1), 1000)
  PR <- ifelse(pars$k2 > 0, BP + pars$dt2 + rexp(1, pars$k2), 1000)
  T <- list(
    BP = BP,
    PR = PR
  )
  return(T)
}

#' Index of a time point
#'
#' @param t Time point in minutes
#' @param tp A list of start, stop and step for the timeline
#'
#' @return Integer index in the time vector corresponding to t;
timeIndex <- function(t, tp, maxn=10000) {
  i <- round((t - tp$start) / tp$step) + 1
  if(i > maxn) i <- maxn
  return(i)
}



#' Create a cell using transition method
#'
#' @param pars Model parameters
#' @param timepars Time parameters
#'
#' @return A cell
tcTransition <- function(pars, timepars) {
  T <- transitionTimes(pars)

  maxn <- 10000
  iBP <- timeIndex(T$BP, timepars, maxn)
  iPR <- timeIndex(T$PR, timepars, maxn)

  # we want to avoid boundary effects
  cell <- rep("B", maxn)
  if(iPR > iBP) cell[(iBP+1):iPR] <- "P"
  cell[(iPR+1):maxn] <- "R"
  cell <- cell[1:timepars$n]

  return(cell)
}


stateTransition <- function(state, i, i2, i3, p1, p2, p3) {
  R <- runif(1)
  if(state == "B") {
    if(R < p1) state <- "P"
  }
  else if(state == "P") {
    if(i >= i3 && R < p3)         state <- "B"
    else if (i >= i2 && R < p2)   state <- "R"
  }
  return(state)
}


#' Create a cell using simulations method
#'
#' @param pars Model parameters
#' @param timepars Time parameters
#'
#' @return A cell
tcSimulation <- function(pars, timepars) {
  cell <- rep("-", timepars$n)
  i1 <- timeIndex(pars$t1, timepars)
  i2 <- timeIndex(pars$t1 + pars$dt2, timepars)
  i3 <- timeIndex(pars$t1 + pars$dt3, timepars)
  p1 <- 1 - exp(-timepars$step * pars$k1)
  p2 <- 1 - exp(-timepars$step * pars$k2)
  p3 <- 1 - exp(-timepars$step * pars$k3)

  cell[1:(i1-1)] <- "B"
  state <- "B"
  for(i in i1:timepars$n) {
    state <- stateTransition(state, i, i2, i3, p1, p2, p3)
    cell[i] <- state
  }

  return(cell)
}





#' Generate timeline for one cell
#'
#' @param pars Model parameters
#' @param timepars A list of start, stop and step for the timeline
#'
#' @return A character vector of colours
#' @export
timelineCell <- function(pars, timepars, method=c("transition", "simulation")) {
  method <- match.arg(method)

  if(method == "transition") {
    cell <- tcTransition(pars, timepars)
  } else if(method == "simulation") {
    cell <- tcSimulation(pars, timepars)
  }

  return(cell)
}



#' Count colours across cells
#'
#' @param cells A matrix with cell colours
#' @param time Time vector
#' @param colours Factor with colours
#'
#' @return Marginal count
#' @export
cellCount <- function(cells, time, colours) {
  cnt <- apply(cells, 2, function(x) table(factor(x, levels=colours)))
  cnt <- as.data.frame(t(cnt))
  cnt$total <- colSums(!is.na(cells))
  cnt$Time <- time
  cnt
}


#' Perform simulation
#'
#' @param chr Initial \code{ChrCom3} object with model parameters.
#' @param nsim Number of simulations.
#' @param method Either "transition" or "simulation"
#'
#' @return A \code{ChrCom3} object with simulation results.
#' @export
generateCells <- function(chr, nsim=1000, method="transition") {
  cells <- t(replicate(nsim, timelineCell(chr$pars, chr$timepars, method)))
  cnt <- cellCount(cells, chr$time, chr$colours)

  chr$cells <- cells
  chr$cnt <- cnt
  chr$nsim <- nsim

  return(chr)
}


smoothColours <- function(cnt, colours, k) {
  cnt[,colours] <- apply(cnt[,colours], 2, function(x) caTools::runmean(x, k))
  return(cnt)
}

#' Melt timelines for plotting
#'
#' @param chr A \code{ChrCom3} object with data
#' @param label1 Name of the first variable in the grid
#' @param label2 Name of the second variable in the grid
#' @param smooth If TRUE, smoothing will be applied
#' @param k Smoothing window size
#'
#' @return Melted data frame
#' @export
meltTimelines <- function(chr, label1="L1", label2="L2", smooth=FALSE, k=5) {
  cnt <- chr$cnt
  colours <- chr$colours
  cnt[,colours] <- cnt[,colours] / cnt$total
  if(smooth) cnt <- smoothColours(cnt, colours, k)
  m <- reshape2::melt(cnt, id.vars="Time", measure.vars=colours, variable.name="Colour", value.name="Count")
  m$X <- label1
  m$Y <- label2
  return(m)
}

#' Plot panel with time lines
#'
#' @param m Melted data frame
#' @param single If true, no facets are applied
#' @param xmin Lower limit on x-axis
#' @param xmax Upper limit on x-axis
#'
#' @export
timelinePanel <- function(m, single=FALSE, xmin=as.numeric(NA), xmax=as.numeric(NA)) {
  cPalette <- c("blue", "pink", "red")
  ggplot(m, aes(x=Time, y=Count)) +
    simple_theme_grid +
    geom_line(aes(colour=Colour), size=1.5) +
    scale_colour_manual(values=cPalette) +
    theme(legend.position="none") +
    xlim(xmin, xmax) +
    labs(x="Time (min)", y="Proportion") +
    geom_vline(xintercept=0, color="grey20", linetype=1) +
    if(!single) facet_grid(X ~ Y)
}

#' Plot time lines
#'
#' @param chr A \code{ChrCom3} object with data
#' @param smooth If TRUE, smoothing will be applied
#' @param k Window size for smoothing
#' @param expdata Experimental data to add to the plot (another \code{ChrCom3} object)
#' @param title Title of the plot
#' @param withpars Title will be replaced with model parameters
#' @param ... Other parameters passed to \code{\link{timelinePanel}}
#'
#' @export
plotTimelines <- function(chr, smooth=FALSE, k=5, expdata=NULL, title='', withpars=FALSE, ...) {
  m <- meltTimelines(chr, smooth=smooth, k=k)
  if(withpars) {
    title <- paste0(lapply(names(chr$pars), function(name) {
      paste0(name, "=", sprintf("%.3g", chr$pars[[name]]))
    }
    ), collapse=", ")
  }
  g <- timelinePanel(m, single=TRUE, ...) +
    geom_vline(xintercept=chr$pars$t1, colour="skyblue", linetype=2) +
    geom_vline(xintercept=(chr$pars$t1+chr$pars$dt2), colour="salmon", linetype=2) +
    geom_vline(xintercept=(chr$pars$t1+chr$pars$dt3), colour="brown", linetype=2)
  if(!is.null(expdata)) {
    exm <- meltTimelines(expdata, smooth=TRUE, k=15)
    rms <- oeError(chr, expdata)
    if(withpars) title <- paste0(title, sprintf(", rms=%.3g", rms))
    g <- g + geom_line(data=exm, aes(colour=Colour), size=0.2)
  }
  g <- g + labs(title=title) + theme(plot.title = element_text(size=10))
  g
}


#' Plot cell colour map
#'
#' @param chr A \code{ChrCom3} object with data
#'
#' @export
plotCells <- function(chr) {
  cells <- chr$cells
  colnames(cells) <- chr$time
  m <- reshape2::melt(cells, varnames=c("Cell", "Time"))
  cPalette <- c("blue", "pink", "red")
  ggplot(m, aes(Time, Cell, fill=value)) +
    simple_theme_grid +
    geom_tile() +
    scale_fill_manual(values=cPalette) +
    theme(legend.position="none") +
    labs(x="Time (min)", y="Cell")
}




#' Read experimental data from CSV file
#'
#' @param file File name
#'
#' @return A \code{ChrCom3} object
#' @export
experimentalData <- function(file) {
  dat <- read.delim(file, header=TRUE, sep=",")
  time <- dat[,1] / 60
  dat <- dat[,2:ncol(dat)]
  # the last non-empty time point
  cut <- max(which(s <- rowSums(dat != '') > 0))
  time <- time[1:cut]
  dat <- dat[1:cut,]
  tdat <- apply(dat, 1, translateVector)
  pars <- c3pars()
  echr <- ChromCom3(pars, time=time, cells=tdat)
}


#' Data-model error
#'
#' Using test of independence between observed and model counts in the limits window.
#'
#' @param chr Model data
#' @param echr Experimental data
#' @param limits Limiting window for calculations.
#'
#' @return RMS over all three curves
#' @export
oeError <- function(chr, echr, limits=c(-90, 30)) {

  getProp <- function(ch, col) {
    x <- ts(ch$cnt[[col]] / ch$cnt$total, start=ch$timepars$start, deltat=ch$timepars$step)
    x[which(is.nan(x) | is.infinite(x))] <- NA
    x <- window(x, start=limits[1], end=limits[2])
  }

  rms <- 0
  for(col in chr$colours) {
    xo <- getProp(echr, col)
    xe <- getProp(chr, col)
    rms <- rms + sqrt(sum((xo - xe)^2))
  }
  return(rms)
}


# Not very good
oeError.chi2 <- function(chr, echr, limits=c(-90, 30)) {
  M <- NULL
  for(col in chr$colours) {
    xo <- ts(echr$cnt[[col]] , start=echr$timepars$start, deltat=echr$timepars$step)
    xe <- ts(chr$cnt[[col]] , start=chr$timepars$start, deltat=chr$timepars$step)
    xo <- window(xo, start=limits[1], end=limits[2])
    xe <- window(xe, start=limits[1], end=limits[2])
    X <- rbind(xo, xe)
    X <- X[,which(colSums(X) > 0)]
    M <- cbind(M, X)
  }
  xi <- chisq.test(M)
  redchi <- xi$statistic / xi$parameter
  return(redchi)
}


# convert a vector of parameters used by optim into c3pars object
vectorPar <- function(p) {
  pars <- c3pars(
    t1 = p[1],
    k1 = p[2],
    k2 = p[3]
  )
  if(length(p) > 3) pars[["dt2"]] <- p[4]
  return(pars)
}

# convert a c3pars object into a vector of parameters used by optim
parVector <- function(pars, n=3) {
  p <- c(pars$t1, pars$k1, pars$k2)
  if(n > 3) p <- c(p, pars$dt2)
  return(p)
}

# error function to minimize; p is a vector of parameters
errorFun <- function(p, echr, nsim) {
  pars <- vectorPar(p)
  chr <- ChromCom3(pars, timepars=list(start=-90, stop=30, step=1))
  chr <- generateCells(chr, nsim=nsim, method="simulation")
  err <- oeError(chr, echr)
  return(err)
}

#' Fit ChromCom data with a model
#'
#' Fits data with a model in which t1, k1, k2 and dt2 are free parameters. If
#' npar=3, only first three of them are used.
#' @param echr A \code{ChromCom3} object with experimental data
#' @param pars A \code{c3pars} object with initial parameters
#' @param npar Number of parameters to use in minimization (3 or 4)
#' @param nsim Number of cells to simulate
#' @param ntry Number of tries in search
#' @param ncores Number of cores
#'
#' @return A \code{ChromCom3} object with the best-fitting model. "rms" field
#'   is added. \code{\link{optim}} function is used for minimization with method "L-BFGS-B".
#' @export
fitChr <- function(echr, pars, npar=3, nsim=1000, ntry=10, ncores=4) {
  stopifnot(is(pars, "c3pars"))

  #chr <- ChromCom3(pars)
  p <- parVector(pars, n=npar)
  lower <- c(-50, 0, 0, 0)
  upper <- c(0, 0.2, 0.2, 30)
  lower <- lower[1:npar]
  upper <- upper[1:npar]

  lopt <- mclapply(1:ntry, function(i) {
    optim(p, errorFun, gr=NULL, echr, nsim, method="L-BFGS-B", lower=lower, upper=upper, control=list(trace=3))
  }, mc.cores = ncores)

  idx.min <- which.min(sapply(1:ntry, function(i) lopt[[i]]$value))
  pars <- vectorPar(lopt[[idx.min]]$p)
  chr <- ChromCom3(pars, timepars=list(start=-90, stop=30, step=1))
  chr <- generateCells(chr, nsim=nsim, method="simulation")
  chr$rms <- lopt[[idx.min]]$value
  return(chr)
}

