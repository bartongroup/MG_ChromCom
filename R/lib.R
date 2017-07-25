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
#' @param colours Colour names
#'
#' @return A \code{ChrCom3} object.
#' @export
ChromCom3 <- function(pars, time=NULL, cells=NULL, colours = c("B", "P", "R")) {
  stopifnot(is(pars, "c3pars"))
  if(is.null(time)) {
    start <- -140
    stop <- 90
    step <- 1
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
#' @param t1
#' @param dt2
#' @param dt3
#' @param k1
#' @param k2
#' @param k3
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
#' @param t Time point in seconds
#' @param tp A list of start, stop and step for the timeline
#'
#' @return Integer index in the time vector corresponding to t;
timeIndex <- function(t, tp, maxn) {
  i <- round((t - tp$start) / tp$step) + 1
  if(i > maxn) i <- maxn
  return(i)
}

#' Generate timeline for one cell
#'
#' @param pars Model parameters
#' @param timepars A list of start, stop and step for the timeline
#'
#' @return A character vector of colours
#' @export
timelineCell <- function(pars, timepars) {
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
#'
#' @return A \code{ChrCom3} object with simulation results.
#' @export
generateCells <- function(chr, nsim=1000) {
  cells <- t(replicate(nsim, timelineCell(chr$pars, chr$timepars)))
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
    geom_vline(xintercept=0, color="grey20", linetype=2) +
    if(!single) facet_grid(X ~ Y)
}

#' Plot time lines
#'
#' @param chr A \code{ChrCom3} object with data
#' @param smooth If TRUE, smoothing will be applied
#' @param k Window size for smoothing
#' @param expdata Experimental data to add to the plot (another \code{ChrCom3} object)
#' @param ... Other parameters passed to \code{\link{timelinePanel}}
#'
#' @export
plotTimelines <- function(chr, smooth=FALSE, k=5, expdata=NULL, title='', ...) {
  m <- meltTimelines(chr, smooth=smooth, k=k)
  g <- timelinePanel(m, single=TRUE, ...)
  if(!is.null(expdata)) {
    exm <- meltTimelines(expdata, smooth=TRUE, k=15)
    g <- g + geom_line(data=exm, aes(colour=Colour), size=0.2) + labs(title=title)
  }
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
#' @param chr Model data
#' @param echr Experimental data
#'
#' @return RMS
#' @export
RMS <- function(chr, echr) {
  rms <- 0
  for(col in chr$colours) {
    E <- ts(chr$cnt[[col]] / chr$cnt$total, start=chr$timepars$start, deltat=chr$timepars$step)
    O <- ts(echr$cnt[[col]] / echr$cnt$total, start=echr$timepars$start, deltat=echr$timepars$step)
    E[which(is.nan(E) | is.infinite(E))] <- NA
    O[which(is.nan(O) | is.infinite(O))] <- NA
    chi2 <- (O - E)**2 / E
    chi2[which(is.infinite(chi2))] <- NA
    rms <- rms + sum(chi2, na.rm=TRUE)
  }
  rms
}
