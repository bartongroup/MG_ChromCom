#' \code{ChrCom3} object constructor
#'
#' @param pars A list of model parameters
#'
#' @return A \code{ChrCom3} object.
ChromCom3 <- function(pars) {
  start <- -80
  stop <- 60
  step <- 1
  time = seq(from=start, to=stop, by=step)
  n <- length(time)

  obj <- list(
    colours = c("BB", "P", "R"),
    timepars = list(
      start = start,
      stop = stop,
      step = step,
      n = n
    ),
    time = time,
    pars = pars,
    cells = NULL,
    cnt = NULL
  )
  class(obj) <- append(class(obj), "ChromCom3")
  return(obj)
}


#' Generate transition times
#'
#' @param pars Parameters of the simulation
#'
#' @return A list with two transition times
transitionTimes <- function(pars) {
  BBP <- pars$t1 + rexp(1, pars$r1)
  if(is.null(pars$t2)) {
    PR <- BBP + rexp(1, pars$r2)
  } else {
    PR <- pars$t2 + rexp(1, pars$r2)
    if(PR < BBP) PR <- BBP
  }
  T <- list(
    BBP = BBP,
    PR = PR
  )
  return(T)
}

#' Index of a time point
#'
#' @param t Time point in seconds
#' @param tp A list of start, stop and step for the timeline
#'
#' @return Integer index in the time vector corresponding to t; returns n if index outside range
timeIndex <- function(t, tp) {
  i <- round((t - tp$start) / tp$step) + 1
  if(i > tp$n) i <- tp$n
  return(i)
}

#' Generate timeline for one cell
#'
#' @param pars Model parameters
#' @param timepars A list of start, stop and step for the timeline
#'
#' @return A character vector of colours
timelineCell <- function(pars, timepars) {
  T <- transitionTimes(pars)

  iBBP <- timeIndex(T$BBP, timepars)
  iPR <- timeIndex(T$PR, timepars)

  cell <- rep('BB', timepars$n)
  if(iBBP <= iPR && iPR < timepars$n) cell[iBBP:iPR] <- 'P'
  if(iPR < timepars$n) cell[(iPR+1):timepars$n] <- 'R'

  return(cell)
}


#' Perform simulation
#'
#' @param chr Initial \code{ChrCom3} object with model parameters.
#' @param nsim Number of simulations.
#'
#' @return A \code{ChrCom3} object with simulation results.
generateCells <- function(chr, nsim=1000) {
  cells <- lapply(1:nsim, function(i) timelineCell(chr$pars, chr$timepars))
  cells <- do.call(rbind, cells)

  cnt <- apply(cells, 2, function(x) table(factor(x, levels=chr$colours)))
  cnt <- as.data.frame(t(cnt))
  cnt$Time <- chr$time
  
  chr$cells <- cells
  chr$cnt <- cnt
  
  return(chr)
}

plotTimelines <- function(chrcom) {
  m <- reshape2::melt(chrcom$cnt, id.vars="Time", variable.name="Colour", value.name="Count")
  cPalette <- c("blue", "pink", "red")
  ggplot(m, aes(x=Time, y=Count)) + 
    geom_line(aes(colour=Colour), size=1.5) + 
    scale_colour_manual(values=cPalette) + 
    theme(legend.position="none")
}


