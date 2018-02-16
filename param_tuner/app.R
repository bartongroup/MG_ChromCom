paths <- c(
  "/cluster/gjb_lab/mgierlinski/R/x86_64-redhat-linux-gnu-library/3.3",
  "/home/mgierlinski/R/x86_64-pc-linux-gnu-library/3.4",
  "/Users/mgierlinski/Library/R/3.4/library"
)
for(path in paths) if(dir.exists(path)) .libPaths(path)

library(shiny)
library(ggplot2)
library(reshape2)
library(caTools)
library(latex2exp)
library(parallel)

source("../R/lib.R")

# main data
D <- lapply(dataFile, experimentalData)

########################################################################


# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
  titlePanel("Parameter tuner"),

  sidebarLayout(
    sidebarPanel(
      radioButtons("dataSelection", "Background data selection", choices=names(dataFile)),
      radioButtons("modelResolution", "Model resolution", choices=list(Low = "low", Medium = "medium", High = "high")),
      radioButtons("modelSwitch", "P->R trigger", choices=c("t1", "t0")),
      sliderInput("t0", "NEB correction", min=-30, max=30, value=0, step=0.1),
      sliderInput("tau1", "Initial timescale", min=0, max=80, value=20, step=0.1),
      sliderInput("k1", "B->P rate", min=0, max=0.5, value=0.05, step=0.001),
      sliderInput("k2", "P->R rate", min=0, max=0.5, value=0.03, step=0.001),
      sliderInput("k3", "P->B rate", min=0, max=0.5, value=0, step=0.001),
      sliderInput("tau2", "P->R delay timescale", min=0, max=80, value=0, step=0.1),
      sliderInput("tau3", "P->B delay timescale", min=0, max=80, value=0, step=0.11)
    ),

    mainPanel(
      p("This app helps finding the best model parameters to match experimental data. The parameters and background data in the figure can be changed in the side panel. There are three levels of model resolution, with increasing accuracy and smoothness but longer computation time. The dashed vertical lines represent time delays. The solid vertical line indicates zero."),
      plotOutput("tPlot"),
      downloadButton('downloadPlot', 'Download Plot')
    )
  )
))


server <- shinyServer(function(input, output) {

  sliderValues <- reactive({
    c3pars(
      t0 = input$t0,
      tau1 = input$tau1,
      k1 = input$k1,
      k2 = input$k2,
      k3 = input$k3,
      tau2 = input$tau2,
      tau3 = input$tau3
    )
  })
  t2ref <- reactive({
    s <- input$modelSwitch
    ifelse(s == "t0", 0, 1)
  })

  getData <- function() {
    dat <- D[[input$dataSelection]]
  }

  getNsim <- function() {
    mr <- input$modelResolution
    if(mr == "low") {
      nsim <- 500
    } else if (mr == "medium") {
      nsim <- 3000
    } else if (mr == "high") {
      nsim <- 30000
    }
  }

  parstr <- function(pars) {
    paste0(lapply(names(pars), function(name) paste0(name, "=", pars[[name]])), collapse=", ")
  }

  parfstr <- function(pars) {
    paste0(lapply(names(pars), function(name) paste0(name, "_", pars[[name]])), collapse=".")
  }

  tPlot <- function(pars) {
    chr <- ChromCom3(pars, timepars=list(start=-150, stop=50, step=1))
    echr <- getData()
    chr <- generateCells(chr, mode="simulation", ncells=getNsim())
    rms <- round(oeError(chr, echr, limits=c(-50, 30)), 1)
    plotTimelines(chr, expdata=echr, xmin=-100, xmax=100, withpars = TRUE)
  }

  output$tPlot <- renderPlot({
    pars <- sliderValues()
    pars$t2ref <- t2ref()
    tPlot(pars)
  }, res=120)

  output$downloadPlot <- downloadHandler(
    filename = function() {
      pars <- sliderValues()
      paste0(parfstr(pars), ".png")
    },
    content = function(file) {
      pars <- sliderValues()
      ggsave(file, plot = tPlot(pars), device = "png")
    }
  )

})

# Run the application
shinyApp(ui = ui, server = server)

