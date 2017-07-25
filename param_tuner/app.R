.libPaths("/cluster/gjb_lab/mgierlinski/R/x86_64-redhat-linux-gnu-library/3.3")
.libPaths("/home/mgierlinski/R/x86_64-pc-linux-gnu-library/3.4")


library(shiny)
library(ggplot2)
library(reshape2)
library(caTools)

source("../R/lib.R")

########################################################################


# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
  titlePanel("Parameter tuner"),

  sidebarLayout(
    sidebarPanel(
      radioButtons("dataSelection", "Background data selection", choices=names(dataFile)),
      radioButtons("modelResolution", "Model resolution", choices=list(Low = "low", Medium = "medium", High = "high")),
      sliderInput("t1", "Start time", min=-80, max=20, value=-27, step=1),
      sliderInput("k1", "B->P rate", min=0, max=0.2, value=0.05, step=0.001),
      sliderInput("k2", "P->R rate", min=0, max=0.2, value=0.03, step=0.001),
      sliderInput("k3", "P->B rate", min=0, max=0.2, value=0, step=0.001),
      sliderInput("dt2", "P->R delay", min=0, max=50, value=0, step=1),
      sliderInput("dt3", "P->B delay", min=0, max=150, value=0, step=1)
    ),

    mainPanel(
      p("This app helps finding the best model parameters to match experimental data. The parameters and background data in the figure can be changed in the side panel. There are three levels of model resolution, with increasing accuracy and smoothness but longer computation time."),
      plotOutput("tPlot"),
      downloadButton('downloadPlot', 'Download Plot')
    )
  )
))


server <- shinyServer(function(input, output) {

  D <- lapply(dataFile, experimentalData)

  sliderValues <- reactive({
    c3pars(
      t1 = input$t1,
      k1 = input$k1,
      k2 = input$k2,
      k3 = input$k3,
      dt2 = input$dt2,
      dt3 = input$dt3
    )
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
    chr <- ChromCom3(pars)
    echr <- getData()
    chr <- generateCells(chr, method="simulation", nsim=getNsim())
    chi2 <- round(oeError(chr, echr), 0)
    plotTimelines(chr, expdata=echr, xmin=-100, xmax=100, title=paste0(parstr(pars), ", chi2=", chi2))
  }

  output$tPlot <- renderPlot({
    pars <- sliderValues()
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

