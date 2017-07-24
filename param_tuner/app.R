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
      sliderInput("t1", "Start time", min=-80, max=20, value=-28, step=1),
      sliderInput("r1", "BB->P rate", min=0, max=0.2, value=0.05, step=0.002),
      sliderInput("r2", "P->R rate", min=0, max=0.2, value=0.03, step=0.002),
      sliderInput("dt", "Time delay", min=0, max=50, value=0, step=1)
    ),

    mainPanel(
      p("This app helps finding the best model parameters to match experimental data. The parameters and background data in the figure can be changed in the side panel."),
      plotOutput("tPlot")
    )
  )
))


server <- shinyServer(function(input, output) {

  D <- lapply(dataFile, experimentalData)

  sliderValues <- reactive({
    list(
      t1 = input$t1,
      r1 = input$r1,
      r2 = input$r2,
      dt = input$dt
    )
  })

  getData <- function() {
    dat <- D[[input$dataSelection]]
  }

  output$tPlot <- renderPlot({
    pars <- sliderValues()
    chr <- ChromCom3(pars)
    echr <- getData()
    #print(chr)
    chr <- generateCells(chr)
    rms <- round(RMS(chr, echr), 0)
    plotTimelines(chr, expdata=echr, xmin=-100, xmax=100, title=paste0("rms = ", rms))
  }, res=120)

})

# Run the application
shinyApp(ui = ui, server = server)

