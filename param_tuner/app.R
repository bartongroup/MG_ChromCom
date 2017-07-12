library(shiny)
library(ggplot2)
library(reshape2)

source("../R/lib.R")


########################################################################


# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
  titlePanel("Parameter tuner"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("t1", "Start time", min=-80, max=0, value=-28, step=1),
      sliderInput("r1", "BB->P rate", min=0, max=0.2, value=0.05, step=0.01),
      sliderInput("r2", "P->R rate", min=0, max=0.2, value=0.03, step=0.01),
      sliderInput("dt", "Time delay", min=0, max=50, value=0, step=1)
    ),
    
    mainPanel(
      plotOutput("tPlot")
    )
  )
))


server <- shinyServer(function(input, output) {

  sliderValues <- reactive({
    list(
      t1 = input$t1,
      r1 = input$r1,
      r2 = input$r2,
      dt = input$dt
    )
  })
  
  
  output$tPlot <- renderPlot({
    pars <- sliderValues()
    chr <- ChromCom3(pars)
    #print(chr)
    chr <- generateCells(chr)
    plotTimelines(chr)
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)

