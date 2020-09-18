library(shiny)
library(shinyWidgets)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Astyanax Metabolomics Study"),
  navlistPanel(
    "Header A",
    tabPanel("Selections",
#       https://stackoverflow.com/questions/27607566/allowing-one-tick-only-in-checkboxgroupinput
      checkboxGroupButtons(
        inputId = "somevalue",
        label = "Select by:",
        choices = c("Metabolites", "Categories"),
        justified = TRUE,
        individual = TRUE,
        status = "primary",
        selected = "Metabolites",
      )
    ),
    tabPanel("Component 2"),
    "Header B",
    tabPanel("Component 3"),
    tabPanel("Component 4"),
    "-----",
    tabPanel("Component 5")
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {

  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({

    x    <- faithful$waiting
    bins <- seq(min(x), max(x), length.out = input$bins + 1)

    hist(x, breaks = bins, col = "#75AADB", border = "white",
         xlab = "Waiting time to next eruption (in mins)",
         main = "Histogram of waiting times")

    })

}

shinyApp(ui = ui, server = server)
