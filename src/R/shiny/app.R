library(shiny)
library(shinyWidgets)
library(readr)

# https://shiny.rstudio.com/articles/dynamic-ui.html
# https://stackoverflow.com/questions/21465411/r-shiny-passing-reactive-to-selectinput-choices/21467399#21467399
# https://stackoverflow.com/questions/33973300/issue-in-dynamic-renderui-in-shiny

primary <- read_csv("out/work/primary/merged-mtic.csv")
compounds <- unique(primary[c("Name","KEGG","HMDB","ChEBI","Category")])
print(compounds)

# https://rdrr.io/cran/shinyWidgets/man/updateCheckboxGroupButtons.html

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Astyanax Metabolomics Study"),
  navlistPanel(
    "Header A",
    tabPanel("Selections",
#       https://stackoverflow.com/questions/27607566/allowing-one-tick-only-in-checkboxgroupinput
      radioButtons(
        inputId = "selection_type",
        label = "Select by:",
        choices = c("Name","KEGG","HMDB","ChEBI","Category"),
        selected = "Name",
#         selected = "Metabolites",
      ),
      uiOutput("selector")
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
  getSelectionType <- reactive({
   input$selection_type
  })
#
#   output$nrows <- reactive({
#     compounds[getSelection()]
#   })
  output$selector <- renderUI({
#     print(input$selector)
    pickerInput(
        inputId = "selector",
        label = "Select/deselect all + format selected",
        choices = compounds[getSelectionType()],
#         selected = input$selector,
        options = list(
  #              https://stackoverflow.com/questions/53609546/how-can-i-have-the-search-option-based-on-typing-letters-in-pickerinput-using-sh
          `actions-box` = TRUE,
          size = 10,
          `selected-text-format` = "count > 3",
          `live-search`=TRUE
        ),
        multiple = TRUE,
      )
  })

}

shinyApp(ui = ui, server = server)
