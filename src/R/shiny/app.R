suppressMessages({
  library(shiny)
  library(shinyWidgets)
  library(readr)
  library(dplyr)
  library(plotly)
  library(heatmaply)
  library(shinycssloaders)
  library(corrr)
  library(tidyr)
})

options(shiny.port = 8080)

# https://shiny.rstudio.com/articles/dynamic-ui.html
# https://stackoverflow.com/questions/21465411/r-shiny-passing-reactive-to-selectinput-choices/21467399#21467399
# https://stackoverflow.com/questions/33973300/issue-in-dynamic-renderui-in-shiny
# https://drsimonj.svbtle.com/exploring-correlations-in-r-with-corrr

primary <- read_csv("out/work/primary/merged-mtic.csv")
primary <- arrange(primary,KEGG)
compounds <- unique(primary[c("Name","KEGG","HMDB","ChEBI","Category")])
# print(compounds)
# print()
# print(as.data.frame(primary))
# print(mtcars)

# https://rdrr.io/cran/shinyWidgets/man/updateCheckboxGroupButtons.html

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Rohner Lab Astyanax Metabolomics Study"),
  navlistPanel(
    "Subset",
    tabPanel("Selections",
#       https://stackoverflow.com/questions/27607566/allowing-one-tick-only-in-checkboxgroupinput
      p("Select which metabolites you want to include in the analysis. You can select via compound name, category, or a number of identifier systems. You can also select from different identifier types and the app will remember your selection for each heading and combine them on the summary page. Aftering making your selection, you can view the \"Summary\" tab to see the metabolites you select or proceed to any of the \"Visualization\" tabs.\""),
      radioButtons(
        inputId = "selection_type",
        label = "Select by:",
        choices = c("Name","KEGG","HMDB","ChEBI","Category"),
        selected = "Name",
#         selected = "Metabolites",
      ),
#       https://shiny.rstudio.com/gallery/creating-a-ui-from-a-loop.html
      lapply(c("Name","KEGG","HMDB","ChEBI","Category"), function(t) {
        conditionalPanel(
          condition = sprintf("input.selection_type == '%s'",t),
          pickerInput(
            inputId = sprintf("selector_%s",t),
            label = "Make selections:",
            choices = unique(compounds[t]),
            options = list(
      #              https://stackoverflow.com/questions/53609546/how-can-i-have-the-search-option-based-on-typing-letters-in-pickerinput-using-sh
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3",
              `live-search`=TRUE
            ),
            multiple = TRUE,
          )
        )
      })
    ),
    tabPanel("Summary",
      dataTableOutput("summary")
    ),
    "Visualization",
    tabPanel("Heatmap",
      tabsetPanel(type="tabs", id="corrPlotTab", selected="featureCorr",
        tabPanel(title = "By Metabolite",
          value = "featureCorr",
          plotlyOutput("featureCorrPlt", height = 600) %>%
            withSpinner(type = 8, color = "#0088cf", size = 1)),
        tabPanel(title = "By Sample",
          value = "sampleCorr",
          plotlyOutput("sampleCorrPlt", height = 600) %>%
            withSpinner(type = 8, color = "#0088cf", size = 1))
      )
    ),
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

  selected_cpds <- reactive(filter(compounds, compounds$Name %in% input$selector_Name | compounds$KEGG %in% input$selector_KEGG | compounds$HMDB %in% input$selector_HMDB | compounds$ChEBI %in% input$selector_ChEBI |  compounds$Category %in% input$selector_Category))

  output$summary <- renderDataTable(selected_cpds)

  ###################################################################
  #           Correlation Plot Tab 1: Sample Correlations           #
  ###################################################################
  output$featureCorrPlt <- renderPlotly({
#     callModule(module = savePlotlyPDF,
#                 id = "download_geneGenePlot",
#                 prefix = "GeneCorrelation_Plot_",
#                 plotlyToSave = reactive({featureCorrPlt[[2]]}))

#     https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
    cpd_data <- filter(primary, primary$Name %in% selected_cpds()$Name)
#     print(cpd_data)
    features <- cpd_data %>% select(Name,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean)
#     print(features)
    features <- as.data.frame(features)
#     features$Name <- as.factor(features$Name)
    rownames(features) <- features$Name
    features <- features[,-1]
    print(head(features))
    print(head(cor(features, use="pairwise.complete.obs")))
#     print("here")
#     stop
#     heatmaply_cor(
#       correlate(features),
#       xlab = "Features",
#       ylab = "Features",
#       k_col = 2,
#       k_row = 2
#     )
  })

}

shinyApp(ui = ui, server = server)
