suppressMessages({
  library(conflicted)
  library(shiny)
  library(shinyWidgets)
  library(readr)
  library(dplyr)
  library(plotly)
  library(heatmaply)
  library(shinycssloaders)
  library(corrr)
  library(tidyr)
  library(htmlwidgets)
  library(webshot)
  library(shinydashboard)
  library(MASS)
})

conflict_prefer("box", "shinydashboard")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

options(shiny.port = 8080)

# https://shiny.rstudio.com/articles/dynamic-ui.html
# https://stackoverflow.com/questions/21465411/r-shiny-passing-reactive-to-selectinput-choices/21467399#21467399
# https://stackoverflow.com/questions/33973300/issue-in-dynamic-renderui-in-shiny
# https://drsimonj.svbtle.com/exploring-correlations-in-r-with-corrr
# https://hssgenomics.shinyapps.io/RNAseq_DRaMA/
# https://stackoverflow.com/questions/28829682/r-shiny-checkboxgroupinput-select-all-checkboxes-by-click
# https://stackoverflow.com/questions/20637248/shiny-4-small-textinput-boxes-side-by-side
# https://rstudio.github.io/shinydashboard/structure.html
# https://rstudio.github.io/shinydashboard/appearance.html
# https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faqhow-are-the-likelihood-ratio-wald-and-lagrange-multiplier-score-tests-different-andor-similar/
# https://stackoverflow.com/questions/37597136/shinydashboard-is-it-not-possible-to-have-nested-menu-sub-items-cant-make-it
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

primary <- read_csv("out/work/primary/merged-mtic.csv") %>% arrange(KEGG)
compounds <- unique(primary[c("Name","KEGG","HMDB","ChEBI","Category")])
lipids <- read_csv("out/work/lipids/merged-lipids.csv", col_types = cols(Saturation = "c", Polarity = "f")) %>% arrange(LMID)
lipid_ids <- unique(lipids[c("LMID","Name","InChIKey","Category","MainClass","Saturation")])

pops = c("Pachon","Tinaja","Surface")
tissues = c("Brain","Muscle","Liver")
conditions = c("30d Starved", "4d Starved", "Refed")

popcolors = c("#b22222","#daa520","#1e90ff")

# https://rdrr.io/cran/shinyWidgets/man/updateCheckboxGroupButtons.html

savePlotlyPDFUI <- function(id, label = "Download PDF File"){
    ns <- NS(id)
    tagList(
        downloadButton(outputId = ns("downloadPDF"), label = label)
    )
}

# From RNADrama, doesn't work
savePlotlyPDF <- function(input, output, session, plotlyToSave, prefix = "",
                          delay = 10, ...){ # these are the default values vwidth = 992, vheight = 744

    namepdf = paste0('Plot_', prefix, Sys.Date(), ".pdf")
    namehtml = paste0('Plot_', prefix, Sys.Date(), ".html")

    output$downloadPDF <- downloadHandler(
        filename = namepdf,
        content =  function(file){
            withProgress(message = 'Saving PDF', style = "notification", value = 0, {
                for (i in 1:5) {
                    incProgress(0.1)
                    Sys.sleep(0.01)
                }
                # check what object is being saved, use webshot fror plotly, ggsave for ggplot
                if(class(plotlyToSave())[1] == "plotly"){
                    filename = namehtml
                    saveWidget(plotlyToSave(), namehtml, selfcontained = TRUE)
                    # vwidth = vwidth, vheight = vheight
                    webshot::webshot(url = namehtml, file = namepdf, delay = delay, ...)
                    file.copy(namepdf, file, overwrite = TRUE)
                    file.remove(namehtml)
                } else if(class(plotlyToSave())[1] == "gg"){
                    ggsave(namepdf, plotlyToSave(), ...)
                    file.copy(namepdf, file, overwrite = TRUE)
                } else if(class(plotlyToSave())[1] == "upset") {
                    filename = namepdf
                    pdf(file = namepdf, ...)
#                     print(plotlyToSave())
                    dev.off()
                    file.copy(namepdf, file, overwrite = TRUE)
                } else if(class(plotlyToSave())[1] == "visNetwork") {
                    filename = namehtml
                    saveWidget(plotlyToSave(), namehtml, selfcontained = TRUE)
                    # visSave(plotlyToSave(), namehtml, selfcontained = TRUE)
                    # vwidth = vwidth, vheight = vheight
                    webshot::webshot(url = namehtml, file = namepdf, delay = delay, ...)
                    file.copy(namepdf, file, overwrite = TRUE)
                    file.remove(namehtml)
                } else{ # supposed to be for base plots, but doesn't work at the moment
                    filename = namepdf
                    pdf(file = namepdf, ...)
                    plotlyToSave()
                    dev.off()
                    file.copy(namepdf, file, overwrite = TRUE)
                }
                for (i in 1:5) {
                    incProgress(0.1)
                    Sys.sleep(0.01)
                }
            })
        }
    )
}

# Define UI for app that draws a histogram ----
ui <- dashboardPage(
  dashboardHeader(title = "Interactive Astyanax Metabolome Project", titleWidth = 500,
                  tags$li(class = "dropdown",
                  tags$p("Rohner Lab, contact: Kyle Medley jmedley@stowers.org ",
                    style = "font-size:12px;color:white"))
  ),
  dashboardSidebar(
#     https://community.rstudio.com/t/how-to-remove-numeric-inputs-spin-button-in-r-shiny/13769/3
    tags$head(tags$style(HTML("
        input[type=number] {
              -moz-appearance:textfield;
        }
        input[type=number]::{
              -moz-appearance:textfield;
        }
        input[type=number]::-webkit-outer-spin-button,
        input[type=number]::-webkit-inner-spin-button {
              -webkit-appearance: none;
              margin: 0;
        }
    "))),
    sidebarMenu(
      id = "theSidebar",
      menuItem("Primary", selected = TRUE, startExpanded = TRUE,
               menuSubItem("Selections", tabName = "primarySelections"),
               menuSubItem("Summary", tabName = "primarySummary"),
               menuItem("PCA", menuSubItem("View Plot", tabName = "primaryPCA"), selectInput("primaryPCAMethod",label="",choices=c("PCA","MDS"))),
               menuSubItem("Correlation", tabName = "primaryCorrelation"),
               menuSubItem("Heatmap", tabName = "primaryHeatmap"),
               menuItem("Quantitative", menuSubItem("View Plot", tabName = "primaryQuantitative"), numericInput(inputId = "primaryQuantPercentileRange", label = "Err bar pct.:", value = 95, min = 1, max = 99, step = 1), checkboxInput(inputId = "primaryQuantShareY", label = "Share Y per row?", value = FALSE), checkboxInput(inputId = "primaryQuantIncludeOutliers", label = "Include Outliers?", value = FALSE))
              ),
      menuItem("Lipids",
               menuSubItem("Selections", tabName = "lipidsSelections"),
               menuSubItem("Summary", tabName = "lipidsSummary"),
               menuItem("PCA", menuSubItem("View Plot", tabName = "lipidsPCA"), selectInput("lipidsPCAMethod",label="",choices=c("PCA","MDS"))),
               menuItem("Correlation", tabName = "lipidsCorrelation"),
               menuItem("Quantitative", menuSubItem("View Plot", tabName = "lipidsQuantitative"), numericInput(inputId = "lipidsQuantPercentileRange", label = "Err bar pct.:", value = 95, min = 1, max = 99, step = 1), checkboxInput(inputId = "lipidsQuantShareY", label = "Share Y per row?", value = FALSE), checkboxInput(inputId = "lipidsQuantIncludeOutliers", label = "Include Outliers?", value = FALSE))
              )
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "primarySelections",
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
      tabItem(tabName = "primarySummary",
        dataTableOutput("primary_summary")
      ),
      tabItem(tabName = "primaryPCA",
          plotlyOutput("primarySamplePCAPlot", height = 600) %>%
            withSpinner(type = 8, color = "#0088cf", size = 1),
          box(title = "Controls",
            width = NULL,
            solidHeader = TRUE,
            status = "primary",
            splitLayout(cellWidths = c("25%","25%","25%","25%"),
              column(6,
                radioButtons("primaryPCAPlotSampleColorBy", label="Color by:", choices=c("Tissue","Population","Condition")),
                checkboxInput(inputId = "primaryPCAPlotSampleIncludeOutliers", label = "Include Outliers?", value = FALSE)
                ),
            checkboxGroupInput("primaryPCAPlotSampleSelectPops", "Populations:", pops, selected = pops),
            checkboxGroupInput("primaryPCAPlotSampleSelectTissues", "Tissues:", tissues, selected = tissues),
            checkboxGroupInput("primaryPCAPlotSampleSelectConditions", "Conditions:", conditions, selected = conditions)
          )
        )
      ),
      tabItem(tabName = "primaryCorrelation",
        tabsetPanel(type="tabs", id="primaryCorrPlotTab", selected="sampleCorr",
          tabPanel(title = "By Sample",
            value = "sampleCorr",
            plotlyOutput("primarySampleCorrPlt", height = 600) %>%
              withSpinner(type = 8, color = "#0088cf", size = 1),
          box(title = "Controls",
            width = NULL,
            solidHeader = TRUE,
            status = "primary",
            splitLayout(cellWidths = c("25%","25%","25%","25%"),
              column(6,
                numericInput(inputId = "primaryCorrPlotSampleNumClusters", label = "Number of Clusters:", value = 3, min = 1, step = 1),
                checkboxInput(inputId = "primaryCorrPlotSampleNormalize", label = "Normalize?", value = FALSE),
                checkboxInput(inputId = "primaryCorrPlotSampleIncludeOutliers", label = "Include Outliers?", value = FALSE),
                ),
              checkboxGroupInput("primaryCorrPlotSampleSelectPops", "Populations:", pops, selected = pops),
              checkboxGroupInput("primaryCorrPlotSampleSelectTissues", "Tissues:", tissues, selected = tissues),
              checkboxGroupInput("primaryCorrPlotSampleSelectConditions", "Conditions:", conditions, selected = conditions)
            )
          )),
          tabPanel(title = "By Metabolite",
            value = "featureCorr",
            plotlyOutput("primaryFeatureCorrPlt", height = 600) %>%
              withSpinner(type = 8, color = "#0088cf", size = 1),
            box(title = "Controls",
              width = NULL,
              solidHeader = TRUE,
              status = "primary",
              splitLayout(cellWidths = c("25%","25%","25%","25%"),
                column(6,
                  numericInput(inputId = "primaryCorrPlotFeatureNumClusters", label = "Number of Clusters:", value = 3, min = 1, step = 1),
                  checkboxInput(inputId = "primaryCorrPlotFeatureNormalize", label = "Normalize?", value = FALSE),
                  checkboxInput(inputId = "primaryCorrPlotFeatureIncludeOutliers", label = "Include Outliers?", value = FALSE),
                  ),
              checkboxGroupInput("primaryCorrPlotFeatureSelectPops", "Populations:", pops, selected = pops),
              checkboxGroupInput("primaryCorrPlotFeatureSelectTissues", "Tissues:", tissues, selected = tissues),
              checkboxGroupInput("primaryCorrPlotFeatureSelectConditions", "Conditions:", conditions, selected = conditions)
            )
          )),
          tabPanel(title = "By Category",
            value = "categoryCorr",
            plotlyOutput("primaryCategoryCorrPlt", height = 600) %>%
              withSpinner(type = 8, color = "#0088cf", size = 1),
            box(title = "Controls",
              width = NULL,
              solidHeader = TRUE,
              status = "primary",
              splitLayout(cellWidths = c("25%","25%","25%","25%"),
                column(6,
                  numericInput(inputId = "primaryCorrPlotCategoryNumClusters", label = "Number of Clusters:", value = 3, min = 1, step = 1),
                  checkboxInput(inputId = "primaryCorrPlotCategoryNormalize", label = "Normalize?", value = FALSE),
                  checkboxInput(inputId = "primaryCorrPlotCategoryIncludeOutliers", label = "Include Outliers?", value = FALSE),
                  ),
              checkboxGroupInput("primaryCorrPlotCategorySelectPops", "Populations:", pops, selected = pops),
              checkboxGroupInput("primaryCorrPlotCategorySelectTissues", "Tissues:", tissues, selected = tissues),
              checkboxGroupInput("primaryCorrPlotCategorySelectConditions", "Conditions:", conditions, selected = conditions)
            )
          ))
        )
#       p(class = 'text-center',
#         savePlotlyPDFUI(id = "download_primarySampleCorrPlot",
#                         label = "Download Plot (PDF)")
#       )
      ),
      tabItem(tabName = "primaryHeatmap",
        plotlyOutput("primaryHeatmapPlt", height = 600) %>%
          withSpinner(type = 8, color = "#0088cf", size = 1),
        box(title = "Controls",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          splitLayout(cellWidths = c("25%","25%","25%","25%"),
            column(6,
              selectInput("primaryHeatmapNormalize",label="",choices=c("row","column")),
              checkboxInput(inputId = "primaryHeatmapIncludeOutliers", label = "Include Outliers?", value = FALSE),
              ),
            checkboxGroupInput("primaryHeatmapSelectPops", "Populations:", pops, selected = pops),
            checkboxGroupInput("primaryHeatmapSelectTissues", "Tissues:", tissues, selected = tissues),
            checkboxGroupInput("primaryHeatmapSelectConditions", "Conditions:", conditions, selected = conditions)
          )
        ),
      ),
      tabItem(tabName = "primaryQuantitative",
#         https://stackoverflow.com/questions/21609436/r-shiny-conditionalpanel-output-value
        conditionalPanel(condition = "output.num_compounds <= 10",
          plotlyOutput("primaryLinePlt") %>%
          withSpinner(type = 8, color = "#0088cf", size = 1)
        ),
        conditionalPanel(condition = "output.num_compounds > 10",
          box(title = "Too many compounds", status = "danger", p("Too many compounds selected (max 10)"))
        )
      ),
      tabItem(tabName = "lipidsSelections",
        p("Select which lipids you want to include in the analysis."),
        radioButtons(
          inputId = "lipid_selection_type",
          label = "Select by:",
          choices = c("LMID","Name","InChIKey","Category","MainClass","Saturation"),
          selected = "LMID",
  #         selected = "Name",
        ),
  #       https://shiny.rstudio.com/gallery/creating-a-ui-from-a-loop.html
        lapply(c("LMID","Name","InChIKey","Category","MainClass","Saturation"), function(t) {
          conditionalPanel(
            condition = sprintf("input.lipid_selection_type == '%s'",t),
            pickerInput(
              inputId = sprintf("lipid_selector_%s",t),
              label = "Make selections:",
              choices = unique(lipid_ids[t]),
              options = list(
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
      tabItem(tabName = "lipidsSummary",
        dataTableOutput("lipids_summary")
      ),
      tabItem(tabName = "lipidsPCA",
          plotlyOutput("lipidsSamplePCAPlot", height = 600) %>%
            withSpinner(type = 8, color = "#0088cf", size = 1),
          box(title = "Controls",
            width = NULL,
            solidHeader = TRUE,
            status = "primary",
            splitLayout(cellWidths = c("25%","25%","25%","25%"),
              column(6,
                radioButtons("lipidsPCAPlotSampleColorBy", label="Color by:", choices=c("Tissue","Population","Condition")),
                checkboxInput(inputId = "lipidsPCAPlotSampleIncludeOutliers", label = "Include Outliers?", value = FALSE)
                ),
            checkboxGroupInput("lipidsPCAPlotSampleSelectPops", "Populations:", pops, selected = pops),
            checkboxGroupInput("lipidsPCAPlotSampleSelectTissues", "Tissues:", tissues, selected = tissues),
            checkboxGroupInput("lipidsPCAPlotSampleSelectConditions", "Conditions:", conditions, selected = conditions)
          )
        )
      ),
      tabItem(tabName = "lipidsCorrelation",
        tabsetPanel(type="tabs", id="lipidsCorrPlotTab", selected="sampleCorr",
          tabPanel(title = "By Sample",
            value = "sampleCorr",
            plotlyOutput("lipidsSampleCorrPlt", height = 600) %>%
              withSpinner(type = 8, color = "#0088cf", size = 1),
          box(title = "Controls",
            width = NULL,
            solidHeader = TRUE,
            status = "primary",
            splitLayout(cellWidths = c("25%","25%","25%","25%"),
              column(6,
                numericInput(inputId = "lipidsCorrPlotSampleNumClusters", label = "Number of Clusters:", value = 3, min = 1, step = 1),
                checkboxInput(inputId = "lipidsCorrPlotSampleNormalize", label = "Normalize?", value = FALSE),
                checkboxInput(inputId = "lipidsCorrPlotSampleIncludeOutliers", label = "Include Outliers?", value = FALSE),
                ),
              checkboxGroupInput("lipidsCorrPlotSampleSelectPops", "Populations:", pops, selected = pops),
              checkboxGroupInput("lipidsCorrPlotSampleSelectTissues", "Tissues:", tissues, selected = tissues),
              checkboxGroupInput("lipidsCorrPlotSampleSelectConditions", "Conditions:", conditions, selected = conditions)
            )
          )),
          tabPanel(title = "By Metabolite",
            value = "featureCorr",
            plotlyOutput("lipidsFeatureCorrPlt", height = 600) %>%
              withSpinner(type = 8, color = "#0088cf", size = 1),
            box(title = "Controls",
              width = NULL,
              solidHeader = TRUE,
              status = "primary",
              splitLayout(cellWidths = c("25%","25%","25%","25%"),
                column(6,
                  numericInput(inputId = "lipidsCorrPlotFeatureNumClusters", label = "Number of Clusters:", value = 3, min = 1, step = 1),
                  checkboxInput(inputId = "lipidsCorrPlotFeatureNormalize", label = "Normalize?", value = FALSE),
                  checkboxInput(inputId = "lipidsCorrPlotFeatureIncludeOutliers", label = "Include Outliers?", value = FALSE),
                  ),
              checkboxGroupInput("lipidsCorrPlotFeatureSelectPops", "Populations:", pops, selected = pops),
              checkboxGroupInput("lipidsCorrPlotFeatureSelectTissues", "Tissues:", tissues, selected = tissues),
              checkboxGroupInput("lipidsCorrPlotFeatureSelectConditions", "Conditions:", conditions, selected = conditions)
            )
          )),
          tabPanel(title = "By Category",
            value = "categoryCorr",
            plotlyOutput("lipidsCategoryCorrPlt", height = 600) %>%
              withSpinner(type = 8, color = "#0088cf", size = 1),
            box(title = "Controls",
              width = NULL,
              solidHeader = TRUE,
              status = "primary",
              splitLayout(cellWidths = c("25%","25%","25%","25%"),
                column(6,
                  numericInput(inputId = "lipidsCorrPlotCategoryNumClusters", label = "Number of Clusters:", value = 3, min = 1, step = 1),
                  checkboxInput(inputId = "lipidsCorrPlotCategoryNormalize", label = "Normalize?", value = FALSE),
                  checkboxInput(inputId = "lipidsCorrPlotCategoryIncludeOutliers", label = "Include Outliers?", value = FALSE),
                  ),
              checkboxGroupInput("lipidsCorrPlotCategorySelectPops", "Populations:", pops, selected = pops),
              checkboxGroupInput("lipidsCorrPlotCategorySelectTissues", "Tissues:", tissues, selected = tissues),
              checkboxGroupInput("lipidsCorrPlotCategorySelectConditions", "Conditions:", conditions, selected = conditions)
            )
          )),
          tabPanel(title = "By Class",
            value = "classCorr",
            plotlyOutput("lipidsClassCorrPlt", height = 600) %>%
              withSpinner(type = 8, color = "#0088cf", size = 1),
            box(title = "Controls",
              width = NULL,
              solidHeader = TRUE,
              status = "primary",
              splitLayout(cellWidths = c("25%","25%","25%","25%"),
                column(3,
                  numericInput(inputId = "lipidsCorrPlotClassNumClusters", label = "Number of Clusters:", value = 3, min = 1, width="50%", step = 1),
                  checkboxInput(inputId = "lipidsCorrPlotClassNormalize", label = "Normalize?", value = FALSE),
                  checkboxInput(inputId = "lipidsCorrPlotClassIncludeOutliers", label = "Include Outliers?", value = FALSE),
                  ),
              checkboxGroupInput("lipidsCorrPlotClassSelectPops", "Populations:", pops, selected = pops),
              checkboxGroupInput("lipidsCorrPlotClassSelectTissues", "Tissues:", tissues, selected = tissues),
              checkboxGroupInput("lipidsCorrPlotClassSelectConditions", "Conditions:", conditions, selected = conditions)
            )
          ))
        )
      ),
      tabItem(tabName = "lipidsQuantitative",
#         https://stackoverflow.com/questions/21609436/r-shiny-conditionalpanel-output-value
        conditionalPanel(condition = "output.num_lipids <= 10",
          plotlyOutput("lipidsLinePlt") %>%
          withSpinner(type = 8, color = "#0088cf", size = 1)
        ),
        conditionalPanel(condition = "output.num_lipids > 10",
          box(title = "Too many lipids", status = "danger", p("Too many lipids selected (max 10)"))
        )
      )
    ),
    conditionalPanel(condition = "0",
      verbatimTextOutput("num_compounds")
    ),
    conditionalPanel(condition = "0",
      verbatimTextOutput("num_lipids")
    )
  )
)


server <- function(input, output) {

  ###################################################################
  #                            Primary                              #
  ###################################################################

  selected_cpds <- reactive(filter(compounds, compounds$Name %in% input$selector_Name | compounds$KEGG %in% input$selector_KEGG | compounds$HMDB %in% input$selector_HMDB | compounds$ChEBI %in% input$selector_ChEBI |  compounds$Category %in% input$selector_Category))

  output$num_compounds <- reactive(nrow(selected_cpds()))
#   https://stackoverflow.com/questions/21609436/r-shiny-conditionalpanel-output-value
  outputOptions(output, "num_compounds", suspendWhenHidden = FALSE)

  output$primary_summary <- renderDataTable(selected_cpds())

  ###################################################################
  #          Primary Sample PCA                                     #
  ###################################################################
  output$primarySamplePCAPlot <- renderPlotly({
    cpd_data <- primary %>% filter(primary$Name %in% selected_cpds()$Name) %>% filter(Population %in% input$primaryPCAPlotSampleSelectPops) %>% filter(Tissue %in% input$primaryPCAPlotSampleSelectTissues) %>% filter(Condition %in% input$primaryPCAPlotSampleSelectConditions)
    if (!input$primaryPCAPlotSampleIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(Name,Population,Tissue,Condition,Replicate,Raw_mTIC)
    spec <- features %>% build_wider_spec(names_from=c("Population","Tissue","Condition","Replicate"),values_from="Raw_mTIC")
    features <- features %>%  pivot_wider_spec(spec, values_fn = mean)
    features <- as.data.frame(features)
    rownames(features) <- features$Name
    features <- features[,-1]
    features <- t(as.matrix(features))
    if (input$primaryPCAMethod == "PCA") {
      pca <- prcomp(features,center=TRUE,scale.=TRUE,rank=2)
      pca <- as.data.frame(pca$x)
    } else {
      pca <- isoMDS(dist(features),k=2)
      pca <- as.data.frame(pca$points)
      colnames(pca) <- c("PC1","PC2")
    }
    pca$Tissue = spec$Tissue
    pca$Population = spec$Population
    pca$Condition = spec$Condition
    colorby <- input$primaryPCAPlotSampleColorBy
    plot_ly(data = as.data.frame(pca), x = ~PC1, y = ~PC2, type = "scatter", mode="markers", color= as.formula(sprintf("~%s",colorby)), height=600)
  })

  ###################################################################
  #          Primary Sample Correlation Plot                        #
  ###################################################################
  output$primarySampleCorrPlt <- renderPlotly({
#     https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
    cpd_data <- primary %>% filter(primary$Name %in% selected_cpds()$Name) %>% filter(Population %in% input$primaryCorrPlotSampleSelectPops) %>% filter(Tissue %in% input$primaryCorrPlotSampleSelectTissues) %>% filter(Condition %in% input$primaryCorrPlotSampleSelectConditions)
    if (!input$primaryCorrPlotSampleIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(Name,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean)
    features <- as.data.frame(features)
    rownames(features) <- features$Name
    features <- features[,-1]
    if (input$primaryCorrPlotSampleNormalize) {
      thecor <- cor(features)
      a <- max(thecor)
      b <- min(thecor)
      thecor <- (thecor - b) / (a-b)
    } else {
      thecor <- cor(features)
    }
    theplt <- heatmaply_cor(
      thecor,
      k_col = input$primaryCorrPlotSampleNumClusters,
      k_row = input$primaryCorrPlotSampleNumClusters
    )
    callModule(module = savePlotlyPDF,
                id = "download_primarySampleCorrPlot",
                prefix = "PrimarySampleCorrPlot_",
                plotlyToSave = reactive(theplt))
    theplt
  })

  ###################################################################
  #           Primary Feature Correlation Plot                      #
  ###################################################################
  output$primaryFeatureCorrPlt <- renderPlotly({
#     https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
    cpd_data <- primary %>% filter(primary$Name %in% selected_cpds()$Name) %>% filter(Population %in% input$primaryCorrPlotFeatureSelectPops) %>% filter(Tissue %in% input$primaryCorrPlotFeatureSelectTissues) %>% filter(Condition %in% input$primaryCorrPlotFeatureSelectConditions)
    if (!input$primaryCorrPlotFeatureIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(Name,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean)
    features <- as.data.frame(features)
    rownames(features) <- features$Name
    features <- features[,-1]
    features <- t(features)
    if (input$primaryCorrPlotFeatureNormalize) {
      thecor <- normalize(cor(features))
    } else {
      thecor <- cor(features)
    }
    theplt <- heatmaply_cor(
      thecor,
      k_col = input$primaryCorrPlotFeatureNumClusters,
      k_row = input$primaryCorrPlotFeatureNumClusters
    )
    callModule(module = savePlotlyPDF,
                id = "download_primaryFeatureCorrPlot",
                prefix = "PrimaryFeatureCorrPlot_",
                plotlyToSave = reactive(theplt))
    theplt
  })

  ###################################################################
  #          Primary Category Correlation Plot                      #
  ###################################################################
  output$primaryCategoryCorrPlt <- renderPlotly({
#     https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
    cpd_data <- primary %>% filter(primary$Name %in% selected_cpds()$Name) %>% filter(Population %in% input$primaryCorrPlotCategorySelectPops) %>% filter(Tissue %in% input$primaryCorrPlotCategorySelectTissues) %>% filter(Condition %in% input$primaryCorrPlotCategorySelectConditions)
    if (!input$primaryCorrPlotCategoryIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(Category,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean)
    features <- as.data.frame(features)
    rownames(features) <- features$Category
    features <- features[,-1]
    features <- t(features)
    if (input$primaryCorrPlotCategoryNormalize) {
      thecor <- normalize(cor(features))
    } else {
      thecor <- cor(features)
    }
    theplt <- heatmaply_cor(
      thecor,
      k_col = input$primaryCorrPlotCategoryNumClusters,
      k_row = input$primaryCorrPlotCategoryNumClusters
    )
    callModule(module = savePlotlyPDF,
                id = "download_primaryCategoryCorrPlot",
                prefix = "PrimaryCategoryCorrPlot_",
                plotlyToSave = reactive(theplt))
    theplt
  })

  ###################################################################
  #          Primary Category Correlation Plot                      #
  ###################################################################
  output$primaryCategoryCorrPlt <- renderPlotly({
#     https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
    cpd_data <- primary %>% filter(primary$Name %in% selected_cpds()$Name) %>% filter(Population %in% input$primaryCorrPlotCategorySelectPops) %>% filter(Tissue %in% input$primaryCorrPlotCategorySelectTissues) %>% filter(Condition %in% input$primaryCorrPlotCategorySelectConditions)
    if (!input$primaryCorrPlotCategoryIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(Category,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean)
    features <- as.data.frame(features)
    rownames(features) <- features$Category
    features <- features[,-1]
    features <- t(features)
    if (input$primaryCorrPlotCategoryNormalize) {
      thecor <- normalize(cor(features))
    } else {
      thecor <- cor(features)
    }
    theplt <- heatmaply_cor(
      thecor,
      k_col = input$primaryCorrPlotCategoryNumClusters,
      k_row = input$primaryCorrPlotCategoryNumClusters
    )
    callModule(module = savePlotlyPDF,
                id = "download_primaryCategoryCorrPlot",
                prefix = "PrimaryCategoryCorrPlot_",
                plotlyToSave = reactive(theplt))
    theplt
  })

  ###################################################################
  #          Primary Heatmap                                        #
  ###################################################################
  output$primaryHeatmapPlt <- renderPlotly({
#     https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
    cpd_data <- primary %>% filter(primary$Name %in% selected_cpds()$Name) %>% filter(Population %in% input$primaryHeatmapSelectPops) %>% filter(Tissue %in% input$primaryHeatmapSelectTissues) %>% filter(Condition %in% input$primaryHeatmapSelectConditions)
    if (!input$primaryHeatmapIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(Name,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean)
    features <- as.data.frame(features)
    rownames(features) <- features$Name
    features <- features[,-1]
    scale <- input$primaryHeatmapNormalize
    theplt <- heatmaply(
      features,
      scale=scale,
      main = "Metabolites vs Samples",
    )
    callModule(module = savePlotlyPDF,
                id = "download_primarySampleCorrPlot",
                prefix = "PrimarySampleCorrPlot_",
                plotlyToSave = reactive(theplt))
    theplt
  })

  ###################################################################
  #                        Primary Line Plot                        #
  ###################################################################
  output$primaryLinePlt <- renderPlotly({
    subplot(lapply(selected_cpds()$Name,
      function(name) {
        sharey <- input$primaryQuantShareY
        subplot(lapply(tissues, function(tissue) {
          cpd_data <- primary %>% filter(primary$Name == name) %>% filter(Tissue == tissue)
          if (!input$primaryQuantIncludeOutliers){
            cpd_data <- cpd_data %>% filter(Outlier == FALSE)
          }
          lower_qt <- (50-input$lipidsQuantPercentileRange/2)/100
          upper_qt <- (50+input$lipidsQuantPercentileRange/2)/100
          cpd_data <- cpd_data %>% group_by(Population,Condition) %>% summarize(.groups="drop_last",Intensity=mean(Raw_mTIC),Std=sd(Raw_mTIC),Lower=mean(Raw_mTIC)-quantile(Raw_mTIC,lower_qt,type=4),Upper=quantile(Raw_mTIC,upper_qt,type=4)-mean(Raw_mTIC)) %>% arrange(factor(Condition, levels = conditions)) %>% arrange(factor(Population, levels = pops))
          cpd_data$Condition <- factor(cpd_data$Condition, levels = conditions)
          cpd_data$Population <- factor(cpd_data$Population, levels=pops)
#           https://stackoverflow.com/questions/37285729/how-to-give-subtitles-for-subplot-in-plot-ly-using-r
          plt <- plot_ly(data = as.data.frame(cpd_data), x = ~Condition, y = ~Intensity, type = "scatter", mode="lines+markers", error_y=~list(symmetric=FALSE,type="data",array=Upper,arrayminus=Lower), color= ~Population, colors=popcolors, legendgroup=~Population, height=250*length(selected_cpds()$Name), showlegend=(tissue == "Brain" && name == selected_cpds()$Name[1])) %>% add_annotations(text = tissue, x = 0.5, y = 1.0, xref = "paper", yref = "paper", xanchor = "middle", yanchor = "top", showarrow = FALSE, font=list(size=15,weight="bold"))
#           https://stackoverflow.com/questions/57253488/how-to-remove-duplicate-legend-entries-w-plotly-subplots/57312776
          if (tissue == "Brain") {
            plt <- plt %>% add_annotations(text = name, x = -0.1, y = 0.5, xref = "paper", yref = "paper", xanchor = "right", yanchor = "middle", showarrow = FALSE, textangle=-90, font=list(size=15,weight="bold"))
          }
          plt
        }), shareY=sharey)
      }), nrows = length(selected_cpds()$Name)
    )
  })

  ###################################################################
  #                            Lipids                               #
  ###################################################################

  selected_lipids <- reactive(filter(lipid_ids, lipid_ids$LMID %in% input$lipid_selector_LMID| lipid_ids$Name %in% input$lipid_selector_Name | lipid_ids$InChIKey %in% input$lipid_selector_InChIKey | lipid_ids$Category %in% input$lipid_selector_Category |  lipid_ids$MainClass %in% input$lipid_selector_MainClass |  lipid_ids$Saturation %in% input$lipid_selector_Saturation))

  output$num_lipids <- reactive(nrow(selected_lipids()))
#   https://stackoverflow.com/questions/21609436/r-shiny-conditionalpanel-output-value
  outputOptions(output, "num_lipids", suspendWhenHidden = FALSE)

  output$lipids_summary <- renderDataTable(selected_lipids())

  ###################################################################
  #          Lipids Sample PCA                                     #
  ###################################################################
  output$lipidsSamplePCAPlot <- renderPlotly({
    cpd_data <- lipids %>% filter(lipids$Name %in% selected_lipids()$Name) %>% filter(Population %in% input$lipidsPCAPlotSampleSelectPops) %>% filter(Tissue %in% input$lipidsPCAPlotSampleSelectTissues) %>% filter(Condition %in% input$lipidsPCAPlotSampleSelectConditions)
    if (!input$lipidsPCAPlotSampleIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(Name,Population,Tissue,Condition,Replicate,Raw_mTIC)
    spec <- features %>% build_wider_spec(names_from=c("Population","Tissue","Condition","Replicate"),values_from="Raw_mTIC")
    features <- features %>%  pivot_wider_spec(spec, values_fn = mean)
    features <- as.data.frame(features)
    features <- features %>% drop_na()
    print(features)
    rownames(features) <- features$Name
    features <- features[,-1]
    features <- t(as.matrix(features))
    if (input$lipidsPCAMethod == "PCA") {
      pca <- prcomp(features,center=TRUE,scale.=TRUE,rank=2)
      pca <- as.data.frame(pca$x)
    } else {
      pca <- isoMDS(dist(features),k=2)
      pca <- as.data.frame(pca$points)
      colnames(pca) <- c("PC1","PC2")
      print(head(pca))
    }
    pca$Tissue = spec$Tissue
    pca$Population = spec$Population
    pca$Condition = spec$Condition
    colorby <- input$lipidsPCAPlotSampleColorBy
    plot_ly(data = as.data.frame(pca), x = ~PC1, y = ~PC2, type = "scatter", mode="markers", color= as.formula(sprintf("~%s",colorby)), height=600)
  })

  ###################################################################
  #           Lipids Sample Correlation Plot                        #
  ###################################################################
  output$lipidsSampleCorrPlt <- renderPlotly({
#     https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
    cpd_data <- lipids %>% filter(lipids$LMID %in% selected_lipids()$LMID) %>% filter(Population %in% input$lipidsCorrPlotSampleSelectPops) %>% filter(Tissue %in% input$lipidsCorrPlotSampleSelectTissues) %>% filter(Condition %in% input$lipidsCorrPlotSampleSelectConditions)
    if (!input$lipidsCorrPlotSampleIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(LMID,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean) %>% drop_na()
    features <- as.data.frame(features)
    rownames(features) <- features$LMID
    features <- features[,-1]
#     print(features)
    if (input$lipidsCorrPlotSampleNormalize) {
      thecor <- cor(features)
      a <- max(thecor)
      b <- min(thecor)
      thecor <- (thecor - b) / (a-b)
    } else {
      thecor <- cor(features)
    }
    theplt <- heatmaply_cor(
      thecor,
      k_col = input$lipidsCorrPlotSampleNumClusters,
      k_row = input$lipidsCorrPlotSampleNumClusters
    )
    callModule(module = savePlotlyPDF,
                id = "download_lipidsSampleCorrPlot",
                prefix = "LipidsSampleCorrPlot_",
                plotlyToSave = reactive(theplt))
    theplt
  })

  ###################################################################
  #            Lipids Feature Correlation Plot                      #
  ###################################################################
  output$lipidsFeatureCorrPlt <- renderPlotly({
#     https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
    cpd_data <- lipids %>% filter(lipids$LMID %in% selected_lipids()$LMID) %>% filter(Population %in% input$lipidsCorrPlotFeatureSelectPops) %>% filter(Tissue %in% input$lipidsCorrPlotFeatureSelectTissues) %>% filter(Condition %in% input$lipidsCorrPlotFeatureSelectConditions)
    if (!input$lipidsCorrPlotFeatureIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(LMID,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean) %>% drop_na()
    features <- as.data.frame(features)
    rownames(features) <- features$LMID
    features <- features[,-1]
    features <- t(features)
    if (input$lipidsCorrPlotFeatureNormalize) {
      thecor <- normalize(cor(features))
    } else {
      thecor <- cor(features)
    }
    theplt <- heatmaply_cor(
      thecor,
      k_col = input$lipidsCorrPlotFeatureNumClusters,
      k_row = input$lipidsCorrPlotFeatureNumClusters
    )
    callModule(module = savePlotlyPDF,
                id = "download_lipidsFeatureCorrPlot",
                prefix = "LipidsFeatureCorrPlot_",
                plotlyToSave = reactive(theplt))
    theplt
  })

  ###################################################################
  #           Lipids Category Correlation Plot                      #
  ###################################################################
  output$lipidsCategoryCorrPlt <- renderPlotly({
#     https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
    cpd_data <- lipids %>% filter(lipids$LMID %in% selected_lipids()$LMID) %>% filter(Population %in% input$lipidsCorrPlotCategorySelectPops) %>% filter(Tissue %in% input$lipidsCorrPlotCategorySelectTissues) %>% filter(Condition %in% input$lipidsCorrPlotCategorySelectConditions)
    if (!input$lipidsCorrPlotCategoryIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(Category,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean) %>% drop_na()
    features <- as.data.frame(features)
    rownames(features) <- features$Category
    features <- features[,-1]
    features <- t(features)
    if (input$lipidsCorrPlotCategoryNormalize) {
      thecor <- normalize(cor(features))
    } else {
      thecor <- cor(features)
    }
    theplt <- heatmaply_cor(
      thecor,
      k_col = input$lipidsCorrPlotCategoryNumClusters,
      k_row = input$lipidsCorrPlotCategoryNumClusters
    )
    callModule(module = savePlotlyPDF,
                id = "download_lipidsCategoryCorrPlot",
                prefix = "LipidsCategoryCorrPlot_",
                plotlyToSave = reactive(theplt))
    theplt
  })

  ###################################################################
  #           Lipids Class Correlation Plot                      #
  ###################################################################
  output$lipidsClassCorrPlt <- renderPlotly({
#     https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
    cpd_data <- lipids %>% filter(lipids$LMID %in% selected_lipids()$LMID) %>% filter(Population %in% input$lipidsCorrPlotClassSelectPops) %>% filter(Tissue %in% input$lipidsCorrPlotClassSelectTissues) %>% filter(Condition %in% input$lipidsCorrPlotClassSelectConditions)
    if (!input$lipidsCorrPlotClassIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(MainClass,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean) %>% drop_na()
    features <- as.data.frame(features)
    rownames(features) <- features$MainClass
    features <- features[,-1]
    features <- t(features)
    if (input$lipidsCorrPlotClassNormalize) {
      thecor <- normalize(cor(features))
    } else {
      thecor <- cor(features)
    }
    theplt <- heatmaply_cor(
      thecor,
      k_col = input$lipidsCorrPlotClassNumClusters,
      k_row = input$lipidsCorrPlotClassNumClusters
    )
    callModule(module = savePlotlyPDF,
                id = "download_lipidsClassCorrPlot",
                prefix = "LipidsClassCorrPlot_",
                plotlyToSave = reactive(theplt))
    theplt
  })

  ###################################################################
  #                         Lipids Line Plot                        #
  ###################################################################
  output$lipidsLinePlt <- renderPlotly({
    subplot(lapply(selected_lipids()$Name,
      function(name) {
        sharey <- input$lipidsQuantShareY
        subplot(lapply(tissues, function(tissue) {
          cpd_data <- lipids %>% filter(lipids$Name == name) %>% filter(Tissue == tissue)
          if (!input$primaryQuantIncludeOutliers){
            cpd_data <- cpd_data %>% filter(Outlier == FALSE)
          }
          lower_qt <- (50-input$primaryQuantPercentileRange/2)/100
          upper_qt <- (50+input$primaryQuantPercentileRange/2)/100
          cpd_data <- cpd_data %>% group_by(Population,Condition) %>% summarize(.groups="drop_last",Intensity=mean(Raw_mTIC),Std=sd(Raw_mTIC),Lower=mean(Raw_mTIC)-quantile(Raw_mTIC,lower_qt,type=4),Upper=quantile(Raw_mTIC,upper_qt,type=4)-mean(Raw_mTIC)) %>% arrange(factor(Condition, levels = conditions)) %>% arrange(factor(Population, levels = pops))
          cpd_data$Condition <- factor(cpd_data$Condition, levels = conditions)
          cpd_data$Population <- factor(cpd_data$Population, levels=pops)
#           https://stackoverflow.com/questions/37285729/how-to-give-subtitles-for-subplot-in-plot-ly-using-r
          plt <- plot_ly(data = cpd_data, x = ~Condition, y = ~Intensity, type = "scatter", mode="lines+markers", error_y=~list(symmetric=FALSE,type="data",array=Upper,arrayminus=Lower), color= ~Population, colors=popcolors, legendgroup=~Population, height=250*length(selected_lipids()$Name), showlegend=(tissue == "Brain" && name == selected_lipids()$Name[1])) %>% add_annotations(text = tissue, x = 0.5, y = 1.0, xref = "paper", yref = "paper", xanchor = "middle", yanchor = "top", showarrow = FALSE, font=list(size=15,weight="bold"))
#           https://stackoverflow.com/questions/57253488/how-to-remove-duplicate-legend-entries-w-plotly-subplots/57312776
          if (tissue == "Brain") {
            plt <- plt %>% add_annotations(text = name, x = -0.2, y = 0.5, xref = "paper", yref = "paper", xanchor = "right", yanchor = "middle", showarrow = FALSE, textangle=-90, font=list(size=15,weight="bold"))
          }
          plt
        }), shareY=sharey)
      }), nrows = length(selected_lipids()$Name)
    )
  })

}

shinyApp(ui = ui, server = server)
