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
  library(htmlwidgets)
  library(webshot)
  library(shinydashboard)
})

options(shiny.port = 8080)

# https://shiny.rstudio.com/articles/dynamic-ui.html
# https://stackoverflow.com/questions/21465411/r-shiny-passing-reactive-to-selectinput-choices/21467399#21467399
# https://stackoverflow.com/questions/33973300/issue-in-dynamic-renderui-in-shiny
# https://drsimonj.svbtle.com/exploring-correlations-in-r-with-corrr
# https://hssgenomics.shinyapps.io/RNAseq_DRaMA/
# https://stackoverflow.com/questions/28829682/r-shiny-checkboxgroupinput-select-all-checkboxes-by-click

primary <- read_csv("out/work/primary/merged-mtic.csv") %>% arrange(KEGG)
compounds <- unique(primary[c("Name","KEGG","HMDB","ChEBI","Category")])
lipids <- read_csv("out/work/lipids/merged-lipids.csv", col_types = cols(Saturation = "c", Polarity = "f")) %>% arrange(LMID)
lipid_ids <- unique(lipids[c("LMID","Name","InChIKey","Category","MainClass","Saturation")])

pops = c("Pachon","Tinaja","Surface")
tissues = c("Muscle","Brain","Liver")
conditions = c("4d Starved", "30d Starved", "Refed")

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
ui <- fluidPage(
  titlePanel("Rohner Lab Astyanax Metabolomics Study"),
  navlistPanel(
    "Primary",
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
      dataTableOutput("primary_summary")
    ),
    tabPanel("Correlation",
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
            column(3,
              numericInput(inputId = "primaryCorrPlotSampleNumClusters", label = "Number of Clusters:", value = 3, min = 1, width="50%", step = 1),
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
              column(3,
                numericInput(inputId = "primaryCorrPlotFeatureNumClusters", label = "Number of Clusters:", value = 3, min = 1, width="50%", step = 1),
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
              column(3,
                numericInput(inputId = "primaryCorrPlotCategoryNumClusters", label = "Number of Clusters:", value = 3, min = 1, width="50%", step = 1),
                checkboxInput(inputId = "primaryCorrPlotCategoryNormalize", label = "Normalize?", value = FALSE),
                checkboxInput(inputId = "primaryCorrPlotCategoryIncludeOutliers", label = "Include Outliers?", value = FALSE),
                ),
            checkboxGroupInput("primaryCorrPlotCategorySelectPops", "Populations:", pops, selected = pops),
            checkboxGroupInput("primaryCorrPlotCategorySelectTissues", "Tissues:", tissues, selected = tissues),
            checkboxGroupInput("primaryCorrPlotCategorySelectConditions", "Conditions:", conditions, selected = conditions)
          )
        ))
      ),
#       p(class = 'text-center',
#         savePlotlyPDFUI(id = "download_primarySampleCorrPlot",
#                         label = "Download Plot (PDF)")
#       )
    ),
    tabPanel("Component 4"),
    "-----",
    "Lipids",
    tabPanel("Selections",
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
    tabPanel("Summary",
      dataTableOutput("lipids_summary")
    ),
    tabPanel("Correlation",
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
            column(3,
              numericInput(inputId = "lipidsCorrPlotSampleNumClusters", label = "Number of Clusters:", value = 3, min = 1, width="50%", step = 1),
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
              column(3,
                numericInput(inputId = "lipidsCorrPlotFeatureNumClusters", label = "Number of Clusters:", value = 3, min = 1, width="50%", step = 1),
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
              column(3,
                numericInput(inputId = "lipidsCorrPlotCategoryNumClusters", label = "Number of Clusters:", value = 3, min = 1, width="50%", step = 1),
                checkboxInput(inputId = "lipidsCorrPlotCategoryNormalize", label = "Normalize?", value = FALSE),
                checkboxInput(inputId = "lipidsCorrPlotCategoryIncludeOutliers", label = "Include Outliers?", value = FALSE),
                ),
            checkboxGroupInput("lipidsCorrPlotCategorySelectPops", "Populations:", pops, selected = pops),
            checkboxGroupInput("lipidsCorrPlotCategorySelectTissues", "Tissues:", tissues, selected = tissues),
            checkboxGroupInput("lipidsCorrPlotCategorySelectConditions", "Conditions:", conditions, selected = conditions)
          )
        ))
      )
    )
  )
)


server <- function(input, output) {

  ###################################################################
  #                            Primary                              #
  ###################################################################

  selected_cpds <- reactive(filter(compounds, compounds$Name %in% input$selector_Name | compounds$KEGG %in% input$selector_KEGG | compounds$HMDB %in% input$selector_HMDB | compounds$ChEBI %in% input$selector_ChEBI |  compounds$Category %in% input$selector_Category))

  output$primary_summary <- renderDataTable(selected_cpds())

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
      thecor <- normalize(cor(features))
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
  #                            Lipids                               #
  ###################################################################

  selected_lipids <- reactive(filter(lipid_ids, lipid_ids$LMID %in% input$lipid_selector_LMID| lipid_ids$Name %in% input$lipid_selector_Name | lipid_ids$InChIKey %in% input$lipid_selector_InChIKey | lipid_ids$Category %in% input$lipid_selector_Category |  lipid_ids$MainClass %in% input$lipid_selector_MainClass |  lipid_ids$Saturation %in% input$lipid_selector_Saturation))

  output$lipids_summary <- renderDataTable(selected_lipids())

  ###################################################################
  #           Lipids Sample Correlation Plot                        #
  ###################################################################
  output$lipidsSampleCorrPlt <- renderPlotly({
#     https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
    cpd_data <- lipids %>% filter(lipids$LMID %in% selected_lipids()$LMID) %>% filter(Population %in% input$lipidsCorrPlotSampleSelectPops) %>% filter(Tissue %in% input$lipidsCorrPlotSampleSelectTissues) %>% filter(Condition %in% input$lipidsCorrPlotSampleSelectConditions)
    if (!input$lipidsCorrPlotSampleIncludeOutliers) {
      cpd_data <- cpd_data %>% filter(Outlier == FALSE)
    }
    features <- cpd_data %>% select(LMID,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean)
    features <- features %>% drop_na()
    print(features)
    features <- as.data.frame(features)
    rownames(features) <- features$LMID
    features <- features[,-1]
    if (input$lipidsCorrPlotSampleNormalize) {
      thecor <- normalize(cor(features))
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
    features <- cpd_data %>% select(LMID,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean)
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
    features <- cpd_data %>% select(Category,Population,Tissue,Condition,Raw_mTIC) %>% pivot_wider(names_from=c("Population","Tissue","Condition"),values_from="Raw_mTIC",values_fn = mean)
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

}

shinyApp(ui = ui, server = server)
