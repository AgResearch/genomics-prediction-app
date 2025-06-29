# R/ui_components.R
# UI components for Genomic Prediction Shiny App

# Custom CSS styles
get_custom_css <- function() {
  tags$head(
    tags$style(HTML("
      .objective-section {
        margin: 15px 0; padding: 15px; background-color: #ffffff;
        border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        border-left: 4px solid #27ae60;
      }
      .objective-text { font-size: 16px; line-height: 1.6; color: #2c3e50; }
      .data-highlight {
        background-color: #e3f2fd; padding: 2px 6px; border-radius: 3px;
        font-weight: 600; color: #1565c0;
      }
      .conclusion-section {
        margin: 15px 0; padding: 15px; background-color: #ffffff;
        border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        border-left: 4px solid #2196F3;
      }
    "))
  )
}

# Study Overview Tab Content
get_study_overview_tab <- function() {
  tabItem(tabName = "studydetails",
          fluidRow(
            box(
              title = "Study Objectives",
              status = "primary", solidHeader = TRUE, width = 12,
              div(class = "objective-section",
                  div(class = "objective-text",
                      tags$strong("Objective 1:"), " Evaluate PCA effectiveness in reducing microbiome dimensionality while retaining essential information."
                  )
              ),
              div(class = "objective-section",
                  div(class = "objective-text",
                      tags$strong("Objective 2:"), " Assess whether incorporating PCA-reduced microbiome data as intermediate traits in NN-GBLUP improves genomic prediction accuracy for methane emissions and feed efficiency traits."
                  )
              )
            )
          ),
          fluidRow(
            box(
              title = "Material & Method",
              status = "info", solidHeader = TRUE, width = 12,
              fluidRow(
                column(width = 6,
                       div(class = "objective-section",
                           h4("🐑 Animals & Genomics"),
                           div(class = "objective-text",
                               tags$ul(
                                 tags$li("Methane Group: ", tags$span("1,051 sheep", class = "data-highlight"), " on grass diet"),
                                 tags$li("RFI Group: ", tags$span("984 sheep", class = "data-highlight"), " on lucerne diet"),
                                 tags$li("Genomic Data: ", tags$span("11,967 SNPs", class = "data-highlight"), " from 15K panel"),
                                 tags$li("Training: 2014-2015 | Test: 2016")
                               )
                           )
                       )
                ),
                column(width = 6,
                       div(class = "objective-section",
                           h4("🦠 Microbiome & Analysis"),
                           div(class = "objective-text",
                               tags$ul(
                                 tags$li("Methane: ", tags$span("240,743", class = "data-highlight"), " microbial sequences"),
                                 tags$li("RFI: ", tags$span("207,393", class = "data-highlight"), " microbial sequences"),
                                 tags$li("PCA: ", tags$span("25-100%", class = "data-highlight"), " variance thresholds"),
                                 tags$li("Models: G, GM, NN-GBLUP")
                               )
                           )
                       )
                )
              )
            )
          ),
          fluidRow(
            box(
              title = "Key Findings & Conclusions",
              status = "success", solidHeader = TRUE, width = 12,
              div(class = "conclusion-section",
                  div(class = "objective-text",
                      tags$strong("Conclusion 1:"), " PCA effectively captured microbial information, with 100% PCA components yielding microbiability estimates comparable to the full dataset."
                  )
              ),
              div(class = "conclusion-section",
                  div(class = "objective-text",
                      tags$strong("Conclusion 2:"), " Improves accuracy for methane emission and residual feed intake (44-91% methane, 26-38% RFI increases). This approach enables genomic selection strategies to reduce methane emissions from ruminant livestock, contributing to sustainable agricultural practices."
                  )
              ),
              div(style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-top: 15px; text-align: center;",
                  tags$em(style = "color: #6c757d; font-size: 14px;", "All detailed results are available in the paper")
              )
            )
          )
  )
}

# Main content tabs function
get_main_analysis_tabs <- function() {
  list(
    # PCA Analysis Tab
    tabItem(tabName = "pca",
            fluidRow(
              box(title = "Figure 1: Principal Component Analysis", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("pcaPlot", height = "500px"))
            ),
            fluidRow(
              box(title = "PCA Summary", status = "info", solidHeader = TRUE, 
                  width = 12, verbatimTextOutput("pcaSummary"))
            )
    ),
    
    # Microbiability Tab
    tabItem(tabName = "microbiability",
            fluidRow(
              box(title = "Figure 2: Microbiability Estimates", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("microbiabilityPlot", height = "500px"))
            ),
            fluidRow(
              box(title = "Microbiability Summary", status = "info", solidHeader = TRUE, 
                  width = 12, verbatimTextOutput("microbiabilitySummary"))
            )
    ),
    
    # Train-Test Results
    tabItem(tabName = "traintest",
            fluidRow(
              box(title = "Figure 3: Methane Group (Train-Test)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("trainTestMethanePlot", height = "600px"))
            ),
            fluidRow(
              box(title = "Figure 4: RFI Group (Train-Test)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("trainTestRFIPlot", height = "600px"))
            )
    ),
    
    # Five-Fold Results
    tabItem(tabName = "fivefold",
            fluidRow(
              box(title = "Figure 5: Methane Group (Five-Fold)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("fiveFoldMethanePlot", height = "600px"))
            ),
            fluidRow(
              box(title = "Figure 6: RFI Group (Five-Fold)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("fiveFoldRFIPlot", height = "600px"))
            )
    )
  )
}
# Supplementary Tables and Figures Tabs
get_supplementary_tabs <- function() {
  list(
    # Tables
    tabItem(tabName = "tables1",
            fluidRow(
              box(title = "Table S1: Train-Test Validation Results", status = "primary", 
                  solidHeader = TRUE, width = 12, DTOutput("tableS1"))
            )
    ),
    
    tabItem(tabName = "tables2",
            fluidRow(
              box(title = "Table S2: Five-Fold Cross-Validation Results", status = "primary", 
                  solidHeader = TRUE, width = 12, DTOutput("tableS2"))
            )
    ),
    
    # Supplementary Figures S1-S4
    tabItem(tabName = "figures1",
            fluidRow(
              box(title = "Figure S1: GM Model Methane Group (Train-Test)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS1Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures2",
            fluidRow(
              box(title = "Figure S2: GM Model RFI Group (Train-Test)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS2Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures3",
            fluidRow(
              box(title = "Figure S3: GM Model Methane Group (Five-Fold)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS3Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures4",
            fluidRow(
              box(title = "Figure S4: GM Model RFI Group (Five-Fold)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS4Plot", height = "600px"))
            )
    ),
    
    # Figures S5-S7
    tabItem(tabName = "figures5",
            fluidRow(
              box(title = "Figure S5: Scatter Plot Comparison (Train-Test)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS5Plot", height = "700px"))
            )
    ),
    
    tabItem(tabName = "figures6",
            fluidRow(
              box(title = "Figure S6: PCA Variance Comparison (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS6Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures7",
            fluidRow(
              box(title = "Figure S7: Microbiability Comparison (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS7Plot", height = "800px"))
            )
    )
  )
}

# Complete tabItems function that combines everything
# REPLACE the get_tab_items() function in R/ui_components.R with this:

# Complete tabItems function that combines everything
get_tab_items <- function() {
  tabItems(
    # Study Overview Tab
    get_study_overview_tab(),
    
    # Main Analysis Tabs  
    tabItem(tabName = "pca",
            fluidRow(
              box(title = "Figure 1: Principal Component Analysis", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("pcaPlot", height = "500px"))
            ),
            fluidRow(
              box(title = "PCA Summary", status = "info", solidHeader = TRUE, 
                  width = 12, verbatimTextOutput("pcaSummary"))
            )
    ),
    
    tabItem(tabName = "microbiability",
            fluidRow(
              box(title = "Figure 2: Microbiability Estimates", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("microbiabilityPlot", height = "500px"))
            ),
            fluidRow(
              box(title = "Microbiability Summary", status = "info", solidHeader = TRUE, 
                  width = 12, verbatimTextOutput("microbiabilitySummary"))
            )
    ),
    
    tabItem(tabName = "traintest",
            fluidRow(
              box(title = "Figure 3: Methane Group (Train-Test)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("trainTestMethanePlot", height = "600px"))
            ),
            fluidRow(
              box(title = "Figure 4: RFI Group (Train-Test)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("trainTestRFIPlot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "fivefold",
            fluidRow(
              box(title = "Figure 5: Methane Group (Five-Fold)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("fiveFoldMethanePlot", height = "600px"))
            ),
            fluidRow(
              box(title = "Figure 6: RFI Group (Five-Fold)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("fiveFoldRFIPlot", height = "600px"))
            )
    ),
    
    # Supplementary Tables and Figures
    tabItem(tabName = "tables1",
            fluidRow(
              box(title = "Table S1: Train-Test Validation Results", status = "primary", 
                  solidHeader = TRUE, width = 12, DTOutput("tableS1"))
            )
    ),
    
    tabItem(tabName = "tables2",
            fluidRow(
              box(title = "Table S2: Five-Fold Cross-Validation Results", status = "primary", 
                  solidHeader = TRUE, width = 12, DTOutput("tableS2"))
            )
    ),
    
    tabItem(tabName = "figures1",
            fluidRow(
              box(title = "Figure S1: GM Model Methane Group (Train-Test)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS1Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures2",
            fluidRow(
              box(title = "Figure S2: GM Model RFI Group (Train-Test)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS2Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures3",
            fluidRow(
              box(title = "Figure S3: GM Model Methane Group (Five-Fold)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS3Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures4",
            fluidRow(
              box(title = "Figure S4: GM Model RFI Group (Five-Fold)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS4Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures5",
            fluidRow(
              box(title = "Figure S5: Scatter Plot Comparison (Train-Test)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS5Plot", height = "700px"))
            )
    ),
    
    tabItem(tabName = "figures6",
            fluidRow(
              box(title = "Figure S6: PCA Variance Comparison (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS6Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures7",
            fluidRow(
              box(title = "Figure S7: Microbiability Comparison (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS7Plot", height = "800px"))
            )
    ),
    
    tabItem(tabName = "figures8",
            fluidRow(
              box(title = "Figure S8: GM Model Methane Group (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS8Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures9",
            fluidRow(
              box(title = "Figure S9: GM Model RFI Group (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS9Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures10",
            fluidRow(
              box(title = "Figure S10: NN-GBLUP Accuracy (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS10Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures11",
            fluidRow(
              box(title = "Figure S11: NN-GBLUP Bias (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS11Plot", height = "600px"))
            )
    )
  )
}
# Remaining Supplementary Figures S8-S11
get_remaining_figures_tabs <- function() {
  list(
    tabItem(tabName = "figures8",
            fluidRow(
              box(title = "Figure S8: GM Model Methane Group (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS8Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures9",
            fluidRow(
              box(title = "Figure S9: GM Model RFI Group (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS9Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures10",
            fluidRow(
              box(title = "Figure S10: NN-GBLUP Accuracy (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS10Plot", height = "600px"))
            )
    ),
    
    tabItem(tabName = "figures11",
            fluidRow(
              box(title = "Figure S11: NN-GBLUP Bias (Log vs CLR)", status = "primary", 
                  solidHeader = TRUE, width = 12, plotOutput("figureS11Plot", height = "600px"))
            )
    )
  )
}

# Additional data loading functions for supplementary figures
load_scatter_data <- function() {
  reactive({
    scatter_path <- "data/supplementarymaterial/Scatterplot"
    gm_data <- read.csv(file.path(scatter_path, "all_predictionsGM.csv"), header = TRUE)
    nngblup_data <- read.csv(file.path(scatter_path, "All_Traits_Predictions_NNGBLUP.csv"), header = TRUE)
    
    names(nngblup_data)[5] <- "Observed"
    names(nngblup_data)[1] <- "Group"
    gm_data <- gm_data[,-5]
    nngblup_data <- nngblup_data[, c(1,3,2,4,5,6)]
    gm_data$Model <- "GM Model"
    nngblup_data$Model <- "NN-GBLUP"
    
    list(gm = gm_data, nngblup = nngblup_data)
  })
}

# Load PCA transformation comparison data for Figure S6
load_pca_transformation_data <- function() {
  reactive({
    pca_path <- "data/PCA"
    methane_pca_log <- read.csv(file.path(pca_path, "combined_pca_variance_explained_methane.csv"), header = TRUE)
    rfi_pca_log <- read.csv(file.path(pca_path, "combined_pca_variance_explained_RFI.csv"), header = TRUE)
    methane_pca_clr <- read.csv(file.path(pca_path, "pca_explained_variancemethanegroup.csv"), header = TRUE)
    rfi_pca_clr <- read.csv(file.path(pca_path, "pca_explained_varianceRFIgroup.csv"), header = TRUE)
    
    process_pc_column <- function(data) {
      if(is.character(data$PC)) {
        data$PC <- as.numeric(gsub("PC", "", data$PC))
      }
      return(data)
    }
    
    methane_pca_log <- process_pc_column(methane_pca_log)
    rfi_pca_log <- process_pc_column(rfi_pca_log)
    methane_pca_clr <- process_pc_column(methane_pca_clr)
    rfi_pca_clr <- process_pc_column(rfi_pca_clr)
    
    fix_variance_scale <- function(data) {
      if(max(data$CumulativeVariance, na.rm = TRUE) <= 1) {
        data$CumulativeVariance <- data$CumulativeVariance * 100
      }
      return(data)
    }
    
    methane_pca_log <- fix_variance_scale(methane_pca_log)
    rfi_pca_log <- fix_variance_scale(rfi_pca_log)
    methane_pca_clr <- fix_variance_scale(methane_pca_clr)
    rfi_pca_clr <- fix_variance_scale(rfi_pca_clr)
    
    methane_pca_log$Transformation <- "Log"
    methane_pca_clr$Transformation <- "CLR"
    rfi_pca_log$Transformation <- "Log"
    rfi_pca_clr$Transformation <- "CLR"
    
    methane_combined <- rbind(methane_pca_log, methane_pca_clr)
    rfi_combined <- rbind(rfi_pca_log, rfi_pca_clr)
    
    list(methane = methane_combined, rfi = rfi_combined)
  })
}