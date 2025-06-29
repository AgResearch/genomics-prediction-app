# app.R
# Main Shiny Application for Genomic Prediction Analysis

# Load required libraries
library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(grid)
library(shinyjs)

# Source helper functions
source("R/data_loading.R")
source("R/plotting_functions.R")
source("R/ui_components.R")

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "Genomic Prediction Analysis"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Study Overview", tabName = "studydetails", icon = icon("info-circle")),
      menuItem("PCA Analysis", tabName = "pca", icon = icon("chart-line")),
      menuItem("Microbiability", tabName = "microbiability", icon = icon("microscope")),
      menuItem("Genomic Prediction", 
               menuSubItem("Train-Test Results", tabName = "traintest", icon = icon("chart-bar")),
               menuSubItem("Five-Fold Results", tabName = "fivefold", icon = icon("chart-area"))),
      menuItem("Additional Files", 
               menuSubItem("Table S1", tabName = "tables1", icon = icon("table")),
               menuSubItem("Table S2", tabName = "tables2", icon = icon("table")),
               menuSubItem("Figure S1", tabName = "figures1", icon = icon("chart-bar")),
               menuSubItem("Figure S2", tabName = "figures2", icon = icon("chart-bar")),
               menuSubItem("Figure S3", tabName = "figures3", icon = icon("chart-area")),
               menuSubItem("Figure S4", tabName = "figures4", icon = icon("chart-area")),
               menuSubItem("Figure S5", tabName = "figures5", icon = icon("chart-line")),
               menuSubItem("Figure S6", tabName = "figures6", icon = icon("line-chart")),
               menuSubItem("Figure S7", tabName = "figures7", icon = icon("chart-column")),
               menuSubItem("Figure S8", tabName = "figures8", icon = icon("chart-simple")),
               menuSubItem("Figure S9", tabName = "figures9", icon = icon("chart-simple")),
               menuSubItem("Figure S10", tabName = "figures10", icon = icon("chart-gantt")),
               menuSubItem("Figure S11", tabName = "figures11", icon = icon("chart-gantt")))
    )
  ),
  dashboardBody(
    useShinyjs(),
    # Custom CSS styles
    tags$head(
      tags$script(async = NA, src = "https://www.googletagmanager.com/gtag/js?id=G-W6HVXEG19E"),
      tags$script(HTML("
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-W6HVXEG19E');
    "))
    ),
    get_custom_css(),
    # Tab items content
    get_tab_items()
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Load all data reactively
  pca_data <- load_pca_data()
  genomic_data <- load_genomic_data()
  supplementary_data <- load_supplementary_data()
  microbiability_data <- load_microbiability_data()
  gm_supplementary_data <- load_gm_supplementary_data()
  scatter_data <- load_scatter_data()
  pca_transformation_data <- load_pca_transformation_data()
  microbiability_transformation_data <- load_microbiability_transformation_data()
  gm_clr_comparison_data <- load_gm_clr_comparison_data()
  nngblup_clr_comparison_data <- load_nngblup_clr_comparison_data()
  
  # PCA Analysis outputs
  output$pcaPlot <- renderPlot({
    create_pca_plot(pca_data())
  })
  
  output$pcaSummary <- renderText({
    create_pca_summary()
  })
  
  # Microbiability outputs
  output$microbiabilityPlot <- renderPlot({
    create_microbiability_plot(microbiability_data())
  })
  
  output$microbiabilitySummary <- renderText({
    create_microbiability_summary(microbiability_data())
  })
  
  # Train-Test plots
  output$trainTestMethanePlot <- renderPlot({
    genomic_list <- genomic_data()
    df_methane <- genomic_list$traintest %>% filter(Trait %in% c("Methane", "Methane ratio", "LWT", "CO2"))
    plot_combined_metrics_separate(df_methane)
  })
  
  output$trainTestRFIPlot <- renderPlot({
    genomic_list <- genomic_data()
    df_rfi <- genomic_list$traintest %>% filter(Trait %in% c("RFI", "Mid intake"))
    plot_combined_metrics_separate_rfi(df_rfi)
  })
  
  # Five-Fold plots
  output$fiveFoldMethanePlot <- renderPlot({
    genomic_list <- genomic_data()
    df_methane <- genomic_list$fivefold %>% filter(Trait %in% c("Methane", "Methane ratio", "LWT"))
    plot_combined_metrics_separate(df_methane)
  })
  
  output$fiveFoldRFIPlot <- renderPlot({
    genomic_list <- genomic_data()
    df_rfi <- genomic_list$fivefold %>% filter(Trait %in% c("RFI", "Mid intake"))
    desired_levels <- c("G", "GM", "PC74", "PC299", "PC600")
    df_rfi_filtered <- df_rfi %>% filter(pc %in% desired_levels)
    plot_combined_metrics_separate_rfi(df_rfi_filtered)
  })
  # Render supplementary tables
  output$tableS1 <- renderDT({
    supp_data <- supplementary_data()
    datatable(supp_data$table_s1,
              options = list(pageLength = 15, autoWidth = TRUE, scrollX = TRUE),
              rownames = FALSE, class = 'cell-border stripe hover compact',
              filter = 'top', style = 'bootstrap4')
  })
  
  output$tableS2 <- renderDT({
    supp_data <- supplementary_data()
    datatable(supp_data$table_s2,
              options = list(pageLength = 15, autoWidth = TRUE, scrollX = TRUE),
              rownames = FALSE, class = 'cell-border stripe hover compact',
              filter = 'top', style = 'bootstrap4')
  })
  
  # Figure S1-S4 outputs (GM Model)
  output$figureS1Plot <- renderPlot({
    gm_data <- gm_supplementary_data()
    df_processed <- gm_data$traintest %>%
      mutate(Trait = dplyr::recode(Trait, "MethaneRatio" = "Methane ratio"))
    plot_gm_metrics_methane(df_processed)
  })
  
  output$figureS2Plot <- renderPlot({
    gm_data <- gm_supplementary_data()
    df_processed <- gm_data$traintest %>%
      mutate(Trait = dplyr::recode(Trait, "AMidIntake" = "Mid intake"))
    plot_gm_metrics_rfi(df_processed)
  })
  
  output$figureS3Plot <- renderPlot({
    gm_data <- gm_supplementary_data()
    df_processed <- gm_data$fivefold %>%
      mutate(Trait = dplyr::recode(Trait, "MethaneRatio" = "Methane ratio"))
    plot_gm_metrics_methane(df_processed)
  })
  
  output$figureS4Plot <- renderPlot({
    gm_data <- gm_supplementary_data()
    plot_gm_metrics_rfi_fivefold(gm_data$fivefold)
  })
  
  # Figure S5: Scatter plots
  output$figureS5Plot <- renderPlot({
    scatter_list <- scatter_data()
    create_scatter_plot(scatter_list)
  })
  
  # Figure S6: PCA Log vs CLR Transformation Comparison
  output$figureS6Plot <- renderPlot({
    transformation_data <- pca_transformation_data()
    create_pca_transformation_plot(transformation_data)
  })
  
  # Figure S7: Microbiability Log vs CLR Comparison  
  output$figureS7Plot <- renderPlot({
    micro_transform_data <- microbiability_transformation_data()
    
    methane_pc_mapping <- c("88", "102", "333", "346", "640", "648", "933", "936", "1030", "1030", "wholedata", "GM")
    rfi_pc_mapping <- c("74", "83", "299", "310", "600", "606", "889", "890", "979", "984", "wholedata", "GM")
    
    methane_variance_data <- create_variance_data(methane_pc_mapping, micro_transform_data$methane_log, micro_transform_data$methane_clr, "Methane")
    methane_ratio_variance_data <- create_variance_data(methane_pc_mapping, micro_transform_data$methane_log, micro_transform_data$methane_clr, "Methane ratio")
    lwt_variance_data <- create_variance_data(methane_pc_mapping, micro_transform_data$methane_log, micro_transform_data$methane_clr, "LWT")
    co2_variance_data <- create_variance_data(methane_pc_mapping, micro_transform_data$methane_log, micro_transform_data$methane_clr, "CO2")
    rfi_variance_data <- create_variance_data(rfi_pc_mapping, micro_transform_data$rfi_log, micro_transform_data$rfi_clr, "RFI")
    mid_intake_variance_data <- create_variance_data(rfi_pc_mapping, micro_transform_data$rfi_log, micro_transform_data$rfi_clr, "Mid intake")
    
    p_methane <- create_variance_plot(methane_variance_data, "Methane")
    p_methane_ratio <- create_variance_plot(methane_ratio_variance_data, "Methane ratio")
    p_lwt <- create_variance_plot(lwt_variance_data, "LWT")
    p_co2 <- create_variance_plot(co2_variance_data, "CO2")
    p_rfi <- create_variance_plot(rfi_variance_data, "RFI")
    p_mid_intake <- create_variance_plot(mid_intake_variance_data, "Mid intake")
    
    grid.arrange(p_methane, p_methane_ratio, p_lwt, p_co2, p_rfi, p_mid_intake, ncol = 2, nrow = 3)
  })
  # Figure S8: GM Model Methane Group Log vs CLR
  output$figureS8Plot <- renderPlot({
    gm_clr_data <- gm_clr_comparison_data()
    methane_variance_mapping <- data.frame(
      PC_Log = c("88", "333", "640", "933", "GM"),
      PC_CLR = c("102", "346", "648", "936", "936"),
      Variance = c("25%", "50%", "75%", "95%", "GM"), stringsAsFactors = FALSE
    )
    methane_variance_data <- prepare_variance_data_gm(gm_clr_data$methane_log, gm_clr_data$methane_clr, methane_variance_mapping)
    create_facet_plots_methane_gm(methane_variance_data)
  })
  
  # Figure S9: GM Model RFI Group Log vs CLR
  output$figureS9Plot <- renderPlot({
    gm_clr_data <- gm_clr_comparison_data()
    rfi_variance_mapping <- data.frame(
      PC_Log = c("74", "299", "600", "889", "GM"),
      PC_CLR = c("83", "310", "606", "890", "890"),
      Variance = c("25%", "50%", "75%", "95%", "GM"), stringsAsFactors = FALSE
    )
    rfi_variance_data <- prepare_variance_data_gm(gm_clr_data$rfi_log, gm_clr_data$rfi_clr, rfi_variance_mapping)
    create_facet_plots_rfi_gm(rfi_variance_data)
  })
  
  # Figure S10: NN-GBLUP All Traits Accuracy Comparison
  output$figureS10Plot <- renderPlot({
    nngblup_data <- nngblup_clr_comparison_data()
    create_all_traits_accuracy_plot(nngblup_data)
  })
  
  # Figure S11: NN-GBLUP All Traits Bias Comparison
  output$figureS11Plot <- renderPlot({
    nngblup_data <- nngblup_clr_comparison_data()
    create_all_traits_bias_plot(nngblup_data)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
