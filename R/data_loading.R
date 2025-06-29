# R/data_loading.R
# Data loading functions for Genomic Prediction Shiny App

# Load PCA data
load_pca_data <- function() {
  reactive({
    pca_path <- "data/PCA"
    methane_pca <- read.csv(file.path(pca_path, "combined_pca_variance_explained_methane.csv"), header = TRUE)
    rfi_pca <- read.csv(file.path(pca_path, "combined_pca_variance_explained_RFI.csv"), header = TRUE)
    list(methane = methane_pca, rfi = rfi_pca)
  })
}

# Load genomic prediction data
load_genomic_data <- function() {
  reactive({
    base_path <- "data/final_combined"
    traintest_data <- read.csv(file.path(base_path, "TrainTest_Combined2.csv"), header = TRUE)
    fivefold_data <- read.csv(file.path(base_path, "FiveFold_Combined2.csv"), header = TRUE)
    
    # Apply trait name recoding
    traintest_data$Trait <- dplyr::recode(traintest_data$Trait, 
                                          "MethaneRatio" = "Methane ratio", "AMidIntake" = "Mid intake")
    fivefold_data$Trait <- dplyr::recode(fivefold_data$Trait, 
                                         "MethaneRatio" = "Methane ratio", "AMidIntake" = "Mid intake")
    
    list(traintest = traintest_data, fivefold = fivefold_data)
  })
}

# Load supplementary tables
# Function to format values with ± standard error (3 significant figures)
format_with_se <- function(value, se) {
  # Format both value and SE to 3 significant figures
  formatted_value <- signif(value, 3)
  formatted_se <- signif(se, 3)
  
  # Handle special formatting for very small numbers
  if (abs(formatted_value) < 0.001 && formatted_value != 0) {
    value_str <- sprintf("%.3e", formatted_value)
  } else if (abs(formatted_value) < 0.1) {
    value_str <- sprintf("%.4f", formatted_value)
  } else if (abs(formatted_value) < 1) {
    value_str <- sprintf("%.3f", formatted_value)
  } else if (abs(formatted_value) < 10) {
    value_str <- sprintf("%.2f", formatted_value)
  } else {
    value_str <- sprintf("%.1f", formatted_value)
  }
  
  # Format SE similarly
  if (abs(formatted_se) < 0.001 && formatted_se != 0) {
    se_str <- sprintf("%.3e", formatted_se)
  } else if (abs(formatted_se) < 0.1) {
    se_str <- sprintf("%.4f", formatted_se)
  } else if (abs(formatted_se) < 1) {
    se_str <- sprintf("%.3f", formatted_se)
  } else if (abs(formatted_se) < 10) {
    se_str <- sprintf("%.2f", formatted_se)
  } else {
    se_str <- sprintf("%.1f", formatted_se)
  }
  
  return(paste0(value_str, " ± ", se_str))
}

# FIXED Load supplementary tables function
load_supplementary_data <- function() {
  reactive({
    supplementary_path <- "data/supplementarymaterial"
    
    # Use the CORRECT filenames that actually exist in your directory
    table_s1_raw <- read.csv(file.path(supplementary_path, "table_s1.csv"), header = TRUE)
    table_s2_raw <- read.csv(file.path(supplementary_path, "table_s2.csv"), header = TRUE)
    
    # Function to process table data into display format
    process_table_data <- function(raw_data) {
      processed_data <- raw_data %>%
        mutate(
          `Genomic Prediction Accuracy` = mapply(format_with_se, Accuracy, Accuracy_SE),
          `Dispersion Bias` = mapply(format_with_se, Bias, Bias_SE)
        ) %>%
        select(
          `Trait Group` = Trait_Group,
          Model,
          Trait,
          `Genomic Prediction Accuracy`,
          `Dispersion Bias`
        )
      
      return(processed_data)
    }
    
    # Process both tables
    table_s1_processed <- process_table_data(table_s1_raw)
    table_s2_processed <- process_table_data(table_s2_raw)
    
    return(list(table_s1 = table_s1_processed, table_s2 = table_s2_processed))
  })
}

# Load microbiability data
load_microbiability_data <- function() {
  reactive({
    micro_path <- "data/Microbiability"
    methane_group <- read.table(file.path(micro_path, "methanegroup.txt"), header = TRUE)
    rfi_group <- read.table(file.path(micro_path, "RFIgroup.txt"), header = TRUE)
    list(methane = methane_group, rfi = rfi_group)
  })
}

# Load GM model supplementary data
load_gm_supplementary_data <- function() {
  reactive({
    gm_path <- "data/Gmmodelonly"
    gm_traintest <- read.csv(file.path(gm_path, "TrainTest.csv"), header = TRUE)
    gm_fivefold <- read.csv(file.path(gm_path, "Fivefold.csv"), header = TRUE)
    gm_fivefold$pc <- ifelse(substring(gm_fivefold$pc, 1, 2) != "GM", 
                             substring(gm_fivefold$pc, 3, nchar(gm_fivefold$pc)), 
                             gm_fivefold$pc)
    gm_fivefold$Trait[gm_fivefold$Trait == "AMidIntake"] <- "Mid intake"
    list(traintest = gm_traintest, fivefold = gm_fivefold)
  })
}

