# Genomic Prediction Analysis Shiny App

An interactive R Shiny application for analyzing genomic prediction results in sheep, focusing on methane emissions and feed efficiency traits using microbiome-enhanced genomic prediction models.

## 🎯 Research Objectives

**Objective 1**: Evaluate PCA effectiveness in reducing rumen microbiome data dimensionality while retaining essential biological information.

**Objective 2**: Assess whether incorporating PCA-reduced microbiome data as intermediate traits in Neural Network GBLUP improves genomic prediction accuracy for methane emissions and feed efficiency traits.

## 🚀 Live Application

Access the interactive application at: https://setalemu.shinyapps.io/genomics_app_deploy/

## 📋 Features

- Interactive visualizations of PCA analysis
- Microbiability estimates across PC thresholds  
- Genomic prediction results (Train-Test and Five-Fold validation)
- Model comparisons (G, GM, NN-GBLUP)
- Downloadable plots and supplementary materials
- Comprehensive data tables

## 🛠️ Local Installation

```r
# Install required packages
packages <- c("shiny", "shinydashboard", "dplyr", "tidyr", 
              "ggplot2", "DT", "RColorBrewer", "gridExtra", 
              "cowplot", "grid", "shinyjs")
install.packages(packages)

# Run the application
shiny::runApp()
