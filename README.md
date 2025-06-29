```markdown
# Genomic Prediction Analysis Shiny App

A comprehensive R Shiny application for analyzing genomic prediction results in livestock, specifically focusing on methane emissions and feed efficiency traits in sheep using microbiome-enhanced genomic prediction models.

## 🎯 Application Overview

This interactive web application visualizes the results of a genomic prediction study that evaluated the effectiveness of Principal Component Analysis (PCA) in reducing microbiome dimensionality while retaining essential information for improving genomic predictions.

## 🚀 Live Application

Access the live application at: https://setalemu.shinyapps.io/genomics_app_deploy/

## 📋 Features

- Interactive visualizations of PCA analysis
- Microbiability estimates across PC thresholds  
- Genomic prediction results (Train-Test and Five-Fold)
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
