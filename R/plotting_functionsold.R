# R/plotting_functions.R
# Plotting functions for Genomic Prediction Shiny App

# PCA plotting function
create_plot <- function(data, group_name, x_breaks) {
  data$CumulativeVariance <- data$CumulativeVariance
  thresholds <- c(25, 50, 75, 95, 100)
  threshold_data <- sapply(thresholds, function(t) {
    which.min(abs(data$CumulativeVariance - t))
  })
  
  p <- ggplot(data, aes(x = PC, y = CumulativeVariance)) +
    geom_line(size = 1, color = "skyblue") +
    geom_point(data = data[threshold_data, ], aes(x = PC, y = CumulativeVariance), color = "skyblue", size = 3) +
    geom_text(data = data[threshold_data, ], 
              aes(x = PC, y = CumulativeVariance, label = paste0(thresholds, "%")), 
              color = "blue", hjust = -0.2, vjust = -0.5, size = 4) +
    scale_x_continuous(breaks = x_breaks, labels = x_breaks, 
                       limits = c(0, max(data$PC) * 1.05)) +
    scale_y_continuous(breaks = seq(0, 100, by = 25), 
                       labels = function(x) paste0(x, "%"),
                       limits = c(0, 105)) +
    labs(title = group_name, x = "Number of Principal Components", y = "Cumulative Variance Explained (%)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(size = 10),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(10, 30, 10, 10))
  return(p)
}

# Create combined PCA plot
create_pca_plot <- function(pca_list) {
  methane_plot <- create_plot(pca_list$methane, "Methane Grass Group", c(88, 333, 640, 932, 1030))
  rfi_plot <- create_plot(pca_list$rfi, "RFI Lucerne Group", c(74, 299, 600, 889, 979))
  grid.arrange(methane_plot, rfi_plot, ncol = 2, widths = c(1.2, 1.2))
}

# Create PCA summary text
create_pca_summary <- function() {
  methane_thresholds <- c(25, 50, 75, 95, 100)
  rfi_thresholds <- c(25, 50, 75, 95, 100)
  
  # Hardcoded PC values for each threshold
  methane_pcs <- c(88, 333, 640, 932, 1030)
  rfi_pcs <- c(74, 299, 600, 889, 979)
  
  paste("PCA Summary:\n",
        "Methane Group - PC counts for variance thresholds:\n",
        paste("  ", methane_thresholds, "%: PC ", methane_pcs, collapse = "\n"), "\n\n",
        "RFI Group - PC counts for variance thresholds:\n",
        paste("  ", rfi_thresholds, "%: PC ", rfi_pcs, collapse = "\n"))
}

# Plotting function for methane group (main analysis)
plot_combined_metrics_separate <- function(df) {
  desired_levels <- c("G", "GM", "PC88", "PC333", "PC640", "PC933")
  df$pc <- factor(df$pc, levels = desired_levels)
  df$Trait <- factor(df$Trait, levels = c("Methane", "Methane ratio", "LWT", "CO2"))
  
  df_long <- df %>%
    pivot_longer(cols = c(Cor, Bias), names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("Cor", "Bias"),
                           labels = c("Genomic prediction accuracy", "Bias")))
  
  accuracy_plot <- df_long %>%
    filter(Metric == "Genomic prediction accuracy") %>%
    ggplot(aes(x = pc, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, nrow = 1) +
    scale_y_continuous(limits = c(0, 0.4)) +
    labs(x = "PC Level", y = "Genomic prediction accuracy") +
    scale_fill_manual(values = c("Genomic prediction accuracy" = "skyblue")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(color = "black", fill = NA))
  
  bias_plot <- df_long %>%
    filter(Metric == "Bias") %>%
    ggplot(aes(x = pc, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, nrow = 1) +
    scale_y_continuous(limits = c(0, 1.9)) +
    labs(x = "PC Level", y = "Bias") +
    scale_fill_manual(values = c("Bias" = "coral")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(color = "black", fill = NA))
  
  cowplot::plot_grid(accuracy_plot, bias_plot, ncol = 1, align = 'v', axis = 'l', rel_heights = c(1, 1))
  
}

# Plotting function for RFI group (main analysis)
plot_combined_metrics_separate_rfi <- function(df) {
  desired_levels <- c("G", "GM", "PC74", "PC299", "PC600", "PC889")
  df_filtered <- df %>% filter(pc %in% desired_levels, Trait %in% c("Mid intake", "RFI"))
  df_filtered$pc <- factor(df_filtered$pc, levels = desired_levels)
  df_filtered$Trait <- factor(df_filtered$Trait, levels = c("RFI", "Mid intake"))
  
  df_long <- df_filtered %>%
    pivot_longer(cols = c(Cor, Bias), names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("Cor", "Bias"),
                           labels = c("Genomic prediction accuracy", "Bias")))
  
  accuracy_plot <- df_long %>%
    filter(Metric == "Genomic prediction accuracy") %>%
    ggplot(aes(x = pc, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, nrow = 1) +
    scale_y_continuous(limits = c(0, 0.5)) +
    labs(x = "PC Level", y = "Genomic prediction accuracy") +
    scale_fill_manual(values = c("Genomic prediction accuracy" = "skyblue")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(color = "black", fill = NA))
  
  bias_plot <- df_long %>%
    filter(Metric == "Bias") %>%
    ggplot(aes(x = pc, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, nrow = 1) +
    scale_y_continuous(limits = c(0, 3.2)) +
    labs(x = "PC Level", y = "Bias") +
    scale_fill_manual(values = c("Bias" = "coral")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(color = "black", fill = NA))
  
  cowplot::plot_grid(accuracy_plot, bias_plot, ncol = 1, align = 'v', axis = 'l', rel_heights = c(1, 1))
}

# Microbiability plotting function
create_microbiability_plot <- function(micro_list) {
  # Process Methane group data
  methane_data <- micro_list$methane
  methane_data$PC_Count <- factor(methane_data$PC_Count,
                                  levels = c("88", "333", "640", "933", "1030", "wholedata"),
                                  labels = c("88 PCs (25%)", "333 PCs (50%)", "640 PCs (75%)", 
                                             "933 PCs (95%)", "1030 PCs (100%)", "Whole Data"))
  methane_data$Trait <- factor(methane_data$Trait, levels = c("Methane", "MethaneRatio", "LWT", "CO2"))  
  methane_data$Trait <- dplyr::recode(methane_data$Trait, "MethaneRatio" = "Methane ratio")
  
  # Create methane plot
  p3 <- ggplot(data = methane_data, aes(x = Trait, y = Microbability, fill = PC_Count)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
    geom_errorbar(aes(ymin = Microbability - SE, ymax = Microbability + SE), 
                  width = 0.2, position = position_dodge(width = 0.9)) +
    labs(x = "Trait", y = "Microbiability", title = "Methane Grass group") +
    theme_minimal() + scale_fill_brewer(palette = "Set3") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 0.8)
  
  # Process RFI group data
  rfi_data <- micro_list$rfi
  rfi_data$PC_Count <- factor(rfi_data$PC_Count,
                              levels = c("74", "299", "600", "889", "979", "wholedata"),
                              labels = c("74 PCs (25%)", "299 PCs (50%)", "600 PCs (75%)", 
                                         "889 PCs (95%)", "979 PCs (100%)", "Whole Data"))
  rfi_data$Trait <- factor(rfi_data$Trait, levels = c("RFI", "AMidIntake"))
  rfi_data$Trait <- dplyr::recode(rfi_data$Trait, "AMidIntake" = "Mid intake")
  
  # Create RFI plot
  p4 <- ggplot(data = rfi_data, aes(x = Trait, y = Microbability, fill = PC_Count)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
    geom_errorbar(aes(ymin = Microbability - SE, ymax = Microbability + SE), 
                  width = 0.2, position = position_dodge(width = 0.9)) +
    labs(x = "Trait", y = "Microbiability", title = "RFI Lucerne group") +
    theme_minimal() + scale_fill_brewer(palette = "Set3") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 0.8)
  
  grid.arrange(p3, p4, ncol = 2)
}

# Microbiability summary function
create_microbiability_summary <- function(micro_list) {
  methane_summary <- micro_list$methane %>%
    group_by(Trait) %>%
    summarise(Min_Microbiability = round(min(Microbability, na.rm = TRUE), 3),
              Max_Microbiability = round(max(Microbability, na.rm = TRUE), 3), .groups = 'drop')
  
  rfi_summary <- micro_list$rfi %>%
    group_by(Trait) %>%
    summarise(Min_Microbiability = round(min(Microbability, na.rm = TRUE), 3),
              Max_Microbiability = round(max(Microbability, na.rm = TRUE), 3), .groups = 'drop')
  
  paste("Microbiability Summary:\nMethane Group:\n",
        paste("  ", methane_summary$Trait, ": ", methane_summary$Min_Microbiability, " - ", methane_summary$Max_Microbiability, collapse = "\n"), 
        "\n\nRFI Group:\n",
        paste("  ", rfi_summary$Trait, ": ", rfi_summary$Min_Microbiability, " - ", rfi_summary$Max_Microbiability, collapse = "\n"))
}
# GM Model plotting function for methane group
plot_gm_metrics_methane <- function(df) {
  desired_levels <- c("GM", "88", "333", "640", "933")
  df_filtered <- df %>%
    filter(pc %in% desired_levels) %>%
    filter(Trait %in% c("Methane", "Methane ratio", "LWT", "CO2")) %>%
    mutate(Trait = factor(Trait, levels = c("Methane", "Methane ratio", "LWT", "CO2")),
           pc = factor(pc, levels = desired_levels))
  
  df_long <- df_filtered %>%
    pivot_longer(cols = c(Cor, Bias), names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("Cor", "Bias"),
                           labels = c("Genomic prediction accuracy", "Bias")))
  
  accuracy_plot <- df_long %>%
    filter(Metric == "Genomic prediction accuracy") %>%
    ggplot(aes(x = pc, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, nrow = 1) +
    scale_y_continuous(limits = c(0, 0.3)) +
    labs(x = "PC Level", y = "Genomic prediction accuracy") +
    scale_fill_manual(values = c("Genomic prediction accuracy" = "skyblue")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  bias_plot <- df_long %>%
    filter(Metric == "Bias") %>%
    ggplot(aes(x = pc, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, nrow = 1) +
    scale_y_continuous(limits = c(0, 1.6)) +
    labs(x = "PC Level", y = "Dispersion bias") +
    scale_fill_manual(values = c("Bias" = "coral")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  cowplot::plot_grid(accuracy_plot, bias_plot, ncol = 1, align = 'v', axis = 'l', rel_heights = c(1, 1))
}

# GM Model plotting function for RFI group
plot_gm_metrics_rfi <- function(df) {
  desired_levels <- c("GM", "74", "299", "600", "889")
  df_filtered <- df %>%
    filter(pc %in% desired_levels) %>%
    filter(Trait %in% c("RFI", "Mid intake")) %>%
    mutate(Trait = factor(Trait, levels = c("RFI", "Mid intake")),
           pc = factor(pc, levels = desired_levels))
  
  df_long <- df_filtered %>%
    pivot_longer(cols = c(Cor, Bias), names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("Cor", "Bias"),
                           labels = c("Genomic prediction accuracy", "Bias")))
  
  accuracy_plot <- df_long %>%
    filter(Metric == "Genomic prediction accuracy") %>%
    ggplot(aes(x = pc, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, nrow = 1) +
    scale_y_continuous(limits = c(0, 0.3)) +
    labs(x = "PC Level", y = "Genomic prediction accuracy") +
    scale_fill_manual(values = c("Genomic prediction accuracy" = "skyblue")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  bias_plot <- df_long %>%
    filter(Metric == "Bias") %>%
    ggplot(aes(x = pc, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, nrow = 1) +
    scale_y_continuous(limits = c(0, 3.2)) +
    labs(x = "PC Level", y = "Dispersion bias") +
    scale_fill_manual(values = c("Bias" = "coral")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  cowplot::plot_grid(accuracy_plot, bias_plot, ncol = 1, align = 'v', axis = 'l', rel_heights = c(1, 1))
}

# Special GM Model plotting function for RFI five-fold (different PC levels)
plot_gm_metrics_rfi_fivefold <- function(df) {
  desired_levels <- c("GM", "74", "299", "600")
  df_filtered <- df %>%
    filter(pc %in% desired_levels) %>%
    filter(Trait %in% c("RFI", "Mid intake")) %>%
    mutate(Trait = factor(Trait, levels = c("RFI", "Mid intake")),
           pc = factor(pc, levels = desired_levels))
  
  df_long <- df_filtered %>%
    pivot_longer(cols = c(Cor, Bias), names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("Cor", "Bias"),
                           labels = c("Genomic prediction accuracy", "Bias")))
  
  accuracy_plot <- df_long %>%
    filter(Metric == "Genomic prediction accuracy") %>%
    ggplot(aes(x = pc, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, nrow = 1) +
    scale_y_continuous(limits = c(0, 0.3)) +
    labs(x = "PC Level", y = "Genomic prediction accuracy") +
    scale_fill_manual(values = c("Genomic prediction accuracy" = "skyblue")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  bias_plot <- df_long %>%
    filter(Metric == "Bias") %>%
    ggplot(aes(x = pc, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, nrow = 1) +
    scale_y_continuous(limits = c(0, 1.6)) +
    labs(x = "PC Level", y = "Dispersion bias") +
    scale_fill_manual(values = c("Bias" = "coral")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  cowplot::plot_grid(accuracy_plot, bias_plot, ncol = 1, align = 'v', axis = 'l', rel_heights = c(1, 1))
}
# Function to create scatter plot for a single trait (Figure S5)
create_trait_plot <- function(gm_data, nngblup_data, trait_name, pc_value, display_name) {
  trait_data <- rbind(
    gm_data %>% filter(Trait == trait_name, PC == pc_value),
    nngblup_data %>% filter(Trait == trait_name, PC == pc_value)
  )
  
  correlations <- trait_data %>%
    group_by(Model) %>%
    summarise(correlation = round(cor(Observed, Predicted, use = "complete.obs"), 3), .groups = "drop") %>%
    mutate(label = paste("r =", correlation))
  
  p <- ggplot(trait_data, aes(x = Observed, y = Predicted)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    facet_wrap(~Model, scales = "fixed") +
    geom_text(data = correlations, aes(x = min(trait_data$Observed, na.rm = TRUE),
                                       y = max(trait_data$Predicted, na.rm = TRUE), label = label),
              hjust = 0, vjust = 1) +
    theme_bw() + labs(title = display_name) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 10, face = "bold"))
  return(p)
}

# Create complete scatter plot (Figure S5)
create_scatter_plot <- function(scatter_list) {
  p1 <- create_trait_plot(scatter_list$gm, scatter_list$nngblup, "AScCH4", 88, "Methane")
  p2 <- create_trait_plot(scatter_list$gm, scatter_list$nngblup, "AScCH4Ratio", 88, "Methane ratio")
  p3 <- create_trait_plot(scatter_list$gm, scatter_list$nngblup, "ARFI", 74, "RFI")
  grid.arrange(p1, p2, p3, ncol = 1)
}

# Function to create overlay plots for PCA transformation comparison (Figure S6)
create_overlay_plot <- function(data, group_name) {
  log_data <- subset(data, Transformation == "Log")
  clr_data <- subset(data, Transformation == "CLR")
  thresholds <- c(25, 50, 75, 95, 100)
  
  log_threshold_data <- sapply(thresholds, function(t) {
    which.min(abs(log_data$CumulativeVariance - t))
  })
  
  clr_threshold_data <- sapply(thresholds, function(t) {
    which.min(abs(clr_data$CumulativeVariance - t))
  })
  
  log_threshold_df <- data.frame(
    PC = log_data$PC[log_threshold_data],
    CumulativeVariance = log_data$CumulativeVariance[log_threshold_data],
    Label = paste0(thresholds, "%"), Transformation = "Log", PositionOffset = 2
  )
  
  clr_threshold_df <- data.frame(
    PC = clr_data$PC[clr_threshold_data],
    CumulativeVariance = clr_data$CumulativeVariance[clr_threshold_data],
    Label = paste0(thresholds, "%"), Transformation = "CLR", PositionOffset = -2
  )
  
  threshold_df <- rbind(log_threshold_df, clr_threshold_df)
  transform_colors <- c("Log" = "blue", "CLR" = "red")
  
  p <- ggplot(data, aes(x = PC, y = CumulativeVariance, color = Transformation, group = Transformation)) +
    geom_line(size = 1) +
    geom_point(data = threshold_df, aes(x = PC, y = CumulativeVariance + PositionOffset, color = Transformation), 
               size = 4, shape = 19) +
    geom_segment(data = threshold_df, aes(x = PC, y = CumulativeVariance, 
                                          xend = PC, yend = CumulativeVariance + PositionOffset, color = Transformation),
                 linetype = "dotted", size = 0.7) +
    geom_text(data = threshold_df, aes(x = PC, y = CumulativeVariance + PositionOffset + ifelse(Transformation == "Log", 4, -4), 
                                       label = Label, color = Transformation), hjust = 0.5, size = 3.3) +
    scale_color_manual(values = transform_colors) +
    scale_x_continuous(limits = c(0, max(data$PC) * 1.05)) +
    scale_y_continuous(breaks = seq(0, 100, by = 25), labels = function(x) paste0(x, "%"), limits = c(0, 105)) +
    labs(title = paste(group_name, "- Log vs CLR Transformation Comparison"),
         x = "Number of Principal Components", y = "Cumulative Variance Explained (%)", color = "Transformation") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.title = element_text(face = "bold", size = 11), axis.text = element_text(size = 9),
          panel.grid.minor = element_blank(), legend.position = "bottom")
  
  for (i in 1:length(thresholds)) {
    log_pc <- log_data$PC[which.min(abs(log_data$CumulativeVariance - thresholds[i]))]
    clr_pc <- clr_data$PC[which.min(abs(clr_data$CumulativeVariance - thresholds[i]))]
    p <- p + geom_vline(xintercept = log_pc, linetype = "dashed", color = "blue", alpha = 0.5) +
      geom_vline(xintercept = clr_pc, linetype = "dashed", color = "red", alpha = 0.5)
  }
  return(p)
}

# Create PCA transformation comparison plot (Figure S6)
create_pca_transformation_plot <- function(transformation_data) {
  methane_overlay <- create_overlay_plot(transformation_data$methane, "Methane Group")
  rfi_overlay <- create_overlay_plot(transformation_data$rfi, "RFI Group")
  grid.arrange(methane_overlay, rfi_overlay, ncol = 2)
}
# Function to create variance data structure for microbiability comparison
create_variance_data <- function(pc_mapping, log_data, clr_data, trait_name) {
  variance_levels <- c("25%", "50%", "75%", "95%", "100%", "Whole Data")
  variance_data <- data.frame(
    Variance = rep(variance_levels, each = 2),
    Transformation = rep(c("Log", "CLR"), times = 6),
    PC = pc_mapping, Microbability = NA, SE = NA, stringsAsFactors = FALSE
  )
  
  log_trait_data <- log_data %>% filter(Trait == trait_name)
  clr_trait_data <- clr_data %>% filter(Trait == trait_name)
  
  for (i in 1:nrow(variance_data)) {
    pc <- variance_data$PC[i]
    transform <- variance_data$Transformation[i]
    
    if (transform == "Log") {
      matching_row <- log_trait_data %>% filter(PC_Count == pc)
      if (nrow(matching_row) > 0) {
        variance_data$Microbability[i] <- matching_row$Microbability[1]
        variance_data$SE[i] <- matching_row$SE[1]
      }
    } else if (transform == "CLR") {
      matching_row <- clr_trait_data %>% filter(PC_Count == pc)
      if (nrow(matching_row) > 0) {
        variance_data$Microbability[i] <- matching_row$Microbability[1]
        variance_data$SE[i] <- matching_row$SE[1]
      }
    }
  }
  
  variance_data$Variance <- factor(variance_data$Variance, levels = variance_levels)
  variance_data$Microbability[is.na(variance_data$Microbability)] <- 0
  variance_data$SE[is.na(variance_data$SE)] <- 0
  return(variance_data)
}

# Function to create variance-based plot for microbiability
create_variance_plot <- function(data, title) {
  p <- ggplot(data, aes(x = Variance, y = Microbability, fill = Transformation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
    geom_errorbar(aes(ymin = Microbability - SE, ymax = Microbability + SE), 
                  width = 0.2, position = position_dodge(width = 0.9)) +
    labs(x = "Variance Explained", y = "Microbiability", title = title) +
    theme_minimal() + scale_fill_manual(values = c("Log" = "skyblue", "CLR" = "indianred")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          legend.position = "bottom") + ylim(0, 0.8)
  return(p)
}

# Function to prepare variance data with PC mappings for GM models
prepare_variance_data_gm <- function(data_log, data_clr, mapping) {
  result <- data.frame(Trait = character(), Variance = character(), Transformation = character(), 
                       Accuracy = numeric(), Bias = numeric(), stringsAsFactors = FALSE)
  
  log_traits <- unique(data_log$Trait)
  clr_traits <- unique(data_clr$Trait)
  all_traits <- union(log_traits, clr_traits)
  
  for (trait in all_traits) {
    log_trait_data <- data_log %>% filter(Trait == trait)
    clr_trait_data <- data_clr %>% filter(Trait == trait)
    
    for (i in 1:nrow(mapping)) {
      pc_log <- mapping$PC_Log[i]
      pc_clr <- mapping$PC_CLR[i]
      variance <- mapping$Variance[i]
      
      log_row <- log_trait_data %>% filter(PC_Count == pc_log)
      if (nrow(log_row) > 0) {
        result <- rbind(result, data.frame(Trait = trait, Variance = variance, Transformation = "Log",
                                           Accuracy = log_row$Accuracy[1], Bias = log_row$Bias[1], stringsAsFactors = FALSE))
      }
      
      clr_row <- clr_trait_data %>% filter(PC_Count == pc_clr)
      if (nrow(clr_row) > 0) {
        result <- rbind(result, data.frame(Trait = trait, Variance = variance, Transformation = "CLR",
                                           Accuracy = clr_row$Accuracy[1], Bias = clr_row$Bias[1], stringsAsFactors = FALSE))
      }
    }
  }
  
  result$Variance <- factor(result$Variance, levels = c("25%", "50%", "75%", "95%", "GM"))
  return(result)
}

# Function to create faceted plots for methane group GM comparison
create_facet_plots_methane_gm <- function(df) {
  df_long <- df %>%
    pivot_longer(cols = c(Accuracy, Bias), names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("Accuracy", "Bias"),
                           labels = c("Genomic prediction accuracy", "Dispersion bias")))
  
  accuracy_plot <- df_long %>%
    filter(Metric == "Genomic prediction accuracy") %>%
    ggplot(aes(x = Variance, y = Value, fill = Transformation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
    geom_text(aes(label = sprintf("%.2f", Value)), position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
    facet_wrap(~ Trait, nrow = 1) + scale_y_continuous(limits = c(0, 0.3)) +
    labs(x = "Variance Explained", y = "Genomic prediction accuracy") +
    scale_fill_manual(values = c("Log" = "skyblue", "CLR" = "indianred")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  
  bias_plot <- df_long %>%
    filter(Metric == "Dispersion bias") %>%
    ggplot(aes(x = Variance, y = Value, fill = Transformation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
    geom_text(aes(label = sprintf("%.2f", Value)), position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
    facet_wrap(~ Trait, nrow = 1) + scale_y_continuous(limits = c(0, 1.6)) +
    labs(x = "Variance Explained", y = "Dispersion bias") +
    scale_fill_manual(values = c("Log" = "skyblue", "CLR" = "indianred")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  
  cowplot::plot_grid(accuracy_plot, bias_plot, ncol = 1, align = 'v', axis = 'l', rel_heights = c(1, 1))
}

# Function to create faceted plots for RFI group GM comparison
create_facet_plots_rfi_gm <- function(df) {
  df_long <- df %>%
    pivot_longer(cols = c(Accuracy, Bias), names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("Accuracy", "Bias"),
                           labels = c("Genomic prediction accuracy", "Dispersion bias")))
  
  accuracy_plot <- df_long %>%
    filter(Metric == "Genomic prediction accuracy") %>%
    ggplot(aes(x = Variance, y = Value, fill = Transformation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
    geom_text(aes(label = sprintf("%.2f", Value)), position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
    facet_wrap(~ Trait, nrow = 1) + scale_y_continuous(limits = c(0, 0.3)) +
    labs(x = "Variance Explained", y = "Genomic prediction accuracy") +
    scale_fill_manual(values = c("Log" = "skyblue", "CLR" = "indianred")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  
  bias_plot <- df_long %>%
    filter(Metric == "Dispersion bias") %>%
    ggplot(aes(x = Variance, y = Value, fill = Transformation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
    geom_text(aes(label = sprintf("%.2f", Value)), position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
    facet_wrap(~ Trait, nrow = 1) + scale_y_continuous(limits = c(0, 3.2)) +
    labs(x = "Variance Explained", y = "Dispersion bias") +
    scale_fill_manual(values = c("Log" = "skyblue", "CLR" = "indianred")) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  
  cowplot::plot_grid(accuracy_plot, bias_plot, ncol = 1, align = 'v', axis = 'l', rel_heights = c(1, 1))
}
# Function to create combined accuracy plot for all traits (Figure S10)
create_all_traits_accuracy_plot <- function(df) {
  all_traits_order <- c("Methane", "Methane ratio", "LWT", "CO2", "RFI", "Mid intake")
  
  all_traits_accuracy_plot <- df %>%
    mutate(Trait = factor(Trait, levels = all_traits_order)) %>%
    ggplot(aes(x = Variance, y = Accuracy, fill = Transformation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
    geom_text(aes(label = round(Accuracy, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, ncol = 3) +
    labs(x = "Variance Explained", y = "NN-GBLUP Accuracy", 
         title = "NN-GBLUP Accuracy Comparison: Log vs CLR Transformation") +
    theme_minimal() + scale_fill_manual(values = c("Log" = "skyblue", "CLR" = "indianred")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "lightgray"), strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold")) + ylim(0, 0.5)
  
  return(all_traits_accuracy_plot)
}

# Function to create combined bias plot for all traits (Figure S11)
create_all_traits_bias_plot <- function(df) {
  all_traits_order <- c("Methane", "Methane ratio", "LWT", "CO2", "RFI", "Mid intake")
  
  all_traits_bias_plot <- df %>%
    mutate(Trait = factor(Trait, levels = all_traits_order)) %>%
    ggplot(aes(x = Variance, y = Bias, fill = Transformation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
    geom_text(aes(label = round(Bias, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    facet_wrap(~ Trait, ncol = 3, scales = "free_y") +
    labs(x = "Variance Explained", y = "NN-GBLUP Bias", 
         title = "NN-GBLUP Bias Comparison: Log vs CLR Transformation") +
    theme_minimal() + scale_fill_manual(values = c("Log" = "skyblue", "CLR" = "indianred")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "lightgray"), strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
  
  return(all_traits_bias_plot)
}

# Additional data loading functions for advanced figures

# Load microbiability transformation comparison data for Figure S7
load_microbiability_transformation_data <- function() {
  reactive({
    micro_path <- "data/Microbiability"
    methane_group_log <- read.table(file.path(micro_path, "methanegroup.txt"), header = TRUE)
    rfi_group_log <- read.table(file.path(micro_path, "RFIgroup.txt"), header = TRUE)
    methane_group_clr <- read.csv(file.path(micro_path, "microbiability_resultsmethane.csv"), header = TRUE)
    rfi_group_clr <- read.csv(file.path(micro_path, "microbiability_resultsRFI.csv"), header = TRUE)
    
    names(methane_group_clr)[names(methane_group_clr) == "Microbiability"] <- "Microbability"
    names(methane_group_clr)[names(methane_group_clr) == "Standard_Error"] <- "SE"
    names(rfi_group_clr)[names(rfi_group_clr) == "Microbiability"] <- "Microbability"
    names(rfi_group_clr)[names(rfi_group_clr) == "Standard_Error"] <- "SE"
    
    methane_group_log$Trait <- dplyr::recode(methane_group_log$Trait, "MethaneRatio" = "Methane ratio")
    methane_group_clr$Trait <- dplyr::recode(methane_group_clr$Trait, "AScCH4" = "Methane", "AScCH4Ratio" = "Methane ratio", "ALWT" = "LWT", "AScCO2" = "CO2")
    rfi_group_log$Trait <- dplyr::recode(rfi_group_log$Trait, "AMidIntake" = "Mid intake")
    rfi_group_clr$Trait <- dplyr::recode(rfi_group_clr$Trait, "ARFI" = "RFI", "AMidIntake" = "Mid intake")
    
    list(methane_log = methane_group_log, methane_clr = methane_group_clr, rfi_log = rfi_group_log, rfi_clr = rfi_group_clr)
  })
}

# Load GM model CLR comparison data for Figures S8 and S9
load_gm_clr_comparison_data <- function() {
  reactive({
    gm_compare_path <- "data/supplementarymaterial/Gmmodelonly"
    gm_log <- read.csv(file.path(gm_compare_path, "TrainTest.csv"), header = TRUE)
    gm_methane_clr <- read.csv(file.path(gm_compare_path, "TrainTest_Methane_GM.csv"), header = TRUE)
    gm_rfi_clr <- read.csv(file.path(gm_compare_path, "TrainTest_RFI_GM.csv"), header = TRUE)
    
    df_methane_log <- gm_log %>%
      filter(pc %in% c("GM", "88", "333", "640", "933")) %>%
      filter(Trait %in% c("Methane", "Methane ratio", "LWT", "CO2")) %>%
      mutate(Trait = factor(Trait, levels = c("Methane", "Methane ratio", "LWT", "CO2")),
             pc = factor(pc, levels = c("GM", "88", "333", "640", "933")),
             Transformation = "Log") %>%
      rename(PC_Count = pc, Accuracy = Cor)
    
    df_rfi_log <- gm_log %>%
      filter(pc %in% c("GM", "74", "299", "600", "889")) %>%
      filter(Trait %in% c("RFI", "Mid intake")) %>%
      mutate(Trait = factor(Trait, levels = c("RFI", "Mid intake")),
             pc = factor(pc, levels = c("GM", "74", "299", "600", "889")),
             Transformation = "Log") %>%
      rename(PC_Count = pc, Accuracy = Cor)
    
    gm_methane_clr <- gm_methane_clr %>%
      mutate(Trait = dplyr::recode(Trait, "AScCH4" = "Methane", "AScCH4Ratio" = "Methane ratio", "ALWT" = "LWT", "AScCO2" = "CO2"),
             Transformation = "CLR")
    
    gm_rfi_clr <- gm_rfi_clr %>%
      mutate(Trait = dplyr::recode(Trait, "ARFI" = "RFI", "AMidIntake" = "Mid intake"), Transformation = "CLR")
    
    list(methane_log = df_methane_log, rfi_log = df_rfi_log, methane_clr = gm_methane_clr, rfi_clr = gm_rfi_clr)
  })
}

# Load NN-GBLUP CLR comparison data for Figures S10 and S11
load_nngblup_clr_comparison_data <- function() {
  reactive({
    nngblup_compare_path <- "data/supplementarymaterial/nngblup"
    log_data_full <- read.table(file.path(nngblup_compare_path, "PCAusingalldatatraintest.txt"), header = TRUE)
    clr_data_full <- read.csv(file.path(nngblup_compare_path, "TrainTest_AllTraits_CLR_NNGBLUP.csv"), header = TRUE)
    
    log_methane_25pc <- "88"
    log_methane_50pc <- "333"
    log_rfi_25pc <- "74"
    log_rfi_50pc <- "299"
    
    clr_methane_25pc <- "102"
    clr_methane_50pc <- "346"
    clr_rfi_25pc <- "83"
    clr_rfi_50pc <- "310"
    
    log_data <- log_data_full %>%
      filter(pc %in% c(log_methane_25pc, log_methane_50pc, log_rfi_25pc, log_rfi_50pc))
    
    clr_data <- clr_data_full %>%
      filter(PC_Count %in% c(clr_methane_25pc, clr_methane_50pc, clr_rfi_25pc, clr_rfi_50pc))
    
    log_data <- log_data %>%
      mutate(Trait = dplyr::recode(Trait, "MethaneRatio" = "Methane ratio", "AMiDintake" = "Mid intake"))
    
    clr_data <- clr_data %>%
      mutate(Trait = dplyr::recode(Trait, "ALWT" = "LWT", "AMidIntake" = "Mid intake", "ARFI" = "RFI",
                                   "AScCH4" = "Methane", "AScCH4Ratio" = "Methane ratio", "AScCO2" = "CO2"))
    
    log_data <- log_data %>%
      mutate(Transformation = "Log",
             Variance = case_when(pc == log_methane_25pc | pc == log_rfi_25pc ~ "25%",
                                  pc == log_methane_50pc | pc == log_rfi_50pc ~ "50%",
                                  TRUE ~ NA_character_), Accuracy = Cor)
    
    clr_data <- clr_data %>%
      mutate(Transformation = "CLR",
             Variance = case_when(PC_Count == clr_methane_25pc | PC_Count == clr_rfi_25pc ~ "25%",
                                  PC_Count == clr_methane_50pc | PC_Count == clr_rfi_50pc ~ "50%",
                                  TRUE ~ NA_character_), pc = PC_Count)
    
    log_data_selected <- log_data %>% select(Trait, Variance, Transformation, pc, Accuracy, Bias)
    clr_data_selected <- clr_data %>% select(Trait, Variance, Transformation, pc = PC_Count, Accuracy, Bias)
    
    combined_data <- rbind(log_data_selected, clr_data_selected)
    combined_data <- combined_data %>%
      mutate(Variance = factor(Variance, levels = c("25%", "50%")),
             Transformation = factor(Transformation, levels = c("Log", "CLR")))
    
    return(combined_data)
  })
}
