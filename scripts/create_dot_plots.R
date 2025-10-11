library(ggplot2)
library(gridExtra)

# Read the TPM data
tpm_data <- read.table("allSamples_tpms.tsv", header = TRUE, row.names = 1)

# Define colors for each treatment
colors <- c("StHOIH" = "#27B4F5", "StMono" = "#09BBC8", "StNSC" = "#056970")

# Define animals and treatments
animals <- c("F1", "F2", "F3")
treatments <- c("StHOIH", "StMono", "StNSC")

# Function to create correlation plot
create_corr_plot <- function(stomach_col, treatment_col, animal, treatment, color) {
  # Extract data
  x <- tpm_data[[stomach_col]]
  y <- tpm_data[[treatment_col]]
  
  # Calculate correlation
  cor_pearson <- cor(x, y, method = "pearson")
  
  # Create data frame
  df <- data.frame(Stomach = x, Treatment = y)
  
  # Get axis limits for annotation placement (in log10 space)
  x_range <- range(log10(x[x > 0]), na.rm = TRUE)
  y_range <- range(log10(y[y > 0]), na.rm = TRUE)
  
  # Create plot with log10 scales
  p <- ggplot(df, aes(x = Stomach, y = Treatment)) +
    geom_point(alpha = 0.3, size = 1, color = color) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_x_log10() +
    scale_y_log10() +
    # Add correlation in top-left corner using actual coordinates
    annotate("text", 
             x = 10^(x_range[1] + 0.05 * diff(x_range)), 
             y = 10^(y_range[2] - 0.05 * diff(y_range)), 
             label = sprintf("r = %.3f", cor_pearson),
             hjust = 0, vjust = 1, size = 10, fontface = "bold",
             color = "black") +
    labs(
      x = paste(animal, "Stomach (TPM)"),
      y = paste(animal, treatment, "(TPM)")
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18)
    )
  
  return(p)
}

# Generate all plots
all_plots <- list()
plot_idx <- 1

for (animal in animals) {
  stomach_col <- paste0(animal, "Stomach")
  
  for (treatment in treatments) {
    treatment_col <- paste0(animal, treatment)
    color <- colors[treatment]
    
    p <- create_corr_plot(stomach_col, treatment_col, animal, treatment, color)
    all_plots[[plot_idx]] <- p
    plot_idx <- plot_idx + 1
  }
}

# Arrange plots in a 3x3 grid (3 animals x 3 treatments)
pdf("stomach_treatment_correlations.pdf", width = 15, height = 15)
grid.arrange(grobs = all_plots, ncol = 3)
dev.off()

# Alternative: Save individual plots per animal
for (i in 1:length(animals)) {
  animal <- animals[i]
  
  # Get the 3 plots for this animal
  animal_plots <- all_plots[((i-1)*3 + 1):(i*3)]
  
  pdf(paste0(animal, "_stomach_treatment_correlations.pdf"), width = 15, height = 5)
  grid.arrange(grobs = animal_plots, ncol = 3)
  dev.off()
}

# Print correlation summary
cat("\n=== Correlation Summary ===\n")
for (animal in animals) {
  cat("\n", animal, ":\n", sep = "")
  stomach_col <- paste0(animal, "Stomach")
  
  for (treatment in treatments) {
    treatment_col <- paste0(animal, treatment)
    cor_val <- cor(tpm_data[[stomach_col]], tpm_data[[treatment_col]], method = "pearson")
    cat("  Stomach vs", treatment, ": r =", sprintf("%.4f", cor_val), "\n")
  }
}

cat("\nPlots saved:\n")
cat("  - stomach_treatment_correlations.pdf (all plots in 3x3 grid)\n")
cat("  - F1_stomach_treatment_correlations.pdf\n")
cat("  - F2_stomach_treatment_correlations.pdf\n")
cat("  - F3_stomach_treatment_correlations.pdf\n")
