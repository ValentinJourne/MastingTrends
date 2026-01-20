# Example Analysis Workflow
# This script demonstrates how to use the MastingTrends package
# to analyze Polish harvest dataset for mast seeding trends

# Clear workspace
rm(list = ls())

# Set working directory (adjust as needed)
# setwd("/path/to/MastingTrends")

# Load the main analysis script
source("scripts/main_analysis.R")

# ===================================================================
# EXAMPLE 1: Complete Analysis Workflow
# ===================================================================

cat("\n=== EXAMPLE 1: Running Complete Analysis ===\n\n")

# Path to your data file
# Replace this with the actual path to your Polish harvest dataset
data_path <- "data/polish_harvest_data.csv"

# Check if example data exists
if (!file.exists(data_path)) {
  cat("Note: Example data file not found at:", data_path, "\n")
  cat("Please place your Polish harvest dataset in the data/ directory\n")
  cat("Expected format: CSV with columns for year, species, and harvest\n\n")
  
  # Create example synthetic data for demonstration
  cat("Creating synthetic example data for demonstration...\n")
  set.seed(42)
  
  example_data <- data.frame(
    year = rep(2000:2020, 3),
    species = rep(c("Oak", "Beech", "Spruce"), each = 21),
    harvest = c(
      # Oak - strong masting with increasing trend
      rnorm(21, mean = seq(50, 80, length.out = 21), sd = 30),
      # Beech - moderate masting with no trend
      rnorm(21, mean = 60, sd = 20),
      # Spruce - weak masting with decreasing trend
      rnorm(21, mean = seq(70, 50, length.out = 21), sd = 10)
    )
  )
  
  # Ensure no negative values
  example_data$harvest <- pmax(example_data$harvest, 0)
  
  # Save example data
  write.csv(example_data, data_path, row.names = FALSE)
  cat("Example data created at:", data_path, "\n\n")
}

# Run the complete analysis
results <- analyze_masting_trends(
  data_path = data_path,
  output_dir = "results",
  year_col = "year",
  species_col = "species",
  harvest_col = "harvest",
  min_year = NULL,  # Optional: set minimum year
  max_year = NULL   # Optional: set maximum year
)

cat("\n=== Analysis results are stored in 'results/' directory ===\n")

# ===================================================================
# EXAMPLE 2: Step-by-step Analysis (for custom workflows)
# ===================================================================

cat("\n=== EXAMPLE 2: Step-by-step Custom Analysis ===\n\n")

# Load utilities
source("scripts/data_utils.R")
source("scripts/trend_analysis.R")
source("scripts/visualization.R")

# Load data
cat("Loading data...\n")
data <- load_harvest_data(data_path)

# Preprocess
cat("Preprocessing data...\n")
data <- preprocess_data(data, remove_na = TRUE)

# Get summary for a specific species
cat("\nAnalyzing Oak species:\n")
oak_data <- data[data$Species == "Oak", ]

# Calculate CV
oak_cv <- calculate_cv(oak_data$Harvest)
cat("  Coefficient of Variation:", round(oak_cv, 3), "\n")

# Detect mast years
oak_mast_years <- detect_mast_years(oak_data)
cat("  Mast years:", paste(oak_mast_years, collapse = ", "), "\n")

# Test for trend
oak_trend <- mann_kendall_trend(oak_data)
cat("  Trend:", oak_trend$trend, "\n")
cat("  P-value:", round(oak_trend$p_value, 4), "\n")

# Calculate autocorrelation
oak_autocor <- calculate_autocorrelation(oak_data, lag = 1)
cat("  Lag-1 autocorrelation:", round(oak_autocor, 3), "\n")

# ===================================================================
# EXAMPLE 3: Custom Visualizations
# ===================================================================

cat("\n=== EXAMPLE 3: Creating Custom Visualizations ===\n\n")

# Create figures directory
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Plot time series for all species
cat("Creating time series plot...\n")
plot_harvest_timeseries(data, save_path = "figures/example_timeseries.png")

# Plot CV trends
cat("Creating CV trends plot...\n")
plot_cv_trends(data, window_size = 5, save_path = "figures/example_cv_trends.png")

# Plot masting summary
cat("Creating masting summary...\n")
metrics <- calculate_masting_metrics(data)
plot_masting_summary(metrics, save_path = "figures/example_summary.png")

cat("\n=== Custom visualizations saved to 'figures/' directory ===\n")

# ===================================================================
# EXAMPLE 4: Comparing Multiple Species
# ===================================================================

cat("\n=== EXAMPLE 4: Comparing Species ===\n\n")

# Calculate metrics for all species
all_metrics <- calculate_masting_metrics(data)

# Sort by CV (masting intensity)
all_metrics <- all_metrics[order(all_metrics$CV, decreasing = TRUE), ]

cat("Species ranked by masting intensity (CV):\n")
print(all_metrics[, c("Species", "CV", "Trend", "P_Value")])

# Identify species with significant trends
significant_trends <- all_metrics[all_metrics$P_Value < 0.05, ]
if (nrow(significant_trends) > 0) {
  cat("\nSpecies with significant trends (p < 0.05):\n")
  print(significant_trends[, c("Species", "Trend", "Slope", "P_Value")])
} else {
  cat("\nNo species showed significant trends.\n")
}

cat("\n=== EXAMPLES COMPLETE ===\n")
cat("Check the 'results/' and 'figures/' directories for outputs\n\n")
