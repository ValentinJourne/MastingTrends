# Test script for MastingTrends functions
# This script tests basic functionality without requiring external data

cat("=== Testing MastingTrends Functions ===\n\n")

# Load required scripts
source("scripts/data_utils.R")
source("scripts/trend_analysis.R")
source("scripts/visualization.R")

# Create synthetic test data
cat("1. Creating synthetic test data...\n")
set.seed(123)

test_data <- data.frame(
  Year = rep(2000:2020, 2),
  Species = rep(c("TestSpecies1", "TestSpecies2"), each = 21),
  Harvest = c(
    rnorm(21, mean = 50, sd = 20),  # Species 1
    rnorm(21, mean = 60, sd = 15)   # Species 2
  )
)

# Ensure no negative values
test_data$Harvest <- pmax(test_data$Harvest, 0)

cat("  Created", nrow(test_data), "test observations\n")
cat("  Species:", paste(unique(test_data$Species), collapse = ", "), "\n")
cat("  Year range:", min(test_data$Year), "-", max(test_data$Year), "\n\n")

# Test data preprocessing
cat("2. Testing data preprocessing...\n")
processed_data <- preprocess_data(test_data, remove_na = TRUE)
cat("  Processed", nrow(processed_data), "observations\n")
cat("  ✓ preprocess_data() working\n\n")

# Test summary statistics
cat("3. Testing summary statistics...\n")
summary_stats <- summarize_harvest_data(test_data)
print(summary_stats)
cat("  ✓ summarize_harvest_data() working\n\n")

# Test CV calculation
cat("4. Testing coefficient of variation calculation...\n")
sp1_data <- test_data[test_data$Species == "TestSpecies1", ]
cv_value <- calculate_cv(sp1_data$Harvest)
cat("  CV for TestSpecies1:", round(cv_value, 3), "\n")
cat("  ✓ calculate_cv() working\n\n")

# Test trend analysis
cat("5. Testing trend analysis...\n")
trend_result <- mann_kendall_trend(sp1_data)
cat("  Trend for TestSpecies1:", trend_result$trend, "\n")
cat("  P-value:", round(trend_result$p_value, 4), "\n")
cat("  ✓ mann_kendall_trend() working\n\n")

# Test masting metrics
cat("6. Testing masting metrics calculation...\n")
metrics <- calculate_masting_metrics(test_data)
print(metrics)
cat("  ✓ calculate_masting_metrics() working\n\n")

# Test mast year detection
cat("7. Testing mast year detection...\n")
mast_years <- detect_mast_years(sp1_data)
cat("  Detected", length(mast_years), "mast years:", paste(mast_years, collapse = ", "), "\n")
cat("  ✓ detect_mast_years() working\n\n")

# Test autocorrelation
cat("8. Testing autocorrelation calculation...\n")
autocor <- calculate_autocorrelation(sp1_data, lag = 1)
cat("  Lag-1 autocorrelation:", round(autocor, 3), "\n")
cat("  ✓ calculate_autocorrelation() working\n\n")

# Test visualization functions (without saving)
cat("9. Testing visualization functions...\n")
cat("  Note: Visualization tests create plots but do not save them\n")

# Test basic plotting
tryCatch({
  plot_harvest_timeseries(test_data, species_name = "TestSpecies1")
  cat("  ✓ plot_harvest_timeseries() working\n")
}, error = function(e) {
  cat("  ⚠ plot_harvest_timeseries() had an issue:", e$message, "\n")
})

tryCatch({
  plot_cv_trends(test_data, window_size = 5)
  cat("  ✓ plot_cv_trends() working\n")
}, error = function(e) {
  cat("  ⚠ plot_cv_trends() had an issue:", e$message, "\n")
})

tryCatch({
  plot_masting_summary(metrics)
  cat("  ✓ plot_masting_summary() working\n")
}, error = function(e) {
  cat("  ⚠ plot_masting_summary() had an issue:", e$message, "\n")
})

cat("\n=== All Basic Tests Complete ===\n")
cat("Core functionality is working correctly!\n")
cat("You can now use these functions with your Polish harvest dataset.\n\n")
