# Main Analysis Script for Polish Harvest Dataset
# Complete workflow for analyzing mast seeding trends

# Load required functions
source("scripts/data_utils.R")
source("scripts/trend_analysis.R")
source("scripts/visualization.R")

#' Run complete masting trend analysis
#'
#' This is the main function that orchestrates the entire analysis workflow
#'
#' @param data_path Path to the Polish harvest dataset file
#' @param output_dir Directory to save results and figures (default: "results")
#' @param year_col Name of the year column (default: "year")
#' @param species_col Name of the species column (default: "species")
#' @param harvest_col Name of the harvest column (default: "harvest")
#' @param min_year Minimum year to include in analysis (optional)
#' @param max_year Maximum year to include in analysis (optional)
#' @return List with analysis results
#' @export
analyze_masting_trends <- function(data_path,
                                   output_dir = "results",
                                   year_col = "year",
                                   species_col = "species",
                                   harvest_col = "harvest",
                                   min_year = NULL,
                                   max_year = NULL) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create figures directory if it doesn't exist
  fig_dir <- file.path(output_dir, "figures")
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE)
  }
  
  cat("=== Masting Trends Analysis ===\n\n")
  
  # Step 1: Load data
  cat("Step 1: Loading data from", data_path, "\n")
  data <- load_harvest_data(data_path, year_col, species_col, harvest_col)
  cat("  Loaded", nrow(data), "observations\n")
  cat("  Species found:", paste(unique(data$Species), collapse = ", "), "\n\n")
  
  # Step 2: Preprocess data
  cat("Step 2: Preprocessing data\n")
  data <- preprocess_data(data, remove_na = TRUE, min_year = min_year, max_year = max_year)
  cat("  After preprocessing:", nrow(data), "observations\n")
  cat("  Year range:", min(data$Year), "-", max(data$Year), "\n\n")
  
  # Step 3: Calculate summary statistics
  cat("Step 3: Calculating summary statistics\n")
  summary_stats <- summarize_harvest_data(data)
  print(summary_stats)
  cat("\n")
  
  # Save summary statistics
  write.csv(summary_stats, 
            file.path(output_dir, "summary_statistics.csv"), 
            row.names = FALSE)
  
  # Step 4: Calculate masting metrics
  cat("Step 4: Calculating masting metrics and trends\n")
  masting_metrics <- calculate_masting_metrics(data)
  print(masting_metrics)
  cat("\n")
  
  # Save masting metrics
  write.csv(masting_metrics, 
            file.path(output_dir, "masting_metrics.csv"), 
            row.names = FALSE)
  
  # Step 5: Detect mast years for each species
  cat("Step 5: Detecting mast years\n")
  species_list <- split(data, data$Species)
  mast_years_list <- lapply(names(species_list), function(sp) {
    mast_yrs <- detect_mast_years(species_list[[sp]])
    cat("  ", sp, ":", length(mast_yrs), "mast years -", 
        paste(mast_yrs, collapse = ", "), "\n")
    data.frame(
      Species = sp,
      Mast_Years = paste(mast_yrs, collapse = ", "),
      N_Mast_Years = length(mast_yrs),
      stringsAsFactors = FALSE
    )
  })
  mast_years_df <- do.call(rbind, mast_years_list)
  cat("\n")
  
  # Save mast years
  write.csv(mast_years_df, 
            file.path(output_dir, "mast_years.csv"), 
            row.names = FALSE)
  
  # Step 6: Create visualizations
  cat("Step 6: Creating visualizations\n")
  
  # Plot time series for all species
  cat("  - Creating time series plot (all species)\n")
  plot_harvest_timeseries(data, 
                          save_path = file.path(fig_dir, "timeseries_all_species.png"))
  
  # Plot time series for each species individually
  for (sp in unique(data$Species)) {
    cat("  - Creating time series plot for", sp, "\n")
    sp_filename <- gsub("[^A-Za-z0-9]", "_", sp)
    plot_harvest_timeseries(data, 
                            species_name = sp,
                            save_path = file.path(fig_dir, paste0("timeseries_", sp_filename, ".png")))
  }
  
  # Plot CV trends
  cat("  - Creating CV trends plot\n")
  plot_cv_trends(data, 
                 save_path = file.path(fig_dir, "cv_trends.png"))
  
  # Plot masting summary
  cat("  - Creating masting summary plot\n")
  plot_masting_summary(masting_metrics,
                       save_path = file.path(fig_dir, "masting_summary.png"))
  
  cat("\n")
  
  # Step 7: Generate summary report
  cat("Step 7: Generating summary report\n")
  report <- generate_summary_report(data, summary_stats, masting_metrics, mast_years_df)
  
  # Save report
  writeLines(report, file.path(output_dir, "analysis_report.txt"))
  cat("  Report saved to:", file.path(output_dir, "analysis_report.txt"), "\n\n")
  
  cat("=== Analysis Complete ===\n")
  cat("Results saved to:", output_dir, "\n")
  
  # Return results
  return(list(
    data = data,
    summary_statistics = summary_stats,
    masting_metrics = masting_metrics,
    mast_years = mast_years_df,
    output_directory = output_dir
  ))
}

#' Generate summary report
#'
#' @param data Original data
#' @param summary_stats Summary statistics
#' @param masting_metrics Masting metrics
#' @param mast_years Mast years data
#' @return Character vector with report text
#' @export
generate_summary_report <- function(data, summary_stats, masting_metrics, mast_years) {
  
  report <- c(
    "===============================================",
    "   MASTING TRENDS ANALYSIS - SUMMARY REPORT   ",
    "===============================================",
    "",
    paste("Analysis Date:", Sys.Date()),
    paste("Total Observations:", nrow(data)),
    paste("Number of Species:", length(unique(data$Species))),
    paste("Year Range:", min(data$Year), "-", max(data$Year)),
    "",
    "-----------------------------------------------",
    "SPECIES OVERVIEW",
    "-----------------------------------------------",
    ""
  )
  
  for (i in 1:nrow(summary_stats)) {
    sp <- summary_stats$species[i]
    report <- c(report,
                paste("Species:", sp),
                paste("  Years of data:", summary_stats$n_years[i]),
                paste("  Year range:", summary_stats$year_range[i]),
                paste("  Mean harvest:", round(summary_stats$mean_harvest[i], 2)),
                paste("  SD harvest:", round(summary_stats$sd_harvest[i], 2)),
                "")
  }
  
  report <- c(report,
              "-----------------------------------------------",
              "MASTING METRICS",
              "-----------------------------------------------",
              "")
  
  for (i in 1:nrow(masting_metrics)) {
    sp <- masting_metrics$Species[i]
    report <- c(report,
                paste("Species:", sp),
                paste("  Coefficient of Variation (CV):", round(masting_metrics$CV[i], 3)),
                paste("  Trend:", masting_metrics$Trend[i]),
                paste("  P-value:", round(masting_metrics$P_Value[i], 4)),
                paste("  Slope:", round(masting_metrics$Slope[i], 4)),
                "")
  }
  
  report <- c(report,
              "-----------------------------------------------",
              "MAST YEARS DETECTED",
              "-----------------------------------------------",
              "")
  
  for (i in 1:nrow(mast_years)) {
    report <- c(report,
                paste("Species:", mast_years$Species[i]),
                paste("  Number of mast years:", mast_years$N_Mast_Years[i]),
                paste("  Mast years:", mast_years$Mast_Years[i]),
                "")
  }
  
  report <- c(report,
              "-----------------------------------------------",
              "INTERPRETATION NOTES",
              "-----------------------------------------------",
              "",
              "Coefficient of Variation (CV):",
              "  - High CV (>1.0) indicates strong masting behavior",
              "  - Moderate CV (0.5-1.0) indicates moderate variability",
              "  - Low CV (<0.5) indicates consistent production",
              "",
              "Trend Analysis:",
              "  - 'increasing': Statistically significant increase over time",
              "  - 'decreasing': Statistically significant decrease over time",
              "  - 'no_trend': No significant temporal trend detected",
              "",
              "Mast Years:",
              "  - Defined as years with harvest > mean + 1 SD",
              "  - Represents years of exceptionally high seed production",
              "",
              "===============================================")
  
  return(report)
}

# Example usage (commented out - uncomment and modify for your data):
# results <- analyze_masting_trends(
#   data_path = "data/polish_harvest_data.csv",
#   output_dir = "results",
#   year_col = "year",
#   species_col = "species",
#   harvest_col = "harvest"
# )
