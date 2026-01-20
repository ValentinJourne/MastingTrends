# Quick Start Guide for MastingTrends

This guide will help you get started with analyzing mast seeding trends in your Polish harvest dataset.

## Step 1: Prepare Your Data

Your dataset should be in CSV format with at least these three columns:
- Year of observation
- Species name
- Harvest/seed production value

Example data structure:
```
year,species,harvest
2000,Oak,45.2
2001,Oak,52.1
2002,Oak,38.9
2000,Beech,55.3
2001,Beech,61.2
2002,Beech,42.8
```

Place your data file in the `data/` directory, for example:
- `data/polish_harvest_data.csv`

## Step 2: Run the Analysis

### Option A: Use the Example Script (Recommended for beginners)

1. Open R or RStudio
2. Set your working directory to the MastingTrends folder:
   ```r
   setwd("/path/to/MastingTrends")
   ```

3. Run the example analysis script:
   ```r
   source("scripts/example_analysis.R")
   ```

This will:
- Load your data (or create example data if your file doesn't exist)
- Run the complete analysis
- Generate all plots and statistics
- Save results to `results/` and `figures/` directories

### Option B: Use the Main Analysis Function (For custom workflows)

```r
# Load the main analysis script
source("scripts/main_analysis.R")

# Run analysis on your data
results <- analyze_masting_trends(
  data_path = "data/polish_harvest_data.csv",
  output_dir = "results",
  year_col = "year",      # Name of year column in your data
  species_col = "species", # Name of species column in your data
  harvest_col = "harvest"  # Name of harvest column in your data
)

# View results
print(results$masting_metrics)
```

## Step 3: Interpret Your Results

After running the analysis, check these output files in the `results/` directory:

### 1. `summary_statistics.csv`
Basic descriptive statistics for each species:
- Number of years with data
- Mean and standard deviation of harvest
- Min and max values

### 2. `masting_metrics.csv`
Key masting metrics:
- **CV (Coefficient of Variation)**: Measure of masting intensity
  - High CV (>1.0) = Strong masting behavior
  - Moderate CV (0.5-1.0) = Moderate variability
  - Low CV (<0.5) = Consistent production
- **Trend**: Direction of temporal trend (increasing/decreasing/no_trend)
- **P_Value**: Statistical significance of trend
- **Slope**: Rate of change over time

### 3. `mast_years.csv`
Years identified as mast years (exceptionally high production) for each species

### 4. `analysis_report.txt`
Human-readable summary of all results with interpretation notes

### 5. `figures/` directory
Visual representations:
- Time series plots for each species
- CV trends over time
- Summary comparison of all species

## Step 4: Common Customizations

### Filter by year range
```r
results <- analyze_masting_trends(
  data_path = "data/polish_harvest_data.csv",
  min_year = 2000,
  max_year = 2020
)
```

### Analyze specific species
```r
source("scripts/data_utils.R")
source("scripts/trend_analysis.R")

data <- load_harvest_data("data/polish_harvest_data.csv")
oak_data <- data[data$Species == "Oak", ]

# Calculate metrics
oak_cv <- calculate_cv(oak_data$Harvest)
oak_trend <- mann_kendall_trend(oak_data)
oak_mast_years <- detect_mast_years(oak_data)

print(paste("Oak CV:", round(oak_cv, 3)))
print(paste("Oak trend:", oak_trend$trend))
print(paste("Mast years:", paste(oak_mast_years, collapse=", ")))
```

### Create custom visualizations
```r
source("scripts/visualization.R")

# Plot specific species
plot_harvest_timeseries(data, species_name = "Oak", add_trend = TRUE)

# Adjust CV window size
plot_cv_trends(data, window_size = 10)
```

## Troubleshooting

### If you get "file not found" errors:
- Make sure you're in the correct working directory: `getwd()`
- Check that your data file exists: `file.exists("data/polish_harvest_data.csv")`

### If you get "column not found" errors:
- Check your column names: `names(read.csv("data/polish_harvest_data.csv"))`
- Specify the correct column names in the function call

### If plots don't display:
- Install ggplot2 for better plots: `install.packages("ggplot2")`
- Or use base R plots (automatic fallback)

### For enhanced functionality:
Install optional packages:
```r
install.packages(c("ggplot2", "readxl", "Kendall", "trend"))
```

## Getting Help

- Check the main README.md for detailed documentation
- Look at example_analysis.R for usage examples
- Review test_functions.R to see how functions work

## Next Steps

After running the basic analysis:

1. Examine the time series plots to identify masting patterns
2. Compare CV values across species
3. Look for species with significant trends
4. Investigate detected mast years
5. Consider ecological factors that might explain patterns

## Citation

When using this toolkit in your research, please cite:
- The original data source (Polish harvest dataset)
- This repository: https://github.com/ValentinJourne/MastingTrends
