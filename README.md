# MastingTrends

A comprehensive R toolkit for analyzing mast seeding trends using the Polish harvest dataset and other seed production data.

## Overview

Mast seeding (or masting) is the synchronous, highly variable production of seeds by a plant population. This repository provides tools to:

- Load and preprocess harvest/seed production datasets
- Detect temporal trends in masting behavior
- Calculate masting metrics (CV, autocorrelation, etc.)
- Identify mast years
- Visualize trends and patterns
- Generate comprehensive analysis reports

## Features

- **Data Loading**: Support for CSV, Excel, and RData formats
- **Trend Analysis**: Mann-Kendall tests, linear trends, and Sen's slope estimation
- **Masting Metrics**: Coefficient of variation (CV), temporal autocorrelation, mast year detection
- **Visualizations**: Time series plots, CV trends, and summary figures
- **Automated Workflow**: Complete end-to-end analysis with a single function call

## Installation

Clone this repository:

```bash
git clone https://github.com/ValentinJourne/MastingTrends.git
cd MastingTrends
```

### Required R Packages

The toolkit uses base R functions but can leverage additional packages for enhanced functionality:

**Optional packages** (install as needed):
```r
install.packages("ggplot2")     # For enhanced visualizations
install.packages("readxl")      # For Excel file support
install.packages("Kendall")     # For Mann-Kendall trend tests
install.packages("trend")       # For Sen's slope estimation
```

The scripts will work without these packages using fallback methods, but installing them is recommended for full functionality.

## Quick Start

### 1. Prepare Your Data

Place your Polish harvest dataset in the `data/` directory. The dataset should be in CSV, Excel, or RData format with the following columns:

- `year`: Year of observation
- `species`: Tree species name
- `harvest`: Seed production or harvest value

Example CSV format:
```csv
year,species,harvest
2000,Oak,45.2
2001,Oak,52.1
2000,Beech,38.5
2001,Beech,41.2
```

### 2. Run the Analysis

Open R and run the example analysis:

```r
# Set working directory to the repository
setwd("/path/to/MastingTrends")

# Run the example analysis
source("scripts/example_analysis.R")
```

Or use the main analysis function directly:

```r
source("scripts/main_analysis.R")

results <- analyze_masting_trends(
  data_path = "data/polish_harvest_data.csv",
  output_dir = "results",
  year_col = "year",
  species_col = "species",
  harvest_col = "harvest"
)
```

### 3. View Results

Results are saved to the `results/` directory:
- `summary_statistics.csv`: Descriptive statistics by species
- `masting_metrics.csv`: Masting metrics and trend analysis
- `mast_years.csv`: Detected mast years for each species
- `analysis_report.txt`: Human-readable summary report
- `figures/`: Visualizations (time series, CV trends, summaries)

## Project Structure

```
MastingTrends/
├── data/                       # Data directory
│   ├── README.md              # Data format documentation
│   └── .gitkeep
├── scripts/                    # Analysis scripts
│   ├── data_utils.R           # Data loading and preprocessing
│   ├── trend_analysis.R       # Statistical analysis functions
│   ├── visualization.R        # Plotting functions
│   ├── main_analysis.R        # Main analysis workflow
│   └── example_analysis.R     # Example usage and tutorials
├── results/                    # Analysis outputs (created automatically)
│   └── .gitkeep
├── figures/                    # Generated figures (created automatically)
│   └── .gitkeep
├── .gitignore                 # Git ignore rules
├── LICENSE                     # License file
└── README.md                   # This file
```

## Usage Examples

### Example 1: Complete Analysis

```r
source("scripts/main_analysis.R")

# Run complete analysis workflow
results <- analyze_masting_trends(
  data_path = "data/polish_harvest_data.csv",
  output_dir = "results"
)

# Access results
print(results$masting_metrics)
print(results$mast_years)
```

### Example 2: Custom Analysis

```r
source("scripts/data_utils.R")
source("scripts/trend_analysis.R")

# Load and preprocess data
data <- load_harvest_data("data/polish_harvest_data.csv")
data <- preprocess_data(data, min_year = 2000, max_year = 2020)

# Calculate metrics for specific species
oak_data <- data[data$Species == "Oak", ]
oak_cv <- calculate_cv(oak_data$Harvest)
oak_trend <- mann_kendall_trend(oak_data)
oak_mast_years <- detect_mast_years(oak_data)

print(paste("Oak CV:", round(oak_cv, 3)))
print(paste("Oak trend:", oak_trend$trend))
print(paste("Oak mast years:", paste(oak_mast_years, collapse=", ")))
```

### Example 3: Visualizations

```r
source("scripts/visualization.R")

# Plot time series
plot_harvest_timeseries(data, species_name = "Oak", add_trend = TRUE)

# Plot CV trends
plot_cv_trends(data, window_size = 5)

# Plot masting summary
metrics <- calculate_masting_metrics(data)
plot_masting_summary(metrics)
```

## Analysis Methods

### Masting Metrics

- **Coefficient of Variation (CV)**: Measures inter-annual variability in seed production. Higher CV indicates stronger masting behavior.
  - High CV (>1.0): Strong masting
  - Moderate CV (0.5-1.0): Moderate variability
  - Low CV (<0.5): Consistent production

### Trend Detection

- **Mann-Kendall Test**: Non-parametric test for monotonic trends
- **Sen's Slope**: Robust estimate of trend magnitude
- **Linear Regression**: Alternative trend estimation

### Mast Year Detection

Mast years are identified as years where seed production exceeds mean + 1 standard deviation.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

See the [LICENSE](LICENSE) file for details.

## References

For more information on mast seeding:
- Kelly, D., & Sork, V. L. (2002). Mast seeding in perennial plants: Why, how, where? Annual Review of Ecology and Systematics, 33(1), 427-447.
- Koenig, W. D., & Knops, J. M. (2005). The mystery of masting in trees. American Scientist, 93(4), 340-347.

## Contact

For questions or issues, please open an issue on GitHub.
