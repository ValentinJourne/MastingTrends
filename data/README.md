# Data Directory

## Polish Harvest Dataset

Place your Polish harvest dataset files in this directory.

### Expected Format

The dataset should contain the following columns:
- `year`: Year of observation
- `species`: Tree species name
- `harvest`: Harvest/seed production value (or similar measure)
- `location`: Geographic location (optional)

### Supported Formats
- CSV files (.csv)
- Excel files (.xlsx, .xls)
- RData files (.RData, .rda)

### Example File Structure
```
data/
├── polish_harvest_data.csv
├── metadata.txt
└── README.md
```

### Data Sources

Add information about your data sources here, including:
- Data collection methods
- Time period covered
- Geographic regions
- Species included
