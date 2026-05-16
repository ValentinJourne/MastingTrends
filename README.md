# MastingTrends

A comprehensive R toolkit for analyzing mast seeding trends using the Polish harvest dataset and other seed production data.
THIS IS ONGOING WORK, and Valentin Journe need to clean the code (sorry)

## Overview

Mast seeding (or masting) is the synchronous, highly variable production of seeds by a plant population. This repository provides tools to:

- Load and preprocess harvest/seed production datasets
- Detect temporal trends in masting behavior
- Calculate rolling CV
- Relate masting metrics / seed production to weather

## Features

- **main_analysisvFinal**: main code for analysis 
- **Trend Analysis**: simple glmm regression (using glmmTMB)
- **Masting Metrics**: Coefficient of variation (CV)
- **Visualizations**: Time series plots, CV trends, and summary figures

## Installation

Clone this repository:

```bash
git clone https://github.com/ValentinJourne/MastingTrends.git
cd MastingTrends
```

### Required R Packages

The toolkit uses base R functions but can leverage additional packages for enhanced functionality:

## Contributing

Contributions are welcome! Please feel free to contact us.

## License

See the [LICENSE](LICENSE) file for details.

## References

## Contact

For questions or issues, please contact Valentin Journe (journe.valentin@gmail.com)
