# MastingTrends

A comprehensive R toolkit for analyzing mast seeding trends using the Polish harvest dataset and other seed production data.

## Overview

Mast seeding (or masting) is the synchronous, highly variable production of seeds by a plant population. This repository provides tools to:

- Load and preprocess harvest/seed production datasets
- Detect temporal trends in masting behavior
- Calculate rolling CV
- Relate masting metrics / seed production to weather

## Features

- **Git_version_main**: main code for analysis 
- **Trend Analysis**: simple glmm regression (using glmmTMB)
- **Data utils**: Calculate coefficient of variation (CV), etc. 
- **Visualizations**: Time series plots, CV trends, and summary figures.

## Installation

Clone this repository:

```bash
git clone https://github.com/ValentinJourne/MastingTrends.git
cd MastingTrends
```
Then how to proceed? once you have this new folder? 
First download data from here https://zenodo.org/records/20658205 and put it in the folder "MastingTrends/data" that contain harvest seed data and temperature (files were too big).
Then go into scripts/Git_version_main.R and you can run the code. 

### Required R Packages

The toolkit uses base R functions but can leverage additional packages for enhanced functionality (ggplot2, glmmTMB, etc)

## References
This analysis is related to the submitted article; Oak masting remains stable despite climate warming available at https://ecoevorxiv.org/repository/view/13095/ 

## Contact

For questions or issues, please contact Valentin Journe (journe.valentin@gmail.com) and Michal Bogdziewicz (michalbogdziewicz@gmail.com)
