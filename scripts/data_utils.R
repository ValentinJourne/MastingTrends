# Data Loading Utilities for Masting Trends Analysis
# Functions to load and preprocess Polish harvest dataset

#' Load Polish harvest data from CSV file
#'
#' @param filepath Path to the CSV file containing harvest data
#' @param year_col Name of the year column (default: "year")
#' @param species_col Name of the species column (default: "species")
#' @param harvest_col Name of the harvest/seed production column (default: "harvest")
#' @return Data frame with standardized column names
#' @export
load_harvest_data <- function(filepath, 
                              year_col = "year", 
                              species_col = "species", 
                              harvest_col = "harvest") {
  
  # Check if file exists
  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }
  
  # Determine file type and load accordingly
  file_ext <- tolower(tools::file_ext(filepath))
  
  data <- switch(file_ext,
                 "csv" = read.csv(filepath, stringsAsFactors = FALSE),
                 "xlsx" = {
                   if (!requireNamespace("readxl", quietly = TRUE)) {
                     stop("Package 'readxl' is required to read Excel files")
                   }
                   readxl::read_excel(filepath)
                 },
                 "xls" = {
                   if (!requireNamespace("readxl", quietly = TRUE)) {
                     stop("Package 'readxl' is required to read Excel files")
                   }
                   readxl::read_excel(filepath)
                 },
                 "rdata" = ,
                 "rda" = {
                   load(filepath)
                   get(ls()[1])  # Return first object
                 },
                 stop(paste("Unsupported file format:", file_ext))
  )
  
  # Convert to data frame if tibble
  data <- as.data.frame(data)
  
  # Check for required columns
  if (!year_col %in% names(data)) {
    stop(paste("Year column not found:", year_col))
  }
  if (!species_col %in% names(data)) {
    stop(paste("Species column not found:", species_col))
  }
  if (!harvest_col %in% names(data)) {
    stop(paste("Harvest column not found:", harvest_col))
  }
  
  # Standardize column names
  data$Year <- data[[year_col]]
  data$Species <- data[[species_col]]
  data$Harvest <- data[[harvest_col]]
  
  # Convert Year to numeric if not already
  data$Year <- as.numeric(data$Year)
  
  # Convert Harvest to numeric if not already
  data$Harvest <- as.numeric(data$Harvest)
  
  return(data)
}

#' Preprocess harvest data
#'
#' @param data Data frame with harvest data
#' @param remove_na Remove rows with NA values (default: TRUE)
#' @param min_year Minimum year to include (optional)
#' @param max_year Maximum year to include (optional)
#' @return Preprocessed data frame
#' @export
preprocess_data <- function(data, remove_na = TRUE, min_year = NULL, max_year = NULL) {
  
  # Remove NA values if requested
  if (remove_na) {
    data <- data[complete.cases(data[, c("Year", "Species", "Harvest")]), ]
  }
  
  # Filter by year range if specified
  if (!is.null(min_year)) {
    data <- data[data$Year >= min_year, ]
  }
  if (!is.null(max_year)) {
    data <- data[data$Year <= max_year, ]
  }
  
  # Sort by species and year
  data <- data[order(data$Species, data$Year), ]
  
  return(data)
}

#' Get summary statistics for harvest data
#'
#' @param data Data frame with harvest data
#' @return Summary statistics by species
#' @export
summarize_harvest_data <- function(data) {
  
  # Split by species
  species_list <- split(data, data$Species)
  
  # Calculate summary statistics for each species
  summaries <- lapply(species_list, function(species_data) {
    list(
      species = unique(species_data$Species),
      n_years = length(unique(species_data$Year)),
      year_range = paste(min(species_data$Year), "-", max(species_data$Year)),
      mean_harvest = mean(species_data$Harvest, na.rm = TRUE),
      sd_harvest = sd(species_data$Harvest, na.rm = TRUE),
      min_harvest = min(species_data$Harvest, na.rm = TRUE),
      max_harvest = max(species_data$Harvest, na.rm = TRUE)
    )
  })
  
  # Convert to data frame
  summary_df <- do.call(rbind, lapply(summaries, as.data.frame))
  rownames(summary_df) <- NULL
  
  return(summary_df)
}
