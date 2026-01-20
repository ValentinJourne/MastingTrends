# Trend Analysis Functions for Mast Seeding
# Statistical methods to detect and analyze masting trends

#' Calculate coefficient of variation (CV)
#'
#' CV is a common metric for measuring masting synchrony
#' Higher CV indicates more variable production (stronger masting)
#'
#' @param x Numeric vector of harvest values
#' @return Coefficient of variation
#' @export
calculate_cv <- function(x) {
  if (length(x) < 2 || all(is.na(x))) {
    return(NA)
  }
  mean_val <- mean(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  if (mean_val == 0) {
    return(NA)
  }
  return(sd_val / mean_val)
}

#' Perform Mann-Kendall trend test
#'
#' Non-parametric test for monotonic trends in time series
#'
#' @param data Data frame with Year and Harvest columns
#' @param species_name Optional species name for filtering
#' @return List with trend test results
#' @export
mann_kendall_trend <- function(data, species_name = NULL) {
  
  # Filter by species if specified
  if (!is.null(species_name)) {
    data <- data[data$Species == species_name, ]
  }
  
  # Sort by year
  data <- data[order(data$Year), ]
  
  # Check if Kendall package is available
  if (!requireNamespace("Kendall", quietly = TRUE)) {
    warning("Package 'Kendall' not available. Installing basic trend calculation.")
    # Simple linear trend as fallback
    if (nrow(data) < 3) {
      return(list(
        trend = "insufficient_data",
        tau = NA,
        p_value = NA,
        slope = NA
      ))
    }
    
    model <- lm(Harvest ~ Year, data = data)
    slope <- coef(model)[2]
    p_value <- summary(model)$coefficients[2, 4]
    
    trend_direction <- ifelse(p_value < 0.05,
                               ifelse(slope > 0, "increasing", "decreasing"),
                               "no_trend")
    
    return(list(
      trend = trend_direction,
      tau = NA,
      p_value = p_value,
      slope = slope
    ))
  }
  
  # Perform Mann-Kendall test
  mk_test <- Kendall::MannKendall(data$Harvest)
  
  # Calculate Sen's slope
  if (requireNamespace("trend", quietly = TRUE)) {
    sen_slope <- trend::sens.slope(data$Harvest)$estimates
  } else {
    # Simple slope estimate
    model <- lm(Harvest ~ Year, data = data)
    sen_slope <- coef(model)[2]
  }
  
  # Determine trend direction
  trend_direction <- ifelse(mk_test$sl < 0.05,
                             ifelse(mk_test$tau > 0, "increasing", "decreasing"),
                             "no_trend")
  
  return(list(
    trend = trend_direction,
    tau = mk_test$tau,
    p_value = mk_test$sl,
    slope = sen_slope
  ))
}

#' Calculate masting metrics for each species
#'
#' @param data Data frame with Year, Species, and Harvest columns
#' @return Data frame with masting metrics by species
#' @export
calculate_masting_metrics <- function(data) {
  
  # Split by species
  species_list <- split(data, data$Species)
  
  # Calculate metrics for each species
  metrics <- lapply(species_list, function(species_data) {
    
    # Basic statistics
    mean_harvest <- mean(species_data$Harvest, na.rm = TRUE)
    sd_harvest <- sd(species_data$Harvest, na.rm = TRUE)
    cv <- calculate_cv(species_data$Harvest)
    
    # Trend analysis
    trend_result <- mann_kendall_trend(species_data)
    
    # Create result
    data.frame(
      Species = unique(species_data$Species),
      N_Years = nrow(species_data),
      Mean_Harvest = mean_harvest,
      SD_Harvest = sd_harvest,
      CV = cv,
      Trend = trend_result$trend,
      Tau = trend_result$tau,
      P_Value = trend_result$p_value,
      Slope = trend_result$slope,
      stringsAsFactors = FALSE
    )
  })
  
  # Combine into single data frame
  metrics_df <- do.call(rbind, metrics)
  rownames(metrics_df) <- NULL
  
  return(metrics_df)
}

#' Detect mast years
#'
#' Identifies years with exceptionally high seed production
#'
#' @param data Data frame with Year and Harvest columns for a single species
#' @param threshold Threshold for defining mast year (default: mean + 1 SD)
#' @return Vector of mast years
#' @export
detect_mast_years <- function(data, threshold = NULL) {
  
  if (is.null(threshold)) {
    threshold <- mean(data$Harvest, na.rm = TRUE) + sd(data$Harvest, na.rm = TRUE)
  }
  
  mast_years <- data$Year[data$Harvest >= threshold]
  
  return(mast_years)
}

#' Calculate temporal autocorrelation
#'
#' Measures synchrony in masting between consecutive years
#'
#' @param data Data frame with Year and Harvest columns
#' @param lag Lag for autocorrelation (default: 1)
#' @return Autocorrelation coefficient
#' @export
calculate_autocorrelation <- function(data, lag = 1) {
  
  # Sort by year
  data <- data[order(data$Year), ]
  
  # Calculate autocorrelation
  if (nrow(data) < lag + 2) {
    return(NA)
  }
  
  acf_result <- acf(data$Harvest, lag.max = lag, plot = FALSE)
  
  return(acf_result$acf[lag + 1])
}
