# Visualization Functions for Mast Seeding Trends
# Create plots and figures for analyzing masting patterns

# Source required utility if not already loaded
if (!exists("calculate_cv")) {
  source("scripts/trend_analysis.R")
}

#' Plot harvest time series
#'
#' @param data Data frame with Year, Species, and Harvest columns
#' @param species_name Optional species name for filtering (plots all if NULL)
#' @param add_trend Add linear trend line (default: TRUE)
#' @param save_path Optional path to save the plot
#' @return ggplot object (if ggplot2 available) or base plot
#' @export
plot_harvest_timeseries <- function(data, species_name = NULL, add_trend = TRUE, save_path = NULL) {
  
  # Filter by species if specified
  if (!is.null(species_name)) {
    data <- data[data$Species == species_name, ]
    plot_title <- paste("Harvest Time Series -", species_name)
  } else {
    plot_title <- "Harvest Time Series - All Species"
  }
  
  # Check if ggplot2 is available
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    
    p <- ggplot2::ggplot(data, ggplot2::aes(x = Year, y = Harvest, color = Species)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = plot_title,
        x = "Year",
        y = "Harvest / Seed Production",
        color = "Species"
      )
    
    if (add_trend && !is.null(species_name)) {
      p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE, linetype = "dashed")
    }
    
    if (!is.null(save_path)) {
      ggplot2::ggsave(save_path, p, width = 10, height = 6)
    }
    
    return(p)
    
  } else {
    # Fallback to base R plotting
    
    # Set up plot
    if (!is.null(save_path)) {
      png(save_path, width = 800, height = 600)
    }
    
    # Get unique species
    species_unique <- unique(data$Species)
    colors <- rainbow(length(species_unique))
    
    # Create plot
    plot(data$Year, data$Harvest, 
         type = "n",
         main = plot_title,
         xlab = "Year",
         ylab = "Harvest / Seed Production")
    
    # Add lines for each species
    for (i in seq_along(species_unique)) {
      species_data <- data[data$Species == species_unique[i], ]
      species_data <- species_data[order(species_data$Year), ]
      lines(species_data$Year, species_data$Harvest, col = colors[i], lwd = 2)
      points(species_data$Year, species_data$Harvest, col = colors[i], pch = 16)
      
      # Add trend line if requested
      if (add_trend && !is.null(species_name)) {
        model <- lm(Harvest ~ Year, data = species_data)
        abline(model, col = colors[i], lty = 2)
      }
    }
    
    # Add legend
    legend("topright", legend = species_unique, col = colors, lwd = 2, pch = 16)
    
    if (!is.null(save_path)) {
      dev.off()
    }
  }
}

#' Plot coefficient of variation over time
#'
#' @param data Data frame with Year, Species, and Harvest columns
#' @param window_size Window size for rolling CV calculation (default: 5 years)
#' @param save_path Optional path to save the plot
#' @return Plot object
#' @export
plot_cv_trends <- function(data, window_size = 5, save_path = NULL) {
  
  # Calculate rolling CV for each species
  species_list <- split(data, data$Species)
  
  cv_data_list <- lapply(species_list, function(species_data) {
    species_data <- species_data[order(species_data$Year), ]
    
    n <- nrow(species_data)
    if (n < window_size) {
      return(NULL)
    }
    
    # Calculate rolling CV
    cv_values <- sapply(window_size:n, function(i) {
      window_data <- species_data[(i - window_size + 1):i, ]
      calculate_cv(window_data$Harvest)
    })
    
    data.frame(
      Species = unique(species_data$Species),
      Year = species_data$Year[window_size:n],
      CV = cv_values,
      stringsAsFactors = FALSE
    )
  })
  
  # Combine data
  cv_data <- do.call(rbind, cv_data_list[!sapply(cv_data_list, is.null)])
  
  if (nrow(cv_data) == 0) {
    warning("Insufficient data to calculate rolling CV")
    return(NULL)
  }
  
  # Plot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    
    p <- ggplot2::ggplot(cv_data, ggplot2::aes(x = Year, y = CV, color = Species)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Coefficient of Variation Over Time (", window_size, "-year window)", sep = ""),
        x = "Year",
        y = "Coefficient of Variation (CV)",
        color = "Species"
      )
    
    if (!is.null(save_path)) {
      ggplot2::ggsave(save_path, p, width = 10, height = 6)
    }
    
    return(p)
    
  } else {
    # Base R plotting
    if (!is.null(save_path)) {
      png(save_path, width = 800, height = 600)
    }
    
    species_unique <- unique(cv_data$Species)
    colors <- rainbow(length(species_unique))
    
    plot(cv_data$Year, cv_data$CV,
         type = "n",
         main = paste("Coefficient of Variation Over Time (", window_size, "-year window)", sep = ""),
         xlab = "Year",
         ylab = "Coefficient of Variation (CV)")
    
    for (i in seq_along(species_unique)) {
      species_cv <- cv_data[cv_data$Species == species_unique[i], ]
      lines(species_cv$Year, species_cv$CV, col = colors[i], lwd = 2)
      points(species_cv$Year, species_cv$CV, col = colors[i], pch = 16)
    }
    
    legend("topright", legend = species_unique, col = colors, lwd = 2, pch = 16)
    
    if (!is.null(save_path)) {
      dev.off()
    }
  }
}

#' Create summary plot of masting metrics
#'
#' @param metrics Data frame from calculate_masting_metrics()
#' @param save_path Optional path to save the plot
#' @return Plot object
#' @export
plot_masting_summary <- function(metrics, save_path = NULL) {
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    
    # Create bar plot of CV by species
    p <- ggplot2::ggplot(metrics, ggplot2::aes(x = Species, y = CV, fill = Trend)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(
        title = "Masting Intensity and Trends by Species",
        x = "Species",
        y = "Coefficient of Variation (CV)",
        fill = "Trend"
      ) +
      ggplot2::scale_fill_manual(
        values = c("increasing" = "red", "decreasing" = "blue", "no_trend" = "gray")
      )
    
    if (!is.null(save_path)) {
      ggplot2::ggsave(save_path, p, width = 10, height = 6)
    }
    
    return(p)
    
  } else {
    # Base R plotting
    if (!is.null(save_path)) {
      png(save_path, width = 800, height = 600)
    }
    
    # Create bar plot
    colors <- ifelse(metrics$Trend == "increasing", "red",
                     ifelse(metrics$Trend == "decreasing", "blue", "gray"))
    
    barplot(metrics$CV,
            names.arg = metrics$Species,
            col = colors,
            main = "Masting Intensity and Trends by Species",
            xlab = "Species",
            ylab = "Coefficient of Variation (CV)",
            las = 2)
    
    legend("topright",
           legend = c("Increasing", "Decreasing", "No Trend"),
           fill = c("red", "blue", "gray"))
    
    if (!is.null(save_path)) {
      dev.off()
    }
  }
}
