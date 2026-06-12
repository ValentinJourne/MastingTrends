#' Reformat and preprocess Polish seed harvest data
#' NOTE that here I finnaly managed dfifferently and did not used it
#'
#' Loads and reformats Polish seed harvest data for masting trend analyses. The
#' function creates a site-by-species identifier, joins site metadata, fills
#' missing demand values, assigns broad functional groups and geographic regions,
#' and calculates previous-year seed production.
#'
#' @param initial.data.harvest.load A data frame containing the raw harvest data.
#'   It must include at least `Species`, `NADL`, `Year`, `Demand`, and `seeds`.
#' @param metatable.site.species.characteristics A metadata table containing site
#'   and species characteristics. It should include variables such as `NADL`,
#'   `Species`, `Latitude`, and `Longitude`.
#' @param replace0 Logical. If `TRUE`, zero seed production values are replaced
#'   by `0.001`. This can be useful before log-transforming seed production.
#'   Default is `TRUE`.
#'
#' @return A preprocessed data frame containing:
#' \describe{
#'   \item{sitenewname}{Combined species-site identifier.}
#'   \item{Year}{Observation year.}
#'   \item{Demand}{Seed demand, with missing values filled where possible.}
#'   \item{seeds}{Seed production, optionally with zeros replaced by `0.001`.}
#'   \item{NADL}{Site identifier.}
#'   \item{Species}{Species name.}
#'   \item{Latitude, Longitude}{Site coordinates.}
#'   \item{FunctionalGroup}{Broad functional group: `"Deciduous"` or `"Coniferous"`.}
#'   \item{NorthSouth}{Geographic grouping based on latitude.}
#'   \item{lag1_seeds}{Previous-year seed production within each site.}
#' }
reformat.harvest.pl.data = function(
  initial.data.harvest.load,
  metatable.site.species.characteristics,
  replace0 = TRUE
) {
  Extend.Harvest.Filled = initial.data.harvest.load %>%
    mutate(sitenewname = paste0(Species, "_", NADL)) %>%
    dplyr::select(sitenewname, Year, Demand, seeds) %>%
    left_join(metatable.site.species.characteristics) %>%
    group_by(sitenewname, NADL, Species, Latitude, Longitude) %>%
    mutate(
      Demand = case_when(
        is.na(Demand) & seeds == 0 ~ 0,
        is.na(Demand) ~ median(Demand, na.rm = TRUE),
        TRUE ~ Demand
      )
    ) %>%
    ungroup() %>%
    mutate(
      FunctionalGroup = if_else(
        Species %in% c("Beech", "Oaks_both"),
        "Deciduous",
        "Coniferous"
      ),
      NorthSouth = if_else(
        Latitude > 52,
        "North",
        "South"
      )
    ) %>%
    arrange(sitenewname, Year) %>%
    group_by(sitenewname) %>%
    mutate(
      lag1_seeds = lag(seeds, n = 1)
    ) %>%
    ungroup()

  # zero_summary <- Extend.Harvest.Filled %>%
  #   group_by(sitenewname) %>%
  #   summarise(
  #     n_obs = n(),
  #     n_zero = sum(seeds == 0, na.rm = TRUE),
  #     prop_zero = n_zero / n_obs
  #   ) %>%
  #   ungroup()

  if (replace0 == T) {
    Extend.Harvest.Filled %>%
      mutate(seeds = if_else(seeds == 0, 0.001, seeds))
  } else {
    Extend.Harvest.Filled
  }
}

#' Aggregate temperature data over selected months
#'
#' Calculates annual mean temperatures for a specified set of months, then
#' computes site-specific temperature anomalies and standardized temperature
#' indices. Output variables are automatically renamed according to the supplied
#' label.
#'
#' @param df A data frame containing monthly temperature data in wide format,
#'   with years as rows and sites as columns.
#' @param months Numeric or character vector specifying the months to include.
#'   For example, `c(3)` for March or `c(4, 5)` for April–May.
#' @param label Character string used to rename output variables (e.g.,
#'   `"Mar"`, `"AprMay"`, `"Dec"`).
#' @param year_col Character. Name of the year column. Default is `"Year"`.
#' @param month_col Character. Name of the month column. Default is `"Month"`.
#' @param drop_cols Character vector of columns to remove before aggregation.
#'   Default is `c("Date", "Day", "Month")`.
#'
#' @return A long-format data frame containing:
#' \describe{
#'   \item{Year}{Year of observation.}
#'   \item{NADL}{Site identifier.}
#'   \item{tmax.harvest_<label>}{Mean temperature across the selected months.}
#'   \item{T_anom_<label>}{Site-specific temperature anomaly relative to the
#'   long-term mean temperature at that site.}
#'   \item{T_z_<label>}{Site-specific standardized temperature (z-score).}
#' }

make_tmax_window <- function(
  df,
  months,
  label,
  year_col = "Year",
  month_col = "Month",
  drop_cols = c("Date", "Day", "Month")
) {
  months_chr <- as.character(months)

  out <- df %>%
    filter(.data[[month_col]] %in% months_chr) %>%
    group_by(.data[[year_col]]) %>%
    dplyr::select(-any_of(drop_cols)) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(
      -all_of(year_col),
      names_to = "NADL",
      values_to = "tmax.harvest"
    ) %>%
    group_by(NADL) %>%
    mutate(
      T_anom = tmax.harvest - mean(tmax.harvest, na.rm = TRUE),
      T_z = (tmax.harvest - mean(tmax.harvest, na.rm = TRUE)) /
        sd(tmax.harvest, na.rm = TRUE)
    ) %>%
    ungroup()

  #rename my var s
  out %>%
    rename(
      !!paste0("tmax.harvest_", label) := tmax.harvest,
      !!paste0("T_anom_", label) := T_anom,
      !!paste0("T_z_", label) := T_z
    )
}

#' Estimate the bootstrap standard error of a trimmed mean
#'
#' Computes a bootstrap estimate of the standard error for a 20% trimmed mean
#' using nonparametric resampling with replacement.
#'
#' @param x A numeric vector.
#' @param nboot Numeric. Number of bootstrap replicates. Default is 2000.
#'
#' @return A single numeric value corresponding to the bootstrap estimate of the
#' standard error of the 20% trimmed mean.
boot_se <- function(x, nboot = 2000) {
  vals <- replicate(nboot, {
    sample_x <- sample(x, replace = TRUE)
    mean(sample_x, trim = 0.2, na.rm = TRUE)
  })
  sd(vals)
}

#' Calculate a 20% trimmed mean
#'
#' Computes the arithmetic mean after removing the lowest and highest 20% of
#' observations. This provides a robust measure of central tendency that is less
#' sensitive to extreme values than the ordinary mean.
#'
#' @param x A numeric vector.
#'
#' @return A single numeric value corresponding to the 20% trimmed mean.
trim_mean <- function(x) mean(x, trim = 0.2, na.rm = TRUE)

#for the theme, but not satisfied finally
#' Apply a custom ggplot2 theme
#'
#' Sets a customized black-and-white ggplot2 theme with simplified panel
#' elements and enlarged axis text. The theme is intended for publication-style
#' figures with minimal visual clutter.
#'
#' @return Invisibly returns the result of `theme_set()`, while changing the
#'   global ggplot2 theme for the current R session.
theme_te <- function()
  theme_set(
    theme_bw(base_size = 16) %+replace%
      theme(
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.background = element_blank(),
        plot.background = element_blank(),
        axis.text = element_text(size = 16)
      )
  )


#' Calculate pairwise synchrony between sites and species
#'
#' Computes pairwise temporal synchrony in seed production between all site pairs
#' for either within-species or cross-species comparisons. Synchrony is measured
#' as the Spearman correlation between seed production time series over
#' overlapping years.
#'
#' @param dat A data frame containing `Species`, `NADL`, `Year`, and `Seeds_cur`.
#' @param coords A data frame containing site coordinates with columns `Species`,
#'   `NADL`, `Longitude`, and `Latitude`.
#' @param sp1 Character. Name of the first species.
#' @param sp2 Character. Name of the second species.
#' @param min_overlap Numeric. Minimum number of overlapping non-missing years
#'   required to calculate a correlation. Default is 5.
#'
#' @return A data frame with one row per site pair, containing:
#' \describe{
#'   \item{site1, site2}{Site identifiers for the compared time series.}
#'   \item{n_overlap}{Number of overlapping years with non-missing seed data.}
#'   \item{correlation}{Spearman correlation between seed production time series.}
#'   \item{lon1, lat1, lon2, lat2}{Coordinates of the two sites.}
#'   \item{distance_km}{Geographic distance between the two sites, in kilometres.}
#'   \item{comparison}{Species comparison label, e.g. `"Pedunculate oak - Sessile oak"`.}
#' }
pairwise_sync <- function(dat, coords, sp1, sp2, min_overlap = 5) {
  #based on Kuba code beforre
  #and cleanded with Claude,ai and then double checked
  dat1 <- dat %>%
    filter(Species == sp1) %>%
    select(NADL, Year, Seeds_cur) %>%
    rename(site1 = NADL, seed1 = Seeds_cur)

  dat2 <- dat %>%
    filter(Species == sp2) %>%
    select(NADL, Year, Seeds_cur) %>%
    rename(site2 = NADL, seed2 = Seeds_cur)

  pairs <- expand_grid(
    site1 = unique(dat1$site1),
    site2 = unique(dat2$site2)
  )

  if (sp1 == sp2) {
    pairs <- pairs %>%
      filter(site1 < site2)
  }

  out <- pairs %>%
    mutate(
      result = map2(site1, site2, function(s1, s2) {
        x <- dat1 %>% filter(site1 == s1)
        y <- dat2 %>% filter(site2 == s2)

        xy <- inner_join(x, y, by = "Year")

        n_overlap <- sum(!is.na(xy$seed1) & !is.na(xy$seed2))

        r <- if (n_overlap >= min_overlap) {
          cor(
            xy$seed1,
            xy$seed2,
            method = "spearman",
            use = "complete.obs"
          )
        } else {
          NA_real_
        }

        tibble(n_overlap = n_overlap, correlation = r)
      })
    ) %>%
    unnest(result)

  coords1 <- coords %>%
    filter(Species == sp1) %>%
    select(site1 = NADL, lon1 = Longitude, lat1 = Latitude)

  coords2 <- coords %>%
    filter(Species == sp2) %>%
    select(site2 = NADL, lon2 = Longitude, lat2 = Latitude)

  out %>%
    left_join(coords1, by = "site1") %>%
    left_join(coords2, by = "site2") %>%
    mutate(
      distance_km = geosphere::distGeo(
        cbind(lon1, lat1),
        cbind(lon2, lat2)
      ) /
        1000,
      comparison = paste(sp1, sp2, sep = " - ")
    )
}
