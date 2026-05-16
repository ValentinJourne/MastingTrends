# Data Loading Utilities for Masting Trends Analysis
# Functions to load and preprocess Polish harvest dataset

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

  #to rename vars
  out %>%
    rename(
      !!paste0("tmax.harvest_", label) := tmax.harvest,
      !!paste0("T_anom_", label) := T_anom,
      !!paste0("T_z_", label) := T_z
    )
}

#calcualted booted se
boot_se <- function(x, nboot = 2000) {
  vals <- replicate(nboot, {
    sample_x <- sample(x, replace = TRUE)
    mean(sample_x, trim = 0.2, na.rm = TRUE)
  })
  sd(vals)
}

trim_mean <- function(x) mean(x, trim = 0.2, na.rm = TRUE)

#for the theme, but not satisfied finally
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
