#' Compute rolling metrics for one variable and rename outputs
#'
#' This helper applies `compute_rolling_metrics()` separately within each site
#' and appends a suffix to the rolling metric columns. It is useful when the same
#' rolling-window summaries are needed for multiple variables, such as seed
#' production, demand, lagged seed production, or climate cues.
#'
#' @param data A data frame containing at least `sitenewname`, `Year`, and the
#'   variable to summarize.
#' @param var Unquoted variable name to summarize.
#' @param window Numeric. Length of the rolling window, in years.
#' @param step Numeric. Step size between rolling windows.
#' @param suffix Character string appended to the names of rolling metric columns.
#'   For example, `"seeds"` gives `mean_value_seeds`, `CVp_seeds`, etc.
#' @param keep Character vector of rolling metric columns to rename. Only columns
#'   that exist in the output of `compute_rolling_metrics()` are renamed.
#'
#' @return A data frame with rolling metrics calculated within each `sitenewname`,
#'   with selected metric columns renamed using the provided suffix.
roll_one <- function(
  data,
  var,
  window,
  step,
  suffix,
  keep = c("mean_value", "sd_value", "CVp", "q25", "q90", "p_zero", "p_nonzero")
) {
  library(dplyr)
  library(rlang)

  varq <- enquo(var)

  out <- data %>%
    group_by(sitenewname) %>%
    group_modify(
      ~ compute_rolling_metrics(.x, !!varq, window = window, step = step)
    ) %>%
    ungroup()

  # Keep only existing rolling columns requested
  keep_exist <- intersect(keep, names(out))

  out %>%
    rename_with(~ paste0(.x, "_", suffix), all_of(keep_exist))
}

#' Fit cue-sensitivity models across rolling-window settings
#'
#' This function computes rolling-window summaries for seed production, demand,
#' lagged seed production, and one or more climate cue variables. It then fits
#' one mixed-effects model per cue to test the association between rolling CV of
#' seed production and mean cue temperature within the same rolling window.
#'
#' @param data A data frame containing seed production, demand, lagged seed
#'   production, climate cues, site identity, and year.
#' @param window Numeric. Length of the rolling window, in years. Default is 10.
#' @param step Numeric. Step size between rolling windows. Default is 5.
#' @param cues Character vector giving the names of climate cue variables to test.
#'   Defaults to `c("tmax.harvest_Dec", "tmax.harvest_Mar",
#'   "tmax.harvest_AprMay")`.
#' @param demand_var Character. Name of the demand variable. Default is `"Demand"`.
#' @param seeds_var Character. Name of the current seed production variable.
#'   Default is `"Seeds_cur"`.
#' @param seeds_lag_var Character. Name of the lagged seed production variable.
#'   Default is `"Seeds_prev"`.
#' @param site_var Character. Name of the site identity variable. Default is
#'   `"sitenewname"`.
#' @param year_var Character. Name of the year variable. Currently retained for
#'   consistency, but the function assumes the rolling output contains `Year`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{rolling_data}{A data frame containing merged rolling-window summaries
#'   for seeds, demand, lagged seeds, and climate cues.}
#'   \item{cue_results}{A tibble containing model results for each cue, including
#'   window size, step size, cue name, estimate, standard error, z statistic,
#'   p-value, and AIC.}
#' }
#'
fit_cue_models_over_window <- function(
  data,
  window = 10,
  step = 5,
  cues = c("tmax.harvest_Dec", "tmax.harvest_Mar", "tmax.harvest_AprMay"),
  demand_var = "Demand",
  seeds_var = "Seeds_cur",
  seeds_lag_var = "Seeds_prev",
  site_var = "sitenewname",
  year_var = "Year"
) {
  library(dplyr)
  library(glmmTMB)

  # rolling for seeds, demand and the lag seeds
  rolling_seeds <- roll_one(data, !!sym(seeds_var), window, step, "seeds")
  rolling_demand <- roll_one(data, !!sym(demand_var), window, step, "demand")
  rolling_lag <- roll_one(
    data,
    !!sym(seeds_lag_var),
    window,
    step,
    "lag1_seeds"
  )

  # rolling for each cue
  rolling_cues <- lapply(cues, function(cn) {
    suf <- sub("tmax\\.harvest_", "T_", cn)
    roll_one(
      data,
      !!sym(cn),
      window,
      step,
      suf,
      keep = c("mean_value", "sd_value")
    )
  })

  # merge all rolling datasets (by site + Year)
  rolling_all <- rolling_seeds %>%
    left_join(
      rolling_demand %>%
        select(
          sitenewname,
          Year,
          starts_with("mean_value_demand"),
          starts_with("sd_value_demand"),
          starts_with("p_zero_demand"),
          starts_with("Year_center")
        ),
      by = c("sitenewname", "Year", "Year_center")
    ) %>%
    left_join(
      rolling_lag %>%
        select(
          sitenewname,
          Year,
          starts_with("mean_value_lag1_seeds"),
          starts_with("CVp_lag1_seeds"),
          starts_with("p_zero_lag1_seeds"),
          starts_with("Year_center")
        ),
      by = c("sitenewname", "Year", "Year_center")
    )

  for (rc in rolling_cues) {
    rolling_all <- rolling_all %>%
      left_join(
        rc %>%
          select(
            sitenewname,
            Year,
            starts_with("mean_value_"),
            starts_with("sd_value_"),
            starts_with("Year_center")
          ),
        by = c("sitenewname", "Year", "Year_center")
      )
  }

  # get my var of interest here
  rolling_all <- rolling_all %>%
    mutate(
      log_CVp_seeds = log(CVp_seeds),
      log_mean_seeds = log1p(mean_value_seeds),
      logit_prop_zero_seeds = car::logit(p_zero_seeds),
    ) %>%
    mutate(Year_c = Year_center - mean(Year_center, na.rm = TRUE))

  # Fit one model per cue
  results <- lapply(cues, function(cn) {
    suf <- sub("tmax\\.harvest_", "T_", cn)
    cue_col <- paste0("mean_value_", suf)

    # model formula original log_CVp_seeds ~ cue + (1|site) + log1p(mean_value_demand)

    fml <- as.formula(paste0(
      "log_CVp_seeds ~ ",
      cue_col,
      " + (1 | ",
      site_var,
      ") + log1p(mean_value_demand)"
    ))

    m <- glmmTMB(fml, data = rolling_all)

    sm <- summary(m)$coefficients$cond
    tibble::tibble(
      window = window,
      step = step,
      cue = cn,
      cue_col = cue_col,
      estimate = sm[cue_col, "Estimate"],
      se = sm[cue_col, "Std. Error"],
      z = sm[cue_col, "z value"],
      p = sm[cue_col, "Pr(>|z|)"],
      AIC = AIC(m)
    )
  })
  #export my lists
  list(
    rolling_data = rolling_all,
    cue_results = dplyr::bind_rows(results)
  )
}

#' Fit temporal trend models across rolling-window settings
#'
#' This function computes rolling-window summaries for seed production, demand,
#' lagged seed production, and optional climate cue variables. It then fits a
#' mixed-effects model testing whether the rolling coefficient of variation of
#' seed production changes through time.
#'
#' @param data A data frame containing seed production, demand, lagged seed
#'   production, site identity, year, and optional climate cue variables.
#' @param window Numeric. Length of the rolling window, in years. Default is 10.
#' @param step Numeric. Step size between rolling windows. Default is 5.
#' @param cues Character vector giving the names of climate cue variables to
#'   summarize. These are included in the returned rolling dataset but are not
#'   used in the trend model. Default is `c("tmax.harvest_Dec",
#'   "tmax.harvest_Mar", "tmax.harvest_AprMay")`.
#' @param demand_var Character. Name of the demand variable. Default is `"Demand"`.
#' @param seeds_var Character. Name of the current seed production variable.
#'   Default is `"Seeds_cur"`.
#' @param seeds_lag_var Character. Name of the lagged seed production variable.
#'   Default is `"Seeds_prev"`.
#' @param site_var Character. Name of the site identity variable. Default is
#'   `"sitenewname"`.
#' @param year_var Character. Name of the year variable. Currently retained for
#'   consistency, but the function assumes the rolling output contains `Year`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{rolling_data}{A data frame containing merged rolling-window summaries.}
#'   \item{cue_results}{A tibble of fixed-effect model estimates, including
#'   window size, step size, AIC, and number of observations.}
#' }
fit_Trends_models_over_window <- function(
  data,
  window = 10,
  step = 5,
  cues = c("tmax.harvest_Dec", "tmax.harvest_Mar", "tmax.harvest_AprMay"),
  demand_var = "Demand",
  seeds_var = "Seeds_cur",
  seeds_lag_var = "Seeds_prev",
  site_var = "sitenewname",
  year_var = "Year"
) {
  library(dplyr)
  library(glmmTMB)

  # rolling for seeds, demand, lag
  rolling_seeds <- roll_one(data, !!sym(seeds_var), window, step, "seeds")
  rolling_demand <- roll_one(data, !!sym(demand_var), window, step, "demand")
  rolling_lag <- roll_one(
    data,
    !!sym(seeds_lag_var),
    window,
    step,
    "lag1_seeds"
  )

  # rolling for each cue (only mean/sd needed)
  rolling_cues <- lapply(cues, function(cn) {
    suf <- sub("tmax\\.harvest_", "T_", cn)
    roll_one(
      data,
      !!sym(cn),
      window,
      step,
      suf,
      keep = c("mean_value", "sd_value")
    )
  })

  # merge all rolling datasets (by site + Year)
  rolling_all <- rolling_seeds %>%
    left_join(
      rolling_demand %>%
        select(
          sitenewname,
          Year,
          starts_with("mean_value_demand"),
          starts_with("sd_value_demand"),
          starts_with("p_zero_demand"),
          starts_with("Year_center")
        ),
      by = c("sitenewname", "Year", "Year_center")
    ) %>%
    left_join(
      rolling_lag %>%
        select(
          sitenewname,
          Year,
          starts_with("mean_value_lag1_seeds"),
          starts_with("CVp_lag1_seeds"),
          starts_with("p_zero_lag1_seeds"),
          starts_with("Year_center")
        ),
      by = c("sitenewname", "Year", "Year_center")
    )

  for (rc in rolling_cues) {
    rolling_all <- rolling_all %>%
      left_join(
        rc %>%
          select(
            sitenewname,
            Year,
            starts_with("mean_value_"),
            starts_with("sd_value_"),
            starts_with("Year_center")
          ),
        by = c("sitenewname", "Year", "Year_center")
      )
  }

  # derived vars
  rolling_all <- rolling_all %>%
    mutate(
      log_CVp_seeds = log(CVp_seeds),
      log_mean_seeds = log1p(mean_value_seeds),
      logit_prop_zero_seeds = car::logit(p_zero_seeds),
    ) %>%
    mutate(Year_c = Year_center - mean(Year_center, na.rm = TRUE))

  # Fit one model
  fml <- as.formula(paste0(
    "log_CVp_seeds ~ ",
    "Year_c + (Year_c | ",
    site_var,
    ") + log1p(mean_value_demand)"
  ))

  m <- glmmTMB(fml, data = rolling_all)

  results = broom.mixed::tidy(m, effects = "fixed") %>%
    mutate(
      window = window,
      step = step,
      AIC = AIC(m),
      nobs = broom.mixed::glance(m)$nobs
    ) %>%
    dplyr::select(-effect, -component, -statistic)

  list(
    rolling_data = rolling_all,
    cue_results = dplyr::bind_rows(results)
  )
}

#' Compute rolling-window reproductive metrics
#'
#' This function calculates rolling-window summaries for a selected variable
#' within a time series. It is designed for masting analyses where reproductive
#' metrics such as mean seed production, standard deviation, coefficient of
#' variation, quantiles, and the proportion of zero years are calculated over
#' moving time windows.
#'
#' @param df A data frame containing a `Year` column and the variable to be
#'   summarized.
#' @param variable Unquoted name of the variable to summarize.
#' @param window Numeric. Length of the rolling window, in years. Default is 10.
#' @param step Numeric. Step size used to retain rolling-window estimates.
#'   Default is 1, meaning that a rolling estimate is retained for every year.
#'   If `step = 5`, only every fifth complete rolling window is retained.
#'
#' @return A data frame containing the original data plus rolling-window
#'   metrics:
#' \describe{
#'   \item{mean_value}{Rolling-window mean.}
#'   \item{sd_value}{Rolling-window standard deviation.}
#'   \item{CVp}{Rolling-window coefficient of variation, calculated as
#'   `sd_value / mean_value`. Values are set to `NA` when the mean is zero.}
#'   \item{p_zero}{Proportion of observations equal to zero in the rolling window.}
#'   \item{p_nonzero}{Proportion of observations greater than zero in the rolling window.}
#'   \item{q25}{25th percentile within the rolling window.}
#'   \item{q90}{90th percentile within the rolling window.}
#'   \item{Year_center}{Approximate center year of the rolling window.}
#'   \item{Year_right_align}{Right-aligned year of the rolling window.}
#' }
#'
compute_rolling_metrics <- function(df, variable, window = 10, step = 1) {
  library(slider) #for the rolling
  library(rlang) #for enqouoing the variable name

  var <- enquo(variable)

  df <- df %>% arrange(Year)

  results = df %>%
    mutate(
      mean_value = slide_dbl(
        !!var,
        mean,
        .before = window - 1,
        .complete = TRUE,
        na.rm = TRUE
      ),

      sd_value = slide_dbl(
        !!var,
        sd,
        .before = window - 1,
        .complete = TRUE,
        na.rm = TRUE
      ),

      CVp = if_else(mean_value == 0, NA_real_, sd_value / mean_value),
      #kCV = sqrt((CVp)^2 / (1 + (CVp)^2)),
      #for another paper they used another method
      #https://pmc.ncbi.nlm.nih.gov/articles/PMC9196089/#CIT0044
      #They found similar "reboustness" if they replace IQR by MAD
      #But I think I might have done some shit because values look super weird
      # median_value = slide_dbl(
      #   !!var,
      #   median,
      #   .before = window - 1,
      #   .complete = TRUE,
      #   na.rm = TRUE
      # ),
      # iqr_value = slide_dbl(
      #   !!var,
      #   EnvStats::iqr,
      #   .before = window - 1,
      #   .complete = TRUE,
      #   na.rm = TRUE
      # ),
      # RobustCV = 0.75 * iqr_value / median_value,
      p_zero = slide_dbl(
        !!var,
        ~ mean(.x == 0, na.rm = TRUE), #mean(c(TRUE, FALSE, TRUE, FALSE, FALSE)) would give proportion here for example 0.4 (2/5)
        .before = window - 1,
        .complete = TRUE
      ),
      p_nonzero = slide_dbl(
        !!var,
        ~ mean(.x > 0, na.rm = TRUE),
        .before = window - 1,
        .complete = TRUE
      ),
      q25 = slide_dbl(
        !!var,
        ~ quantile(.x, 0.25, na.rm = TRUE),
        .before = window - 1,
        .complete = TRUE
      ),

      q90 = slide_dbl(
        !!var,
        ~ quantile(.x, 0.90, na.rm = TRUE),
        .before = window - 1,
        .complete = TRUE
      ),

      Year_center = Year - floor((window - 1) / 2),
      Year_right_align = Year
    ) %>%
    filter(!is.na(mean_value))

  # I did same as Jessie to remove step years with filterstep function
  if (step > 1) {
    results <- results %>%
      mutate(.idx = row_number()) %>%
      filter((.idx - 1) %% step == 0) %>%
      select(-.idx)
  }

  return(results)
}


#' Compute a rolling mean
#' CODE FROM JESSIE FOEST , Ecology Letters
#'
#' Calculates a right-aligned rolling mean for a selected variable using
#' `zoo::rollapplyr()`. Each value represents the mean of the current
#' observation and the previous `win - 1` observations.
#'
#' @param data A data frame containing the variable to summarize.
#' @param val Character string giving the name of the variable for which the
#'   rolling mean should be calculated. Default is `"collection_current"`.
#' @param win Numeric. Length of the rolling window. Default is 10.
#'
#' @return A numeric vector of the same length as the input variable containing
#'   the rolling means. The first `win - 1` values are returned as `NA` because
#'   a complete window is required.
rollingmean <- function(data, val = "collection_current", win = 10) {
  x <- data[[val]]
  zoo::rollapplyr(
    # right
    x,
    width = win,
    fill = NA_real_,
    FUN = function(v) {
      mean(v, na.rm = TRUE)
    }
  )
}


#' Compute a rolling standard deviation
#'Code obtained from Jessie Foest, Ecology Letters
#' Calculates a right-aligned rolling standard deviation for a selected variable
#' using `zoo::rollapplyr()`. Each value represents the standard deviation of
#' the current observation and the previous `win - 1` observations.
#'
#' @param data A data frame containing the variable to summarize.
#' @param val Character string giving the name of the variable for which the
#'   rolling standard deviation should be calculated. Default is
#'   `"collection_current"`.
#' @param win Numeric. Length of the rolling window. Default is 10.
#'
#' @return A numeric vector of the same length as the input variable containing
#'   the rolling standard deviations. The first `win - 1` values are returned
#'   as `NA` because a complete window is required.

rollingsd <- function(data, val = "collection_current", win = 10) {
  x <- data[[val]]
  zoo::rollapplyr(
    # right
    x,
    width = win,
    fill = NA_real_,
    FUN = function(v) {
      sd(v, na.rm = TRUE)
    }
  )
}

#' Compute a rolling quantile
#'Code from Jessie Foest, Ecology Letters
#' Calculates a right-aligned rolling quantile for a selected variable using
#' `zoo::rollapplyr()`. Each value represents the specified quantile of the
#' current observation and the previous `win - 1` observations.
#'
#' @param data A data frame containing the variable to summarize.
#' @param val Character string giving the name of the variable for which the
#'   rolling quantile should be calculated. Default is `"collection_current"`.
#' @param prob Numeric. Quantile probability to compute, typically between 0 and
#'   1 (e.g., `0.25` for the 25th percentile, `0.90` for the 90th percentile).
#' @param win Numeric. Length of the rolling window. Default is 10.
#'
#' @return A numeric vector of the same length as the input variable containing
#'   the rolling quantile values. The first `win - 1` values are returned as
#'   `NA` because a complete window is required.
rollingquant <- function(data, val = "collection_current", prob, win = 10) {
  x <- data[[val]]
  zoo::rollapplyr(
    x,
    width = win,
    fill = NA_real_,
    FUN = function(v) {
      as.numeric(stats::quantile(v, probs = prob, na.rm = TRUE))
    }
  )
}

#' Compute a rolling coefficient of variation
#'
#' Calculates the rolling coefficient of variation (CV) for a selected variable
#' using rolling means and rolling standard deviations computed over a moving
#' time window.
#'
#' @param data A data frame containing the variable to summarize.
#' @param val Character string giving the name of the variable for which the
#'   rolling coefficient of variation should be calculated. Default is
#'   `"collection_current"`.
#' @param win Numeric. Length of the rolling window. Default is 10.
#'
#' @return A numeric vector of the same length as the input variable containing
#'   rolling coefficients of variation, calculated as:
#'
#'   \deqn{CV = \frac{\sigma}{\mu}}
#'
#'   where \eqn{\sigma} is the rolling standard deviation and \eqn{\mu} is the
#'   rolling mean. The first `win - 1` values are returned as `NA` because a
#'   complete window is required.
rollingCVp <- function(data, val = "collection_current", win = 10) {
  mu <- rollingmean(data, val = val, win = win)
  sig <- rollingsd(data, val = val, win = win)
  sig / mu
}


#' Compute pairwise correlations with a minimum number of overlapping years
#'Code adaptaed from Jakub Szymkowiak here and below
#' Calculates a pairwise Spearman correlation matrix among sites while requiring
#' a minimum number of years with non-missing observations for each site pair.
#' Site pairs with fewer overlapping observations than `min.no.years` are
#' assigned `NA`.
#'
#' @param x A wide-format data frame or matrix where columns are sites and rows
#'   are years or time steps.
#' @param sites Character vector giving the site names to include. These should
#'   match column names in `x`.
#' @param min.no.years Numeric. Minimum number of overlapping non-missing
#'   observations required to calculate a pairwise correlation.
#'
#' @return A square correlation matrix with site names as row and column names.
#'   Values are Spearman correlation coefficients. Pairs with insufficient
#'   overlap are returned as `NA`.
corNyrs <- function(x, sites, min.no.years) {
  cor.mat <- matrix(NA, ncol = ncol(x), nrow = ncol(x)) ## correlation matrix template to be filled by a loop
  rownames(cor.mat) <- sites
  colnames(cor.mat) <- sites

  ## loop calculating pairwise correlations taking into account no. of overlapping years (as specified in min.no.years argument)
  for (i in 1:length(sites)) {
    for (j in 1:length(sites)) {
      site1.tmp <- sites[i]
      site2.tmp <- sites[j]
      cor.data.tmp <- subset(x, select = c(site1.tmp, site2.tmp))
      cor.data.tmp <- na.omit(cor.data.tmp)

      {
        if (nrow(cor.data.tmp) >= min.no.years) {
          cor.tmp <- cor(
            cor.data.tmp,
            use = "pairwise.complete.obs",
            method = "spearman"
          )
          cor.mat[i, j] <- cor.tmp[1, 2]
          cor.mat[j, i] <- cor.tmp[1, 2]
        } else {
          cor.mat[i, j] <- NA
          cor.mat[j, i] <- NA
        }
      }
    }
  }

  return(cor.mat) ## correlation matrix
}

#' Convert a symmetric matrix into a pairwise vector format
#'
#' Transforms a square matrix (e.g., a correlation, similarity, or distance
#' matrix) into a long-format data frame containing all unique pairwise
#' combinations of rows/columns and their associated matrix values.
#'
#' @param x A square matrix with row and column names representing the entities
#'   being compared (e.g., sites).
#' @param var.name Character string giving the name of the output column
#'   containing matrix values.
#' @param id.vars Character vector of length two specifying the names of the
#'   columns identifying the paired entities.
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{<id.vars[1]>}{First entity in the pair.}
#'   \item{<id.vars[2]>}{Second entity in the pair.}
#'   \item{<var.name>}{Value extracted from the matrix for that pair.}
#' }
matvec <- function(x, var.name, id.vars) {
  data.sim <- t(combn(colnames(x), 2))
  data.sim <- data.frame(data.sim, sim = x[data.sim])
  colnames(data.sim)[1:2] <- id.vars
  colnames(data.sim)[3] <- var.name
  return(data.sim)
}

## function used to calculate distance-decay of within-tail synchrony -----
## this function is copy-pasted here for convenience but has not been developed by the authors of the present study
## the original source of the code: Walter, J.A., Castorani, M.C., Bell, T.W., Sheppard, L., Cavanaugh, K.C. & Reuman, D.C. (2022). Tail-dependent spatial synchrony arises from nonlinear driver–response relationships. Ecology Letters, 25, 1189–1201
splineFit <- function(
  distmat,
  zmat,
  nresamp = 1000,
  quantiles = c(0, 0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99, 1)
) {
  triang <- lower.tri(distmat)
  distmat <- distmat
  xemp <- distmat[triang]
  yemp <- zmat[triang]
  drop.NaNs <- !is.na(yemp)
  dfs = sqrt(nrow(distmat))
  out <- list()

  emp.spline <- smooth.spline(xemp[drop.NaNs], yemp[drop.NaNs], df = dfs)
  out$emp.spline <- emp.spline

  resamp.splines <- matrix(NA, nrow = nresamp, ncol = length(emp.spline$y))
  for (ii in 1:nresamp) {
    shuffle <- sample(1:nrow(distmat), size = nrow(distmat), replace = TRUE)
    xres <- distmat[shuffle, shuffle][triang]
    yres <- zmat[shuffle, shuffle][triang]
    drop.NaNs <- !is.na(yres)
    xres <- xres[drop.NaNs]
    yres <- yres[drop.NaNs]
    yres <- yres[!(xres == 0)]
    xres <- xres[!(xres == 0)]
    res.spline <- smooth.spline(xres, yres, df = dfs)
    resamp.splines[ii, ] <- predict(res.spline, x = emp.spline$x)$y
  }
  out$resamp.splines <- resamp.splines
  out$spline.quantiles <- apply(resamp.splines, 2, quantile, probs = quantiles)
  return(out)
}

#' Rescale values to a specified range
#'
#' Linearly rescales a numeric vector, matrix, or data frame to a user-defined
#' range. The minimum value of the input is mapped to `first` and the maximum
#' value is mapped to `last`, with all intermediate values transformed
#' proportionally.
#'
#' @param x A numeric vector, matrix, or data frame to be rescaled.
#' @param first Numeric. The lower bound of the desired output range.
#' @param last Numeric. The upper bound of the desired output range.
#'
#' @return An object of the same dimensions as `x`, with values rescaled to the
#'   interval [`first`, `last`].
ReScale <- function(x, first, last) {
  (last - first) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) *
    (x - min(x, na.rm = TRUE)) +
    first
}

#' Calculate rolling-window spatial synchrony in seed production
#'
#' This function estimates spatial synchrony in seed production across sites
#' using moving time windows. Within each window, seed production time series are
#' log-transformed, pairwise Spearman correlations are calculated between sites,
#' and correlations are averaged at the site level to obtain a site-specific
#' synchrony estimate.
#'
#' @param dataset A data frame containing seed production data. It must include
#'   the columns `sitenewname`, `NADL`, `Year`, `Seeds_cur`, `Longitude`,
#'   `Latitude`, and `Species`.
#' @param n.years Numeric. Length of the moving window in years. Default is 10.
#' @param step Numeric. Step size used to move the rolling window forward.
#'   Default is 10.
#'
#' @return A data frame with one row per site per time window, containing:
#' \describe{
#'   \item{NADL}{Site identifier.}
#'   \item{mean.synch}{Mean pairwise synchrony of the site with all other sites
#'   in the same window. Synchrony is rescaled to the interval [0, 1].}
#'   \item{mean.dist}{Mean geographic distance from the site to all other sites.}
#'   \item{sd.synch}{Standard deviation of pairwise synchrony values for the site.}
#'   \item{sd.dist}{Standard deviation of pairwise geographic distances.}
#'   \item{Species}{Species name.}
#'   \item{window.no}{Window number.}
#'   \item{start.year}{First year of the moving window.}
#'   \item{end.year}{Last year of the moving window.}
#' }
calculate.sync = function(dataset, n.years = 10, step = 10) {
  dataset = dataset %>%
    group_by(sitenewname, Longitude, Latitude) %>%
    mutate(has_0 = ifelse(any(Seeds_cur == 0), "True", "False")) %>%
    mutate(
      log.seed = ifelse(has_0 == 'True', log(1 + Seeds_cur), log(Seeds_cur))
    ) %>%
    ungroup()

  min.year <- min(dataset$Year)
  max.year <- max(dataset$Year) - n.years + 1
  start.year <- seq(min.year, max.year, step)
  end.year <- start.year + n.years - 1
  years <- data.frame(start.year, end.year)
  # years[6,] <- c(2013, 2021) ## manually add the last 9-yr window
  years$window.no <- 1:nrow(years)
  spatial.cor.time <- data.frame()

  ## rolling window synchrony loop
  for (j in 1:nrow(years)) {
    ## seeds
    ## select window
    start.year <- years$start.year[j]
    end.year <- years$end.year[j]
    ## filter dataset
    data.tmp <- subset(
      dataset,
      dataset$Year >= start.year & dataset$Year <= end.year
    )

    ## coordinates and distance matrix
    xy.tmp <- data.tmp %>%
      group_by(NADL) %>%
      summarise(
        Longitude = mean(Longitude, na.rm = TRUE),
        Latitude = mean(Latitude, na.rm = TRUE)
      )
    xy.tmp <- na.omit(xy.tmp)

    coords <- as.matrix(xy.tmp[, 2:3])
    geog.dist <- geosphere::distm(coords)
    colnames(geog.dist) <- xy.tmp$NADL
    rownames(geog.dist) <- xy.tmp$NADL
    geog.dist <- geog.dist / 1000 ## so that the distance is in km

    ## seed production data into wide format
    data.tmp <- data.tmp %>% dplyr::select(NADL, Year, log.seed)
    data.tmp <- data.tmp[data.tmp$NADL %in% xy.tmp$NADL, ] ## match data frames
    # data.tmp$Species <- NULL
    data.tmp <- data.tmp %>%
      pivot_wider(names_from = NADL, values_from = log.seed)
    data.tmp <- data.tmp[, 2:ncol(data.tmp)]
    sites <- colnames(data.tmp)
    # data.tmp[is.na(data.tmp)] <- 0
    cormat <- corNyrs(data.tmp, sites = sites, min.no.years = 5)
    #cormat <- cor(data.tmp, method = "spearman", use = "pairwise.complete.obs")

    ## spatial spline
    cormat <- ReScale(cormat, 0, 1)
    #sncf <- splineFit(geog.dist, cormat, quantiles = c(0.025, 0.5, 0.975))
    diag(cormat) <- NA
    seeds.synch <- matvec(
      cormat,
      var.name = "seeds.synch",
      id.vars = c("site1", "site2")
    )
    diag(geog.dist) <- NA ## matrix of pairwise distances
    spat.dist <- matvec(
      geog.dist,
      var.name = "spat.dist",
      id.vars = c("site1", "site2")
    )

    plotdata <- merge(
      seeds.synch,
      spat.dist,
      by = c("site1", "site2"),
      all = TRUE
    )
    plotdata$spat.dist <- plotdata$spat.dist / 1000 #convert km

    sncf.plot = plotdata %>%
      group_by(site1) %>%
      summarise(
        mean.synch = mean(seeds.synch, na.rm = TRUE),
        mean.dist = mean(spat.dist, na.rm = TRUE),
        sd.synch = sd(seeds.synch, na.rm = TRUE),
        sd.dist = sd(spat.dist, na.rm = TRUE)
      ) %>%
      rename(NADL = site1) %>%
      mutate(
        Species = unique(dataset$Species),
        window.no = years$window.no[j],
        start.year = years$start.year[j],
        end.year = years$end.year[j]
      )

    ## extract results into data frame for easy plotting
    #sncf.plot <- data.frame(
    #Distance = sncf$emp.spline$x,
    #Correlation = sncf$emp.spline$y,
    #LCL = sncf$spline.quantiles[1, ],
    #UCL = sncf$spline.quantiles[3, ],
    # Species = unique(dataset$Species),
    # window.no = years$window.no[j],
    # start.year = years$start.year[j],
    # end.year = years$end.year[j]
    #)
    spatial.cor.time <- rbind(spatial.cor.time, sncf.plot)
  }
  ## export
  return(spatial.cor.time)
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


#' Generate prediction curves from a glmmTMB model
#'
#' Creates population-level predictions from a fitted `glmmTMB` model across the
#' observed range of a focal predictor. The function can optionally
#' back-transform predictions from the log scale and also returns partial
#' residuals generated with `ggeffects`.
#'
#' @param model A fitted `glmmTMB` model.
#' @param data A data frame containing the variables used in the model.
#' @param x_var Character. Name of the focal predictor to vary along the
#'   prediction curve.
#' @param demand_var Character. Name of the demand variable to hold constant.
#'   Default is `"mean_value_demand"`.
#' @param site_var Character. Name of the site grouping variable. Default is
#'   `"sitenewname"`.
#' @param backtransform Character. Either `"none"` or `"exp"`. Use `"exp"` when
#'   the model response was fitted on the log scale and predictions should be
#'   returned on the original response scale.
#' @param n_points Numeric. Number of points used to construct the prediction
#'   sequence. Default is 200.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{newdat}{A data frame containing the prediction grid, fitted values,
#'   standard errors, and approximate 95% confidence intervals.}
#'   \item{part_resid}{A tibble containing partial residuals produced by
#'   `ggeffects::residualize_over_grid()`.}
#' }
predict_glmmtmb_curve <- function(
  model,
  data,
  x_var,
  demand_var = "mean_value_demand",
  site_var = "sitenewname",
  backtransform = c("none", "exp"),
  n_points = 200
) {
  backtransform <- match.arg(backtransform)

  # sequence of x values on original scale
  x_seq <- seq(
    min(data[[x_var]], na.rm = TRUE),
    max(data[[x_var]], na.rm = TRUE),
    length.out = n_points
  )

  # add log varible
  newdat <- data.frame(
    x_tmp = x_seq,
    demand_tmp = median(data[[demand_var]], na.rm = TRUE)
  )

  names(newdat)[1] <- x_var
  names(newdat)[2] <- demand_var

  # add site column for population-level prediction
  newdat[[site_var]] <- NA

  # prediction on model scale
  pred <- predict(
    model,
    newdata = newdat,
    re.form = NA, # fixed effects only
    se.fit = TRUE,
    type = "response"
  )

  newdat$fit <- pred$fit
  newdat$se <- pred$se.fit
  newdat$lwr <- newdat$fit - 1.96 * newdat$se
  newdat$upr <- newdat$fit + 1.96 * newdat$se

  # optional back-transformation
  if (backtransform == "exp") {
    newdat <- newdat %>%
      mutate(
        fit = exp(fit),
        lwr = exp(lwr),
        upr = exp(upr)
      )
  }

  dat = ggeffects::ggpredict(
    model,
    terms = paste0(x_var, " [all]")
  )

  part_resid <- tibble(
    ggeffects::residualize_over_grid(
      dat,
      model
    )
  )

  list(
    newdat = newdat,
    part_resid = part_resid
  )
}


#' Generate prediction data for a Tweedie cue-response model
#'
#' Creates a prediction data frame from a fitted `glmmTMB` Tweedie model across
#' the observed range of a climate cue. The function uses the original cue
#' variable for the prediction grid, converts it to the scaled variable used in
#' the model, and returns fitted values with approximate 95% confidence
#' intervals.
#'
#' @param model A fitted `glmmTMB` model, typically with `family = tweedie()`.
#' @param data A data frame containing the original cue variable, the scaled cue
#'   variable, `Seeds_prev`, `Demand`, and `sitenewname`.
#' @param x_var Character. Name of the original, unscaled cue variable.
#' @param x_scaled_var Character. Name of the scaled cue variable used in the
#'   fitted model.
#'
#' @return A data frame containing:
#' \describe{
#'   \item{x_var}{Prediction values on the original cue scale.}
#'   \item{x_scaled_var}{Prediction values on the scaled cue scale used in the model.}
#'   \item{Seeds_prev}{Previous seed production, fixed at its median.}
#'   \item{Demand}{Demand, fixed at its median.}
#'   \item{sitenewname}{Set to `NA` for population-level predictions.}
#'   \item{fit}{Predicted seed production on the response scale.}
#'   \item{se}{Standard error of the prediction.}
#'   \item{lwr}{Lower approximate 95% confidence interval, truncated at zero.}
#'   \item{upr}{Upper approximate 95% confidence interval.}
#' }
plot_tweedie_cue_data <- function(
  model,
  data,
  x_var,
  x_scaled_var
) {
  # coefficient from the fitted model
  sm <- summary(model)$coefficients$cond
  beta <- sm[x_scaled_var, "Estimate"]
  se <- sm[x_scaled_var, "Std. Error"]
  pval <- sm[x_scaled_var, "Pr(>|z|)"]

  # original-scale range for predictions
  #with min max
  x_seq <- seq(
    min(data[[x_var]], na.rm = TRUE),
    max(data[[x_var]], na.rm = TRUE),
    length.out = 200
  )

  # scaling values used in the data
  x_mean <- mean(data[[x_var]], na.rm = TRUE)
  x_sd <- sd(data[[x_var]], na.rm = TRUE)

  newdat <- data.frame(
    #Seeds_prev = median(data$Seeds_prev, na.rm = TRUE),
    Seeds_prev = median(data$Seeds_prev, na.rm = T),
    Demand = median(data$Demand, na.rm = TRUE),
    sitenewname = NA
  )

  newdat <- newdat[rep(1, length(x_seq)), , drop = FALSE]
  newdat[[x_var]] <- x_seq
  newdat[[x_scaled_var]] <- (x_seq - x_mean) / x_sd

  pred <- predict(
    model,
    newdata = newdat,
    type = "response",
    se.fit = TRUE,
    re.form = NA
  )

  newdat$fit <- pred$fit
  newdat$se <- pred$se.fit
  newdat$lwr <- pmax(0, newdat$fit - 1.96 * newdat$se)
  newdat$upr <- newdat$fit + 1.96 * newdat$se
  return(as.data.frame(newdat))
}
