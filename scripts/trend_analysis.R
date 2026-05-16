# Trend Analysis Functions for Mast Seeding

#check sensivity slopes
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

  # Fit one model per cue
  results <- lapply(cues, function(cn) {
    suf <- sub("tmax\\.harvest_", "T_", cn)
    cue_col <- paste0("mean_value_", suf)

    # model formula: log_CVp_seeds ~ cue + (1|site) + log1p(mean_value_demand)

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

  list(
    rolling_data = rolling_all,
    cue_results = dplyr::bind_rows(results)
  )
}


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

# Statistical methods to detect and analyze masting trends

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


# Get standard devation for each 10 year window

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

## Get quantiles running for each time-series ====

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


rollingCVp <- function(data, val = "collection_current", win = 10) {
  mu <- rollingmean(data, val = val, win = win)
  sig <- rollingsd(data, val = val, win = win)
  sig / mu
}


#from kuba to get spatil sync
## corNyrs function
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

##
## matvec
## Function transforming matrices into vectors.
##

## matvec function
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

## ReScale function -----
ReScale <- function(x, first, last) {
  (last - first) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) *
    (x - min(x, na.rm = TRUE)) +
    first
}

## windows setup
#step for windows shift
#and window size
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

#plot tweedie reg
plot_tweedie_cue <- function(
  model,
  data,
  x_var,
  x_scaled_var,
  panel_title,
  x_lab,
  ylab.text = 60000,
  col = "#E64B35FF"
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

  # x_seq <- seq(
  #   quantile(data[[x_var]], 0.05, na.rm = TRUE),
  #   quantile(data[[x_var]], 0.95, na.rm = TRUE),
  #   length.out = 200
  # )

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

  lab <- paste0(
    "slope = ",
    round(beta, 3),
    "\nSE = ",
    round(se, 3),
    "\np = ",
    signif(pval, 2)
  )

  bin_width <- 0.3

  dat = ggeffects::ggpredict(
    model,
    terms = paste0(x_scaled_var, " [all]")
  )

  part_resid <- tibble(
    ggeffects::residualize_over_grid(
      dat,
      model
    )
  )

  part_resid$x_original <- part_resid$x * x_sd + x_mean

  part_resid$bin_center <- round(part_resid$x_original / bin_width) * bin_width

  summary_df <- part_resid %>%
    group_by(bin_center) %>%
    summarise(
      med_y = median(predicted, na.rm = TRUE),
      q25 = quantile(predicted, 0.25, na.rm = TRUE),
      q75 = quantile(predicted, 0.75, na.rm = TRUE),
      n = n(),
      mad = mad(predicted, na.rm = TRUE),
      se_robust = mad / sqrt(n)
      #tm = trim_mean(predicted),
      #se = boot_se(predicted),
    ) %>%
    filter(n > 1)

  #ggplot(data, aes(x = .data[[x_var]], y = Seeds_cur)) +
  ggplot(summary_df, aes(x = bin_center, y = med_y)) +
    #geom_point(alpha = 0.2, size = 1, col = "grey10") +
    #geom_errorbar(
    #  aes(ymin = q25, ymax = q75),
    #  width = 0,
    #  alpha = .2,
    #  col = "grey30"
    #) +
    # geom_errorbar(
    #   aes(ymin = mad - se_robust, ymax = mad + se_robust),
    #   width = 0,
    #   alpha = .2
    # ) +

    #geom_hex() +
    #geom_rug()+
    geom_ribbon(
      data = newdat,
      aes(x = .data[[x_var]], ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      fill = col,
      alpha = 0.5
    ) +
    geom_line(
      data = newdat,
      aes(x = .data[[x_var]], y = fit),
      inherit.aes = FALSE,
      color = col,
      linewidth = 1.1
    ) +
    #annotate(
    #  "text",
    #  x = min(data[[x_var]], na.rm = TRUE),
    #  y = ylab.text, #max(data$Seeds_cur, na.rm = TRUE),
    #  label = lab,
    #  hjust = 0,
    #  vjust = 1,
    #  size = 2.5
    #) +
    #theme_minimal() +
    scale_fill_distiller(palette = "Blues", direction = -1) +
    labs(
      title = panel_title,
      x = x_lab,
      y = "Seed production"
    )
}

#predict glmmTMB
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


#predict specific value
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
