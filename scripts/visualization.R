#' Plot Tweedie model predictions for a climate cue
#'
#' Creates a cue-response plot from a fitted `glmmTMB` Tweedie model. The
#' function predicts seed production across the observed range of a climate cue,
#' back-converts the cue to its original temperature scale if the model used a
#' scaled predictor, and overlays the fitted regression line with a 95% confidence
#' ribbon.
#'
#' @param model A fitted `glmmTMB` model, typically with `family = tweedie()`.
#' @param data A data frame containing the original cue variable, the scaled cue
#'   variable, `Seeds_prev`, `Demand`, and `sitenewname`.
#' @param x_var Character. Name of the original, unscaled cue variable to plot on
#'   the x-axis.
#' @param x_scaled_var Character. Name of the scaled cue variable used in the
#'   fitted model.
#' @param panel_title Character. Plot title.
#' @param x_lab Character. Label for the x-axis.
#' @param ylab.text Numeric. y-axis position for optional coefficient annotation.
#'   Currently not used because the annotation code is commented out.
#' @param col Character. Colour used for the fitted line and confidence ribbon.
#'   Default is `"#E64B35FF"`.
#'
#' @return A `ggplot` object showing predicted seed production across the cue
#'   gradient with a fitted line and 95% confidence ribbon.
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
