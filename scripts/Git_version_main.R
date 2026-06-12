# Main Analysis Script for Polish Harvest Dataset
# Complete workflow for analyzing mast seeding trends
#system("git rm --cached scripts/main_analysis.R")
#system("git rm --cached scripts/main_analysisvFinal.R")
#system("git rm --cached scripts/main_analysis_may2026.R")
#system("git rm --cached functions.R")

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(GGally)
library(patchwork)
# Load required functions
source("scripts/data_utils.R")
source("scripts/trend_analysis.R")
source("scripts/visualization.R")
theme_set(cowplot::theme_cowplot())
#install.packages("showtext")
library(showtext)
library(stringi)
# Add Roboto from Google Fonts
font_add_google("Roboto", "roboto")
# Enable automatic use
showtext_auto()

############################################################################
#Data seed production harvest load
############################################################################
#updated file from Jessie Foest - Nature Climate Change 2026
polish_harvest <- data.table::fread(
  file = here::here("data", "Harvest_Poland.csv"),
  sep = ","
) %>%
  as_tibble()

# Prepare initial data
polish_harvest_df_simpl.init <- polish_harvest %>%
  # rm/correct wrong values
  group_by(Species, NADL, Year) %>%
  mutate(n = n()) %>%
  filter(
    !(n > 1 &
      is.na(Prev_year_harvest))
  ) %>%
  unique() %>%
  ungroup() %>%
  mutate(
    NADL = case_when(grepl("Koście", NADL) ~ "Kościerzyna", TRUE ~ NADL)
  ) %>%
  dplyr::select(
    Species,
    NADL,
    Year,
    Seeds_prev = Prev_year_harvest,
    Demand,
    Latitude = lat,
    Longitude = long,
    Seed_source_area
  ) %>%
  group_by(Species, NADL, Latitude, Longitude) %>%
  # When seed source area was recorded (from 2007 onwards), there should be at least some area from which was harvested.
  filter(!all(is.na(Seed_source_area) & Year > 2007)) %>%
  ungroup()

# Fix Oak issue: add "Oaks_both" from separate spp
oak_append <- polish_harvest_df_simpl.init %>%
  filter(Species %in% c("Pedunculate oak", "Sessile oak")) %>%
  group_by(NADL, Latitude, Longitude, Year) %>%
  summarise(
    Species = "Oaks_both",
    Seeds_prev = sum(Seeds_prev, na.rm = TRUE),
    Demand = sum(Demand, na.rm = TRUE),
    Seed_source_area = mean(Seed_source_area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  anti_join(
    polish_harvest_df_simpl.init,
    by = c("Species", "NADL", "Latitude", "Longitude", "Year")
  )

# Combine back into the main data
polish_harvest_df_simpl <- bind_rows(
  polish_harvest_df_simpl.init,
  oak_append
) %>%
  arrange(Species, NADL, Year)


# Continue processing (complete, and impplicit zeroes)
#did the same as Jessie
polish_harvest_df_simpl_compl <- polish_harvest_df_simpl %>%
  group_by(Species, NADL, Latitude, Longitude) %>%
  complete(Year = full_seq(c(Year, min(Year) - 1), 1)) %>%
  # make implicit explicit as per data owner instructions
  mutate(Seeds_prev = replace_na(Seeds_prev, 0)) %>%
  # ts with limited info (convergence)
  filter(sum(is.na(Seeds_prev) | Seeds_prev == 0) / n() < 0.9) %>%
  arrange(Species, NADL, Year) %>%
  # Seeds_prev was actually the 'realised' harvest column in this dataset, so we use that to get Seeds_cur.
  mutate(Seeds_cur = lead(Seeds_prev, 1)) %>%
  slice(-1) %>% # Remove first incomplete year after adding lead
  # Imputing demand
  mutate(
    Demand = case_when(
      #is.na(Demand) & Seeds_cur == 0 ~ 0,
      is.na(Demand) ~ mean(Demand, na.rm = TRUE),
      TRUE ~ Demand
    )
  ) %>%
  dplyr::select(-Seed_source_area) %>%
  ungroup() %>%
  filter(!is.na(Seeds_cur)) %>%
  distinct()


# Long time series
#ok here I have the last dataset from Jessie code
#the idea is to combine oaks together because they have different sample (aggregated together before 2008 pedunculate and sessile)
polish_decade <- polish_harvest_df_simpl_compl %>%
  group_by(Species, NADL) %>%
  filter(n_distinct(Year) > 10) %>%
  ungroup() %>%
  mutate(sitenewname = paste0(Species, "_", NADL))

#check proporition of 0
prop_zero_by_site <- polish_decade %>%
  group_by(Species, NADL) %>%
  summarise(
    n_obs = n(),
    n_zero_cur = sum(Seeds_cur == 0, na.rm = TRUE),
    prop_zero_cur = n_zero_cur / n_obs,
    .groups = "drop"
  )

#initial.data.harvest.load.kubav %>% filter(Year == 1988)
year_gaps <- polish_decade %>%
  group_by(NADL) %>%
  summarise(
    min_year = min(Year, na.rm = TRUE),
    max_year = max(Year, na.rm = TRUE),
    n_obs = n_distinct(Year),
    expected_n = max_year - min_year + 1,
    n_missing = expected_n - n_obs,
    continuous = n_missing == 0,
    .groups = "drop"
  )


#ok now I want those with less than 80 percent 0 and Oaks both (because the data contain also other oaks separately)
Extend.Harvest.Filled = polish_decade %>%
  left_join(
    prop_zero_by_site %>% dplyr::select(Species, NADL, prop_zero_cur),
    by = c("Species", "NADL")
  ) %>%
  filter(Species == "Oaks_both" & prop_zero_cur < 0.8) #replace here by 0.4 or 0.7, does not affect final figures


#check correlation between demand and current year seeds
plot(
  log1p(Extend.Harvest.Filled$Demand),
  log1p(Extend.Harvest.Filled$Seeds_cur)
)
cor(log1p(Extend.Harvest.Filled$Demand), log1p(Extend.Harvest.Filled$Seeds_cur))
#basic plot test
Extend.Harvest.Filled %>%
  filter(Species == "Oaks_both") %>%
  group_by(Year) %>%
  mutate(zero = sum(Demand == 0, na.rm = T), n = n(), prop_zero = zero / n) %>%
  ungroup() %>%
  ggplot(aes(Year, prop_zero)) +
  geom_point() +
  ylab("Porpotion 0 for demand (replacement for demand)")


############################################################################
#Data daily temperature from EOEBS
############################################################################
#extraction done by Jakub Szymkowiak in 2025, PNAS
initial.tmax.harvest.v1 <- readxl::read_excel(
  here::here("data", "tmax.xlsx")
) %>%
  #filter(Year >= (1987 - 1) & Year < 2023) %>%
  rename(
    #have to rename some because of polish letters
    Karniszewice = Karnieszewice,
    Kwidzyń = Kwidzyn,
    `Połczyn Zdrój` = Połczyn,
    `Międzyrzec Podlaski` = Międzyrzec,
    `Krosno Odrzańskie` = Krosno,
    `Buda stalowska` = `Buda Stalowska`
  )

############################################################################
#Get cues in oaks but not only
############################################################################
#format specific window
initial.tmax.harvest.aprmay <- make_tmax_window(
  initial.tmax.harvest.v1,
  months = c(4, 5),
  label = "AprMay"
)

initial.tmax.harvest.may <- make_tmax_window(
  initial.tmax.harvest.v1,
  months = c(5),
  label = "May"
)

initial.tmax.harvest.march <- make_tmax_window(
  initial.tmax.harvest.v1,
  months = c(3),
  label = "Mar"
)

initial.tmax.harvest.december <- make_tmax_window(
  initial.tmax.harvest.v1,
  months = c(12),
  label = "Dec"
) %>%
  mutate(Year = Year + 1) # Shift December to the following year
#(I add 1 to climate because I merge this file to seed prod later) for indexing
#(the new year climate for example 1951, is related to seed production 1951 and also
#because I added 1 year to the old year, make the new year 1951)

initial.tmax.harvest.JuneJulyAugust1 <- make_tmax_window(
  initial.tmax.harvest.v1,
  months = c(6, 7, 8),
  label = "JunJulAug1"
) %>%
  mutate(Year = Year + 1)

initial.tmax.harvest.JuneJuly1 <- make_tmax_window(
  initial.tmax.harvest.v1,
  months = c(6, 7),
  label = "JunJul1"
) %>%
  mutate(Year = Year + 1)

initial.tmax.harvest.JuneJulyAugust2 <- make_tmax_window(
  initial.tmax.harvest.v1,
  months = c(6, 7, 8),
  label = "JunJulAug2"
) %>%
  mutate(Year = Year + 2)

#combination
initial.tmax.harvest.full.year = initial.tmax.harvest.v1 %>%
  group_by(Year) %>%
  dplyr::select(-Date, -Day, -Month) %>%
  summarise_all(mean) %>%
  pivot_longer(-Year) %>%
  rename(NADL = name, tmax.fullyear.harvest = value)

############################################################################
#Combined dataset for oaks seed and cues
############################################################################
#focus on oaks first
Extend.Harvest.Filled %>%
  left_join(initial.tmax.harvest.aprmay) %>%
  left_join(initial.tmax.harvest.march) %>%
  left_join(initial.tmax.harvest.december) %>%
  left_join(initial.tmax.harvest.may) %>%
  left_join(initial.tmax.harvest.JuneJulyAugust1) %>%
  left_join(initial.tmax.harvest.JuneJulyAugust2) %>%
  left_join(initial.tmax.harvest.full.year) %>%
  group_by(sitenewname) %>%
  mutate(
    scaled.tmax.harvest_Mar = scale(tmax.harvest_Mar),
    scaled.tmax.harvest_AprMay = scale(tmax.harvest_AprMay),
    scaled.tmax.harvest_JunJulAug1 = scale(tmax.harvest_JunJulAug1) #
    #og_previous_seeds = log1p(Seeds_prev)
  ) %>%
  ungroup() %>%
  mutate(Year_stand = scale(Year)) -> Oaks_Harvest #%>% filter(Seeds_cur > 0)


ggplot(
  Oaks_Harvest,
  aes(x = Year, y = Seeds_cur, group = sitenewname, color = sitenewname)
) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Year",
    y = "Seed"
  ) +
  theme(legend.position = "none")
############################################################################
#Cues analysis
############################################################################
#because I have to standardised year, because some model did not converge
sd_year <- sd(Oaks_Harvest$Year, na.rm = TRUE)

#not reported finally
# cues.december = glmmTMB(
#   tmax.harvest_Dec ~ Year + (Year | sitenewname),
#   data = Oaks_Harvest
# )
# summary(cues.december)

#cues in march
cues.march = glmmTMB(
  tmax.harvest_Mar ~ Year_stand + (Year_stand | sitenewname),
  data = Oaks_Harvest
)
summary(cues.march)
slope_stand_march <- fixef(cues.march)$cond["Year_stand"]
slope_year_march <- slope_stand_march / sd_year
se_stand_march <- summary(cues.march)$coefficient$cond[
  "Year_stand",
  "Std. Error"
]
se_year_march <- se_stand_march / sd_year

plot(ggeffects::ggpredict(cues.march))
#8.32/6.65 #for the fold increase

#cues in April May
cues.aprmay = glmmTMB(
  tmax.harvest_AprMay ~ Year_stand + (Year_stand | sitenewname),
  data = Oaks_Harvest
)
summary(cues.aprmay)
plot(ggeffects::ggpredict(cues.aprmay))
#17.01/16.16#for the fold increase

slope_stand_aprmay <- fixef(cues.aprmay)$cond["Year_stand"]
se_stand_aprmay <- summary(cues.aprmay)$coefficient$cond[
  "Year_stand",
  "Std. Error"
]
slope_year_aprmay <- slope_stand_aprmay / sd_year
se_year_aprmay <- se_stand_aprmay / sd_year

plot(ggeffects::ggpredict(cues.aprmay))

#cues in JJA 1
cues.JJA1 = glmmTMB(
  tmax.harvest_JunJulAug1 ~ Year_stand + (Year_stand | sitenewname),
  data = Oaks_Harvest
)
summary(cues.JJA1)
plot(ggeffects::ggpredict(cues.JJA1))
#24.90 / 22.27#for the fold increase
slope_stand_JJA1 <- fixef(cues.JJA1)$cond["Year_stand"]
se_stand_JJA1 <- summary(cues.JJA1)$coefficient$cond["Year_stand", "Std. Error"]
slope_year_JJA1 <- slope_stand_JJA1 / sd_year
se_year_JJA1 <- se_stand_JJA1 / sd_year

#cues in JJA2, but not reported
cues.JJ2 = glmmTMB(
  tmax.harvest_JunJulAug2 ~ Year + (1 | sitenewname),
  data = Oaks_Harvest
)
summary(cues.JJA1)
summary(cues.JJ2)

############################################################################
#Seed direct count analysis
############################################################################

#March
model.tweedie.year.march <- glmmTMB(
  Seeds_cur ~
    scaled.tmax.harvest_Mar +
      log1p(Seeds_prev) +
      log1p(Demand) +
      (1 | sitenewname), #Seeds_prev  + (1 | sitenewname),
  data = Oaks_Harvest,
  family = tweedie() #default log link
)
summary(model.tweedie.year.march)

#April May
model.tweedie.year.aprmay <- glmmTMB(
  Seeds_cur ~
    scaled.tmax.harvest_AprMay +
      log1p(Seeds_prev) +
      log1p(Demand) +
      (1 | sitenewname),
  data = Oaks_Harvest,
  family = tweedie()
)
summary(model.tweedie.year.aprmay)

#JJA1
model.tweedie.year.JJA1 <- glmmTMB(
  Seeds_cur ~
    scaled.tmax.harvest_JunJulAug1 +
      log1p(Seeds_prev) +
      log1p(Demand) +
      (1 | sitenewname),
  data = Oaks_Harvest,
  family = tweedie()
)
summary(model.tweedie.year.JJA1)


model.tweedie.year.JJ2 <- glmmTMB(
  Seeds_cur ~
    tmax.harvest_JunJulAug2 +
      log1p(Seeds_prev) +
      log1p(Demand) +
      (1 | sitenewname),
  data = Oaks_Harvest,
  family = tweedie()
)
summary(model.tweedie.year.JJ2)

#model check
#take some time
#check the qqplot
#performance::model_performance(model.tweedie.year.march)
#performance::check_model(model.tweedie.year.march)
#performance::model_performance(model.tweedie.year.aprmay)
#performance::check_model(model.tweedie.year.aprmay)
#performance::model_performance(model.tweedie.year.JJA1)
#performance::check_model(model.tweedie.year.JJA1)

#AIC model
AIC(
  model.tweedie.year.march,
  model.tweedie.year.aprmay,
  model.tweedie.year.JJA1,
  model.tweedie.year.JJ2
)

#ok we decided to report one model per cue,
#because later analysis with rolling cues, the correlation where too high
#see the ggpairs plot below

#but we still tried model with 3 cues together
#first check correlation between cues. It is not thaht high
#but as said, once using rolling cues it was too high
#so we wanted to keep method consistent
cor(
  Oaks_Harvest[, c(
    "tmax.harvest_Dec",
    "tmax.harvest_Mar",
    "tmax.harvest_AprMay"
  )],
  use = "complete.obs"
)
model.tweedie.year.3cues = glmmTMB(
  Seeds_cur ~
    tmax.harvest_Mar +
      tmax.harvest_Mar +
      tmax.harvest_AprMay +
      log1p(Seeds_prev) +
      log1p(Demand) +
      (1 | sitenewname),
  data = Oaks_Harvest,
  family = tweedie()
)

summary(model.tweedie.year.3cues)


#model prediction
value510JJA = plot_tweedie_cue_data(
  model = model.tweedie.year.JJA1,
  data = Oaks_Harvest,
  x_var = "tmax.harvest_JunJulAug1",
  x_scaled_var = "scaled.tmax.harvest_JunJulAug1"
)
value510 = plot_tweedie_cue_data(
  model = model.tweedie.year.march,
  data = Oaks_Harvest,
  x_var = "tmax.harvest_Mar",
  x_scaled_var = "scaled.tmax.harvest_Mar"
)
#5 degree = 2245.592
#10 degree = 1626.629

#finally we reported coef plot
# p_JJA1 <- plot_tweedie_cue(
#   model = model.tweedie.year.JJA1,
#   data = Oaks_Harvest,
#   x_var = "tmax.harvest_JunJulAug1",
#   x_scaled_var = "scaled.tmax.harvest_JunJulAug1",
#   panel_title = "Previous summer (JJA1)",
#   x_lab = "Temperature (°C)",
#   ylab.text = 4000
# ) +
#   hrbrthemes::theme_ipsum_rc(
#     base_size = 14,
#     axis_title_size = 16,
#     plot_title_face = NULL,
#     base_family = "roboto"
#   )
# p_march <- plot_tweedie_cue(
#   model = model.tweedie.year.march,
#   data = Oaks_Harvest,
#   x_var = "tmax.harvest_Mar",
#   x_scaled_var = "scaled.tmax.harvest_Mar",
#   panel_title = "Early spring (Mar)",
#   x_lab = "Temperature (°C)",
#   ylab.text = 4000
# ) +
#   hrbrthemes::theme_ipsum_rc(
#     base_size = 14,
#     axis_title_size = 16,
#     plot_title_face = NULL,
#     base_family = "roboto"
#   )
# p_aprmay <- plot_tweedie_cue(
#   model = model.tweedie.year.aprmay,
#   data = Oaks_Harvest,
#   x_var = "tmax.harvest_AprMay",
#   x_scaled_var = "scaled.tmax.harvest_AprMay",
#   panel_title = "Late spring (Apr–May)",
#   x_lab = "Temperature (°C)"
# ) +
#   hrbrthemes::theme_ipsum_rc(
#     base_size = 14,
#     axis_title_size = 16,
#     plot_title_face = NULL,
#     base_family = "roboto"
#   )
# Cues.seeds = (p_JJA1 +
#   theme(plot.margin = margin(5, 5, 5, 5)) +
#   ylim(0, 4100) |
#   p_march +
#     ylim(0, 4100) +
#     (theme(
#       axis.text.y = element_blank(),
#       axis.title.y = element_blank(),
#       plot.margin = margin(5, 5, 5, 5)
#     )) |
#   p_aprmay +
#     ylim(0, 4100) +
#     (theme(
#       axis.text.y = element_blank(),
#       axis.title.y = element_blank(),
#       plot.margin = margin(5, 5, 5, 5)
#     )))
#
# cowplot::save_plot(
#   "figures/cues.pdf",
#   Cues.seeds,
#   nrow = 1.2,
#   ncol = 1.8
# )

############################################################################
#Rolling window of seed production
############################################################################
rolling_seeds <- Oaks_Harvest %>%
  group_by(sitenewname) %>%
  group_modify(
    ~ compute_rolling_metrics(.x, Seeds_cur, window = 10, step = 5)
  ) %>%
  ungroup() %>%
  rename_with(
    ~ paste0(.x, "_seeds"),
    c(mean_value, sd_value, CVp, q25, q90, p_zero, p_nonzero)
  )

# ggplot(rolling_seeds, aes(CVp_seeds, mean_value_seeds)) + geom_point()
# ggplot(rolling_seeds, aes(CVp_seeds, p_nonzero_seeds)) +
#   geom_point()
# ggplot(rolling_seeds, aes(CVp_seeds, p_zero_seeds)) +
#   geom_point()

rolling_demand <- Oaks_Harvest %>%
  group_by(sitenewname) %>%
  group_modify(~ compute_rolling_metrics(.x, Demand, window = 10, step = 5)) %>%
  ungroup() %>%
  rename_with(
    ~ paste0(.x, "_demand"),
    c(mean_value, sd_value, CVp, q25, q90, p_zero, p_nonzero)
  )

rolling_seeds_lag <- Oaks_Harvest %>%
  group_by(sitenewname) %>%
  group_modify(
    ~ compute_rolling_metrics(.x, Seeds_prev, window = 10, step = 5)
  ) %>%
  ungroup() %>%
  rename_with(
    ~ paste0(.x, "_lag1_seeds"),
    c(mean_value, sd_value, CVp, q25, q90, p_zero, p_nonzero)
  )

#here I am rolling cues, but i am interested by the average cue value (matching the seed value obs)
rolling_cues_AprMay <- Oaks_Harvest %>%
  group_by(sitenewname) %>%
  group_modify(
    ~ compute_rolling_metrics(.x, tmax.harvest_AprMay, window = 10, step = 5)
  ) %>%
  ungroup() %>%
  rename_with(
    ~ paste0(.x, "_T_AprMay"),
    c(mean_value, sd_value, CVp, q25, q90)
  )

rolling_cues_Dec <- Oaks_Harvest %>%
  group_by(sitenewname) %>%
  group_modify(
    ~ compute_rolling_metrics(.x, tmax.harvest_Dec, window = 10, step = 5)
  ) %>%
  ungroup() %>%
  rename_with(
    ~ paste0(.x, "_T_Dec"),
    c(mean_value, sd_value, CVp)
  )

rolling_cues_Mar <- Oaks_Harvest %>%
  group_by(sitenewname) %>%
  group_modify(
    ~ compute_rolling_metrics(.x, tmax.harvest_Mar, window = 10, step = 5)
  ) %>%
  ungroup() %>%
  rename_with(
    ~ paste0(.x, "_T_Mar"),
    c(mean_value, sd_value, CVp)
  )

rolling_cues_JunJulAug1 <- Oaks_Harvest %>%
  group_by(sitenewname) %>%
  group_modify(
    ~ compute_rolling_metrics(
      .x,
      tmax.harvest_JunJulAug1,
      window = 10,
      step = 5
    )
  ) %>%
  ungroup() %>%
  rename_with(
    ~ paste0(.x, "_T_JunJulAug1"),
    c(mean_value, sd_value, CVp)
  )

rolling_tmax_fullyear <- Oaks_Harvest %>%
  group_by(sitenewname) %>%
  group_modify(
    ~ compute_rolling_metrics(.x, tmax.fullyear.harvest, window = 10, step = 5)
  ) %>%
  dplyr::select(-p_zero, -p_nonzero, -q25, -q90) %>%
  ungroup() %>%
  rename_with(~ paste0(.x, "_Tfull"), c(mean_value, sd_value, CVp))


rolling_metrics_oaks <- rolling_seeds %>%
  left_join(
    rolling_demand %>%
      dplyr::select(
        sitenewname,
        Year,
        mean_value_demand,
        sd_value_demand,
        CVp_demand,
        p_zero_demand,
        Year_center
      )
  ) %>%
  left_join(
    rolling_seeds_lag %>%
      dplyr::select(
        sitenewname,
        Year,
        mean_value_lag1_seeds,
        CVp_lag1_seeds,
        p_zero_lag1_seeds,
        Year_center
      )
  ) %>%
  left_join(
    rolling_cues_Mar %>%
      dplyr::select(
        sitenewname,
        Year,
        mean_value_T_Mar,
        sd_value_T_Mar,
        Year_center
      )
  ) %>%
  left_join(
    rolling_cues_AprMay %>%
      dplyr::select(
        sitenewname,
        Year,
        mean_value_T_AprMay,
        sd_value_T_AprMay,
        Year_center
      )
  ) %>%
  left_join(
    rolling_cues_Dec %>%
      dplyr::select(
        sitenewname,
        Year,
        mean_value_T_Dec,
        Year_center
      )
  ) %>%
  left_join(
    rolling_cues_JunJulAug1 %>%
      dplyr::select(
        sitenewname,
        Year,
        mean_value_T_JunJulAug1,
        sd_value_T_JunJulAug1,
        Year_center
      )
  ) %>%
  left_join(
    rolling_tmax_fullyear %>%
      dplyr::select(
        sitenewname,
        Year,
        mean_value_Tfull,
        sd_value_Tfull,
        Year_center
      )
  ) %>%
  mutate(
    log_CVp_seeds = log(CVp_seeds),
    log_mean_seeds = log1p(mean_value_seeds),
    logit_prop_zero_seeds = car::logit(p_zero_seeds),
    log_mean_value_demand = log1p(mean_value_demand)
  )


zero_year <- rolling_metrics_oaks %>%
  group_by(Year) %>%
  summarise(
    prop_zero = mean(Seeds_cur == 0, na.rm = TRUE),
    .groups = "drop"
  )

############################################################################
#Ananlysis olling window of seed production
############################################################################

rolling_metrics_oaks %>%
  select(mean_value_T_Mar, mean_value_T_AprMay) %>%
  ggpairs(
    lower = list(
      continuous = wrap("smooth", method = "lm", se = FALSE, alpha = 0.4)
    ),
    upper = list(
      continuous = wrap("cor", size = 4)
    ),
    diag = list(
      continuous = wrap("densityDiag")
    )
  ) +
  theme_minimal()


ggplot(zero_year, aes(Year, prop_zero)) +
  geom_point() +
  geom_vline(xintercept = 2008, linetype = "dashed") +
  theme_minimal() +
  labs(y = "Proportion of zero seed years")
rolling_metrics_oaks %>%
  mutate(period = if_else(Year < 2008, "pre", "post")) %>%
  filter(Seeds_cur > 0) %>%
  ggplot(aes(log1p(Seeds_cur), fill = period)) +
  geom_density(alpha = 0.4) +
  theme_minimal()

############################################################################
#rolling cues trends
############################################################################

rolling_metrics_oaks <- rolling_metrics_oaks %>%
  mutate(
    Year_bin = cut(
      Year_center,
      breaks = seq(
        floor(min(Year_center, na.rm = TRUE) / 5) * 5,
        ceiling(max(Year_center, na.rm = TRUE) / 5) * 5,
        by = 5
      ),
      include.lowest = TRUE
    )
  ) %>%
  mutate(Year_c = Year_center - mean(Year_center, na.rm = TRUE))

plot(rolling_metrics_oaks$Year_c, rolling_metrics_oaks$Year_center)

#cjeck corelation
cor(
  rolling_metrics_oaks$mean_value_demand,
  rolling_metrics_oaks$mean_value_seeds,
  use = "complete.obs"
)

cor(
  rolling_metrics_oaks$CVp_seeds,
  rolling_metrics_oaks$CVp_demand,
  use = "complete.obs"
)

#I checked JJF paper, and she tool into account the deamnde as offset or covaraite, but
#more seed you have more demand you will have, so I will just check correlation first
#but I am not sure I would include both of them
#because I just increase the risk to inflate my stat models
#also the trends where too different between sites
#so better to include random slope and itnercept

#random slope for each sites
Temporal.Trend.CVp.model_rs <- glmmTMB(
  log_CVp_seeds ~
    Year_c +
      log1p(mean_value_demand) +
      (Year_c | sitenewname),
  data = rolling_metrics_oaks
)
summary(Temporal.Trend.CVp.model_rs)
ranef_slopes <- ranef(Temporal.Trend.CVp.model_rs)$cond$sitenewname
head(ranef_slopes)

#for the plot option below
#new option
breaks5 <- seq(
  floor(min(rolling_metrics_oaks$Year_center, na.rm = TRUE)),
  ceiling(max(rolling_metrics_oaks$Year_center, na.rm = TRUE)) + 5,
  by = 5
)

bin_df <- rolling_metrics_oaks %>%
  mutate(
    Year_bin = cut(
      Year_center,
      breaks = breaks5,
      include.lowest = TRUE,
      right = FALSE
    )
  )

bin_counts <- bin_df %>%
  group_by(Year_bin) %>%
  summarise(
    x = mean(Year_center, na.rm = TRUE), #Year_center
    n = sum(!is.na(CVp_seeds)),
    .groups = "drop"
  )


#make the main figures
newdat <- data.frame(
  Year_c = seq(
    min(rolling_metrics_oaks$Year_c),
    max(rolling_metrics_oaks$Year_c),
    length.out = 100
  ),
  mean_value_demand = mean(
    rolling_metrics_oaks$mean_value_demand,
    na.rm = TRUE
  ),
  sitenewname = NA # important for population-level prediction
)

pred <- predict(
  Temporal.Trend.CVp.model_rs,
  newdata = newdat,
  re.form = NA, # remove random effects
  se.fit = TRUE
)

newdat$fit_log <- pred$fit
newdat$se_log <- pred$se.fit

# 95% CI
newdat$lwr_log <- newdat$fit_log - 1.96 * newdat$se_log
newdat$upr_log <- newdat$fit_log + 1.96 * newdat$se_log

newdat <- newdat %>%
  mutate(
    fit_CV = exp(fit_log),
    lwr_CV = exp(lwr_log),
    upr_CV = exp(upr_log),
    Year_center = Year_c + mean(rolling_metrics_oaks$Year_center, na.rm = TRUE)
  )

checking.bins = bin_df %>%
  dplyr::select(sitenewname, Year_bin, Year_center) %>%
  group_by(sitenewname, Year_bin) %>%
  tally()
View(checking.bins) #should be 1 for all

#now the plto
boxplot.cv.trend <- ggplot(bin_df, aes(x = Year_bin, y = CVp_seeds)) +
  gghalves::geom_half_violin(
    aes(group = Year_bin),
    side = "l",
    fill = "#E64B35FF",
    color = NA,
    alpha = 0.2
  ) +
  gghalves::geom_half_boxplot(
    aes(group = Year_bin),
    side = "r",
    alpha = .1,
    lwd = 0.5,
    fill = "#F39B7FFF",
    errorbar.length = 0.4,
    color = "grey20",
    outlier.shape = NA
  ) +
  geom_text(
    data = bin_counts,
    aes(x = Year_bin, y = 3.2, label = paste0("n=", n)),
    size = 4,
    family = "sans",
    col = "black"
  ) +
  ylim(0., 3.2) +
  #coord_cartesian(ylim = c(0.15, 3.3)) +
  labs(x = "Year", y = "CVp") +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  ) +
  scale_x_discrete(
    breaks = c(bin_counts$Year_bin),
    label = c(1995, 2000, 2005, 2010, 2015, 2020)
  )

boxplot.cv.trend

cowplot::save_plot(
  "figures/cv.rolling.decline.boxplot.pdf",
  boxplot.cv.trend,
  nrow = 1.1,
  ncol = 1
)

#hisrogram of the slopes
#finnaly not reproted
# fixed_slope <- fixef(Temporal.Trend.CVp.model_rs)$cond["Year_c"]
# site_slopes <- fixed_slope + ranef_slopes$Year_c
# site_slopes_df <- data.frame(
#   sitenewname = rownames(ranef_slopes),
#   slope = site_slopes
# ) %>%
#   mutate(
#     percent_change = (exp(slope) - 1) * 100,
#     backtransf_tempslope = exp(slope),
#     group_tempslope = case_when(slope <= 0 ~ "<=0", TRUE ~ ">0")
#   )
#
# mean(site_slopes_df$percent_change, na.rm = TRUE)
#
# VarCorr(Temporal.Trend.CVp.model_rs)
# col.density = ggplot(site_slopes_df, aes(x = slope, fill = group_tempslope)) +
#   geom_histogram(binwidth = 0.01, color = "black", alpha = 0.7) +
#   theme_minimal() +
#   labs(
#     x = "Slope (log scale)",
#     y = "Count of Sites"
#   ) +
#   theme(legend.title = element_blank(), legend.position = c(.8, .8))
#
# col.density
#
# newdat_pop <- data.frame(
#   Year_c = seq(
#     min(rolling_metrics_oaks$Year_c, na.rm = TRUE),
#     max(rolling_metrics_oaks$Year_c, na.rm = TRUE),
#     length.out = 100
#   ),
#   mean_value_demand = mean(
#     rolling_metrics_oaks$mean_value_demand,
#     na.rm = TRUE
#   ),
#   sitenewname = NA
# )
#
# pred_pop <- predict(
#   Temporal.Trend.CVp.model_rs,
#   newdata = newdat_pop,
#   re.form = NA,
#   se.fit = TRUE
# )
#
# year_mean <- mean(rolling_metrics_oaks$Year_center, na.rm = TRUE)
#
# newdat_pop <- newdat_pop %>%
#   mutate(
#     fit_log = pred_pop$fit,
#     se_log = pred_pop$se.fit,
#     lwr_log = fit_log - 1.96 * se_log,
#     upr_log = fit_log + 1.96 * se_log,
#     fit_CV = exp(fit_log),
#     lwr_CV = exp(lwr_log),
#     upr_CV = exp(upr_log),
#     Year_center = Year_c + year_mean
#   )
#
#
# newdat_site <- expand.grid(
#   Year_c = seq(
#     min(rolling_metrics_oaks$Year_c, na.rm = TRUE),
#     max(rolling_metrics_oaks$Year_c, na.rm = TRUE),
#     length.out = 100
#   ),
#   sitenewname = unique(rolling_metrics_oaks$sitenewname)
# )
#
# newdat_site$mean_value_demand <- mean(
#   rolling_metrics_oaks$mean_value_demand,
#   na.rm = TRUE
# )
#
# newdat_site$fit_log <- predict(
#   Temporal.Trend.CVp.model_rs,
#   newdata = newdat_site,
#   re.form = NULL # include random effects
# )
#
# newdat_site <- newdat_site %>%
#   mutate(
#     fit_CV = exp(fit_log),
#     Year_center = Year_c + year_mean
#   )
#
#
# decline.pop = ggplot() +
#   geom_point(
#     data = rolling_metrics_oaks,
#     aes(x = Year_center, y = CVp_seeds),
#     alpha = 0.12,
#     size = 0.1
#   ) +
#   geom_line(
#     data = newdat_site %>% left_join(site_slopes_df),
#     aes(
#       x = Year_center,
#       y = fit_CV,
#       group = sitenewname,
#       col = group_tempslope
#     ),
#     alpha = 0.25,
#     linewidth = 0.4
#   ) +
#   # global ribbon
#   geom_ribbon(
#     data = newdat_pop,
#     aes(x = Year_center, ymin = lwr_CV, ymax = upr_CV),
#     fill = "grey20",
#     alpha = 0.2
#   ) +
#   # global line
#   geom_line(
#     data = newdat_pop,
#     aes(x = Year_center, y = fit_CV),
#     color = "black",
#     linewidth = 1.2
#   ) +
#   labs(
#     x = "Year",
#     y = "Rolling CV of seed production"
#   ) +
#   theme(legend.position = "none")
#
#
# col.density <- col.density +
#   labs(x = "Slope", y = "Density") +
#   theme(
#     axis.title = element_text(size = 8),
#     axis.text = element_text(size = 7)
#   )
# cowplot::save_plot(
#   "figures/cv.rolling.decline.pdf",
#   decline.pop +
#     patchwork::inset_element(
#       col.density + theme(legend.position = "none"),
#       left = 0.62,
#       bottom = 0.55,
#       right = 0.98,
#       top = 0.98
#     ),
#   nrow = 1,
#   ncol = .8
# )

############################################################################
#SRolling cues of seed in relation to cues
############################################################################
Temporal.Trend.CVp.model.cues1 <- glmmTMB::glmmTMB(
  log(CVp_seeds) ~
    mean_value_T_Dec +
      (1 | sitenewname) +
      log_mean_value_demand,
  rolling_metrics_oaks
)
summary(Temporal.Trend.CVp.model.cues1)
insight::find_transformation(Temporal.Trend.CVp.model.cues1)
#need this, just to make sure later ggpredict would adjust transofrmation
Temporal.Trend.CVp.model.cues2 <- glmmTMB::glmmTMB(
  log(CVp_seeds) ~
    mean_value_T_Mar +
      (1 | sitenewname) +
      log_mean_value_demand,
  rolling_metrics_oaks
)
summary(Temporal.Trend.CVp.model.cues2)

Temporal.Trend.CVp.model.cues3 <- glmmTMB::glmmTMB(
  log(CVp_seeds) ~
    mean_value_T_AprMay +
      (1 | sitenewname) +
      log_mean_value_demand,
  rolling_metrics_oaks
)
summary(Temporal.Trend.CVp.model.cues3)

Temporal.Trend.CVp.model.cues4 <- glmmTMB::glmmTMB(
  log(CVp_seeds) ~
    mean_value_T_JunJulAug1 +
      (1 | sitenewname) +
      log_mean_value_demand,
  rolling_metrics_oaks
)
summary(Temporal.Trend.CVp.model.cues4)

AIC(
  Temporal.Trend.CVp.model.cues1,
  Temporal.Trend.CVp.model.cues2,
  Temporal.Trend.CVp.model.cues3,
  Temporal.Trend.CVp.model.cues4
)

Temporal.Trend.CVp.model.all.cues <- glmmTMB::glmmTMB(
  log_CVp_seeds ~
    mean_value_T_Dec +
      mean_value_T_Mar +
      mean_value_T_AprMay +
      (1 | sitenewname) +
      log_mean_value_demand,
  rolling_metrics_oaks
)
summary(Temporal.Trend.CVp.model.all.cues)

#Figure 3
pred_JJ1 <- predict_glmmtmb_curve(
  model = Temporal.Trend.CVp.model.cues4,
  data = rolling_metrics_oaks,
  x_var = "mean_value_T_JunJulAug1",
  demand_var = "log_mean_value_demand",
  site_var = "sitenewname",
  backtransform = "exp"
)

JJ1.rolCV = ggplot() +
  geom_point(
    data = pred_JJ1$part_resid,
    aes(x = x, y = predicted),
    alpha = 0.1,
    col = "black",
    size = .3
  ) +
  geom_line(
    data = pred_JJ1$newdat,
    aes(x = mean_value_T_JunJulAug1, y = fit),
    color = "#E64B35FF",
    linetype = "dashed",
    linewidth = .8
  ) +
  scale_fill_viridis_c() +
  labs(
    x = "Temperature in JJAT1",
    y = "CVp (partial residuals)"
  ) +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  )

pred_Mar <- predict_glmmtmb_curve(
  model = Temporal.Trend.CVp.model.cues2,
  data = rolling_metrics_oaks,
  x_var = "mean_value_T_Mar",
  demand_var = "log_mean_value_demand",
  site_var = "sitenewname",
  backtransform = "exp"
)
Mar.rolCV = ggplot() +
  geom_point(
    data = pred_Mar$part_resid,
    aes(x = x, y = predicted),
    alpha = 0.1,
    col = "black",
    size = .3
  ) +
  geom_ribbon(
    data = pred_Mar$newdat,
    aes(x = mean_value_T_Mar, ymin = lwr, ymax = upr),
    fill = "#E64B35FF",
    alpha = 0.2
  ) +
  geom_line(
    data = pred_Mar$newdat,
    aes(x = mean_value_T_Mar, y = fit),
    color = "#E64B35FF",
    linewidth = .8
  ) +
  scale_fill_viridis_c() +
  labs(
    x = "Temperature in March",
    y = "Rolling CV of seed production"
  ) +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  )

pred_AprMay <- predict_glmmtmb_curve(
  model = Temporal.Trend.CVp.model.cues3,
  data = rolling_metrics_oaks,
  x_var = "mean_value_T_AprMay",
  demand_var = "log_mean_value_demand",
  site_var = "sitenewname",
  backtransform = "exp"
)
AprMay.rolCV = ggplot() +
  geom_point(
    data = pred_AprMay$part_resid,
    aes(x = x, y = predicted),
    alpha = 0.1,
    col = "black",
    size = .3
  ) +
  geom_line(
    data = pred_AprMay$newdat,
    aes(x = mean_value_T_AprMay, y = fit),
    color = "#E64B35FF",
    linewidth = .8,
    linetype = "dashed"
  ) +
  scale_fill_viridis_c() +
  labs(
    x = "Temperature in April-May",
    y = "Rolling CV of seed production"
  ) +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  )
#good now
#those figure will ocme back below

############################################################################
#Check synchronys as well
############################################################################
#check for correlation
# sync.oaks.rolling = calculate.sync(Oaks_Harvest, n.years = 10, step = 5) %>%
#   mutate(Year_c = end.year - mean(end.year, na.rm = TRUE))
#
# sync.oaks.rolling %>%
#   ggplot(aes(x = end.year, y = mean.synch)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_classic()
#
#
# Temporal.Trend.Synchrony <- glmmTMB::glmmTMB(
#   mean.synch ~
#     end.year +
#       (1 | NADL),
#   family = glmmTMB::beta_family(link = "logit"),
#   sync.oaks.rolling
# )
# summary(Temporal.Trend.Synchrony)
# #did the same as well here
# Temporal.Trend.Synchrony.model_rs <- glmmTMB(
#   mean.synch ~
#     Year_c +
#       (Year_c | NADL),
#   family = glmmTMB::beta_family(link = "logit"),
#   data = sync.oaks.rolling
# )
# summary(Temporal.Trend.Synchrony.model_rs)
# ranef_slopes_sync <- ranef(Temporal.Trend.Synchrony.model_rs)$cond$NADL
# fixed_slope_sync <- fixef(Temporal.Trend.Synchrony.model_rs)$cond["Year_c"]
# site_slopes_sync <- fixed_slope_sync + ranef_slopes_sync$Year_c
# site_slopes_df_sync <- data.frame(
#   sitenewname = rownames(ranef_slopes_sync),
#   slope = site_slopes_sync
# ) %>%
#   mutate(
#     group_tempslope = case_when(slope <= 0 ~ "<=0", TRUE ~ ">0")
#   )
# ggplot(site_slopes_df_sync, aes(x = slope, fill = group_tempslope)) +
#   geom_histogram(binwidth = 0.001, color = "black", alpha = 0.7) +
#   theme_minimal() +
#   labs(
#     title = "Distribution of Site-Specific Slopes \nfor Year Effect on sync (logit)",
#     x = "Slope (logit scale)",
#     y = "Count of Sites"
#   )

############################################################################
#Check other step sensitvity only for oaks here
############################################################################
#that takes some time
library(purrr)
out10 <- fit_cue_models_over_window(
  data = Oaks_Harvest,
  window = 10,
  step = 5,
  #cues = c("tmax.harvest_Dec", "tmax.harvest_Mar", "tmax.harvest_AprMay"),
  cues = c(
    "tmax.harvest_JunJulAug1",
    "tmax.harvest_Mar",
    "tmax.harvest_AprMay"
  )
)

out10$cue_results


wins <- seq(6, 20, by = 2)

res_all <- map_dfr(wins, function(w) {
  fit_cue_models_over_window(Oaks_Harvest, window = w, step = 5)$cue_results
})

res_all

ggplot(res_all, aes(x = window, y = estimate, color = cue)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(y = "Cue effect on log(CV)", x = "Window length (years)")


step <- seq(1:10)
res_all_step <- map_dfr(step, function(s) {
  fit_cue_models_over_window(
    Oaks_Harvest,
    window = 10,
    step = s,
    cues = c(
      "tmax.harvest_JunJulAug1",
      "tmax.harvest_Mar",
      "tmax.harvest_AprMay"
    )
  )$cue_results
})


ggplot(res_all_step, aes(x = step, y = estimate, color = cue)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(y = "Cue effect on log(CV)", x = "Steps (years)")

ggplot(
  res_all_step,
  aes(x = step, y = estimate, shape = cue_col, col = log10(p))
) +
  geom_point() +
  facet_grid(. ~ cue_col) +
  scale_color_viridis_c()

res_all_step %>%
  ggplot(aes(x = estimate, color = cue)) +
  geom_density() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  )

##################################
#trend sensivity slope
step <- seq(1:10)
res_all_step_trends_LMM <- map_dfr(step, function(s) {
  fit_Trends_models_over_window(
    Oaks_Harvest,
    window = 10,
    step = s
  )$cue_results
})

table_latex = res_all_step_trends_LMM %>%
  mutate(
    term = recode(
      term,
      "(Intercept)" = "Intercept",
      "Year_c" = "Year",
      "log1p(mean_value_demand)" = "Demand"
    ),
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    p.value = signif(p.value, 3),
    AIC = round(AIC, 1)
  ) %>%
  dplyr::select(step, term, estimate, std.error, p.value, AIC, nobs)

table_latex %>%
  kableExtra::kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    caption = "Mixed-effects model estimates across rolling-window and step-size combinations",
    align = "c l c c c c c",
    col.names = c(
      "Step",
      "Term",
      "Estimate",
      "SE",
      "$p$",
      "AIC",
      "nobs"
    ),
    escape = FALSE
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("hold_position", "repeat_header"),
    font_size = 10
  )

############################################################################
#OK so here going back to cue sensivity plot
############################################################################
#check cue trends
#basiclaly I am doing the coefficient plot instead of plotign regression
#i combined those coefficient value to the beech models
pred_march <- ggeffects::ggpredict(cues.march, terms = "Year_stand [all]") %>%
  as.data.frame() %>%
  mutate(cue = "March")

pred_aprmay <- ggeffects::ggpredict(cues.aprmay, terms = "Year_stand [all]") %>%
  as.data.frame() %>%
  mutate(cue = "AprMay")

pred_JJ1 <- ggeffects::ggpredict(cues.JJA1, terms = "Year_stand [all]") %>%
  as.data.frame() %>%
  mutate(cue = "JunJulAug1")

# pred_JJ2 <- ggeffects::ggpredict(cues.JJ2, terms = "Year [all]") %>%
#   as.data.frame() %>%
#   mutate(cue = "JunJul2")

pred_all <- bind_rows(pred_march, pred_aprmay, pred_JJ1)

# extract fixed-effect slopes (direction)
slope_march <- fixef(cues.march)$cond["Year"]
slope_aprmay <- fixef(cues.aprmay)$cond["Year"]
slope_JJ1 <- fixef(cues.JJA1)$cond["Year"]
#slope_JJ2 <- fixef(cues.JJ2)$cond["Year"]

slope_df <- tibble(
  cue = c("March", "AprMay", "JunJulAug1"),
  slope = c(
    as.numeric(slope_march),
    as.numeric(slope_aprmay),
    as.numeric(slope_JJ1)
  ),
  slope_dir = if_else(slope > 0, "positive", "negative")
)

pred_all <- pred_all %>%
  left_join(slope_df, by = "cue") %>%
  mutate(
    Year = x *
      attr(Oaks_Harvest$Year_stand, "scaled:scale") +
      attr(Oaks_Harvest$Year_stand, "scaled:center")
  )


cue.JJ1 = ggplot(
  pred_all %>% filter(cue == "JunJulAug1"),
  aes(x = Year, y = predicted)
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2,
    fill = "red"
  ) +
  geom_line(linewidth = 1.2, col = "red") +
  labs(
    x = "Year",
    y = "JJA1 cue"
  )

cue.march = ggplot(
  pred_all %>% filter(cue == "March"),
  aes(x = Year, y = predicted)
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2,
    fill = "red"
  ) +
  geom_line(linewidth = 1.2, col = "red") +
  labs(
    x = "Year",
    y = "March cue"
  )

cue.AprMay = ggplot(
  pred_all %>% filter(cue == "AprMay"),
  aes(x = Year, y = predicted)
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2,
    fill = "red"
  ) +
  geom_line(linewidth = 1.2, col = "red") +
  labs(
    x = "Year",
    y = "April-May cue"
  )

cue.JJ1 + cue.march + cue.AprMay

###################################
#now the same with beech

Extend.Harvest.Filled.Beech = polish_decade %>%
  left_join(
    prop_zero_by_site %>% dplyr::select(Species, NADL, prop_zero_cur),
    by = c("Species", "NADL")
  ) %>%
  filter(Species == "Beech" & prop_zero_cur < 0.8)

Extend.Harvest.Filled.Beech %>%
  left_join(initial.tmax.harvest.aprmay) %>%
  left_join(initial.tmax.harvest.march) %>%
  left_join(initial.tmax.harvest.december) %>%
  left_join(initial.tmax.harvest.may) %>%
  left_join(initial.tmax.harvest.JuneJuly1) %>%
  left_join(initial.tmax.harvest.full.year) %>%
  group_by(sitenewname) %>%
  mutate(
    scaled.tmax.harvest_Mar = scale(tmax.harvest_Mar),
    scaled.tmax.harvest_AprMay = scale(tmax.harvest_AprMay),
    scaled.tmax.harvest_JunJul1 = scale(tmax.harvest_JunJul1)
  ) %>%
  ungroup() %>%
  mutate(Year_stand = scale(Year)) -> Beech_Harvest


cues.JJA.beech = glmmTMB(
  tmax.harvest_JunJul1 ~ Year_stand + (Year_stand | sitenewname),
  data = Beech_Harvest
)
summary(cues.JJA.beech)

ggeffects::ggpredict(
  cues.JJA.beech,
  terms = "Year_stand [all]"
)
24.31 / 21.73

beech.cue.JJA1 = ggeffects::ggpredict(
  cues.JJA.beech,
  terms = "Year_stand [all]"
) %>%
  as.data.frame() %>%
  mutate(cue = "JJA.beech") %>%
  mutate(
    Year = x *
      attr(Beech_Harvest$Year_stand, "scaled:scale") +
      attr(Beech_Harvest$Year_stand, "scaled:center")
  ) %>%
  ggplot(
    aes(x = Year, y = predicted)
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2,
    fill = "#3C5488FF"
  ) +
  geom_line(linewidth = 1.2, col = "#3C5488FF") +
  labs(
    x = "Year",
    y = "JJ(t-1) cue"
  ) +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  )

#now with beech seeds
model.tweedie.year.JJ1.beech <- glmmTMB(
  Seeds_cur ~
    scaled.tmax.harvest_JunJul1 +
      log1p(Seeds_prev) +
      log1p(Demand) +
      (1 | sitenewname),
  data = Beech_Harvest,
  family = tweedie()
)
summary(model.tweedie.year.JJ1.beech)
p_JJ1_beech <- plot_tweedie_cue(
  model = model.tweedie.year.JJ1.beech,
  data = Beech_Harvest,
  x_var = "tmax.harvest_JunJul1",
  x_scaled_var = "scaled.tmax.harvest_JunJul1",
  panel_title = "",
  x_lab = "Temperature (°C)",
  ylab.text = 20000,
  col = "#3C5488FF"
) +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  )

beech.tweedie.pred.da = plot_tweedie_cue_data(
  model = model.tweedie.year.JJ1.beech,
  data = Beech_Harvest,
  x_var = "tmax.harvest_JunJul1",
  x_scaled_var = "scaled.tmax.harvest_JunJul1"
)
beech.tweedie.pred.da
497 / 119
#now lets do rolling beech
list.beech.model = fit_Trends_models_over_window(
  Beech_Harvest,
  window = 10,
  cues = c(
    "tmax.harvest_JunJul1",
    "tmax.harvest_Mar",
    "tmax.harvest_AprMay"
  ),
  step = 5
)

rolling.beech.data = list.beech.model$rolling_data %>%
  mutate(
    log_mean_value_demand = log1p(mean_value_demand),
    log_seed_previous = log1p(Seeds_prev)
  )

Temporal.Trend.CVp.model.cuesJJA.beech <- glmmTMB::glmmTMB(
  log(CVp_seeds) ~
    mean_value_T_JunJul1 +
      (1 | sitenewname) +
      log_mean_value_demand,
  rolling.beech.data
)
summary(Temporal.Trend.CVp.model.cuesJJA.beech)


#make a tidy broom plot
figure.coef = bind_rows(
  broom::tidy(
    model.tweedie.year.JJ1.beech,
    conf.level = 0.95,
    conf.int = T
  ) %>%
    mutate(model = "Beech", cue = "June-JulyT1"),
  broom::tidy(
    model.tweedie.year.JJA1,
    conf.level = 0.95,
    conf.int = T
  ) %>%
    mutate(model = "Oaks", cue = "JJAT1"),
  broom::tidy(
    model.tweedie.year.march,
    conf.level = 0.95,
    conf.int = T
  ) %>%
    mutate(model = "Oaks", cue = "March"),
  broom::tidy(
    model.tweedie.year.aprmay,
    conf.level = 0.95,
    conf.int = T
  ) %>%
    mutate(model = "Oaks", cue = "April-May")
) %>%
  filter(
    term == "scaled.tmax.harvest_JunJul1" |
      term == "scaled.tmax.harvest_JunJulAug1" |
      term == "scaled.tmax.harvest_Mar" |
      term == "scaled.tmax.harvest_AprMay"
  ) %>%
  mutate(
    cue = factor(
      cue,
      levels = rev(c("JJAT1", "March", "April-May", "June-JulyT1"))
    ),
    model = factor(model, levels = c("Oaks", "Beech"))
  ) %>%
  ggplot(aes(x = cue, y = estimate, color = model)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.1,
    position = position_dodge(width = 0.5)
  ) +
  labs(
    x = "Species",
    y = "Estimated Effect of Cue on log(CVp)"
  ) +
  scale_color_manual(
    values = c("Oaks" = "#E64B35FF", "Beech" = "#3C5488FF")
  ) +
  theme(legend.title = element_blank()) +
  coord_flip() +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  xlab("") +
  ylab("Sensitivity to weather cues") +
  theme(legend.title = element_blank())
figure.coef

cowplot::save_plot(
  "figures/cues.pdf",
  figure.coef,
  nrow = 1.1,
  ncol = 1.2
)


#Figure 3
pred_JJA1_beech <- predict_glmmtmb_curve(
  model = Temporal.Trend.CVp.model.cuesJJA.beech,
  data = rolling.beech.data,
  x_var = "mean_value_T_JunJul1",
  demand_var = "log_mean_value_demand",
  site_var = "sitenewname",
  backtransform = "exp"
)


# plot(ggeffects::ggpredict(
#   Temporal.Trend.CVp.model.cuesJJA.beech,
#   terms = "mean_value_T_JunJul1"
# ))

View(pred_JJA1_beech$newdat)
#2.21/1.50
JJ1.rolCV.beech = ggplot() +
  geom_point(
    data = pred_JJA1_beech$part_resid,
    aes(x = x, y = predicted),
    alpha = 0.1,
    size = .3,
    col = "black"
  ) +
  geom_line(
    data = pred_JJA1_beech$newdat,
    aes(x = mean_value_T_JunJul1, y = fit),
    color = "#3C5488FF",
    linewidth = .8
  ) +
  geom_ribbon(
    data = pred_JJA1_beech$newdat,
    aes(x = mean_value_T_JunJul1, ymin = lwr, ymax = upr),
    fill = "#3C5488FF",
    alpha = 0.2
  ) +
  scale_fill_viridis_c() +
  labs(
    x = "Temperature in June-JulyT1",
    y = "CVp (partial residuals)"
  ) +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  )
JJ1.rolCV.beech

# beech.fig.5 = p_JJ1_beech +
#   scale_x_continuous(
#     breaks = c(16, 18, 20, 22, 24, 26, 28),
#     labels = c(16, 18, 20, 22, 24, 26, 28)
#   ) +
#   theme(
#     plot.margin = margin(7, 7, 7, 7)
#   ) +
#   theme(legend.position = c(.8, .8))
#
#
# beech.fig.5
# cowplot::save_plot(
#   "figures/beech.fig.5.pdf",
#   beech.fig.5,
#   nrow = 1,
#   ncol = 1.1
# )

ylim.def = c(-.2, 2.8)
rolling.cv.cues = JJ1.rolCV +
  ylim(ylim.def) +
  theme(
    plot.margin = margin(7, 7, 7, 7)
  ) +
  (Mar.rolCV +
    ylim(ylim.def) +
    ylab("") +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(7, 7, 7, 7)
    )) +
  (AprMay.rolCV +
    ylim(ylim.def) +
    ylab("") +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(7, 7, 7, 7)
    )) +
  (JJ1.rolCV.beech +
    ylim(ylim.def) +
    scale_x_continuous(
      breaks = c(18, 20, 22, 24, 26),
      labels = c(18, 20, 22, 24, 26)
    ) +
    theme(
      plot.margin = margin(7, 7, 7, 7),
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )) +
  plot_layout(nrow = 1, ncol = 4) +
  plot_annotation(tag_levels = 'A', tag_suffix = ')')
rolling.cv.cues
cowplot::save_plot(
  "figures/rolling.cv.cues.pdf",
  rolling.cv.cues,
  nrow = 1.1,
  ncol = 2
)

#####################################################################################
#MAPS figure 1
#####################################################################################
library(sf)
data.points.maps = bind_rows(Beech_Harvest, Oaks_Harvest) %>%
  dplyr::select(
    Longitude,
    Latitude,
    Species
  ) %>%
  distinct() %>%
  mutate(Species = if_else(Species == "Oaks_both", "Oaks", Species))

plot_locations <- st_as_sf(
  data.points.maps,
  coords = c("Longitude", "Latitude")
) %>%
  st_set_crs(4326) |>
  st_transform(3035) |> # ETRS89 / LAEA Europe (meters)
  st_jitter(factor = .004) |> # ~500 m jitter
  st_transform(4326)


poland <- rnaturalearth::ne_countries(
  country = c("Poland"),
  scale = "large",
  returnclass = "sf"
) |>
  st_make_valid()

around.poland <- rnaturalearth::ne_countries(
  country = c(
    "Germany",
    "Czechia",
    "Slovakia",
    "Ukraine",
    "Belarus",
    "Lithuania"
  ),
  scale = "large",
  returnclass = "sf"
) |>
  st_make_valid()

plot_locations$Species = factor(
  plot_locations$Species,
  levels = c("Oaks", "Beech")
)


g_pol = ggplot() +
  geom_sf(
    data = around.poland,
    fill = "grey90",
    colour = "black",
    linewidth = .05
  ) +
  geom_sf(data = poland, fill = "grey100", colour = "black", linewidth = .6) +
  geom_sf(
    data = plot_locations,
    aes(color = Species, fill = Species),
    alpha = .8,
    size = 2,
    stroke = .5
  ) +
  facet_grid(. ~ Species) +
  scale_y_continuous(
    breaks = c(50, 55),
    labels = function(b) paste0(b, "°N"),
    position = "right"
  ) +
  coord_sf(xlim = c(12, 26), ylim = c(49, 56), expand = FALSE) +
  scale_x_continuous(
    breaks = c(15, 20, 25),
    labels = function(b) paste0(b, "°E")
  ) +
  labs(x = NULL, y = NULL) +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  ) +
  #theme_bw(base_size = 16) +
  theme(plot.margin = margin(10, 90, 10, 10)) +
  scale_fill_manual(
    values = c("#E64B35FF", "#3C5488FF")
  ) +
  scale_color_manual(
    values = c("#E64B35FF", "#3C5488FF")
  ) +
  theme(legend.position = "none") +
  theme(panel.spacing = unit(0.1, "cm"))
g_pol


cowplot::save_plot(
  here::here('figures/maps.oaks.harvest.pdf'),
  g_pol,
  #dpi = 300,
  nrow = 1,
  ncol = 1.5
)

############################################################################
#correlation between oaks
############################################################################
comparison.oaks = polish_decade %>%
  left_join(
    prop_zero_by_site %>% dplyr::select(Species, NADL, prop_zero_cur),
    by = c("Species", "NADL")
  ) %>%
  filter(prop_zero_cur < 0.8) %>%
  filter(Species %in% c("Sessile oak", "Pedunculate oak")) %>%
  dplyr::select(Species, NADL, Year, Seeds_cur) %>%
  group_by(NADL, Year) %>%
  pivot_wider(names_from = Species, values_from = Seeds_cur) %>%
  group_by(NADL) %>%
  summarise(
    n_years = sum(!is.na(`Pedunculate oak`) & !is.na(`Sessile oak`)),
    cor = ifelse(
      n_years >= 5, # minimum overlap
      cor(
        `Pedunculate oak`,
        `Sessile oak`,
        method = "spearman",
        use = "complete.obs"
      ),
      NA_real_
    ),
    .groups = "drop"
  )

summary(comparison.oaks)
mean(comparison.oaks$cor, na.rm = TRUE)
sd(comparison.oaks$cor, na.rm = TRUE)
hist(comparison.oaks$cor)

############################################################################
#distance decay
############################################################################
#response to Andrew comments
library(tidyverse)
library(geosphere)

# data cleaning
make_site_series <- function(data, sp) {
  data %>%
    filter(Species == sp, prop_zero_cur < 0.8) %>%
    dplyr::select(NADL, Year, Seeds_cur, Latitude, Longitude)
}

robur <- make_site_series(
  polish_decade %>%
    left_join(
      prop_zero_by_site %>%
        dplyr::select(Species, NADL, prop_zero_cur),
      by = c("Species", "NADL")
    ),
  "Pedunculate oak"
)

petraea <- make_site_series(
  polish_decade %>%
    left_join(
      prop_zero_by_site %>%
        dplyr::select(Species, NADL, prop_zero_cur),
      by = c("Species", "NADL")
    ),
  "Sessile oak"
)

#function basde on Kuba code cleaned with Claude stuff
pairwise_cor_dist <- function(df1, df2, pair_label, min_overlap = 5) {
  sites1 <- unique(df1$NADL)
  sites2 <- unique(df2$NADL)

  cross <- expand.grid(NADL1 = sites1, NADL2 = sites2, stringsAsFactors = FALSE)

  if (pair_label != "robur–petraea") {
    cross <- cross %>% filter(NADL1 < NADL2)
  }

  map_dfr(seq_len(nrow(cross)), function(i) {
    s1 <- cross$NADL1[i]
    s2 <- cross$NADL2[i]

    ts1 <- df1 %>% filter(NADL == s1) %>% dplyr::select(Year, Seeds_cur)
    ts2 <- df2 %>% filter(NADL == s2) %>% dplyr::select(Year, Seeds_cur)

    shared <- inner_join(ts1, ts2, by = "Year", suffix = c(".1", ".2"))
    if (nrow(shared) < min_overlap) return(NULL)

    r <- cor(
      shared$Seeds_cur.1,
      shared$Seeds_cur.2,
      method = "spearman",
      use = "complete.obs"
    )

    coords1 <- df1 %>%
      filter(NADL == s1) %>%
      slice(1) %>%
      dplyr::select(Longitude, Latitude)
    coords2 <- df2 %>%
      filter(NADL == s2) %>%
      slice(1) %>%
      dplyr::select(Longitude, Latitude)
    dist_km <- geosphere::distHaversine(
      c(coords1$Longitude, coords1$Latitude),
      c(coords2$Longitude, coords2$Latitude)
    ) /
      1000

    tibble(NADL1 = s1, NADL2 = s2, r = r, dist_km = dist_km, pair = pair_label)
  })
}

#merign
dd <- bind_rows(
  pairwise_cor_dist(robur, robur, "Robur–Robur"),
  pairwise_cor_dist(petraea, petraea, "Petraea–Petraea"),
  pairwise_cor_dist(robur, petraea, "Robur–Petraea")
)


spearman.sync.oaks = ggplot(
  dd %>%
    mutate(
      paired.name = if_else(
        pair == "Robur–Robur",
        "Q. robur–Q. robur",
        if_else(
          pair == "Petraea–Petraea",
          "Q. petraea–Q. petraea",
          "Q. robur–Q. petraea"
        )
      )
    ),
  aes(x = dist_km, y = r)
) +
  #geom_point(alpha = 0.25, size = 1.2, shape = 16) +
  geom_hex() +
  geom_smooth(
    method = "loess",
    span = 0.6,
    se = F,
    col = "red",
    alpha = 0.12,
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.4
  ) +
  facet_grid(. ~ paired.name) +
  scale_color_distiller(palette = "Blues", direction = -1) +
  scale_fill_distiller(palette = "Blues", direction = -1) +
  scale_x_continuous(labels = scales::comma_format(suffix = " km")) +
  labs(
    x = "Distance between sites (km)",
    y = "Spearman correlation",
    colour = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")
spearman.sync.oaks
cowplot::save_plot(
  here::here('figures/correlation.sync.pdf'),
  spearman.sync.oaks,
  nrow = 1,
  ncol = 1.2
)
