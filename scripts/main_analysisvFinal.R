# Main Analysis Script for Polish Harvest Dataset
# Complete workflow for analyzing mast seeding trends
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
install.packages("showtext")
library(showtext)

# Add Roboto from Google Fonts
font_add_google("Roboto", "roboto")
# Enable automatic use
showtext_auto()

#HERE IS THE UPDATED FROM JESSIE
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

polish_harvest_df_simpl.init %>%
  filter(str_detect(Species, "Oak|oak")) %>%
  arrange(Species, NADL, Year) %>%
  ggplot(aes(x = Year, y = Seeds_prev, group = NADL)) +
  geom_line() +
  facet_grid(. ~ Species)

polish_harvest_df_simpl.init %>%
  filter(str_detect(Species, "Oak|oak")) %>%
  filter(str_detect(Species, "Oak|oak")) %>%
  group_by(NADL, Year, Species) %>%
  summarise(
    prop_na = mean(is.na(Demand)),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = Year, y = NADL, fill = prop_na)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Proportion NA",
    option = "plasma",
    limit = c(0, 1)
  ) +
  theme_minimal() +
  labs(
    x = "Year",
    y = "Site"
  ) +
  facet_grid(. ~ Species)

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


oak_dat <- polish_decade %>%
  left_join(
    prop_zero_by_site %>% select(Species, NADL, prop_zero_cur),
    by = c("Species", "NADL")
  ) %>%
  filter(prop_zero_cur < 0.8) %>%
  filter(Species %in% c("Sessile oak", "Pedunculate oak")) %>%
  filter(Year >= 2007, Year <= 2022) %>%
  select(Species, NADL, Year, Seeds_cur, Latitude, Longitude)

coords <- oak_dat %>%
  group_by(Species, NADL) %>%
  summarise(
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE),
    .groups = "drop"
  )

pairwise_sync <- function(dat, coords, sp1, sp2, min_overlap = 5) {
  #based on Kuba code beforre
  #and cleanded with GPT
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

sync_robur <- pairwise_sync(
  oak_dat,
  coords,
  sp1 = "Pedunculate oak",
  sp2 = "Pedunculate oak",
  min_overlap = 5
)

sync_petraea <- pairwise_sync(
  oak_dat,
  coords,
  sp1 = "Sessile oak",
  sp2 = "Sessile oak",
  min_overlap = 5
)

sync_cross <- pairwise_sync(
  oak_dat,
  coords,
  sp1 = "Pedunculate oak",
  sp2 = "Sessile oak",
  min_overlap = 5
)

sync_all <- bind_rows(sync_robur, sync_petraea, sync_cross) %>%
  filter(!is.na(correlation), !is.na(distance_km))

ggplot(sync_all, aes(x = distance_km, y = correlation, color = comparison)) +
  #geom_point(alpha = 0.25, size = 1) +
  geom_smooth(method = "loess", se = F, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_minimal() +
  labs(
    x = "Distance between sites (km)",
    y = "Spearman synchrony",
    color = "Comparison"
  )

summary(comparison.oaks)
mean(comparison.oaks$cor, na.rm = TRUE)
sd(comparison.oaks$cor, na.rm = TRUE)
hist(comparison.oaks$cor)
plotrix::std.error(comparison.oaks$cor)

#check correlation between demand and current year seeds
plot(
  log1p(Extend.Harvest.Filled$Demand),
  log1p(Extend.Harvest.Filled$Seeds_cur)
)
cor(log1p(Extend.Harvest.Filled$Demand), log1p(Extend.Harvest.Filled$Seeds_cur))


Extend.Harvest.Filled %>%
  group_by(Year) %>%
  mutate(zero = sum(Demand == 0, na.rm = T), n = n(), prop_zero = zero / n) %>%
  ungroup() %>%
  ggplot(aes(Year, prop_zero)) +
  geom_point() +
  ylab("Porpotion 0 for demand (no data replacement for demand)")

Extend.Harvest.Filled %>%
  filter(Species == "Oaks_both") %>%
  group_by(NADL, Year) %>%
  summarise(
    n_total = n(),
    n_na = sum(is.na(Demand)),
    prop_na = n_na / n_total,
    .groups = "drop"
  )

#should not have NA for demand because I replaced them with mean, but I want to check if there is any pattern in the NA before replacement
Extend.Harvest.Filled %>%
  filter(Species == "Oaks_both") %>%
  group_by(NADL, Year) %>%
  summarise(
    prop_na = mean(is.na(Demand)),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = Year, y = NADL, fill = prop_na)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Proportion NA",
    option = "plasma",
    limit = c(0, 1)
  ) +
  theme_minimal() +
  labs(
    x = "Year",
    y = "Site"
  )
Extend.Harvest.Filled %>%
  filter(Species == "Oaks_both") %>%
  group_by(Year) %>%
  mutate(zero = sum(Demand == 0, na.rm = T), n = n(), prop_zero = zero / n) %>%
  ungroup() %>%
  ggplot(aes(Year, prop_zero)) +
  geom_point() +
  ylab("Porpotion 0 for demand (replacement for demand)")


# versionkuba. = initial.data.harvest.load.kubav %>%
#   filter(Species == "Oaks_both" & Year < 1990) %>%
#   dplyr::select(seeds, NADL, Year) %>%
#   rename(seeds_kuba = seeds)
#
# versionval. = Extend.Harvest.Filled %>%
#   filter(Species == "Oaks_both" & Year < 1990) %>%
#   dplyr::select(seeds, NADL, Year) %>%
#   rename(seeds_val = seeds)
#
# test. = full_join(versionkuba., versionval.) %>%
#   filter(is.na(seeds_kuba) | is.na(seeds_val))
#
# test. %>% filter(is.na(seeds_kuba) | is.na(seeds_val)) %>% distinct(NADL)

#plot(log1p(seeds)~log1p(lag1_seeds), data = Extend.Harvest.Filled)

#analysis with tmzx spring now
initial.tmax.harvest.v1 <- readxl::read_excel(
  here::here("data", "tmax.xlsx")
) %>%
  #filter(Year >= (1987 - 1) & Year < 2023) %>%
  rename(
    Karniszewice = Karnieszewice,
    Kwidzyń = Kwidzyn,
    `Połczyn Zdrój` = Połczyn,
    `Międzyrzec Podlaski` = Międzyrzec,
    `Krosno Odrzańskie` = Krosno,
    `Buda stalowska` = `Buda Stalowska`
  )

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
  mutate(Year = Year + 1) # Shift December to the following year (I add 1 to climate because I merge this file to seed prod later) for indexing (the new year climate for example 1951, is related to seed production 1951 and also because I added 1 year to the old year, make the new year 1951)

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

#cues models
ggplot(Oaks_Harvest, aes(x = tmax.harvest_Dec)) + geom_histogram()
ggplot(Oaks_Harvest, aes(x = tmax.harvest_Mar)) + geom_histogram()
ggplot(Oaks_Harvest, aes(x = tmax.harvest_AprMay)) + geom_histogram()
ggplot(Oaks_Harvest, aes(x = tmax.harvest_May)) + geom_histogram()
ggplot(Oaks_Harvest, aes(x = tmax.harvest_JunJulAug1)) + geom_histogram()

#because I have to standardised year, because some model did not converge
sd_year <- sd(Oaks_Harvest$Year, na.rm = TRUE)


cues.december = glmmTMB(
  tmax.harvest_Dec ~ Year + (Year | sitenewname),
  data = Oaks_Harvest
)
summary(cues.december)

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


cues.JJ2 = glmmTMB(
  tmax.harvest_JunJulAug2 ~ Year + (1 | sitenewname),
  data = Oaks_Harvest
)
summary(cues.JJA1)
summary(cues.JJ2)

Oaks_Harvest %>%
  dplyr::select(
    Year,
    #tmax.harvest_May,
    #tmax.harvest_Dec,
    #tmax.harvest_MarApr,
    #T_anom_Dec,
    T_anom_Mar,
    T_anom_AprMay,
    T_anom_JunJulAug1,
    T_anom_JunJulAug2,
    #tmax.fullyear.harvest,
    NADL
  ) %>%
  pivot_longer(-c(Year, NADL), names_to = "cue", values_to = "tmax") %>%
  ggplot(aes(x = Year, y = tmax, col = cue)) +
  geom_point(alpha = .01) +
  geom_smooth(method = "lm") +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_color_viridis_d(option = "turbo") +
  ylab("Temperature anomalies")

##############################################################################
##############################################################################
#############SEED PRODUCTION MODELS###########################################
##############################################################################
#####################################
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
hist(log1p(Oaks_Harvest$Seeds_prev))
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

# plot(
#   ggeffects::ggpredict(
#     model.tweedie.year.march,
#     terms = "scaled.tmax.harvest_Mar [all]"
#   ),
#   show_data = T
# )

AIC(
  # model.tweedie.year.december,
  model.tweedie.year.march,
  model.tweedie.year.aprmay
)

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
AIC(model.tweedie.year.JJA1)
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
AIC(model.tweedie.year.JJ2)
#check correlation between cues
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

performance::model_performance(model.tweedie.year.JJA1)
performance::check_model(model.tweedie.year.JJA1)
value510JJA = plot_tweedie_cue_data(
  model = model.tweedie.year.JJA1,
  data = Oaks_Harvest,
  x_var = "tmax.harvest_JunJulAug1",
  x_scaled_var = "scaled.tmax.harvest_JunJulAug1"
)
p_JJA1 <- plot_tweedie_cue(
  model = model.tweedie.year.JJA1,
  data = Oaks_Harvest,
  x_var = "tmax.harvest_JunJulAug1",
  x_scaled_var = "scaled.tmax.harvest_JunJulAug1",
  panel_title = "Previous summer (JJA1)",
  x_lab = "Temperature (°C)",
  ylab.text = 4000
) +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  )

# me <- ggeffects::predict_response(
#   model.tweedie.year.march,
#   terms = "scaled.tmax.harvest_Mar [all]"
# )
# part_resid_anom_c_worldclim <- tibble(
#   ggeffects::residualize_over_grid(me, model.tweedie.year.march)
# )
#
# sim_res <- simulateResiduals(model.tweedie.year.march)
# plot(sim_res)

# dat <- model.frame(model.tweedie.year.march)
#
# dat$.fitted_resp <- predict(model.tweedie.year.march, type = "response")
# dat$.fitted_link <- predict(model.tweedie.year.march, type = "link")
# dat$.resid_resp <- residuals(model.tweedie.year.march, type = "response")
# dat$.resid_pears <- residuals(model.tweedie.year.march, type = "pearson")
# range(dat$.fitted_resp, na.rm = TRUE)
# range(dat$.fitted_link, na.rm = TRUE)
# range(dat$.resid_resp, na.rm = TRUE)
#pred_term <- predict(model.tweedie.year.march, type = "link")
#dat$partial_resid <- residuals(model.tweedie.year.march, type = "pearson") +
#  pred_term

## overall prediction (main effect) ====
# pred_CV_anom_c_worldclim <- ggeffects::ggpredict(
#   model.tweedie.year.march,
#   terms = "scaled.tmax.harvest_Mar [all]"
# )
#
# ## Decline from baseline values  for different MAT ====
# dat = ggeffects::ggpredict(
#   model.tweedie.year.march,
#   terms = "scaled.tmax.harvest_Mar [all]"
# ) %>%
#   as_tibble()
#
#
# ## partial residuals
# part_resid_anom_c_worldclim <- tibble(
#   ggeffects::residualize_over_grid(
#     pred_CV_anom_c_worldclim,
#     CVclimtrend_anom_worldclim_int
#   )
# ) %>%
#   bind_cols(CVclimtrend_anom_worldclim_int$frame) %>%
#   mutate(`scale(rol_demand)` = `scale(rol_demand)`[, 1]) %>%
#   left_join(
#     buki_metr_step5_w_WorldClim %>%
#       mutate(
#         `scale(rol_demand)` = scale(rol_demand)[, 1],
#         `log(rol_CV)` = log(rol_CV)
#       )
#   )

# plot(me, show_residuals = T, alpha = .1) + ylim(0, 30000)
# plot(me, show_data = T, alpha = .1)
# visreg::visreg(model.tweedie.year.march)
# gb <- ggplot_build(p)
# length(gb$data)
# gb$data[[1]]
# gb$data[[2]]
# gb$data[[3]]
# resid_points <- gb$data[[1]]
# head(resid_points)

# dat <- model.frame(model.tweedie.year.march)
#
# dat <- model.frame(model.tweedie.year.march)
#
# dat$.fitted_resp <- predict(model.tweedie.year.march, type = "response")
# dat$.resid_resp <- residuals(model.tweedie.year.march, type = "response")
# dat$.resid_pearson <- residuals(model.tweedie.year.march, type = "pearson")
# ggplot(dat, aes(x = scaled.tmax.harvest_Mar, y = .fitted_resp)) +
#   geom_point(alpha = 0.5) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_classic()
#
# ggplot(dat, aes(x = scaled.tmax.harvest_Mar, y = .resid_resp)) +
#   geom_point(alpha = 0.5) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_classic()
# ggplot(dat, aes(x = scaled.tmax.harvest_Mar, y = .resid_pearson)) +
#   geom_point(alpha = 0.5) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_classic()
performance::check_model(model.tweedie.year.march)

predict(model.tweedie.year.march, newdata = c(5, 10))
value510 = plot_tweedie_cue_data(
  model = model.tweedie.year.march,
  data = Oaks_Harvest,
  x_var = "tmax.harvest_Mar",
  x_scaled_var = "scaled.tmax.harvest_Mar"
)
#5 degree = 2245.592
#10 degree = 1626.629

p_march <- plot_tweedie_cue(
  model = model.tweedie.year.march,
  data = Oaks_Harvest,
  x_var = "tmax.harvest_Mar",
  x_scaled_var = "scaled.tmax.harvest_Mar",
  panel_title = "Early spring (Mar)",
  x_lab = "Temperature (°C)",
  ylab.text = 4000
) +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  )
performance::check_model(model.tweedie.year.aprmay)
p_aprmay <- plot_tweedie_cue(
  model = model.tweedie.year.aprmay,
  data = Oaks_Harvest,
  x_var = "tmax.harvest_AprMay",
  x_scaled_var = "scaled.tmax.harvest_AprMay",
  panel_title = "Late spring (Apr–May)",
  x_lab = "Temperature (°C)"
) +
  hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  )

Cues.seeds = (p_JJA1 +
  theme(plot.margin = margin(5, 5, 5, 5)) +
  ylim(0, 4100) |
  p_march +
    ylim(0, 4100) +
    (theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )) |
  p_aprmay +
    ylim(0, 4100) +
    (theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )))

cowplot::save_plot(
  "figures/cues.pdf",
  Cues.seeds,
  nrow = 1.2,
  ncol = 1.8
)

# m_zi <- glmmTMB(
#   Seeds_cur ~ tmax.harvest + log1p(Demand) + (1 | sitenewname),
#   ziformula = ~ tmax.harvest + Seeds_prev, # + log1p(Demand),
#   family = tweedie(), #detualt is log llink
#   data = Oaks_Harvest
# )
# summary(m_zi)
#
# res.zi = simulateResiduals(m_zi, quantreg = T)
# plot(res.zi) #those test are also shitty when too much observation

#let try to use climwin here
library(climwin)
initial.tmax.harvest.v1 #weather
Oaks_Harvest #seed production
# n1 <- unique(Oaks_Harvest_climwin_data$NADL)
# n2 <- unique(initial.tmax.harvest_climwin_data$NADL)
# matching_NADL <- intersect(n1, n2)
# matching_NADL
# only_in_oaks <- setdiff(n1, n2)
# only_in_oaks
# only_in_climate <- setdiff(n2, n1)
# only_in_climate #should have names because some sites dont have seed data for oaks
# data.frame(
#   dataset = c(
#     rep("Both", length(matching_NADL)),
#     rep("Only Oaks_Harvest", length(only_in_oaks)),
#     rep("Only Climate", length(only_in_climate))
#   ),
#   NADL = c(matching_NADL, only_in_oaks, only_in_climate)
# )
# Oaks_Harvest_climwin_data = Oaks_Harvest %>%
#   mutate(Date = as.Date(paste0(Year, '-10-01'))) %>%
#   mutate(NADL = as_factor(NADL), log.seed = log1p(Seeds_cur))
#
# #plot(log.seed ~ Seeds_cur, Oaks_Harvest_climwin_data)
# initial.tmax.harvest_climwin_data = initial.tmax.harvest.v1 %>%
#   pivot_longer(
#     -c(Date, Day, Month, Year),
#     names_to = "NADL",
#     values_to = "tmax"
#   ) %>%
#   mutate(NADL = as_factor(NADL))
#
#
# climwin_output <- slidingwin(
#   xvar = list(
#     Temp = initial.tmax.harvest_climwin_data$tmax
#   ),
#   cdate = initial.tmax.harvest_climwin_data$Date,
#   bdate = Oaks_Harvest_climwin_data$Date,
#   baseline = lm(log.seed ~ 1, data = Oaks_Harvest_climwin_data),
#   cinterval = "month",
#   range = c(9, 0), #number months before
#   type = "absolute",
#   refday = c(1, 10), #absolute based on 01-10 1st october
#   stat = "mean",
#   func = "lin",
#   cmissing = "method1",
#   spatial = list(
#     Oaks_Harvest_climwin_data$NADL,
#     initial.tmax.harvest_climwin_data$NADL
#   )
# )
# check.cw = climwin_output[[1]]$BestModelData
# hist(check.cw$climate)
# check = initial.tmax.harvest_climwin_data %>%
#   right_join(Oaks_Harvest_climwin_data %>% dplyr::select(NADL)) %>%
#   filter(Month %in% c(5)) %>%
#   group_by(NADL, Year) %>%
#   summarise(mean = mean(tmax, na.rm = TRUE))
# hist(check$mean)

# newdat <- Oaks_Harvest %>%
#   summarise(
#     tmin = min(tmax.harvest, na.rm = TRUE),
#     tmax = max(tmax.harvest, na.rm = TRUE)
#   ) %>%
#   tidyr::expand_grid(
#     tmax.harvest = seq(tmin, tmax, length.out = 200),
#     sitenewname = NA
#   )
#
# #adpat code from previous Jims code
# invlogit <- function(x) 1 / (1 + exp(-x))
# b_cond <- fixef(m_zi)$cond
# V_cond <- vcov(m_zi)$cond
#
# b_zi <- fixef(m_zi)$zi
# V_zi <- vcov(m_zi)$zi
#
# X_cond <- model.matrix(delete.response(terms(m_zi, component = "cond")), newdat)
# X_zi <- model.matrix(delete.response(terms(m_zi, component = "zi")), newdat)
#
# R <- 1000
# B_cond <- MASS::mvrnorm(R, mu = b_cond, Sigma = V_cond)
# B_zi <- MASS::mvrnorm(R, mu = b_zi, Sigma = V_zi)
# eta_cond <- X_cond %*% t(B_cond)
# eta_zi <- X_zi %*% t(B_zi)
# mu_cond <- exp(eta_cond) # conditional mean
# p_zi <- invlogit(eta_zi)
# Ey <- (1 - p_zi) * mu_cond #Overall expected mean including ZI: E[y] = (1 - p_zi) * mu_cond
# pred_df <- newdat %>%
#   mutate(
#     pred = apply(Ey, 1, median),
#     lo = apply(Ey, 1, quantile, probs = 0.025),
#     hi = apply(Ey, 1, quantile, probs = 0.975)
#   )
#
# #now aggregate obs for the original data
# bin_width <- 0.5
#
# obs_agg <- Oaks_Harvest %>%
#   mutate(
#     t_bin = floor(tmax.harvest / bin_width) * bin_width #to agg bins
#   ) %>%
#   group_by(sitenewname, t_bin) %>%
#   summarise(
#     tmax_h = mean(tmax.harvest, na.rm = TRUE),
#     seeds_mean = mean(seeds, na.rm = TRUE),
#     n = n(),
#     .groups = "drop"
#   )
#
# ggplot() +
#   # ribbon (95% interval)
#   #geom_ribbon(
#   #  data = pred_df,
#   #  aes(x = tmax.harvest, ymin = lo, ymax = hi),
#   #  alpha = 0.2
#   #) +
#   # median prediction
#   geom_line(
#     data = pred_df,
#     aes(x = tmax.harvest, y = log1p(pred)),
#     linewidth = 1
#   ) +
#   # aggregated points (per site per temperature bin)
#   geom_point(
#     data = obs_agg %>% filter(seeds_mean > 0),
#     aes(x = tmax_h, y = log1p(seeds_mean)),
#     alpha = 0.1,
#     size = 1.5
#   ) +
#   labs(
#     x = "Spring temperature (tmax.harvest)",
#     y = "Seeds (observed and predicted mean) (log trans for visualisation)"
#   ) +
#   theme_classic()

# newdat <- Oaks_Harvest %>%
#   summarise(
#     tmin = min(tmax.harvest, na.rm = TRUE),
#     tmax = max(tmax.harvest, na.rm = TRUE)
#   ) %>%
#   tidyr::expand_grid(
#     tmax.harvest = seq(tmin, tmax, length.out = 200),
#     sitenewname = NA
#   )
# newdat$Ey <- predict(
#   m_zi,
#   newdata = newdat,
#   type = "response",
#   re.form = NA
# )
# ggplot(newdat, aes(tmax.harvest, Ey)) +
#   geom_line(linewidth = 1) +
#   labs(x = "Spring temperature (tmax.harvest)",
#        y = "Predicted mean seeds (E[y])") +
#   theme_classic()
# eta_zi <- predict(
#   m_zi,
#   newdata = newdat,
#   type = "link",
#   zitype = "zlink",   # gives zi linear predictor
#   re.form = NA
# )
#
# newdat$p_zi <- plogis(eta_zi)
# ggplot(newdat, aes(tmax.harvest, p_zi)) +
#   geom_line(linewidth = 1) +
#   labs(x = "Spring temperature (tmax.harvest)",
#        y = "Predicted P(structural zero)") +
#   ylim(0, 1) +
#   theme_classic()
# #conditional prediction of being nn zero
# newdat$Ey_pos <- newdat$Ey / pmax(1 - newdat$p_zi, 1e-6)
# ggplot(newdat, aes(tmax.harvest, Ey_pos)) +
#   geom_line(linewidth = 1) +
#   labs(x = "Spring temperature (tmax.harvest)",
#        y = "Predicted mean seeds given non-zero (approx.)") +
#   theme_classic()

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

# rolling_seeds <- Oaks_Harvest %>%
#   group_by(sitenewname) %>%
#   group_modify(~ compute_rolling_metrics(.x, Seeds_cur, window = 10, step = 1)) %>%
#   ungroup() %>%
#   rename_with(
#     ~ paste0(.x, "_seeds"),
#     c(mean_value, sd_value, CVp, q25, q90, p_zero, p_nonzero)
#   )
ggplot(rolling_seeds, aes(CVp_seeds, mean_value_seeds)) + geom_point()
ggplot(rolling_seeds, aes(CVp_seeds, p_nonzero_seeds)) +
  geom_point()
ggplot(rolling_seeds, aes(CVp_seeds, p_zero_seeds)) +
  geom_point()

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

# rolling_tmax_anomalies <- Oaks_Harvest %>%
#   group_by(sitenewname) %>%
#   group_modify(
#     ~ compute_rolling_metrics(.x, T_spring_z_May, window = 10)
#   ) %>%
#   dplyr::select(-p_zero, -p_nonzero, -q25, -q90) %>%
#   ungroup() %>%
#   rename_with(~ paste0(.x, "_TMayzanomalies"), c(mean_value, sd_value, CVp))

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

zero_year <- rolling_metrics_oaks %>%
  group_by(Year) %>%
  summarise(
    prop_zero = mean(Seeds_cur == 0, na.rm = TRUE),
    .groups = "drop"
  )

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
##################################################
##################################################
##################################################
#different plot
ggplot(
  rolling_metrics_oaks,
  aes(x = Year_center, y = CVp_seeds, group = sitenewname, color = sitenewname)
) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Rolling Coefficient of Variation (CV) for Oak Seed Harvest",
    x = "Year",
    y = "Coefficient of Variation (CV)"
  ) +
  theme(legend.position = "none")

ggplot(
  rolling_metrics_oaks %>%
    mutate(Year = as_factor(Year)),
  aes(
    x = p_zero_seeds,
    y = CVp_seeds,
    color = Year,
    fill = Year
  )
) +
  geom_point(alpha = .4, shape = 21, size = .6) +
  geom_smooth(method = "lm", se = F) +

  theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

cv_summary <- rolling_metrics_oaks %>%
  group_by(Year_center) %>%
  summarise(
    med = median(CVp_seeds, na.rm = TRUE),
    q25 = quantile(CVp_seeds, 0.25, na.rm = TRUE),
    q75 = quantile(CVp_seeds, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(cv_summary, aes(x = Year_center)) +
  geom_ribbon(
    aes(ymin = q25, ymax = q75),
    fill = "darkred",
    col = "darkred",
    alpha = 0.3
  ) +
  geom_line(
    aes(y = med),
    linewidth = 1.2,
    color = "red"
  ) +
  theme_minimal() +
  ylab("CV")

cv_summary <- rolling_metrics_oaks %>%
  group_by(Year_center) %>%
  summarise(
    med = median(log_CVp_seeds, na.rm = TRUE),
    q25 = quantile(log_CVp_seeds, 0.25, na.rm = TRUE),
    q75 = quantile(log_CVp_seeds, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
ggplot(cv_summary, aes(x = Year_center, y = med)) +
  geom_pointrange(
    aes(ymin = q25, ymax = q75),
    color = "darkred",
    linewidth = 0.6
  ) +
  geom_point(
    color = "darkred",
    size = 2
  ) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  ylab("log CVp")

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


##################################################
##################################################
##################################################
ggplot(
  rolling_metrics_oaks,
  aes(
    x = mean_value_T_Mar,
    y = mean_value_Tfull,
    col = as_factor(Year_center)
  )
) +
  geom_point(alpha = .4) +
  scale_color_viridis_d() +
  geom_smooth(method = "lm", se = F) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_flip() +
  facet_wrap(. ~ Year_center) +
  theme(legend.position = "none")


ggplot(
  rolling_metrics_oaks,
  aes(
    x = mean_value_demand,
    y = CVp_seeds,
    group = sitenewname,
    color = sitenewname
  )
) +
  geom_point() +
  theme(legend.position = "none")


ggplot(
  rolling_metrics_oaks,
  aes(y = p_zero_seeds, x = CVp_seeds)
) +
  geom_point()

hist(rolling_metrics_oaks$log_CVp_seeds)
hist(rolling_metrics_oaks$CVp_seeds)
hist(rolling_metrics_oaks$mean_value_demand)
#test = Oaks_Harvest %>% filter(sitenewname=="Oaks_both_Andrychów")
#just double check with Jesie F functions, and it is working the same
#rollingCVp(test, val = "seeds", win = 10)

#I checked JJF paper, and she tool into account the deamnde as offset or covaraite, but
#more seed you have more demand you will have, so I will just check correlation first
#but I am not sure I would include both of them
#because I just increase the risk to inflate my stat models
Temporal.Trend.Mean <- glmmTMB::glmmTMB(
  mean_value_seeds ~
    Year_center +
      (1 | sitenewname),
  family = glmmTMB::tweedie(),
  rolling_metrics_oaks
)
summary(Temporal.Trend.Mean)
visreg::visreg(Temporal.Trend.Mean)
plot((mean_value_demand) ~ mean_value_seeds, data = rolling_metrics_oaks)

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

ggplot(rolling_metrics_oaks, aes(x = CVp_seeds, y = log1p(mean_value_seeds))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()

Temporal.Trend.CVp.model1 <- glmmTMB::glmmTMB(
  log_CVp_seeds ~
    Year +
      (1 | sitenewname) +
      log_mean_seeds +
      log1p(mean_value_demand),
  rolling_metrics_oaks
)
summary(Temporal.Trend.CVp.model1)
Temporal.Trend.CVp.model2 <- glmmTMB::glmmTMB(
  log_CVp_seeds ~
    Year +
      (1 | sitenewname) +
      log1p(mean_value_demand),
  rolling_metrics_oaks
)
summary(Temporal.Trend.CVp.model2)


MuMIn::AICc(Temporal.Trend.CVp.model1, Temporal.Trend.CVp.model2)
visreg::visreg(Temporal.Trend.CVp.model2)

#make bnetter prediction
newdat <- data.frame(
  Year = seq(
    min(rolling_metrics_oaks$Year),
    max(rolling_metrics_oaks$Year),
    length.out = 100
  ),
  mean_value_demand = mean(
    rolling_metrics_oaks$mean_value_demand,
    na.rm = TRUE
  ),
  sitenewname = NA # important for population-level prediction
)

pred <- predict(
  Temporal.Trend.CVp.model2,
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
    upr_CV = exp(upr_log)
  )

ggplot(newdat, aes(x = Year, y = fit_CV)) +
  geom_ribbon(aes(ymin = lwr_CV, ymax = upr_CV), fill = "red", alpha = 0.2) +
  geom_line(color = "red", linewidth = 1.2) +
  theme_minimal() +
  labs(
    y = "Predicted CV",
    title = "Temporal trend in CV (population-level)"
  )

prediction.cv = ggplot() +
  geom_point(
    data = rolling_metrics_oaks,
    aes(x = Year, y = CVp_seeds),
    alpha = 0.1,
    size = .8
  ) +
  geom_ribbon(
    data = newdat,
    aes(x = Year, ymin = lwr_CV, ymax = upr_CV),
    fill = "red",
    alpha = 0.2
  ) +
  geom_line(
    data = newdat,
    aes(x = Year, y = fit_CV),
    color = "red",
    linewidth = 1.2
  ) +
  theme_minimal() +
  xlab("Year") +
  ylab("Rolling CV of seed production") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  ylim(0.15, 3.3)
prediction.cv

ggplot(
  newdat %>% filter(Year %in% c(1997, 2022)) %>% mutate(Year = as_factor(Year)),
  aes(x = Year, y = fit_CV)
) +
  geom_pointrange(aes(ymin = lwr_CV, ymax = upr_CV)) +
  theme_minimal() +
  labs(
    y = "Predicted CV"
  ) +
  ylim(0.15, 3.3)


ggplot() +
  geom_point(
    data = rolling_metrics_oaks,
    aes(x = Year, y = CVp_seeds),
    alpha = 0.1,
    size = .8
  ) +
  geom_bin2d(bins = 70) +
  geom_ribbon(
    data = newdat,
    aes(x = Year, ymin = lwr_CV, ymax = upr_CV),
    fill = "red",
    alpha = 0.2
  ) +
  geom_line(
    data = newdat,
    aes(x = Year, y = fit_CV),
    color = "red",
    linewidth = 1.2
  ) +
  theme_minimal() +
  xlab("Year") +
  ylab("Rolling CV of seed production") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  ylim(0.15, 3.3)

cowplot::save_plot(
  "figures/cv.rolling.decline.pdf",
  prediction.cv,
  nrow = 1,
  ncol = 1
)


#boxplot option
bin_df <- rolling_metrics_oaks %>%
  mutate(
    Year_bin = cut_width(
      Year_center, #Year_center
      width = 5,
      boundary = min(Year_center, na.rm = TRUE) #Year_center
    )
  )

bin_counts <- bin_df %>%
  group_by(Year_bin) %>%
  summarise(
    x = mean(Year_center, na.rm = TRUE), #Year_center
    n = sum(!is.na(CVp_seeds)),
    .groups = "drop"
  )
RColorBrewer::brewer.pal(9, "OrRd")

# boxplot.cv.trend = ggplot() +
#   geom_boxplot(
#     data = bin_df,
#     aes(
#       x = Year,
#       y = CVp_seeds,
#       group = Year_bin
#     ),
#     width = 3.5,
#     fill = "#FFF7EC",
#     color = "grey20",
#     #outlier.shape = NA,
#     alpha = 0.8,
#     size = .5
#   ) +
#   geom_ribbon(
#     data = newdat,
#     aes(x = Year, ymin = lwr_CV, ymax = upr_CV),
#     fill = "#D7301F",
#     alpha = 0.2
#   ) +
#   geom_line(
#     data = newdat,
#     aes(x = Year, y = fit_CV),
#     color = "#D7301F",
#     linewidth = 1.2
#   ) +
#   geom_text(
#     data = bin_counts,
#     aes(x = x, y = 3.2, label = paste0("n=", n)),
#     size = 4,
#     family = "Helvetica"
#   ) +
#   coord_cartesian(ylim = c(0.15, 3.3)) +
#   #hrbrthemes::theme_ipsum(base_size = 12, plot_title_size = 20) +
#   labs(x = "Year (start)", y = "Rolling CV of seed production")

# cowplot::save_plot(
#   "figures/cv.rolling.decline.boxplot.pdf",
#   boxplot.cv.trend,
#   nrow = 1,
#   ncol = .8
# )
##########
rolling_metrics_oaks %>%
  group_by(Year) %>%
  summarise(mean.cv = mean(CVp_seeds, na.rm = T)) %>%
  ggplot(
    aes(x = Year, y = mean.cv)
  ) +
  geom_point() +
  geom_smooth(method = "lm")

rolling_metrics_oaks %>%
  group_by(Year) %>%
  summarise(
    mean.cv = mean(CVp_seeds, na.rm = T),
    mean.pzero = mean(p_zero_seeds)
  ) %>%
  ggplot(
    aes(x = mean.cv, y = mean.pzero, col = Year)
  ) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_viridis_c()

ggplot(rolling_metrics_oaks, aes(x = CVp_seeds, y = p_zero_seeds)) +
  geom_point()


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

boxplot.cv.trend = ggplot() +
  geom_boxplot(
    data = bin_df,
    aes(
      x = Year_bin, #Year_center,
      y = CVp_seeds
    ),
    fill = "#FFF7EC",
    color = "grey20",
    alpha = 0.8,
    size = .5
  ) +
  geom_jitter(
    data = bin_df,
    aes(
      x = Year_bin, #Year_center,
      y = CVp_seeds,
      group = Year_bin
    ),
    color = "black",
    size = 0.1,
    alpha = 0.1
  ) +
  #coord_cartesian(ylim = c(0.15, 3.3)) +
  labs(x = "Year", y = "Rolling CV of seed production") +
  theme(axis.text.x = element_text(size = 12, angle = 45))
boxplot.cv.trend


boxplot.cv.trend <- ggplot(bin_df, aes(x = Year_center, y = CVp_seeds)) +
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
    aes(x = x, y = 3.2, label = paste0("n=", n)),
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
  scale_x_continuous(
    breaks = c(1995, 2000, 2005, 2010, 2015, 2020),
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
fixed_slope <- fixef(Temporal.Trend.CVp.model_rs)$cond["Year_c"]
site_slopes <- fixed_slope + ranef_slopes$Year_c
site_slopes_df <- data.frame(
  sitenewname = rownames(ranef_slopes),
  slope = site_slopes
) %>%
  mutate(
    percent_change = (exp(slope) - 1) * 100,
    backtransf_tempslope = exp(slope),
    group_tempslope = case_when(slope <= 0 ~ "<=0", TRUE ~ ">0")
  )

mean(site_slopes_df$percent_change, na.rm = TRUE)

VarCorr(Temporal.Trend.CVp.model_rs)
col.density = ggplot(site_slopes_df, aes(x = slope, fill = group_tempslope)) +
  geom_histogram(binwidth = 0.01, color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    x = "Slope (log scale)",
    y = "Count of Sites"
  ) +
  theme(legend.title = element_blank(), legend.position = c(.8, .8))

# 1. population-level prediction
# ---------------------------
newdat_pop <- data.frame(
  Year_c = seq(
    min(rolling_metrics_oaks$Year_c, na.rm = TRUE),
    max(rolling_metrics_oaks$Year_c, na.rm = TRUE),
    length.out = 100
  ),
  mean_value_demand = mean(
    rolling_metrics_oaks$mean_value_demand,
    na.rm = TRUE
  ),
  sitenewname = NA
)

pred_pop <- predict(
  Temporal.Trend.CVp.model_rs,
  newdata = newdat_pop,
  re.form = NA,
  se.fit = TRUE
)

year_mean <- mean(rolling_metrics_oaks$Year_center, na.rm = TRUE)

newdat_pop <- newdat_pop %>%
  mutate(
    fit_log = pred_pop$fit,
    se_log = pred_pop$se.fit,
    lwr_log = fit_log - 1.96 * se_log,
    upr_log = fit_log + 1.96 * se_log,
    fit_CV = exp(fit_log),
    lwr_CV = exp(lwr_log),
    upr_CV = exp(upr_log),
    Year_center = Year_c + year_mean
  )

# ---------------------------
# 2. site-specific prediction
# ---------------------------
newdat_site <- expand.grid(
  Year_c = seq(
    min(rolling_metrics_oaks$Year_c, na.rm = TRUE),
    max(rolling_metrics_oaks$Year_c, na.rm = TRUE),
    length.out = 100
  ),
  sitenewname = unique(rolling_metrics_oaks$sitenewname)
)

newdat_site$mean_value_demand <- mean(
  rolling_metrics_oaks$mean_value_demand,
  na.rm = TRUE
)

newdat_site$fit_log <- predict(
  Temporal.Trend.CVp.model_rs,
  newdata = newdat_site,
  re.form = NULL # include random effects
)

newdat_site <- newdat_site %>%
  mutate(
    fit_CV = exp(fit_log),
    Year_center = Year_c + year_mean
  )


decline.pop = ggplot() +
  # raw data
  geom_point(
    data = rolling_metrics_oaks,
    aes(x = Year_center, y = CVp_seeds),
    alpha = 0.12,
    size = 0.1
  ) +
  geom_line(
    data = newdat_site %>% left_join(site_slopes_df),
    aes(
      x = Year_center,
      y = fit_CV,
      group = sitenewname,
      col = group_tempslope
    ),
    alpha = 0.25,
    linewidth = 0.4
  ) +
  # global ribbon
  geom_ribbon(
    data = newdat_pop,
    aes(x = Year_center, ymin = lwr_CV, ymax = upr_CV),
    fill = "grey20",
    alpha = 0.2
  ) +
  # global line
  geom_line(
    data = newdat_pop,
    aes(x = Year_center, y = fit_CV),
    color = "black",
    linewidth = 1.2
  ) +
  labs(
    x = "Year",
    y = "Rolling CV of seed production"
  ) +
  theme(legend.position = "none")


col.density <- col.density +
  labs(x = "Slope", y = "Density") +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )
cowplot::save_plot(
  "figures/cv.rolling.decline.pdf",
  decline.pop +
    patchwork::inset_element(
      col.density + theme(legend.position = "none"),
      left = 0.62,
      bottom = 0.55,
      right = 0.98,
      top = 0.98
    ),
  nrow = 1,
  ncol = .8
)

m_temp <- glmmTMB(
  mean_value_T_Mar ~ Year_c + (Year_c | sitenewname),
  data = rolling_metrics_oaks
)
summary(m_temp)
VarCorr(m_temp) #ok
fixed_slope_temp <- fixef(m_temp)$cond["Year_c"]
ranef_temp <- ranef(m_temp)$cond$sitenewname

temp_slopes <- data.frame(
  sitenewname = rownames(ranef_temp),
  slope_temp = fixed_slope_temp + ranef_temp$Year_c
)

slope_compare <- site_slopes_df %>%
  left_join(temp_slopes, by = "sitenewname")
ggplot(slope_compare, aes(x = slope_temp, y = slope)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(
    x = "Temperature trend (°C per year)",
    y = "CV trend (log-scale per year)"
  )
summary(lm(slope ~ slope_temp, data = slope_compare))
site_temp_slopes <- rolling_metrics_oaks %>%
  group_by(sitenewname) %>%
  do({
    fit <- lm(mean_value_T_AprMay ~ Year_c, data = .)
    data.frame(slope = coef(fit)["Year_c"])
  }) %>%
  ungroup()

ggplot(site_temp_slopes, aes(x = slope)) +
  geom_density(fill = "steelblue", alpha = .4) +
  theme_minimal()

sd(site_slopes_df$slope)
sd(temp_slopes$slope)
sd(site_temp_slopes$slope, na.rm = T)

Temporal.Trend.CVp.cue <- glmmTMB(
  log_CVp_seeds ~
    Year_c +
      mean_value_T_AprMay +
      Year_c:mean_value_T_AprMay +
      (1 | sitenewname),
  data = rolling_metrics_oaks
)
summary(Temporal.Trend.CVp.cue)

# rolling_metrics_oaks.v2 = rolling_metrics_oaks %>%
#   group_by(sitenewname) %>%
#   mutate(
#     T_within = mean_value_T_AprMay - mean(mean_value_T_AprMay, na.rm = TRUE),
#     T_between = mean(mean_value_T_AprMay, na.rm = TRUE)
#   ) %>%
#   ungroup()
#
# ggplot(rolling_metrics_oaks.v2, aes(x = T_within, y = T_between)) + geom_point()
# cue.cv.modv2 = glmmTMB(
#   log_CVp_seeds ~ Year_c + T_within + T_between + (1 | sitenewname),
#   data = rolling_metrics_oaks.v2
# )
# summary(cue.cv.modv2)

library(stringi)
formaps.cv = site_slopes_df %>%
  left_join(
    Oaks_Harvest %>%
      dplyr::select(sitenewname, NADL, Latitude, Longitude) %>%
      distinct()
  ) %>%
  mutate(
    nadl_key = NADL %>%
      stri_trans_general("Latin-ASCII") %>% # AUGUSTÓW -> AUGUSTOW
      tolower() %>%
      trimws()
  ) %>%
  mutate(
    nadl_key = dplyr::recode(
      #still missing nowa deba ?
      nadl_key,
      'bielsk podlaski' = 'bielsk',
      "karniszewice" = "karnieszewice", #not exactly same area?
      'krosno odrzanskie' = 'krosno',
      'miedzyrzec podlaski' = 'miedzyrzec',
      "minsk mazowiecki" = "minsk",
      "piotrkow trybunalski" = "piotrkow",
      "polczyn zdroj" = "polczyn",
      "starogard gdanski" = "starogard"
    )
  )

plot_locations <- sf::st_as_sf(
  formaps.cv,
  coords = c("Longitude", "Latitude")
) %>%
  sf::st_set_crs(4326) |>
  sf::st_transform(3035) |> # ETRS89 / LAEA Europe (meters)
  sf::st_jitter(factor = .004) |> # ~500 m jitter
  sf::st_transform(4326)


poland <- rnaturalearth::ne_countries(
  country = c("Poland"),
  scale = "large",
  returnclass = "sf"
) |>
  sf::st_make_valid()

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
  sf::st_make_valid()

ggplot() +
  geom_sf(
    data = around.poland,
    fill = "grey95",
    colour = "black",
    linewidth = .05
  ) +
  geom_sf(data = poland, fill = "grey99", colour = "black", linewidth = .6) +
  geom_sf(
    data = plot_locations,
    aes(color = slope, fill = slope),
    alpha = .8,
    size = 1,
    stroke = 1
  ) +
  scale_y_continuous(
    breaks = c(50, 55),
    labels = function(b) paste0(b, "°N"),
    position = "right"
  ) +
  coord_sf(xlim = c(14, 26), ylim = c(49, 56), expand = FALSE) +
  scale_x_continuous(
    breaks = c(15, 20, 25),
    labels = function(b) paste0(b, "°E")
  ) +
  labs(x = NULL, y = NULL) +
  hrbrthemes::theme_ft_rc(base_size = 14, axis_title_size = 16) +
  theme(plot.margin = margin(10, 90, 10, 10)) +
  theme(legend.position = "right") +
  theme(panel.spacing = unit(0.1, "cm")) +
  scale_fill_distiller(palette = "BrBG") +
  scale_color_distiller(palette = "BrBG")
#scale_color_viridis_c() +
#scale_fill_viridis_c()

poland_admin1 <- rnaturalearth::ne_states(
  country = "Poland",
  returnclass = "sf"
) |>
  sf::st_make_valid()

sites_sf <- sf::st_as_sf(
  formaps.cv,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

sites_with_region <- sf::st_join(
  sites_sf,
  poland_admin1,
  join = sf::st_within
)

region_slopes <- sites_with_region %>%
  sf::st_drop_geometry() %>%
  group_by(name) %>%
  summarise(
    mean_slope = mean(slope, na.rm = TRUE),
    n_sites = n()
  )
poland_map <- poland_admin1 %>%
  left_join(region_slopes, by = "name")

ggplot() +
  geom_sf(
    data = around.poland,
    fill = "grey95",
    colour = "black",
    linewidth = .05
  ) +
  geom_sf(
    data = poland_map,
    aes(fill = mean_slope),
    colour = "black",
    linewidth = .3
  ) +
  scale_fill_distiller(
    palette = "BrBG",
    na.value = "grey90",
    name = "Mean slope"
  ) +
  coord_sf(xlim = c(14, 26), ylim = c(49, 56), expand = FALSE) +
  theme_minimal() +
  theme(legend.position = "right")

#obtained from https://gis.openforestdata.pl/layers/geonode%3Anadlesnictwa_wgs84/metadata_detail
nadl_sf <- sf::st_read(
  "/Users/valentinjourne/MyGitProject/MastingTrends/nadlesnictwa_wgs84/nadlesnictwa_wgs84.shp"
) # or whatever the .shp is named

nadl_sf2 <- nadl_sf %>%
  mutate(
    nadl_key = NAZWA %>%
      stri_trans_general("Latin-ASCII") %>%
      tolower() %>%
      trimws()
  )

#combine both
nadl_map <- nadl_sf2 %>%
  left_join(
    formaps.cv,
    by = "nadl_key"
  )

ggplot(nadl_map) +
  geom_sf(aes(fill = slope), color = "black", linewidth = 0.1) +
  scale_fill_distiller(palette = "BrBG", na.value = "black") +
  theme_minimal() +
  geom_sf(
    data = around.poland,
    fill = "grey95",
    colour = "black",
    linewidth = .05
  ) +
  scale_y_continuous(
    breaks = c(50, 55),
    labels = function(b) paste0(b, "°N"),
    position = "right"
  ) +
  coord_sf(xlim = c(14, 26), ylim = c(49, 56), expand = FALSE) +
  scale_x_continuous(
    breaks = c(15, 20, 25),
    labels = function(b) paste0(b, "°E")
  ) +
  labs(x = NULL, y = NULL) +
  hrbrthemes::theme_ipsum(base_size = 14, axis_title_size = 16) +
  theme(plot.margin = margin(10, 90, 10, 10)) +
  theme(legend.position = "right") +
  theme(panel.spacing = unit(0.1, "cm")) +
  scale_fill_distiller(palette = "BrBG") +
  scale_color_distiller(palette = "BrBG")

Temporal.Trend.ZeroProp.quad <- glmmTMB::glmmTMB(
  logit_prop_zero_seeds ~
    poly(Year, 2) + #, raw = TRUE
      (1 | sitenewname) +
      log1p(mean_value_demand),
  rolling_metrics_oaks
)
Temporal.Trend.ZeroProp.sing <- glmmTMB::glmmTMB(
  logit_prop_zero_seeds ~
    Year +
      (1 | sitenewname) +
      log1p(mean_value_demand),
  rolling_metrics_oaks
)
MuMIn::AICc(Temporal.Trend.ZeroProp.quad, Temporal.Trend.ZeroProp.sing)
visreg::visreg(Temporal.Trend.ZeroProp.quad)
summary(Temporal.Trend.ZeroProp.quad)
summary(Temporal.Trend.ZeroProp.sing)
hist(rolling_metrics_oaks$p_zero_seeds)
co <- fixef(Temporal.Trend.ZeroProp.quad)$cond
b <- co["poly(Year, 2, raw = TRUE)1"]
c <- co["poly(Year, 2, raw = TRUE)2"]

Year_tip <- -b / (2 * c)
Year_tip

#library(segmented)

# m_lm <- lm(
#   logit_prop_zero_seeds ~ Year + log1p(mean_value_demand),
#   data = rolling_metrics_oaks
# )
# m_lm <- lm(
#   log_CVp_seeds ~ Year + log1p(mean_value_demand),
#   data = rolling_metrics_oaks
# )
# rolling_metrics_oaks$log_CVp_seeds
# m_seg <- segmented(
#   m_lm,
#   seg.Z = ~Year,
#   psi = list(Year = 2000) # initial guess
# )
# summary(m_seg)

ggplot(
  rolling_metrics_oaks,
  aes(x = logit_prop_zero_seeds, y = log_CVp_seeds)
) +
  geom_point()
summary(lm(log_CVp_seeds ~ logit_prop_zero_seeds, data = rolling_metrics_oaks))
summary(lm(log_CVp_seeds ~ p_zero_seeds, data = rolling_metrics_oaks))

hist(rolling_metrics_oaks$logit_prop_zero_seeds)
hist(rolling_metrics_oaks$log_CVp_seeds)
ggplot(
  rolling_metrics_oaks,
  aes(x = logit_prop_zero_seeds, y = log_CVp_seeds)
) +
  geom_point()
summary(Temporal.Trend.CVp.model1)
summary(Temporal.Trend.CVp.model2)

visreg::visreg(Temporal.Trend.CVp.model1)
visreg::visreg(Temporal.Trend.CVp.model2)
library(DHARMa)
res.mod1 = simulateResiduals(Temporal.Trend.CVp.model1)
res.mod2 = simulateResiduals(Temporal.Trend.CVp.model2)
plot(res.mod1)
plot(res.mod2)

################
#model check cues
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

visreg::visreg(Temporal.Trend.CVp.model.cues1)
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
  # geom_point(
  #   data = rolling_metrics_oaks,
  #   aes(x = mean_value_T_JunJulAug1, y = CVp_seeds),
  #   alpha = 0.1,
  #   size = 1
  # ) +
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
  # geom_point(
  #   data = rolling_metrics_oaks,
  #   aes(x = mean_value_T_Mar, y = CVp_seeds),
  #   alpha = 0.1,
  #   size = 1
  # ) +
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
  # geom_point(
  #   data = rolling_metrics_oaks,
  #   aes(x = mean_value_T_AprMay, y = CVp_seeds),
  #   alpha = 0.1,
  #   size = 1
  # ) +
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


# Temporal.Trend.Mean.Spring <- glmmTMB::glmmTMB(
#   mean_value_seeds ~
#     #Year_center + #if we include spring and year, it would flip the relation because of strong cor between the var
#     mean_value_Tspring +
#       log1p(mean_value_demand) +
#       (1 | sitenewname),
#   family = glmmTMB::tweedie(),
#   rolling_metrics_oaks
# )
# summary(Temporal.Trend.Mean.Spring)
#try to do a small SEM
#library(glmmTMB)
m_cue <- glmmTMB(
  mean_value_T_AprMay ~ Year_center + (1 | sitenewname),
  data = rolling_metrics_oaks
)
visreg::visreg(m_cue, ylab = "Partial Effect Spring Temp (Apr-May)")
m_cue_variab <- glmmTMB(
  sd_value_T_AprMay ~ Year_center + (1 | sitenewname),
  data = rolling_metrics_oaks
)
summary(m_cue_variab)
visreg::visreg(m_cue_variab)

m_mean_seed <- glmmTMB(
  log_mean_seeds ~ sd_value_T_AprMay + Year_center + (1 | sitenewname),
  data = rolling_metrics_oaks
)
visreg::visreg(m_mean_seed, ylab = "Partial Effect mean seed production (log)")
performance::model_performance(m_mean_seed)

m_mean_seed_year <- glmmTMB(
  log_mean_seeds ~ Year_center + (1 | sitenewname),
  data = rolling_metrics_oaks
)
summary(m_mean_seed_year)
library(MuMIn)
AICc(m_mean_seed_year, m_mean_seed)
visreg::visreg(m_mean_seed_year, ylab = "Partial Effect seed prod (log)")
m_cv <- glmmTMB(
  log_CVp_seeds ~
    log_mean_seeds +
      mean_value_T_AprMay +
      Year_center +
      (1 | sitenewname),
  data = rolling_metrics_oaks
)
summary(m_cv)
performance::model_performance(m_cv)

library(piecewiseSEM)

sem_model <- psem(
  m_cue,
  m_mean_seed,
  m_cv
)
plot(sem_model)
summary(sem_model)
summary(sem_model[[1]])
summary(sem_model[[2]])
summary(sem_model[[3]])
coefs <- coefs(sem_model, standardize = "scale")
coefs_clean <- coefs[, nzchar(names(coefs))]


library(igraph)
library(ggraph)

edges <- coefs_clean %>%
  dplyr::transmute(
    from = Predictor,
    to = Response,
    weight = Estimate
  )

nodes <- data.frame(
  name = unique(c(edges$from, edges$to))
)

g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
ggraph(g, layout = "tree") +
  geom_edge_link(
    aes(label = round(weight, 2), width = abs(weight)),
    arrow = arrow(length = unit(4, "mm")),
    end_cap = circle(4, "mm"),
    color = "grey"
  ) +
  geom_node_label(
    aes(label = name),
    size = 4,
    fill = "white"
  ) +
  theme_void()


#check for correlation
sync.oaks.rolling = calculate.sync(Oaks_Harvest, n.years = 10, step = 5) %>%
  mutate(Year_c = end.year - mean(end.year, na.rm = TRUE))

sync.oaks.rolling %>%
  ggplot(aes(x = end.year, y = mean.synch)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()


Temporal.Trend.Synchrony <- glmmTMB::glmmTMB(
  mean.synch ~
    end.year +
      (1 | NADL),
  family = glmmTMB::beta_family(link = "logit"),
  sync.oaks.rolling
)
summary(Temporal.Trend.Synchrony)


Temporal.Trend.Synchrony.model_rs <- glmmTMB(
  mean.synch ~
    Year_c +
      (Year_c | NADL),
  family = glmmTMB::beta_family(link = "logit"),
  data = sync.oaks.rolling
)
ranef_slopes_sync <- ranef(Temporal.Trend.Synchrony.model_rs)$cond$NADL
fixed_slope_sync <- fixef(Temporal.Trend.Synchrony.model_rs)$cond["Year_c"]
site_slopes_sync <- fixed_slope_sync + ranef_slopes_sync$Year_c
site_slopes_df_sync <- data.frame(
  sitenewname = rownames(ranef_slopes_sync),
  slope = site_slopes_sync
) %>%
  mutate(
    group_tempslope = case_when(slope <= 0 ~ "<=0", TRUE ~ ">0")
  )
ggplot(site_slopes_df_sync, aes(x = slope, fill = group_tempslope)) +
  geom_histogram(binwidth = 0.001, color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Distribution of Site-Specific Slopes \nfor Year Effect on sync (logit)",
    x = "Slope (logit scale)",
    y = "Count of Sites"
  )

formaps.sync = site_slopes_df_sync %>%
  rename(NADL = sitenewname) %>%
  left_join(
    Oaks_Harvest %>%
      dplyr::select(sitenewname, NADL, Latitude, Longitude) %>%
      distinct()
  ) %>%
  mutate(
    nadl_key = NADL %>%
      stri_trans_general("Latin-ASCII") %>% # AUGUSTÓW -> AUGUSTOW
      tolower() %>%
      trimws()
  ) %>%
  mutate(
    nadl_key = dplyr::recode(
      #still missing nowa deba ?
      nadl_key,
      'bielsk podlaski' = 'bielsk',
      "karniszewice" = "karnieszewice", #not exactly same area?
      'krosno odrzanskie' = 'krosno',
      'miedzyrzec podlaski' = 'miedzyrzec',
      "minsk mazowiecki" = "minsk",
      "piotrkow trybunalski" = "piotrkow",
      "polczyn zdroj" = "polczyn",
      "starogard gdanski" = "starogard"
    )
  )

nadl_map_sync <- nadl_sf2 %>%
  full_join(
    formaps.sync,
    by = "nadl_key"
  )

ggplot(nadl_map_sync) +
  geom_sf(aes(fill = slope), color = "black", linewidth = 0.1) +
  scale_fill_distiller(palette = "BrBG", na.value = "black") +
  theme_minimal() +
  geom_sf(
    data = around.poland,
    fill = "grey95",
    colour = "black",
    linewidth = .05
  ) +
  scale_y_continuous(
    breaks = c(50, 55),
    labels = function(b) paste0(b, "°N"),
    position = "right"
  ) +
  coord_sf(xlim = c(14, 26), ylim = c(49, 56), expand = FALSE) +
  scale_x_continuous(
    breaks = c(15, 20, 25),
    labels = function(b) paste0(b, "°E")
  ) +
  labs(x = NULL, y = NULL) +
  +hrbrthemes::theme_ipsum_rc(
    base_size = 14,
    axis_title_size = 16,
    plot_title_face = NULL,
    base_family = "roboto"
  ) +
  theme(plot.margin = margin(10, 90, 10, 10)) +
  theme(legend.position = "right") +
  theme(panel.spacing = unit(0.1, "cm")) +
  scale_fill_distiller(palette = "BrBG") +
  scale_color_distiller(palette = "BrBG")


test.rs.slope = formaps.sync %>%
  dplyr::select(NADL, slope) %>%
  rename(slope.sync = slope) %>%
  left_join(
    formaps.cv %>% dplyr::select(NADL, slope) %>% rename(slope.cv = slope),
    by = "NADL"
  )
cor(test.rs.slope$slope.sync, test.rs.slope$slope.cv, use = "complete.obs")
ggplot(test.rs.slope, aes(x = slope.sync, y = slope.cv)) +
  geom_point() +
  geom_smooth(method = 'lm')

#library(semEff)
#library(lme4)
#can not supprot glmmTMB
# sem_model_2 <- psem(
#   lmer(
#     mean_value_Tspring ~
#       Year_center +
#         (1 | sitenewname),
#     data = rolling_metrics_oaks,
#     REML = T
#   ),
#   lmer(
#     log_mean_seeds ~
#       mean_value_Tspring +
#         Year_center +
#         (1 | sitenewname),
#     data = rolling_metrics_oaks,
#     REML = T
#   ),
#   lmer(
#     log_CVp_seeds ~
#       log_mean_seeds +
#         mean_value_Tspring +
#         Year_center +
#         (1 | sitenewname),
#     data = rolling_metrics_oaks,
#     REML = T
#   )
# )
#
# piecewiseSEM:::plot.psem(
#   piecewiseSEM::as.psem(sem_model_2),
#   node_attrs = data.frame(
#     shape = "rectangle",
#     color = "black",
#     fillcolor = "grey"
#   ),
#   layout = "tree"
# )
#
# system.time(
#   keeley.sem.boot <- bootEff(
#     sem_model_2,
#     R = 100,
#     seed = 13,
#     ran.eff = "sitenewname"
#   )
# )
#
# (keeley.sem.eff <- semEff(keeley.sem.boot))
# summary(keeley.sem.eff)
#
# ggplot(rolling_metrics_oaks, aes(x = log_CVp_seeds, y = log_mean_seeds)) +
#   geom_point()
# #think about some detreding to remove the decline seed production
# detrended.rolling_metrics_oaks <- rolling_metrics_oaks %>%
#   group_by(sitenewname) %>%
#   mutate(
#     CV_resid_log = resid(lm(log_CVp_seeds ~ log_mean_seeds)),
#     CV_resid = resid(lm(CVp_seeds ~ mean_value_seeds))
#   ) %>%
#   ungroup()

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
library(purrr)

wins <- seq(6, 20, by = 2)

res_all <- map_dfr(wins, function(w) {
  fit_cue_models_over_window(Oaks_Harvest, window = w, step = 5)$cue_results
})

res_all
library(ggplot2)

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

#check cue trends
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
  mutate(log_mean_value_demand = log1p(mean_value_demand))

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


plot(ggeffects::ggpredict(
  Temporal.Trend.CVp.model.cuesJJA.beech,
  terms = "mean_value_T_JunJul1"
))

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

beech.fig.5 = p_JJ1_beech +
  scale_x_continuous(
    breaks = c(16, 18, 20, 22, 24, 26, 28),
    labels = c(16, 18, 20, 22, 24, 26, 28)
  ) +
  theme(
    plot.margin = margin(7, 7, 7, 7)
  ) +
  theme(legend.position = c(.8, .8))
#beech.cue.JJA1 +
#theme(
#  plot.margin = margin(7, 7, 7, 7)
#) +

beech.fig.5
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
