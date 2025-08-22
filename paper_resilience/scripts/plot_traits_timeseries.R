library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(readr)
library(purrr)
library(tidyr)
library(zoo)
library(moments)


# Read the CSV file
table_2y <- read.csv("~/Desktop/CAETE-DVM-alloc-allom-including_alloc2_Cm2/paper_resilience/tables_results/PLS_alive_traits_MAN_30prec_2y.csv") %>%
  mutate(Source = "2y")

table_8y <- read.csv("~/Desktop/CAETE-DVM-alloc-allom-including_alloc2_Cm2/paper_resilience/tables_results/PLS_alive_traits_MAN_30prec_8y.csv") %>%
  mutate(Source = "8y")

table_regclim <- read.csv("~/Desktop/CAETE-DVM-alloc-allom-including_alloc2_Cm2/paper_resilience/tables_results/PLS_alive_traits_MAN_regularclimate.csv") %>%
  mutate(Source = "regclim")


# ---- Combine scenarios and set a consistent order for plotting ----
df <- bind_rows(table_regclim, table_2y, table_8y) %>%
  mutate(Source = factor(Source, levels = c("regclim", "2y", "8y")))

# ---- Transform the year in run year, instead of the date itself ----
# ---- Build a YEAR -> run_year index PER SCENARIO (1,2,3,...) ----
year_map <- df %>%
  distinct(Source, YEAR) %>%
  arrange(Source, YEAR) %>%
  group_by(Source) %>%
  mutate(run_year = row_number()) %>%
  ungroup()

# ---- COLORS: set colors for each scenario (for plotting) ----
cols <- c("regclim" = "steelblue3", "2y" = "lightsalmon2", "8y" = "seagreen4")

# --- MEAN TEMPORAL SERIES--- #
# Here the mean value is calculated using the weight of occupation. Otherwise,
# the mean would be "puxada" to become bigger or lower depending on the rarity
# of strategies
# =====================================================================================
# TRAIT 1 — Wood density
# =====================================================================================

# Compute weighted mean per year and scenario (returns NA if sum(OC)==0)
wd_means <- df %>%
  group_by(Source, YEAR) %>%
  summarise(
    sum_OC = sum(OC, na.rm = TRUE),
    weighted_mean_wd = ifelse(
      is.finite(sum_OC) & sum_OC > 0,
      sum(wd_random * OC, na.rm = TRUE) / sum_OC,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  left_join(year_map, by = c("Source", "YEAR")) %>%
  arrange(Source, run_year)

# Plot WD weighted mean
p_wd <- ggplot(wd_means, aes(x = run_year, y = weighted_mean_wd, color = Source)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  scale_x_continuous(breaks = seq(0, max(wd_means$run_year, na.rm = TRUE), by = 5)) +
  labs(x = "Year", y = expression("WD weighted mean (g/cm"^3*")"), color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_wd

# =====================================================================================
# TRAIT 2 — Specific Leaf Area
# =====================================================================================
# * 1000 is to change the unit to m2/kg
sla_means <- df %>%
  group_by(Source, YEAR) %>%
  summarise(
    sum_OC = sum(OC, na.rm = TRUE),
    weighted_mean_sla = ifelse(
      is.finite(sum_OC) & sum_OC > 0,
      sum(sla_random * OC *1000, na.rm = TRUE) / sum_OC,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  left_join(year_map, by = c("Source", "YEAR")) %>%
  arrange(Source, run_year)

p_sla <- ggplot(sla_means, aes(x = run_year, y = weighted_mean_sla, color = Source)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  scale_x_continuous(breaks = seq(0, max(wd_means$run_year, na.rm = TRUE), by = 5)) +
  labs(x = "Year", y = expression("SLA weighted mean (m"^"2"*"/kg)"), color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_sla

# =====================================================================================
# TRAIT 3 — g1 (g1)
# =====================================================================================
g1_means <- df %>%
  group_by(Source, YEAR) %>%
  summarise(
    sum_OC = sum(OC, na.rm = TRUE),
    weighted_mean_g1 = ifelse(
      is.finite(sum_OC) & sum_OC > 0,
      sum(g1 * OC, na.rm = TRUE) / sum_OC,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  left_join(year_map, by = c("Source", "YEAR")) %>%
  arrange(Source, run_year)

p_g1 <- ggplot(g1_means, aes(x = run_year, y = weighted_mean_g1, color = Source)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  scale_x_continuous(breaks = seq(0, max(wd_means$run_year, na.rm = TRUE), by = 5)) +
  labs(x = "Year", y = expression("g1 weighted mean (kPa"^{0.5}*")"), color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_g1

#grid.arrange(p_wd, p_sla, p_g1, ncol = 2, nrow = 2)


# --- Complementary analysis (variance, amplitude and skewness) --- #
# Variance (spread), amplitude (min–max range), and skewness (shape) 
# of each trait per year, weighted by occupancy (OC).

# --- AMPLITUDE TEMPORAL SERIES--- #
# Here the amplitude is not weighted once it is only the range of values that exist.
# Even if a PLS occupies only 0.1% of the gridcell, it still defines the min and max value

# --- a) Non normalized
# Compute WD amplitude per scenario-year (no weights)
wd_amp <- df %>%
  filter(OC > 0) %>%   # only alive PLS
  group_by(Source, YEAR) %>%
  summarise(
    wd_min = min(wd_random, na.rm = TRUE),
    wd_max = max(wd_random, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(wd_amp = wd_max - wd_min) %>%
  left_join(year_map, by = c("Source", "YEAR")) %>%
  arrange(Source, run_year)

# Plot
p_wd_amp<-ggplot(wd_amp, aes(x = run_year, y = wd_amp, color = Source)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cols) +
  labs(x = "Year", y = "WD amplitude",
       color = "Scenario") +
  theme_minimal(base_size = 16) +
  facet_wrap(~Source, scales = "fixed")

p_wd_amp 

# Compute sla amplitude per scenario-year (no weights)
sla_amp <- df %>%
  filter(OC > 0) %>%   # only alive PLS
  group_by(Source, YEAR) %>%
  summarise(
    sla_min = min(sla_random, na.rm = TRUE),
    sla_max = max(sla_random, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sla_amp = sla_max - sla_min) %>%
  left_join(year_map, by = c("Source", "YEAR")) %>%
  arrange(Source, run_year)

# Plot
p_sla_amp<-ggplot(sla_amp, aes(x = run_year, y = sla_amp, color = Source)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "SLA amplitude",
       color = Scenario) +
  theme_minimal(base_size = 16) +
  facet_wrap(~Source, scales = "fixed")

p_sla_amp 


# Compute g1 amplitude per scenario-year (no weights)
g1_amp <- df %>%
  filter(OC > 0) %>%   # only alive PLS
  group_by(Source, YEAR) %>%
  summarise(
    g1_min = min(g1, na.rm = TRUE),
    g1_max = max(g1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(g1_amp = g1_max - g1_min) %>%
  left_join(year_map, by = c("Source", "YEAR")) %>%
  arrange(Source, run_year)

# Plot
p_g1_amp<-ggplot(g1_amp, aes(x = run_year, y = g1_amp, color = Source)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "g1 amplitude (unweighted)") +
  theme_minimal(base_size = 16) +
  facet_wrap(~Source, scales = "free_y")

p_g1_amp 


# --- b) Normalized

# ---  min/max for 0–1 rescaling ---
wd_min <- min(df$wd_random, na.rm = TRUE)
wd_max <- max(df$wd_random, na.rm = TRUE)
sla_min <- min(df$sla_random, na.rm = TRUE)
sla_max <- max(df$sla_random, na.rm = TRUE)
g1_min <- min(df$g1, na.rm = TRUE)
g1_max <- max(df$g1, na.rm = TRUE)

# 0–1 rescaled versions (clip is optional; usually not needed if mins/maxes come from df)
df <- df %>%
  dplyr::mutate(
    wd_01  = (wd_random  - wd_min)  / (wd_max  - wd_min),
    sla_01 = (sla_random - sla_min) / (sla_max - sla_min),
    g1_01  = (g1         - g1_min)  / (g1_max  - g1_min)
  )


# =====================================================================================
# Variance
# =====================================================================================


# ---- Compute WD stats per year and scenario ----
wd_stats <- df %>%
  group_by(Source, YEAR) %>%
  summarise(
    sum_OC = sum(OC, na.rm = TRUE),
    wd_mean = ifelse(sum_OC > 0, sum(wd_random * OC, na.rm = TRUE) / sum_OC, NA_real_),
    wd_var  = ifelse(sum_OC > 0,
                     sum(OC * (wd_random - sum(wd_random * OC, na.rm = TRUE) / sum_OC)^2, na.rm = TRUE) / sum_OC,
                     NA_real_),
    wd_min  = ifelse(sum_OC > 0, min(wd_random, na.rm = TRUE), NA_real_),
    wd_max  = ifelse(sum_OC > 0, max(wd_random, na.rm = TRUE), NA_real_),
    wd_skew = ifelse(sum_OC > 0,
                     {
                       mu <- sum(wd_random * OC, na.rm = TRUE) / sum_OC
                       m2 <- sum(OC * (wd_random - mu)^2, na.rm = TRUE) / sum_OC
                       m3 <- sum(OC * (wd_random - mu)^3, na.rm = TRUE) / sum_OC
                       if (m2 > 0) m3 / (m2^(3/2)) else NA_real_
                     },
                     NA_real_)
  ) %>%
  ungroup() %>%
  left_join(year_map, by = c("Source", "YEAR")) %>%
  arrange(Source, run_year)

# ---- Add amplitude (max - min) ----
wd_stats <- wd_stats %>%
  mutate(wd_amp = wd_max - wd_min)

p_wd_var <- ggplot(wd_stats, aes(x = run_year, y = wd_var, color = Source)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "WD variance", color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_wd_var

p_wd_amp <- ggplot(wd_stats, aes(x = run_year, y = wd_amp, color = Source)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "WD amplitude", color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_wd_amp

p_wd_skew <- ggplot(wd_stats, aes(x = run_year, y = wd_skew, color = Source)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "WD skewness", color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_wd_skew


###COMPARAR "O QUE EXISTE DO QUE DOMINA" (CHAT GPT)
# write.csv(combined_data_g1, "/home/bianca/bianca/CAETE-DVM-alloc-allom/scripts/traits/g1_weighted_mean.csv")
# write.csv(combined_data_wd, "/home/bianca/bianca/CAETE-DVM-alloc-allom/scripts/traits/wd_weighted_mean.csv")
# write.csv(combined_data_sla, "/home/bianca/bianca/CAETE-DVM-alloc-allom/scripts/traits/sla_weighted_mean.csv")
# 
# combined_data_regclim = cbind(mean_g1_regclim, mean_sla_regclim, mean_wd_regclim)
# combined_data_8y = cbind(mean_g1_8y,mean_sla_8y,mean_wd_8y)
# combined_data_2y = cbind(mean_g1_2y,mean_sla_2y,mean_wd_2y)
# 
# write.csv(combined_data_regclim, "/home/bianca/bianca/CAETE-DVM-alloc-allom/scripts/traits/traits_weightedmean_regclim.csv")
# write.csv(combined_data_8y, "/home/bianca/bianca/CAETE-DVM-alloc-allom/scripts/traits/traits_weightedmean_8y.csv")
# write.csv(combined_data_2y, "/home/bianca/bianca/CAETE-DVM-alloc-allom/scripts/traits/traits_weightedmean_2y.csv")
