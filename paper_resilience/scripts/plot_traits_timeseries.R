library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(readr)
library(purrr)
library(tidyr)
library(zoo)
library(moments)
library(viridis)


# Read the CSV file
table_2y <- read.csv("~/Desktop/CAETE-DVM-alloc-allom-including_alloc2_Cm2/paper_resilience/tables_results/PLS_alive_traits_MAN_30prec_2y.csv") %>%
  mutate(Scenario = "2y")

table_4y <- read.csv("~/Desktop/CAETE-DVM-alloc-allom-including_alloc2_Cm2/paper_resilience/tables_results/PLS_alive_traits_MAN_30prec_4y.csv") %>%
  mutate(Scenario = "4y")

table_6y <- read.csv("~/Desktop/CAETE-DVM-alloc-allom-including_alloc2_Cm2/paper_resilience/tables_results/PLS_alive_traits_MAN_30prec_6y.csv") %>%
  mutate(Scenario = "6y")

table_8y <- read.csv("~/Desktop/CAETE-DVM-alloc-allom-including_alloc2_Cm2/paper_resilience/tables_results/PLS_alive_traits_MAN_30prec_8y.csv") %>%
  mutate(Scenario = "8y")

table_regclim <- read.csv("~/Desktop/CAETE-DVM-alloc-allom-including_alloc2_Cm2/paper_resilience/tables_results/PLS_alive_traits_MAN_regularclimate.csv") %>%
  mutate(Scenario = "regclim")


# ---- Combine scenarios and set a consistent order for plotting ----
df <- bind_rows(table_regclim, table_2y, table_4y, table_6y, table_8y) %>%
  mutate(Scenario = factor(Scenario, levels = c("regclim", "2y", "4y", "6y", "8y")))

#  Change SLA unit  (*1000 to m2/kg):
df <- df %>% mutate(sla_random = sla_random * 1000)

# ---- Transform the year in run year, instead of the date itself ----
# ---- Build a YEAR -> run_year index PER SCENARIO (1,2,3,...) ----
year_map <- df %>%
  distinct(Scenario, YEAR) %>%
  arrange(Scenario, YEAR) %>%
  group_by(Scenario) %>%
  mutate(run_year = row_number()) %>%
  ungroup()


cols <- c(
  "regclim" = "#1b9e77", 
  "2y"      = "#D55E00",
  "4y"      = "#e6ab02",
  "6y"      = "#e7298a", 
  "8y"      = "#7570b3" 
)

# --- MEAN TEMPORAL SERIES--- #
# Here the mean value is calculated using the weight of occupation. Otherwise,
# the mean would be "puxada" to become bigger or lower depending on the rarity
# of strategies

# =====================================================================================
# FUNCTION: compute weighted temporal mean for a given trait (returns NA if sum(OC)==0)
# =====================================================================================
compute_weighted_mean <- function(df, trait_col, year_map) {
  # df        : your data frame (with Scenario, YEAR, OC and traits)
  # trait_col : the name of the trait column as string (e.g. "wd_random")
  # year_map  : table that links YEAR -> run_year per Scenario
  
  df %>%
    group_by(Scenario, YEAR) %>%
    summarise(
      sum_OC = sum(OC, na.rm = TRUE),
      weighted_mean = ifelse(
        is.finite(sum_OC) & sum_OC > 0,
        sum(.data[[trait_col]] * OC, na.rm = TRUE) / sum_OC,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    left_join(year_map, by = c("Scenario", "YEAR")) %>%
    arrange(Scenario, run_year)
}

# Compute means
wd_means  <- compute_weighted_mean(df, "wd_random", year_map)
sla_means <- compute_weighted_mean(df, "sla_random", year_map)
g1_means  <- compute_weighted_mean(df, "g1", year_map)

# =====================================================================================
# FUNCTION: plot weighted temporal mean for a given trait
# =====================================================================================
# =====================================================================================
# FUNCTION: plot weighted temporal mean for a given trait, with optional scenario filter
# =====================================================================================
plot_weighted_mean <- function(means_df, trait_label, y_lab, scenarios = NULL) {
  # means_df   : dataframe gerado pela compute_weighted_mean()
  # trait_label: string para usar como título (ex: "WD", "SLA", "g1")
  # y_lab      : rótulo do eixo Y (expression ou string)
  # scenarios  : vetor opcional com cenários a incluir (ex: c("2y", "8y", "regclim"))
  
  if (!is.null(scenarios)) {
    means_df <- means_df %>% dplyr::filter(Scenario %in% scenarios)
  }
  
  ggplot(means_df, aes(x = run_year, y = weighted_mean, color = Scenario)) +
    geom_line(linewidth = 1.2, na.rm = TRUE) +
    scale_color_manual(values = cols) +
    scale_x_continuous(breaks = seq(0, max(means_df$run_year, na.rm = TRUE), by = 5)) +
    labs(x = "Year", 
         y = y_lab, 
         color = "Scenario",
         title = "") +
    theme_minimal(base_size = 16) +
    theme(legend.position = "bottom")
}

# Plot all scenarios
p_wd_all <- plot_weighted_mean(wd_means, "WD", expression("WD weighted mean (g/cm"^3*")"))
p_sla_all <- plot_weighted_mean(sla_means, "SLA", expression("SLA weighted mean (m"^2*"/kg)"))
p_g1_all  <- plot_weighted_mean(g1_means,  "g1",  expression("g"[1]*" weighted mean (kPa"^{-0.5}*")"))

# Plot subset (only 2y, 8y and regclim)
p_wd_subset <- plot_weighted_mean(
  wd_means, "WD", expression("WD weighted mean (g/cm"^3*")"),
  scenarios = c("2y", "8y", "regclim")
)

p_sla_subset <- plot_weighted_mean(
  sla_means, "SLA", expression("SLA weighted mean (m"^2*"/kg)"),
  scenarios = c("2y", "8y", "regclim")
)

p_g1_subset <- plot_weighted_mean(
  g1_means, "g1",  expression("g"[1]*" weighted mean (kPa"^{-0.5}*")"),
  scenarios = c("2y", "8y", "regclim")
)


# Show plot
p_wd_all
p_sla_all
p_g1_all
p_wd_subset
p_sla_subset
p_g1_subset

#grid.arrange(p_wd, p_sla, p_g1, ncol = 2, nrow = 2)

ggplot(df, aes(x = factor(YEAR), y = wd_random, color = Scenario)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +   # boxplot sem outliers automáticos
  geom_point(aes(size = OC), alpha = 0.6, position = position_jitter(width = 0.2)) +
  scale_size_continuous(range = c(0.1, 4)) +
  facet_wrap(~Scenario, scales = "fixed") +
  labs(x = "Year", y = "Wood density (g/cm³)", size = "Occupancy") +
  theme_minimal(base_size = 14)

library(ggridges)

ggplot(df, aes(x = wd_random, y = factor(YEAR), height = ..density.., fill = Scenario)) +
  geom_density_ridges(aes(weight = OC), scale = 3, alpha = 0.6, rel_min_height = 0.01) +
  facet_wrap(~Scenario, scales = "free_x") +
  labs(x = "Wood density (g/cm³)", y = "Year") +
  theme_minimal(base_size = 14)

library(dplyr)
library(Hmisc)
library(ggplot2)

# ---- Weighted quantiles per Scenario–YEAR (robust to empty groups) ----
wd_quant <- df %>%
  group_by(Scenario, YEAR) %>%
  summarise(
    # vectors for this group
    q05 = {
      x  <- wd_random
      w  <- OC
      ok <- is.finite(x) & is.finite(w) & w > 0
      x  <- x[ok]; w <- w[ok]
      if (length(x) == 0 || sum(w) == 0) NA_real_ else {
        w <- w / sum(w)  # normalize weights
        as.numeric(Hmisc::wtd.quantile(x, weights = w, probs = 0.05, na.rm = TRUE))
      }
    },
    q50 = {
      x  <- wd_random
      w  <- OC
      ok <- is.finite(x) & is.finite(w) & w > 0
      x  <- x[ok]; w <- w[ok]
      if (length(x) == 0 || sum(w) == 0) NA_real_ else {
        w <- w / sum(w)
        as.numeric(Hmisc::wtd.quantile(x, weights = w, probs = 0.50, na.rm = TRUE))
      }
    },
    q95 = {
      x  <- wd_random
      w  <- OC
      ok <- is.finite(x) & is.finite(w) & w > 0
      x  <- x[ok]; w <- w[ok]
      if (length(x) == 0 || sum(w) == 0) NA_real_ else {
        w <- w / sum(w)
        as.numeric(Hmisc::wtd.quantile(x, weights = w, probs = 0.95, na.rm = TRUE))
      }
    },
    .groups = "drop"
  ) %>%
  left_join(year_map, by = c("Scenario", "YEAR")) %>%
  arrange(Scenario, run_year)

# ---- Plot: ribbon (5–95%) + median, faceted by scenario ----
ggplot(wd_quant, aes(x = run_year, color = Scenario, fill = Scenario)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.20, color = NA) +
  geom_line(aes(y = q50), linewidth = 1.2) +
  facet_wrap(~ Scenario, ncol = 1, scales = "fixed") +
  scale_x_continuous(breaks = seq(1, max(wd_quant$run_year, na.rm = TRUE), by = 5)) +
  labs(x = "Running year", y = "Wood density (g/cm³)",
       title = "WD: weighted median (line) and 5–95% span (ribbon)") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

# --- Heatmap of occupancy in Year × WD space (WD) ---
library(ggplot2)

ggplot(df, aes(x = YEAR, y = wd_random)) +
  stat_bin2d(aes(weight = OC), bins = 40) +   # OC used as weight
  facet_wrap(~ Scenario, ncol = 1, scales = "fixed") +
  scale_fill_viridis_c(option = "C", trans = "log10", name = "Weighted count") +
  labs(x = "Running year", y = "Wood density (g/cm³)",
       title = "WD: occupancy-weighted density over time (log scale)") +
  theme_minimal(base_size = 15)



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
  group_by(Scenario, YEAR) %>%
  summarise(
    wd_min = min(wd_random, na.rm = TRUE),
    wd_max = max(wd_random, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(wd_amp = wd_max - wd_min) %>%
  left_join(year_map, by = c("Scenario", "YEAR")) %>%
  arrange(Scenario, run_year)

# Plot
p_wd_amp<-ggplot(wd_amp, aes(x = run_year, y = wd_amp, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cols) +
  labs(x = "Year", y = "WD amplitude",
       color = "Scenario") +
  theme_minimal(base_size = 16) +
  facet_wrap(~Scenario, scales = "fixed")

p_wd_amp 

# Compute sla amplitude per scenario-year (no weights)
sla_amp <- df %>%
  filter(OC > 0) %>%   # only alive PLS
  group_by(Scenario, YEAR) %>%
  summarise(
    sla_min = min(sla_random, na.rm = TRUE),
    sla_max = max(sla_random, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sla_amp = sla_max - sla_min) %>%
  left_join(year_map, by = c("Scenario", "YEAR")) %>%
  arrange(Scenario, run_year)

# Plot
p_sla_amp<-ggplot(sla_amp, aes(x = run_year, y = sla_amp, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "SLA amplitude",
       color = Scenario) +
  theme_minimal(base_size = 16) +
  facet_wrap(~Scenario, scales = "fixed")

p_sla_amp 


# Compute g1 amplitude per scenario-year (no weights)
g1_amp <- df %>%
  filter(OC > 0) %>%   # only alive PLS
  group_by(Scenario, YEAR) %>%
  summarise(
    g1_min = min(g1, na.rm = TRUE),
    g1_max = max(g1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(g1_amp = g1_max - g1_min) %>%
  left_join(year_map, by = c("Scenario", "YEAR")) %>%
  arrange(Scenario, run_year)

# Plot
p_g1_amp<-ggplot(g1_amp, aes(x = run_year, y = g1_amp, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "g1 amplitude (unweighted)") +
  theme_minimal(base_size = 16) +
  facet_wrap(~Scenario, scales = "free_y")

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
  group_by(Scenario, YEAR) %>%
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
  left_join(year_map, by = c("Scenario", "YEAR")) %>%
  arrange(Scenario, run_year)

sla_stats <- df %>%
  group_by(Scenario, YEAR) %>%
  summarise(
    sum_OC = sum(OC, na.rm = TRUE),
    sla_mean = ifelse(sum_OC > 0, sum(sla_random * OC, na.rm = TRUE) / sum_OC, NA_real_),
    sla_var  = ifelse(sum_OC > 0,
                     sum(OC * (sla_random - sum(sla_random * OC, na.rm = TRUE) / sum_OC)^2, na.rm = TRUE) / sum_OC,
                     NA_real_),
    sla_min  = ifelse(sum_OC > 0, min(sla_random, na.rm = TRUE), NA_real_),
    sla_max  = ifelse(sum_OC > 0, max(sla_random, na.rm = TRUE), NA_real_),
    sla_skew = ifelse(sum_OC > 0,
                     {
                       mu <- sum(sla_random * OC, na.rm = TRUE) / sum_OC
                       m2 <- sum(OC * (sla_random - mu)^2, na.rm = TRUE) / sum_OC
                       m3 <- sum(OC * (sla_random - mu)^3, na.rm = TRUE) / sum_OC
                       if (m2 > 0) m3 / (m2^(3/2)) else NA_real_
                     },
                     NA_real_)
  ) %>%
  ungroup() %>%
  left_join(year_map, by = c("Scenario", "YEAR")) %>%
  arrange(Scenario, run_year)


# ---- Add amplitude (max - min) ----
wd_stats <- wd_stats %>%
  mutate(wd_amp = wd_max - wd_min)

sla_stats <- sla_stats %>%
  mutate(sla_amp = sla_max - sla_min)


p_wd_var <- ggplot(wd_stats, aes(x = run_year, y = wd_var, color = Scenario)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "WD variance", color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_wd_var

p_sla_var <- ggplot(sla_stats, aes(x = run_year, y = sla_var, color = Scenario)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "sla variance", color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_sla_var

p_wd_amp <- ggplot(wd_stats, aes(x = run_year, y = wd_amp, color = Scenario)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "WD amplitude", color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_wd_amp

p_sla_amp <- ggplot(sla_stats, aes(x = run_year, y = sla_amp, color = Scenario)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "sla amplitude", color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_sla_amp

p_wd_skew <- ggplot(wd_stats, aes(x = run_year, y = wd_skew, color = Scenario)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "WD skewness", color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_wd_skew

p_sla_skew <- ggplot(sla_stats, aes(x = run_year, y = sla_skew, color = Scenario)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_manual(values = cols) +
  labs(x = "Running year", y = "sla skewness", color = "Scenario") +
  theme_minimal(base_size = 16) + theme(legend.position = "bottom")
p_sla_skew


###COMPARAR "O QUE EXISTE DO QUE DOMINA" (CHAT GPT)


library(TPD)

df_first <- df %>%
  filter(YEAR == 1979, Scenario == "2y")

df_last <- df %>%
  filter(YEAR == 2004, Scenario == "2y")

df2 = rbind(df_first, df_last)

TPD_2y = TPDs(species = df2$YEAR , traits = df2$wd_ocp, alpha = 1)

plotTPD(TPD_2y
        )
library(ggplot2)

plot_box_start_end <- function(df, trait_col, scenario, y_lab) {
  d <- df %>% filter(Scenario == scenario)
  years <- range(d$YEAR, na.rm = TRUE)   # start & end
  d2 <- d %>% filter(YEAR %in% years)
  
  ggplot(d2, aes(x = factor(YEAR), y = .data[[trait_col]])) +
    geom_boxplot(outlier.shape = NA, fill = "grey85", alpha = 0.5) +
    geom_jitter(aes(size = OC), width = 0.15, alpha = 0.7) +
    scale_size_continuous(name = "Occupancy", range = c(1,6)) +
    labs(x = "Year", y = y_lab,
         title = paste0(scenario, " — start vs end")) +
    theme_minimal(base_size = 14)
}

plot_box_start_end(df, "wd_random", "4y", expression("Wood density (g/cm"^3*")"))

library(Hmisc)

summarise_quantiles <- function(df, trait_col, scenario) {
  d <- df %>% filter(Scenario == scenario)
  years <- range(d$YEAR, na.rm = TRUE)
  d %>% filter(YEAR %in% years) %>%
    group_by(YEAR) %>%
    summarise(
      mean_w = sum(.data[[trait_col]] * OC, na.rm = TRUE) / sum(OC, na.rm = TRUE),
      q05 = as.numeric(wtd.quantile(.data[[trait_col]], weights = OC, probs = 0.05, na.rm = TRUE)),
      q95 = as.numeric(wtd.quantile(.data[[trait_col]], weights = OC, probs = 0.95, na.rm = TRUE)),
      .groups = "drop"
    )
}

q_df <- summarise_quantiles(df, "wd_random", "4y")

ggplot(q_df, aes(x = factor(YEAR), y = mean_w)) +
  geom_point(size = 3, color = "darkred") +
  geom_errorbar(aes(ymin = q05, ymax = q95), width = 0.2, color = "darkred") +
  labs(x = "Year", y = expression("Wood density (g/cm"^3*")")) +
  theme_minimal(base_size = 14)

ggplot(d2, aes(x = factor(YEAR), y = .data[[trait_col]], fill = factor(YEAR))) +
  geom_violin(alpha = 0.4) +
  geom_jitter(aes(size = OC), width = 0.1, alpha = 0.7) +
  scale_size_continuous(name = "Occupancy", range = c(1,6)) +
  labs(x = "Year", y = y_lab, fill = "Year") +
  theme_minimal(base_size = 14)

library(dplyr)
library(ggplot2)

plot_violin_start_end <- function(df, trait_col, scenario, y_lab) {
  # Filtra o cenário desejado
  d <- df %>% filter(Scenario == scenario)
  # Pega ano inicial e final
  years <- range(d$YEAR, na.rm = TRUE)
  # Subconjunto com só esses anos
  d2 <- d %>% filter(YEAR %in% years)
  
  ggplot(d2, aes(x = factor(YEAR), y = .data[[trait_col]], fill = factor(YEAR))) +
    geom_violin(alpha = 0.4, trim = FALSE) +  # trim = FALSE para mostrar toda a faixa
    geom_jitter(aes(size = OC), width = 0.1, alpha = 0.7) +
    scale_size_continuous(name = "Occupancy", range = c(1,6)) +
    labs(x = "Year", y = y_lab, fill = "Year",
         title = paste0(scenario, " — start vs end")) +
    theme_minimal(base_size = 14)
}

# Exemplo para WD no cenário 2y
plot_violin_start_end(df, "wd_random", "2y",
                      expression("Wood density (g/cm"^3*")"))


library(Hmisc)

cwm_df <- d2 %>%
  group_by(YEAR) %>%
  summarise(
    mean_w = sum(wd_random * OC, na.rm = TRUE) / sum(OC, na.rm = TRUE)
  )



ggplot(d2, aes(x = factor(YEAR), y = wd_random, fill = factor(YEAR))) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_jitter(aes(size = OC), width = 0.1, alpha = 0.7) +
  geom_point(data = cwm_df, aes(x = factor(YEAR), y = mean_w),
             color = "black", size = 3, shape = 18) +
  scale_size_continuous(name = "Occupancy", range = c(1,6)) +
  labs(x = "Year", y = expression("Wood density (g/cm"^3*")"), fill = "Year") +
  theme_minimal(base_size = 14)

library(dplyr)
library(ggplot2)
library(Hmisc)

plot_violin_start_end <- function(df, trait_col, scenario, y_lab,
                                  start_year = NULL, end_year = NULL) {
  # Filtra o cenário desejado
  d <- df %>% filter(Scenario == scenario)
  
  # Se não especificar, usa mínimo e máximo
  if (is.null(start_year)) start_year <- min(d$YEAR, na.rm = TRUE)
  if (is.null(end_year))   end_year   <- max(d$YEAR, na.rm = TRUE)
  
  # Subconjunto com só os anos escolhidos e apenas PLS com OC > 0
  d2 <- d %>%
    filter(YEAR %in% c(start_year, end_year), OC > 0)
  
  # Calcular médias ponderadas (CWM) por ano
  cwm_df <- d2 %>%
    group_by(YEAR) %>%
    summarise(
      mean_w = sum(.data[[trait_col]] * OC, na.rm = TRUE) / sum(OC, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(d2, aes(x = factor(YEAR), y = .data[[trait_col]], fill = factor(YEAR))) +
    geom_violin(alpha = 0.4, trim = FALSE) +
    geom_jitter(aes(size = OC), width = 0.1, alpha = 0.7) +
    geom_point(data = cwm_df, aes(x = factor(YEAR), y = mean_w),
               color = "black", size = 3, shape = 18) +
    scale_size_continuous(name = "Occupancy", range = c(1,6)) +
    labs(x = "Year", y = y_lab, fill = "Year",
         title = paste0(scenario, " — start (", start_year, ") vs end (", end_year, ")")) +
    theme_minimal(base_size = 14)
}

plot_violin_start_end(
  df, "wd_random", "2y",
  y_lab = expression("Wood density (g/cm"^3*")"),
  start_year = 1979, end_year = 2004
)


library(dplyr)
library(ggplot2)
library(Hmisc)

plot_violin_start_end <- function(df, trait_col, scenario, y_lab,
                                  start_year = NULL, end_year = NULL) {
  # Filtra o cenário desejado
  d <- df %>% filter(Scenario == scenario)
  
  # Se não especificar, usa mínimo e máximo
  if (is.null(start_year)) start_year <- min(d$YEAR, na.rm = TRUE)
  if (is.null(end_year))   end_year   <- max(d$YEAR, na.rm = TRUE)
  
  # Subconjunto com só os anos escolhidos e apenas PLS com OC > 0
  d2 <- d %>%
    filter(YEAR %in% c(start_year, end_year), OC > 0)
  
  # Calcular médias ponderadas (CWM) por ano
  cwm_df <- d2 %>%
    group_by(YEAR) %>%
    summarise(
      mean_w = sum(.data[[trait_col]] * OC, na.rm = TRUE) / sum(OC, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(d2, aes(x = factor(YEAR), y = .data[[trait_col]], fill = factor(YEAR))) +
    geom_violin(alpha = 0.4, trim = FALSE) +
    geom_jitter(aes(size = OC), width = 0.1, alpha = 0.7) +
    # Linha horizontal na média ponderada
    geom_crossbar(data = cwm_df, aes(x = factor(YEAR), y = mean_w, ymin = mean_w, ymax = mean_w),
                  width = 0.4, fatten = 2, color = "black") +
    scale_size_continuous(name = "Occupancy", range = c(1,6)) +
    labs(x = "Year", y = y_lab, fill = "Year",
         title = paste0(scenario, " — start (", start_year, ") vs end (", end_year, ")")) +
    theme_minimal(base_size = 14)
}
plot_violin_start_end(
  df, "wd_random", "2y",
  y_lab = expression("Wood density (g/cm"^3*")"),
  start_year = 1979, end_year = 2004
)

plot_violin_start_end(
  df, "sla_random", "2y",
  y_lab = expression("Wood density (g/cm"^3*")"),
  start_year = 1979, end_year = 2004
)
plot_violin_start_end(
  df, "g1", "2y",
  y_lab = expression("Wood density (g/cm"^3*")"),
  start_year = 1979, end_year = 2004
)


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
