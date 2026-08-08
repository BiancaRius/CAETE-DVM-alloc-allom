# ==============================================================================
# Recovery analysis of biogeochemical processes under recurrent imposed drought
# ==============================================================================

# This script evaluates the recovery dynamics of multiple biogeochemical
# processes simulated by CAETÊ under recurrent drought experiments.

# The drought experiments impose a fixed precipitation reduction (30%) on top of the
# original climate at different recurrence intervals. The analysis is designed to
# compare each drought scenario against a paired control simulation driven by the
# same background climate but without imposed precipitation reduction.

# The main goal is to separate the effect of the imposed drought treatment from
# background interannual climatic variability. For this reason, each process is
# expressed as an anomaly relative to the control simulation:

# relative anomaly = (process_drought - process_control) / process_control

# Values below zero indicate that the process value in the drought scenario is
# lower than in the control simulation for the same year. Values above zero
# indicate that the process value is higher than in the control simulation.

# It follows the rationale analysis of Mitchell et al., 2016 - GCB : An ecoclimatic framework for evaluating the resilience of vegetation to water

# For each process and each imposed drought
# event, the script identifies:
# 1. the maximum deficit relative to the control;
# 2. the year in which this maximum deficit occurs;
# 3. the half-recovery target, defined as halfway between the maximum deficit
# and full recovery relative to the control;
# 4. whether half recovery is reached before the next imposed drought event;
# 5. whether half recovery is reached later in the simulation;
# 6. the remaining deficit before the next drought event or by the end of the
# time series.

# Half recovery does not mean recovery to half of the original process value.
# Instead, it represents recovery of half of the drought-induced deficit relative
# to the control simulation. For example, if the maximum relative anomaly is
# -0.40, the half-recovery target is -0.20.

# This analysis quantifies short-term inter-event recovery, delayed recovery,
# incomplete recovery, and cumulative functional deficits under recurrent
# imposed drought. It is intended to support the evaluation of ecosystem
# resistance and recovery across multiple processes and drought frequencies.
# ==============================================================================


# Load required packages
library(tidyverse)

input_path <- "~/Desktop/CAETE-DVM-alloc-allom/scripts/yearly_mean_tables"

regclim <- read.csv(file.path(input_path, "MAN_regularclimate_yearly.csv"))

freq_8y <- read.csv(file.path(input_path, "MAN_30prec_8y_yearly.csv"))

freq_6y <- read.csv(file.path(input_path, "MAN_30prec_6y_yearly.csv"))

freq_4y <- read.csv(file.path(input_path, "MAN_30prec_4y_yearly.csv"))

freq_2y <- read.csv(file.path(input_path, "MAN_30prec_2y_yearly.csv"))

###############################################################
# Preparing data

# Merge tables and include a column for scenario identification
all_scenarios <- bind_rows(
  regclim = regclim,
  `8y` = freq_8y,
  `6y` = freq_6y,
  `4y` = freq_4y,
  `2y` = freq_2y,
  .id = "scenario"
) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("regclim", "8y", "6y", "4y", "2y")
    )
  )

# Change the name of the column "date" to "year"
names(all_scenarios)[names(all_scenarios) == "date"] <- "year"


###############################################################
# Remove invalid years after collapse in the 2-year scenario

last_valid_year_2y <- 2007

all_scenarios <- all_scenarios %>%
  filter(
    as.character(scenario) != "2y" |
      year <= last_valid_year_2y
  )

###############################################################

# Define the biogeochemical processes to be analyzed
# IMPORTANT: Edit this vector according to the exact column names in your tables.
process_vars <- c(
  "npp",
  "photo",
  "cleaf",
  "croot",
  "cwood",
  "csap",
  "cheart",
  "evapm",
  "csto",
  "ctotal", 
  "ls"
)

process_labels <- c(
  "npp"    = "NPP",
  "ctotal" = "Total carbon",
  "evapm"  = "Evapotranspiration",
  "ls"     = "N. of surviving strategies",
  "cleaf"  = "Leaf carbon",
  "croot"  = "Fine-root carbon",
  "cwood"  = "Aboveground woody carbon",
  "csap"   = "Sapwood carbon",
  "cheart" = "Heartwood carbon"
)

# Keep only variables that are actually present in the data
process_vars <- process_vars[process_vars %in% names(all_scenarios)]

# Select only the control scenario
control_wide <- all_scenarios %>%
  filter(scenario == "regclim")

# Keep only year and selected biogeochemical processes for the control
control_wide <- control_wide %>%
  select(year, all_of(process_vars))

# Convert the control table from wide to long format
control_long <- control_wide %>%
  pivot_longer(
    cols = all_of(process_vars),
    names_to = "process",
    values_to = "control_value"
  )

# Select only drought scenarios
drought_wide <- all_scenarios %>%
  filter(scenario != "regclim")

# Keep only scenario, year, and selected biogeochemical processes
drought_wide <- drought_wide %>%
  select(scenario, year, all_of(process_vars))

# Convert drought scenarios from wide to long format
drought_long <- drought_wide %>%
  pivot_longer(
    cols = all_of(process_vars),
    names_to = "process",
    values_to = "drought_value"
  )

# Join drought scenarios with the paired control by year and process
analysis_data <- drought_long %>%
  left_join(
    control_long,
    by = c("year", "process")
  )

###############################################################


###############################################################
# Calculate relative anomaly and ratio to control

analysis_data <- analysis_data %>%
  mutate(
    # Relative anomaly compared with the paired control
    relative_anomaly = if_else(
      control_value != 0,
      (drought_value - control_value) / control_value,
      NA_real_
    ),
    
    # Ratio to control
    ratio_to_control = if_else(
      control_value != 0,
      drought_value / control_value,
      NA_real_
    )
  ) %>%
  arrange(scenario, process, year)

###############################################################

###############################################################
# Define imposed drought event years

# Define drought event years for each recurrence frequency
drought_events <- bind_rows(
  tibble(
    scenario = "8y",
    frequency = 8,
    event_year = seq(1980, 2016, by = 8)
  ),
  
  tibble(
    scenario = "6y",
    frequency = 6,
    event_year = seq(1980, 2016, by = 6)
  ),
  
  tibble(
    scenario = "4y",
    frequency = 4,
    event_year = seq(1980, 2016, by = 4)
  ),
  
  tibble(
    scenario = "2y",
    frequency = 2,
    event_year = seq(1980, 2016, by = 2)
  )
)

###############################################################
# Calculate the response window for each of the scenarios

# Event-response window
# ---------------------
# For each drought event i, the event-response window extends from the drought
# year (t0_i) to the year immediately preceding the subsequent drought:
#
# W_response,i = [t0_i, t0_i+1 - 1]
#
# For a fixed drought recurrence interval of f years:
#
# W_response,i = [t0_i, t0_i + f - 1]
#
# This window contains all annual observations associated with the current
# drought before the next imposed drought occurs.

# Assign an identification number to each drought event
drought_events <- drought_events %>%
  group_by(scenario) %>%
  arrange(event_year, .by_group = TRUE) %>%
  mutate(
    event_id = row_number()
  ) %>%
  ungroup()

# Identify the last year available in the time series
last_simulation_year <- max(
  analysis_data$year,
  na.rm = TRUE
)

# Calculate the year of the subsequent drought
response_windows <- drought_events %>%
  mutate(
    next_event_year = event_year + frequency
  )

# Define the beginning of each response window
response_windows <- response_windows %>%
  mutate(
    response_start = event_year
  )

# Calculate the end of each response window
response_windows <- response_windows %>%
  mutate(
    response_end_expected = next_event_year - 1
  )

# Limit response windows to the available simulation period (
# because considering only the frequency would surpass the last simulation year)
response_windows <- response_windows %>%
  mutate(
    response_end = pmin(
      response_end_expected,
      last_simulation_year
    )
  )

# Identify whether the complete response window was observed

response_windows <- response_windows %>%
  mutate(
    response_window_complete =
      response_end_expected <= last_simulation_year
  )
#Uncomment if something changed in this table
#write_csv(response_windows, "~/Desktop/CAETE-DVM-alloc-allom/scripts/recovery_rate_analysis/response_windows.csv" )
###############################################################


###############################################################
# Calculate event metrics for all selected processes
###############################################################

# Use every variable previously defined in process_vars
summary_processes <- process_vars

# Numerical tolerance used throughout the calculations
metric_tolerance <- 1e-10


###############################################################
# Combine all processes with their response windows

all_processes_response_data <- analysis_data %>%
  filter(
    process %in% summary_processes
  ) %>%
  mutate(
    scenario = as.character(scenario)
  ) %>%
  inner_join(
    response_windows %>%
      mutate(
        scenario = as.character(scenario)
      ),
    by = "scenario",
    relationship = "many-to-many"
  ) %>%
  filter(
    year >= response_start,
    year <= response_end
  )


###############################################################
# Identify the maximum deficit within each response window

all_processes_maximum_deficit <- all_processes_response_data %>%
  group_by(
    scenario,
    process,
    frequency,
    event_id,
    event_year,
    response_start,
    response_end,
    response_window_complete
  ) %>%
  slice_min(
    order_by = relative_anomaly,
    n = 1,
    with_ties = FALSE,
    na_rm = TRUE
  ) %>%
  ungroup() %>%
  transmute(
    scenario,
    process,
    frequency,
    event_id,
    event_year,
    response_start,
    response_end,
    response_window_complete,
    
    # Year of maximum deficit
    t_min = year,
    
    # Signed minimum anomaly
    anomaly_min = relative_anomaly,
    
    # Positive magnitude of the maximum deficit
    maximum_deficit = if_else(
      anomaly_min < 0,
      -anomaly_min,
      0
    ),
    
    # Time to maximum deficit, defined only when a deficit occurred
    time_to_maximum_deficit = if_else(
      maximum_deficit > metric_tolerance,
      t_min - response_start,
      NA_real_
    )
  )


###############################################################
# Create a lookup table for the pre-event anomalies

all_processes_anomaly_lookup <- analysis_data %>%
  filter(
    process %in% summary_processes
  ) %>%
  transmute(
    scenario = as.character(scenario),
    process,
    pre_event_year = year,
    anomaly_pre = relative_anomaly
  )


###############################################################
# Add the pre-event anomaly

all_processes_event_metrics <- all_processes_maximum_deficit %>%
  mutate(
    pre_event_year = event_year - 1
  ) %>%
  left_join(
    all_processes_anomaly_lookup,
    by = c(
      "scenario",
      "process",
      "pre_event_year"
    )
  )

###############################################################
# Identify the anomaly during the drought-event year

all_processes_event_year_anomaly <- all_processes_response_data %>%
  filter(
    year == event_year
  ) %>%
  transmute(
    scenario,
    process,
    event_id,
    event_year,
    anomaly_event_year = relative_anomaly
  )


# Add the drought-event-year anomaly to the event table

all_processes_event_metrics <- all_processes_event_metrics %>%
  left_join(
    all_processes_event_year_anomaly,
    by = c(
      "scenario",
      "process",
      "event_id",
      "event_year"
    )
  )

###############################################################
# Identify anomaly at the end of each response window

all_processes_end_anomaly <- all_processes_response_data %>%
  filter(
    year == response_end
  ) %>%
  transmute(
    scenario,
    process,
    event_id,
    event_year,
    response_end,
    anomaly_end = relative_anomaly
  )


# Add end-of-window anomaly to the event table

all_processes_event_metrics <- all_processes_event_metrics %>%
  left_join(
    all_processes_end_anomaly,
    by = c(
      "scenario",
      "process",
      "event_id",
      "event_year",
      "response_end"
    )
  )


###############################################################
# Calculate impact and recovery metrics

all_processes_event_metrics <- all_processes_event_metrics %>%
  mutate(
    # Signed change relative to the pre-event condition
    incremental_anomaly_change =
      anomaly_min - anomaly_pre,
    
    # Signed change during the drought-event year
    initial_anomaly_change =
      anomaly_event_year - anomaly_pre,
    
    # Positive loss occurring during the drought-event year
    initial_event_year_loss = pmax(
      0,
      anomaly_pre - anomaly_event_year
    ),
    
    # Positive magnitude of the loss caused by the current event
    event_specific_incremental_loss = pmax(
      0,
      anomaly_pre - anomaly_min
    ),
    
    # Amount recovered after the maximum deficit
    recovered_amount =
      anomaly_end - anomaly_min,
    
    # Fraction of the event-specific loss recovered
    event_specific_recovery_fraction = if_else(
      event_specific_incremental_loss > metric_tolerance,
      recovered_amount / event_specific_incremental_loss,
      NA_real_
    ),
    
    # Remaining deficit relative to the control
    residual_deficit = pmax(
      0,
      -anomaly_end
    ),
    
    # Recovery classification
    event_specific_recovery_status = case_when(
      is.na(event_specific_incremental_loss) |
        is.na(recovered_amount) ~
        "Not calculated: missing data",
      
      event_specific_incremental_loss <= metric_tolerance ~
        "No incremental loss",
      
      event_specific_recovery_fraction >=
        (1 - metric_tolerance) ~
        "Complete recovery",
      
      event_specific_recovery_fraction >
        metric_tolerance ~
        "Partial recovery",
      
      TRUE ~
        "No recovery"
    ),
    
    # Target corresponding to recovery of 50% of the
    # event-specific incremental loss
    half_recovery_target = if_else(
      event_specific_incremental_loss > metric_tolerance,
      
      anomaly_min +
        0.5 * (anomaly_pre - anomaly_min),
      
      NA_real_
    )
  )


###############################################################
# Identify the first year in which half recovery was reached

all_processes_half_recovery <- all_processes_response_data %>%
  left_join(
    all_processes_event_metrics %>%
      select(
        scenario,
        process,
        event_id,
        t_min,
        half_recovery_target
      ),
    by = c(
      "scenario",
      "process",
      "event_id"
    )
  ) %>%
  filter(
    year > t_min,
    year <= response_end,
    !is.na(half_recovery_target),
    relative_anomaly >= half_recovery_target
  ) %>%
  group_by(
    scenario,
    process,
    event_id,
    t_min
  ) %>%
  summarise(
    half_recovery_year = min(year),
    .groups = "drop"
  ) %>%
  mutate(
    half_recovery_time =
      half_recovery_year - t_min
  ) %>%
  select(
    scenario,
    process,
    event_id,
    half_recovery_year,
    half_recovery_time
  )


###############################################################
# Add half-recovery results to the event table

all_processes_event_metrics <- all_processes_event_metrics %>%
  left_join(
    all_processes_half_recovery,
    by = c(
      "scenario",
      "process",
      "event_id"
    )
  ) %>%
  mutate(
    half_recovery_observed =
      !is.na(half_recovery_year),
    
    half_recovery_status = case_when(
      is.na(event_specific_incremental_loss) ~
        "Not calculated: missing data",
      
      event_specific_incremental_loss <= metric_tolerance ~
        "Not applicable: no incremental loss",
      
      half_recovery_observed ~
        "Observed",
      
      !response_window_complete ~
        "Not observed before simulation end",
      
      TRUE ~
        "Not observed before next drought"
    )
  )


###############################################################
# Final summary table containing all processes

all_processes_summary_table <- all_processes_event_metrics %>%
  mutate(
    process = factor(
      process,
      levels = summary_processes
    ),
    
    scenario = factor(
      scenario,
      levels = c("8y", "6y", "4y", "2y")
    )
  ) %>%
  select(
    # Process identification
    process,
    
    # Event identification
    scenario,
    event_id,
    event_year,
    response_start,
    response_end,
    response_window_complete,
    
    # Block A: long-term condition
    anomaly_pre,
    
    # Block B: event impact
    t_min,
    anomaly_min,
    maximum_deficit,
    time_to_maximum_deficit,
    event_specific_incremental_loss,
    
    # Block C: event recovery
    anomaly_end,
    recovered_amount,
    event_specific_recovery_fraction,
    event_specific_recovery_status,
    half_recovery_target,
    half_recovery_year,
    half_recovery_time,
    half_recovery_status,
    
    anomaly_event_year,
    initial_anomaly_change,
    initial_event_year_loss
  ) %>%
  arrange(
    process,
    scenario,
    event_year
  )

all_processes_summary_table %>%
  count(
    process,
    scenario,
    name = "number_of_events"
  ) %>%
  print(
    n = Inf
  )

# write_csv(
#   all_processes_summary_table,
#   "~/Desktop/CAETE-DVM-alloc-allom/scripts/recovery_rate_analysis/all_processes_summary_table.csv"
# )

###############################################################
# Annual anomaly trajectories for the main ecosystem variables

selected_processes <- c(
  "npp",
  "evapm",
  "ctotal",
  "ls"
)

annual_anomaly_plot_data <- analysis_data %>%
  filter(
    process %in% selected_processes
  ) %>%
  mutate(
    # Order and labels of drought-frequency scenarios
    scenario = factor(
      scenario,
      levels = c("8y", "6y", "4y", "2y"),
      labels = c(
        "8y",
        "6y",
        "4y",
        "2y"
      )
    ),
    
    # Order and labels of panels
    process = factor(
      process,
      levels = selected_processes,
      labels = c(
        "NPP",
        "Evapotranspiration",
        "Total carbon",
        "Surviving strategies"
      )
    )
  )


###############################################################
# Colorblind-friendly colors for drought recurrence intervals

scenario_colors <- c(
  "8y" = "#0072B2",
  "6y" = "#009E73",
  "4y" = "#E69F00",
  "2y" = "#D55E00"
)

###############################################################
# Plot annual anomaly trajectories

annual_anomaly_plot <- ggplot(
  annual_anomaly_plot_data,
  aes(
    x = year,
    y = relative_anomaly,
    color = scenario,
    group = scenario
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey60"
  ) +
  
  geom_line(
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~ process,
    ncol = 2
  ) +
  
  scale_color_manual(
    name = "Drought recurrence",
    values = scenario_colors
  ) +
  
  scale_x_continuous(
    breaks = seq(1980, 2015, by = 5)
  ) +
  
  labs(
    x = "Year",
    y = "Relative anomaly"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(
      margin = margin(r = 10)
    ),
    axis.title.x = element_text(
      margin = margin(t = 10)
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      color = "grey92"
    ),
    panel.border = element_rect(
      color = "grey40",
      fill = NA,
      linewidth = 0.8
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 1.4
      )
    )
  )

annual_anomaly_plot

###############################################################
# Save figure in publication quality

output_path <- paste0(
  "~/Desktop/CAETE-DVM-alloc-allom/scripts/",
  "recovery_rate_analysis/figures"
)

dir.create(
  output_path,
  recursive = TRUE,
  showWarnings = FALSE
)


# # Vector format: preferred for line plots
# 
# ggsave(
#   filename = file.path(
#     output_path,
#     "annual_anomaly_trajectories.pdf"
#   ),
#   plot = annual_anomaly_plot,
#   device = cairo_pdf,
#   width = 180,
#   height = 145,
#   units = "mm",
#   bg = "white"
# )
# 
# 
# # High-resolution raster format
# 
# ggsave(
#   filename = file.path(
#     output_path,
#     "annual_anomaly_trajectories.tiff"
#   ),
#   plot = annual_anomaly_plot,
#   device = "tiff",
#   width = 180,
#   height = 145,
#   units = "mm",
#   dpi = 600,
#   compression = "lzw",
#   bg = "white"
# )

###############################################################
# Pre-event anomaly trajectories

pre_event_plot_data <- all_processes_summary_table %>%
  filter(
    process %in% c(
      "npp",
      "evapm",
      "ctotal",
      "ls"
    )
  ) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("8y", "6y", "4y", "2y"),
      labels = c(
        "8y",
        "6y",
        "4y",
        "2y"
      )
    ),
    
    process = factor(
      as.character(process),
      levels = c(
        "npp",
        "evapm",
        "ctotal",
        "ls"
      ),
      labels = c(
        "NPP",
        "Evapotranspiration",
        "Total carbon",
        "Surviving strategies"
      )
    )
  )



###############################################################
# Plot pre-event anomaly trajectories

pre_event_anomaly_plot <- ggplot(
  pre_event_plot_data,
  aes(
    x = event_year,
    y = anomaly_pre,
    color = scenario,
    group = scenario
  )
) +
  
  # Control reference
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey60"
  ) +
  
  # Connect successive drought events
  geom_line(
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  
  # Show the condition immediately before each event
  geom_point(
    size = 1.2,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~ process,
    ncol = 2
  ) +
  
  scale_color_manual(
    name = "Drought recurrence",
    values = scenario_colors
  ) +
  
  scale_x_continuous(
    breaks = seq(1980, 2015, by = 5)
  ) +
  
  labs(
    x = "Drought-event year",
    y = "Pre-event relative anomaly"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(
      margin = margin(r = 10)
    ),
    axis.title.x = element_text(
      margin = margin(t = 10)
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      color = "grey92"
    ),
    panel.border = element_rect(
      color = "grey40",
      fill = NA,
      linewidth = 0.8
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 1.4,
        size = 2.5
      )
    )
  )

pre_event_anomaly_plot

# ggsave(
#   filename = file.path(
#     output_path,
#     "pre_event_anomaly_trajectories.pdf"
#   ),
#   plot = pre_event_anomaly_plot,
#   device = cairo_pdf,
#   width = 180,
#   height = 145,
#   units = "mm",
#   bg = "white"
# )
# 
# ggsave(
#   filename = file.path(
#     output_path,
#     "pre_event_anomaly_trajectories.tiff"
#   ),
#   plot = pre_event_anomaly_plot,
#   device = "tiff",
#   width = 180,
#   height = 145,
#   units = "mm",
#   dpi = 600,
#   compression = "lzw",
#   bg = "white"
# )

###############################################################
# Heatmap of integrated negative deficit: 1980–2016
###############################################################

# Variables included in the analysis
selected_processes <- c(
  "npp",
  "evapm",
  "ctotal",
  "ls"
)

# Common analysis period
analysis_start_year <- 1980
analysis_end_year   <- 2016

expected_number_of_years <-
  analysis_end_year - analysis_start_year + 1


###############################################################
# Calculate the integrated negative deficit
#
# D_v,f = -sum[min(Anomaly_v,f,t, 0)]
#
# Positive anomalies are assigned zero and therefore do not
# compensate for years in which the variable was below control.

integrated_deficit_data <- analysis_data %>%
  filter(
    process %in% selected_processes,
    year >= analysis_start_year,
    year <= analysis_end_year
  ) %>%
  mutate(
    scenario = as.character(scenario)
  ) %>%
  group_by(
    process,
    scenario
  ) %>%
  summarise(
    # Number of years containing valid anomaly values
    number_of_years_observed =
      n_distinct(year[!is.na(relative_anomaly)]),
    
    # Identify missing anomaly values
    has_missing_values =
      any(is.na(relative_anomaly)),
    
    # Integrated negative deficit
    raw_integrated_deficit =
      -sum(
        pmin(relative_anomaly, 0),
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Only retain the metric when all 37 years are available
    integrated_negative_deficit = if_else(
      number_of_years_observed ==
        expected_number_of_years &
        !has_missing_values,
      
      raw_integrated_deficit,
      NA_real_
    )
  )


###############################################################
# Prepare data for the heatmap

integrated_deficit_heatmap_data <-
  integrated_deficit_data %>%
  mutate(
    # Order of drought-recurrence intervals
    scenario = factor(
      scenario,
      levels = c("8y", "6y", "4y", "2y")
    ),
    
    # Reverse factor order so NPP appears at the top
    process = factor(
      as.character(process),
      levels = c(
        "ls",
        "ctotal",
        "evapm",
        "npp"
      ),
      labels = c(
        "Surviving strategies",
        "Total carbon",
        "Evapotranspiration",
        "NPP"
      )
    ),
    
    # Values displayed inside cells
    deficit_label = case_when(
      as.character(scenario) == "2y" &
        number_of_years_observed <
        expected_number_of_years ~
        "Collapsed",
      
      is.na(integrated_negative_deficit) ~
        "NA",
      
      TRUE ~
        sprintf(
          "%.3f",
          integrated_negative_deficit
        )
    )
  )


###############################################################
# Define text colors according to cell background

maximum_integrated_deficit <- max(
  integrated_deficit_heatmap_data$
    integrated_negative_deficit,
  na.rm = TRUE
)

integrated_deficit_heatmap_data <-
  integrated_deficit_heatmap_data %>%
  mutate(
    # With direction = 1, smaller values have dark colors
    # and larger values have light colors
    label_color = case_when(
      is.na(integrated_negative_deficit) ~
        "grey30",
      
      integrated_negative_deficit <=
        0.45 * maximum_integrated_deficit ~
        "white",
      
      TRUE ~
        "black"
    )
  )


###############################################################
# Create the heatmap

integrated_deficit_heatmap <- ggplot(
  integrated_deficit_heatmap_data,
  aes(
    x = scenario,
    y = process,
    fill = integrated_negative_deficit
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 1
  ) +
  
  geom_text(
    aes(
      label = deficit_label,
      color = label_color
    ),
    size = 4.2,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  # Common color scale across all variables
  # Larger deficits are represented by lighter colors
  scale_fill_viridis_c(
    name = paste0(
      "Integrated negative deficit"
    ),
    option = "C",
    direction = 1,
    na.value = "grey85"
  ) +
  
  scale_color_identity() +
  
  # Keep cells square
  coord_equal() +
  
  labs(
    x = "Drought recurrence interval",
    y = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    axis.title.x = element_text(
      face = "bold",
      margin = margin(t = 10)
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text.y = element_text(
      face = "bold",
      size = 12
    ),
    
    panel.grid = element_blank(),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    legend.position = "right"
  )


# Display the heatmap
integrated_deficit_heatmap


###############################################################
# Check temporal coverage before interpreting the heatmap

integrated_deficit_data %>%
  select(
    process,
    scenario,
    number_of_years_observed,
    has_missing_values,
    integrated_negative_deficit
  ) %>%
  arrange(
    process,
    scenario
  ) %>%
  print(
    n = Inf
  )


# Vector PDF

# ggsave(
#   filename = file.path(
#     output_path,
#     "integrated_negative_deficit_heatmap.pdf"
#   ),
#   plot = integrated_deficit_heatmap,
#   device = cairo_pdf,
#   width = 180,
#   height = 125,
#   units = "mm",
#   bg = "white"
# )
# 
# 
# # TIFF at 600 dpi
# 
# ggsave(
#   filename = file.path(
#     output_path,
#     "integrated_negative_deficit_heatmap.tiff"
#   ),
#   plot = integrated_deficit_heatmap,
#   device = "tiff",
#   width = 180,
#   height = 125,
#   units = "mm",
#   dpi = 600,
#   compression = "lzw",
#   bg = "white"
# )

# ==============================================================================
# Event-specific magnitude of ecosystem responses to recurrent drought
# ==============================================================================
# ------------------------------------------------------------------------------
# Prepare the event-level data
# ------------------------------------------------------------------------------

incremental_loss_plot_data <- all_processes_summary_table %>%
  mutate(
    process_key = tolower(as.character(process)),
    
    process_label = case_when(
      process_key == "npp" ~ "NPP",
      
      process_key %in% c(
        "evapm",
        "evapotranspiration",
        "et"
      ) ~ "Evapotranspiration",
      
      process_key %in% c(
        "ctotal",
        "total_carbon",
        "total carbon"
      ) ~ "Total carbon",
      
      process_key %in% c(
        "ls",
        "surviving_strategies",
        "number_of_surviving_strategies",
        "surviving strategies"
      ) ~ "Surviving strategies",
      
      TRUE ~ NA_character_
    ),
    
    scenario = factor(
      as.character(scenario),
      levels = c("8y", "6y", "4y", "2y")
    ),
    
    window_status = if_else(
      response_window_complete,
      "Complete window",
      "Incomplete window"
    ),
    
    window_status = factor(
      window_status,
      levels = c(
        "Complete window",
        "Incomplete window"
      )
    ),
    
    process_label = factor(
      process_label,
      levels = c(
        "NPP",
        "Evapotranspiration",
        "Total carbon",
        "Surviving strategies"
      )
    )
  ) %>%
  filter(
    !is.na(event_specific_incremental_loss),
    !is.na(process_label),
    !is.na(scenario)
  )

# ------------------------------------------------------------------------------
# Create the event-specific incremental-loss figure
# ------------------------------------------------------------------------------

incremental_loss_plot <- ggplot(
  incremental_loss_plot_data,
  aes(
    x = event_year,
    y = event_specific_incremental_loss,
    color = scenario,
    group = scenario
  )
) +
  
  # Connect successive drought events within each scenario
  geom_line(
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  
  # Zero means that the event did not intensify
  # the pre-event deficit
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    linetype = "dashed",
    color = "grey60"
  ) +
  
  # One panel for each ecosystem process
  facet_wrap(
    facets = vars(process_label),
    ncol = 2,
    scales = "fixed"
  ) +
  
  # Same colors used in the previous figures
  scale_color_manual(
    name = "Drought recurrence",
    values = scenario_colors,
    breaks = c("8y", "6y", "4y", "2y"),
    labels = c("8y", "6y", "4y", "2y"),
    drop = FALSE
  ) +
  
  scale_x_continuous(
    breaks = seq(1980, 2015, by = 5),
    limits = c(1979, 2017),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.10)
    )
  ) +
  
  labs(
    x = "Drought-event year",
    y = "Event-specific incremental loss"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title = element_text(face = "bold"),
    
    strip.background = element_rect(
      fill = "grey94",
      color = "grey70",
      linewidth = 0.4
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    
    panel.spacing = grid::unit(
      0.8,
      "lines"
    ),
    
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 1.4
      )
    )
  )


# Display the figure
incremental_loss_plot

ggsave(
  filename = file.path(
    output_path,
    "event_specific_incremental_loss.tiff"
  ),
  plot = incremental_loss_plot,
  device = "tiff",
  width = 180,
  height = 145,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# ==============================================================================
# Time to maximum event-specific deficit
# ==============================================================================

# ------------------------------------------------------------------------------
# Define process order, labels, colors, and shapes
# ------------------------------------------------------------------------------

process_levels <- c(
  "npp",
  "evapm",
  "ctotal",
  "ls"
)

process_labels <- c(
  "NPP",
  "Evapotranspiration",
  "Total carbon",
  "Surviving strategies"
)

# Colorblind-friendly palette
process_colors <- c(
  "NPP" = "#0072B2",                   # Blue
  "Evapotranspiration" = "#009E73",    # Bluish green
  "Total carbon" = "#E69F00",          # Orange
  "Surviving strategies" = "#D55E00"  # Vermillion
)

# Shapes provide an additional distinction between processes
process_shapes <- c(
  "NPP" = 16,
  "Evapotranspiration" = 17,
  "Total carbon" = 15,
  "Surviving strategies" = 18
)


# ------------------------------------------------------------------------------
# Prepare event-level data
# ------------------------------------------------------------------------------

time_to_maximum_deficit_plot_data <-
  all_processes_summary_table %>%
  filter(
    process %in% process_levels,
    
    # The definitive minimum cannot be identified when the
    # response window is incomplete.
    response_window_complete,
    
    # Time is undefined when the event produced no incremental
    # loss relative to the pre-event condition.
    !is.na(time_to_maximum_deficit)
  ) %>%
  mutate(
    # Preserve the drought-recurrence order
    scenario = factor(
      as.character(scenario),
      levels = c(
        "8y",
        "6y",
        "4y",
        "2y"
      )
    ),
    
    # Include the possible time range in each panel label
    scenario_label = factor(
      as.character(scenario),
      levels = c(
        "8y",
        "6y",
        "4y",
        "2y"
      ),
      labels = c(
        "8-year interval (0–7 years)",
        "6-year interval (0–5 years)",
        "4-year interval (0–3 years)",
        "2-year interval (0–1 year)"
      )
    ),
    
    # Define process labels and their order
    process_label = factor(
      as.character(process),
      levels = process_levels,
      labels = process_labels
    )
  )

# ------------------------------------------------------------------------------
# Prepare drought-event years for the vertical reference lines
# ------------------------------------------------------------------------------

drought_event_lines <- all_processes_summary_table %>%
  filter(
    process %in% process_levels,
    !is.na(event_year)
  ) %>%
  mutate(
    # Assign each drought event to its corresponding facet
    scenario_label = factor(
      as.character(scenario),
      levels = c(
        "8y",
        "6y",
        "4y",
        "2y"
      ),
      labels = c(
        "8-year interval (0–7 years)",
        "6-year interval (0–5 years)",
        "4-year interval (0–3 years)",
        "2-year interval (0–1 year)"
      )
    )
  ) %>%
  distinct(
    scenario_label,
    event_year
  )

# ------------------------------------------------------------------------------
# Define the possible vertical range for each recurrence interval
#
# These values ensure that each panel displays the complete range
# allowed by its response-window duration, even when the observed
# values do not reach the theoretical maximum.
# ------------------------------------------------------------------------------

time_axis_limits <- tibble::tibble(
  scenario_label = factor(
    c(
      "8-year interval (0–7 years)",
      "6-year interval (0–5 years)",
      "4-year interval (0–3 years)",
      "2-year interval (0–1 year)"
    ),
    levels = c(
      "8-year interval (0–7 years)",
      "6-year interval (0–5 years)",
      "4-year interval (0–3 years)",
      "2-year interval (0–1 year)"
    )
  ),
  
  event_year = analysis_start_year,
  
  minimum_time = 0,
  
  maximum_time = c(
    7,
    5,
    3,
    1
  )
)


# ------------------------------------------------------------------------------
# Create the point figure
# ------------------------------------------------------------------------------

time_to_maximum_deficit_plot <- ggplot(
  time_to_maximum_deficit_plot_data,
  aes(
    x = event_year,
    y = time_to_maximum_deficit,
    color = process_label,
    shape = process_label
  )
) +
  
  # Force each panel to show the complete range allowed by its
  # corresponding response-window duration.
  geom_blank(
    data = time_axis_limits,
    aes(
      x = event_year,
      y = minimum_time
    ),
    inherit.aes = FALSE
  ) +
  
  geom_blank(
    data = time_axis_limits,
    aes(
      x = event_year,
      y = maximum_time
    ),
    inherit.aes = FALSE
  ) +
  # Mark the calendar year of each imposed drought event
  geom_vline(
    data = drought_event_lines,
    aes(
      xintercept = event_year
    ),
    inherit.aes = FALSE,
    color = "grey82",
    linetype = "dashed",
    linewidth = 0.35
  ) +
  
  # Mark time zero relative to the drought event
  geom_hline(
    yintercept = 0,
    color = "grey45",
    linetype = "solid",
    linewidth = 0.45
  ) +
  
  # Slightly separate processes observed during the same event
  geom_point(
    position = position_dodge(
      width = 0.9
    ),
    size = 2.8,
    stroke = 0.8
  ) +
  
  
  # Slightly separate processes that share the same event year
  geom_point(
    position = position_dodge(
      width = 0.9
    ),
    size = 2.8,
    stroke = 0.8
  ) +
  
  # Each panel represents one drought-recurrence interval.
  # Free vertical scales are necessary because the maximum possible
  # time differs among recurrence intervals.
  facet_wrap(
    facets = vars(scenario_label),
    ncol = 2,
    scales = "free_y"
  ) +
  
  scale_color_manual(
    name = "Ecosystem indicator",
    values = process_colors,
    breaks = process_labels,
    drop = FALSE
  ) +
  
  scale_shape_manual(
    name = "Ecosystem indicator",
    values = process_shapes,
    breaks = process_labels,
    drop = FALSE
  ) +
  
  scale_x_continuous(
    breaks = seq(
      1980,
      2015,
      by = 5
    ),
    limits = c(
      1979,
      2017
    ),
    expand = expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +
  
  # Only integer values are meaningful because the analysis uses
  # annual data.
  scale_y_continuous(
    breaks = 0:7,
    expand = expansion(
      add = c(
        0.15,
        0.25
      )
    )
  ) +
  
  labs(
    x = "Drought-event year",
    y = "Time to maximum deficit (years)"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    # Horizontal grid lines facilitate comparison of discrete times
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank(),
    
    strip.background = element_rect(
      fill = "grey94",
      color = "grey70",
      linewidth = 0.4
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    
    panel.spacing = grid::unit(
      0.9,
      "lines"
    ),
    
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 3
      )
    ),
    
    shape = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )


# Display the figure
time_to_maximum_deficit_plot

