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
# Plot anomalies for selected processes (main figure)

facet_labels <- c(
  "8y" = "8-Year Interval",
  "6y" = "6-Year Interval",
  "4y" = "4-Year Interval",
  "2y" = "2-Year Interval"
)

process_colors <- c(
  "npp" = "#004B8D",      
  "evapm" = "#00874E",    
  "ctotal" = "#D9A000",   
  "ls" = "#BF4A00"        
)

# Select only the processes that will be shown in this plot
selected_processes <- c("npp", "evapm", "ctotal", "ls")

# Prepare data for plotting only the selected processes
plot_data <- analysis_data %>%
  filter(process %in% selected_processes) %>%
  # Group by scenario and process to evaluate the timeline sequentially
  group_by(scenario, process) %>%
  mutate(
    # Create a cumulative flag that turns TRUE once the system collapses
    # in the '2y' scenario (anomaly drops to -1 or reaches 1).
    # Using 0.99 accounts for potential floating-point imprecision.
    has_collapsed = cumany(scenario == "2y" & abs(year) >2007),
    
    # Replace the relative anomaly with NA from the collapse point onwards
    # so ggplot stops drawing the line
    relative_anomaly = if_else(has_collapsed, NA_real_, relative_anomaly)
  ) %>%
  ungroup() %>%
  mutate(
    process = factor(process, levels = selected_processes),
    scenario_label = factor(
      scenario,
      levels = names(facet_labels),
      labels = facet_labels
    )
  )

anomaly_plot <- ggplot(
  plot_data,
  aes(
    x = year,
    y = relative_anomaly,
    color = process,
    group = process
  )
) +
  # Zero reference line
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey60"
  ) +
  geom_line(
    linewidth = 1.0,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ scenario_label,
    ncol = 2,
    strip.position = "top"
  ) +
  scale_color_manual(
    name = NULL,
    values = process_colors[selected_processes],
    labels = process_labels[selected_processes]
  ) +
  scale_x_continuous(
    breaks = seq(1980, 2015, by = 5)
  ) +
  labs(
    x = "Year",
    y = "Relative anomaly"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(face = "bold", size = 12),
    
    # Facet titles
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 13,
      hjust = 0.05,
      margin = margin(t = 10, b = 5)
    ),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    
    # Rotate x-axis labels
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5)))

anomaly_plot

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
# Uncomment if something changed in this table
#write_csv(response_windows, "~/Desktop/CAETE-DVM-alloc-allom/scripts/recovery_rate_analysis/response_windows.csv" )
###############################################################

###############################################################
# Now let's calculate the mximum déficit for each scenario x variable x drought event

# NPP
# Select NPP data from all drought scenarios

npp_all_scenarios <- analysis_data %>%
  filter(
    process == "npp"
  ) %>%
  mutate(
    scenario = as.character(scenario)
  )

# Combine NPP anomalies with the response windows
npp_response_data <- npp_all_scenarios %>%
  left_join(
    response_windows,
    by = "scenario",
    relationship = "many-to-many"
  )

# Keep only years belonging to each response window
npp_response_data <- npp_response_data %>%
  filter(
    year >= response_start,
    year <= response_end
  )

# Calculate the maximum NPP deficit for each drought event
npp_maximum_deficit <- npp_response_data %>%
  
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
  
  # Select the minimum relative anomaly within each response window
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
    
    # Year when the minimum anomaly occurred
    t_min = year,
    
    # Signed minimum anomaly
    anomaly_min = relative_anomaly,
    
    # Positive magnitude of the maximum deficit
    maximum_deficit = if_else(
      anomaly_min < 0,
      -anomaly_min,
      0
    ),
    
    # Number of years from the drought to the maximum deficit
    time_to_maximum_deficit = t_min - response_start
  )

# Define scenario labels and their plotting order

npp_maximum_deficit <- npp_maximum_deficit %>%
  mutate(
    scenario_label = factor(
      scenario,
      levels = c("8y", "6y", "4y", "2y"),
      labels = c(
        "8-Year Interval",
        "6-Year Interval",
        "4-Year Interval",
        "2-Year Interval"
      )
    ),
    
    window_status = if_else(
      response_window_complete,
      "Complete window",
      "Incomplete window"
    )
  )

###############################################################
# Plot the time required to reach the maximum NPP deficit
# for all drought scenarios

npp_time_to_deficit_plot <- ggplot(
  npp_maximum_deficit,
  aes(
    x = event_year,
    y = time_to_maximum_deficit,
    fill = window_status
  )
) +
  geom_col(
    width = 1.3
  ) +
  geom_text(
    aes(label = time_to_maximum_deficit),
    vjust = -0.4,
    size = 3.5
  ) +
  facet_wrap(
    ~ scenario_label,
    ncol = 2
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Complete window" = "#004B8D",
      "Incomplete window" = "grey65"
    )
  ) +
  scale_x_continuous(
    breaks = seq(1980, 2016, by = 4)
  ) +
  scale_y_continuous(
    breaks = 0:7,
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Time to maximum NPP deficit",
    x = "Drought-event year",
    y = "Time to maximum deficit (years)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    )
  )

npp_time_to_deficit_plot

###############################################################
# Create a lookup table for pre-event NPP anomalies
# The pre-event anomaly will be used to calculate event-specific incremental loss
npp_anomaly_lookup <- analysis_data %>%
  filter(
    process == "npp"
  ) %>%
  transmute(
    scenario = as.character(scenario),
    process,
    
    # This year will be matched with the pre-event year
    pre_event_year = year,
    
    # Relative anomaly observed in that year
    anomaly_pre = relative_anomaly
  )

# Define the year immediately preceding each drought event
npp_event_metrics <- npp_maximum_deficit %>%
  mutate(
    scenario = as.character(scenario),
    pre_event_year = event_year - 1
  )

# Add the anomaly observed immediately before each drought event (pre-event anomaly)
npp_event_metrics <- npp_event_metrics %>%
  left_join(
    npp_anomaly_lookup,
    by = c(
      "scenario",
      "process",
      "pre_event_year"
    )
  )

# Calculate the signed anomaly change after each drought event 
# note this metric uses the minimum anomaly, that is the maximum 
# deficit in relation to the control
npp_event_metrics <- npp_event_metrics %>%
  mutate(
    incremental_anomaly_change =
      anomaly_min - anomaly_pre
  )

# Just make the previous variable as a positive magnitude
# It separates the legacy effect from the impact of this specific 
# drought event
npp_event_metrics <- npp_event_metrics %>%
  mutate(
    event_specific_incremental_loss = pmax(
      0,
      anomaly_pre - anomaly_min
    )
  )
###############################################################

###############################################################
# Identify the NPP anomaly at the end of each response window
npp_end_anomaly <- npp_response_data %>%
  filter(
    year == response_end
  ) %>%
  transmute(
    scenario = as.character(scenario),
    process,
    event_id,
    event_year,
    response_end,
    
    # Relative anomaly at the end of the response window
    anomaly_end = relative_anomaly
  )

# Add the end-of-window anomaly to the event metrics table
npp_event_metrics <- npp_event_metrics %>%
  left_join(
    npp_end_anomaly,
    by = c(
      "scenario",
      "process",
      "event_id",
      "event_year",
      "response_end"
    )
  )
###############################################################


###############################################################
# Calculate the TOTAL residual NPP deficit at the end of each
# response window in relation to the baseline

npp_event_metrics <- npp_event_metrics %>%
  mutate(
    residual_deficit = pmax(
      0,
      -anomaly_end
    )
  )

###############################################################

###############################################################
# Calculate the amount recovered after the minimum anomaly

# Numerical tolerance used to identify negligible incremental losses
recovery_tolerance <- 1e-10
# Calculate event-specific recovery metrics

npp_event_metrics <- npp_event_metrics %>%
  mutate(
    # Amount recovered between the minimum anomaly and the end
    # of the response window
    recovered_amount =
      anomaly_end - anomaly_min,
    
    # Fraction of the event-specific incremental loss that was recovered
    #
    # 0   = no recovery
    # 1   = complete recovery to the pre-event condition
    # < 0 = further deterioration after the minimum
    # > 1 = recovery beyond the pre-event condition
    event_specific_recovery_fraction = if_else(
      event_specific_incremental_loss > loss_tolerance,
      recovered_amount / event_specific_incremental_loss,
      NA_real_
    )
  )

# Prepare data for plotting
npp_recovery_plot_data <- npp_event_metrics %>%
  mutate(
    recovery_percentage =
      event_specific_recovery_fraction * 100,
    
    scenario_label = factor(
      scenario,
      levels = c("8y", "6y", "4y", "2y"),
      labels = c(
        "8-Year Interval",
        "6-Year Interval",
        "4-Year Interval",
        "2-Year Interval"
      )
    ),
    
    window_status = if_else(
      response_window_complete,
      "Complete window",
      "Incomplete window"
    )
  )


# Plot event-specific recovery

npp_recovery_plot <- ggplot(
  data = npp_recovery_plot_data,
  aes(
    x = event_year,
    y = recovery_percentage
  )
) +
  
  # Connect only events with complete response windows
  geom_line(
    data = npp_recovery_plot_data %>%
      filter(
        response_window_complete,
        !is.na(recovery_percentage)
      ),
    aes(group = 1),
    color = "grey55",
    linewidth = 0.8
  ) +
  
  # Plot only events for which recovery could be calculated
  geom_point(
    data = npp_recovery_plot_data %>%
      filter(!is.na(recovery_percentage)),
    aes(color = window_status),
    size = 3
  ) +
  
  # Complete recovery to the pre-event condition
  geom_hline(
    yintercept = 100,
    linetype = "dashed",
    color = "grey35"
  ) +
  
  # No recovery after the minimum anomaly
  geom_hline(
    yintercept = 0,
    color = "grey70"
  ) +
  
  facet_wrap(
    ~ scenario_label,
    ncol = 2
  ) +
  
  scale_color_manual(
    name = NULL,
    values = c(
      "Complete window" = "#004B8D",
      "Incomplete window" = "grey65"
    )
  ) +
  
  scale_x_continuous(
    breaks = seq(1980, 2016, by = 4)
  ) +
  
  scale_y_continuous(
    labels = scales::label_number(
      accuracy = 1,
      suffix = "%"
    )
  ) +
  
  labs(
    title = "Event-specific NPP recovery",
    x = "Drought-event year",
    y = "Recovery of event-specific incremental loss"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )


# Display the plot
npp_recovery_plot
###############################################################


###############################################################
# Classify event-specific recovery at the end of each
# response window

npp_event_metrics <- npp_event_metrics %>%
  mutate(
    event_specific_recovery_status = case_when(
      
      # Recovery cannot be calculated when no incremental
      # loss occurred
      is.na(event_specific_recovery_fraction) ~
        "No incremental loss",
      
      # The pre-event condition was reached or exceeded
      event_specific_recovery_fraction >=
        (1 - recovery_tolerance) ~
        "Complete recovery",
      
      # Part of the incremental loss was recovered
      event_specific_recovery_fraction >
        recovery_tolerance ~
        "Partial recovery",
      
      # No meaningful improvement occurred after the minimum
      TRUE ~
        "No recovery"
    )
  )

npp_event_metrics %>%
  select(
    scenario,
    event_year,
    anomaly_pre,
    anomaly_min,
    anomaly_end,
    event_specific_incremental_loss,
    recovered_amount,
    event_specific_recovery_fraction,
    event_specific_recovery_status,
    response_window_complete
  )

npp_event_metrics %>%
  count(
    scenario,
    event_specific_recovery_status
  )


###############################################################
# Calculate the event-specific half-recovery target

npp_event_metrics <- npp_event_metrics %>%
  mutate(
    half_recovery_target = if_else(
      event_specific_incremental_loss > loss_tolerance,
      
      anomaly_min +
        0.5 * (anomaly_pre - anomaly_min),
      
      NA_real_
    )
  )

###############################################################
# Identify the first year in which half recovery was reached

npp_half_recovery_results <- npp_response_data %>%
  left_join(
    npp_event_metrics %>%
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
# Add half-recovery results to the main event table

npp_event_metrics <- npp_event_metrics %>%
  select(
    -any_of(
      c(
        "half_recovery_year",
        "half_recovery_time",
        "half_recovery_observed",
        "half_recovery_status"
      )
    )
  ) %>%
  left_join(
    npp_half_recovery_results,
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
      is.na(half_recovery_target) ~
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


# Prepare half-recovery time data for plotting

npp_half_recovery_plot_data <- npp_event_metrics %>%
  mutate(
    # Time available to observe recovery after the minimum anomaly
    time_available_after_min =
      response_end - t_min,
    
    # Use the observed half-recovery time when available.
    # Otherwise, use the available follow-up time to show censoring.
    half_recovery_plot_time = if_else(
      half_recovery_observed,
      half_recovery_time,
      time_available_after_min
    ),
    
    # Define the observation status
    half_recovery_plot_status = case_when(
      half_recovery_observed ~
        "Half recovery observed",
      
      !response_window_complete ~
        "Not observed: simulation ended",
      
      TRUE ~
        "Not observed before next drought"
    ),
    
    # Define labels for observed and censored events
    half_recovery_label = case_when(
      half_recovery_observed ~
        as.character(half_recovery_time),
      
      TRUE ~
        paste0(">", time_available_after_min)
    )
  ) %>%
  
  # Remove events without an incremental loss
  filter(
    !is.na(half_recovery_target)
  )

###############################################################
# Plot event-specific NPP half-recovery time

npp_half_recovery_plot <- ggplot(
  npp_half_recovery_plot_data,
  aes(
    x = event_year,
    y = half_recovery_plot_time,
    color = half_recovery_plot_status,
    shape = half_recovery_plot_status
  )
) +
  
  geom_point(
    size = 3.2
  ) +
  
  geom_text(
    aes(label = half_recovery_label),
    vjust = -0.8,
    size = 3.5,
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~ scenario_label,
    ncol = 2
  ) +
  
  scale_color_manual(
    name = NULL,
    values = c(
      "Half recovery observed" = "#0072B2",
      "Not observed before next drought" = "#D55E00",
      "Not observed: simulation ended" = "grey60"
    )
  ) +
  
  scale_shape_manual(
    name = NULL,
    values = c(
      "Half recovery observed" = 16,
      "Not observed before next drought" = 17,
      "Not observed: simulation ended" = 15
    )
  ) +
  
  scale_x_continuous(
    breaks = seq(1980, 2016, by = 4)
  ) +
  
  scale_y_continuous(
    breaks = 0:7,
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  
  labs(
    title = "Event-specific NPP half-recovery time",
    x = "Drought-event year",
    y = "Time after maximum deficit (years)"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

npp_half_recovery_plot

###############################################################
# Long term anomaly analysis
# Plot the NPP condition immediately before each drought event

npp_pre_event_plot <- npp_event_metrics %>%
  ggplot(
    aes(
      x = event_year,
      y = anomaly_pre
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_line(
    aes(group = scenario),
    color = "#009E73",
    linewidth = 0.8
  ) +
  geom_point(
    color = "#009E73",
    size = 1.5
  ) +
  facet_wrap(
    ~ scenario_label,
    ncol = 2
  ) +
  scale_x_continuous(
    breaks = seq(1980, 2015, by = 5)
  ) +
  labs(
    title = "Pre-event NPP anomaly across successive droughts",
    x = "Drought-event year",
    y = "NPP anomaly before the drought"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

npp_pre_event_plot

###############################################################
# Organize all NPP event metrics into a single analysis table

npp_summary_table <- npp_event_metrics %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("8y", "6y", "4y", "2y")
    )
  ) %>%
  select(
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
    half_recovery_status
  ) %>%
  arrange(
    scenario,
    event_year
  )

###############################################################
# Inspect all NPP events in the 8-year drought scenario

npp_summary_table %>%
  filter(
    scenario == "8y"
  ) %>%
  print(
    n = Inf,
    width = Inf
  )

#Uncomment if something changed
#write_csv(npp_summary_table, "~/Desktop/CAETE-DVM-alloc-allom/scripts/recovery_rate_analysis/npp_summary_table.csv" )


