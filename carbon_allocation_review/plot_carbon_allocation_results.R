# =============================================================================
# Carbon allocation time-series plots
#
# This script reads the daily and final CSV files produced by the simplified
# Fortran carbon-allocation driver and creates a set of diagnostic figures.
#
# Expected input files:
#   carbon_allocation_daily.csv
#   carbon_allocation_final.csv
#
# Main outputs:
#   carbon_allocation_plots/*.png
#   carbon_allocation_plots/all_time_series_plots.pdf
#
# Run from the directory containing the CSV files:
#   Rscript plot_carbon_allocation_results.R
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load required packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "The following packages are required but not installed: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them with:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# -----------------------------------------------------------------------------
# 2. User settings
# -----------------------------------------------------------------------------

daily_file <- "~/Desktop/CAETE-DVM-alloc-allom/carbon_allocation_review/carbon_allocation_daily.csv"
final_file <- "~/Desktop/CAETE-DVM-alloc-allom/carbon_allocation_review/carbon_allocation_final.csv"

output_directory <- "~/Desktop/CAETE-DVM-alloc-allom/carbon_allocation_review/carbon_allocation_plots"

# Use 1 to plot every day. Increase this value for very long simulations.
plot_every_n_days <- 1L

# Image resolution used for PNG outputs.
png_resolution <- 300

# -----------------------------------------------------------------------------
# 3. Check input files and create the output directory
# -----------------------------------------------------------------------------

if (!file.exists(daily_file)) {
  stop(
    paste0(
      "Daily output file not found: ", daily_file,
      "\nRun this script from the directory containing the Fortran outputs, ",
      "or change 'daily_file' at the top of the script."
    ),
    call. = FALSE
  )
}

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# 4. Read and validate the daily output
# -----------------------------------------------------------------------------

daily <- read_csv(
  daily_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_daily_columns <- c(
  "npp_case",
  "npp_rate",
  "day",
  "year",
  "npp_daily",
  "leaf",
  "root",
  "sapwood",
  "heartwood",
  "height",
  "storage",
  "total_demand",
  "carbon_allocated",
  "leaf_demand",
  "root_demand",
  "sapwood_demand",
  "delta_leaf",
  "delta_root",
  "delta_sapwood",
  "unmet_storage_deficit",
  "unpaid_carbon_deficit",
  "leaf_turnover",
  "root_turnover",
  "sapwood_turnover",
  "storage_turnover",
  "heartwood_turnover",
  "leaf_root_residual",
  "pipe_model_residual"
)

missing_daily_columns <- setdiff(required_daily_columns, names(daily))

if (length(missing_daily_columns) > 0) {
  stop(
    paste0(
      "The daily file is missing the following required columns: ",
      paste(missing_daily_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

# Create a readable scenario label.
# The autocorrelated-NPP driver includes target_mean_npp, whereas the constant-NPP
# driver may contain only npp_rate.
if ("target_mean_npp" %in% names(daily)) {
  daily <- daily %>%
    mutate(
      scenario = paste0(
        "Case ", npp_case,
        " | target mean NPP = ",
        sprintf("%.2f", target_mean_npp)
      )
    )
} else {
  daily <- daily %>%
    mutate(
      scenario = paste0(
        "Case ", npp_case,
        " | NPP = ",
        sprintf("%.2f", npp_rate)
      )
    )
}

# Calculate aggregate carbon pools that are useful for interpretation.
daily <- daily %>%
  mutate(
    living_carbon = leaf + root + sapwood,
    stem_carbon = sapwood + heartwood,
    structural_carbon = leaf + root + sapwood + heartwood,
    total_plant_carbon = structural_carbon + storage
  )

# Reduce plotting density when requested, while retaining the first and final day.
daily_plot <- daily %>%
  group_by(npp_case) %>%
  filter(
    day == min(day) |
      day == max(day) |
      day %% plot_every_n_days == 0L
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 5. Read the final output when available
# -----------------------------------------------------------------------------

final_results <- NULL

if (file.exists(final_file)) {
  final_results <- read_csv(
    final_file,
    show_col_types = FALSE,
    progress = FALSE
  )

  if ("target_mean_npp" %in% names(final_results)) {
    final_results <- final_results %>%
      mutate(
        scenario = paste0(
          "Case ", npp_case,
          " | target mean NPP = ",
          sprintf("%.2f", target_mean_npp)
        )
      )
  } else if ("npp_rate" %in% names(final_results)) {
    final_results <- final_results %>%
      mutate(
        scenario = paste0(
          "Case ", npp_case,
          " | NPP = ",
          sprintf("%.2f", npp_rate)
        )
      )
  } else {
    final_results <- final_results %>%
      mutate(scenario = paste0("Case ", npp_case))
  }
} else {
  warning(
    paste0(
      "Final output file not found: ", final_file,
      "\nTime-series plots will still be produced, but final-state plots will be skipped."
    ),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 6. Common plot theme and saving function
# -----------------------------------------------------------------------------

plot_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    plot.title.position = "plot"
  )

save_plot <- function(plot_object, filename, width = 11, height = 7) {
  ggsave(
    filename = file.path(output_directory, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = png_resolution
  )
}

# -----------------------------------------------------------------------------
# 7. NPP time series
# -----------------------------------------------------------------------------

plot_npp <- ggplot(
  daily_plot,
  aes(x = year, y = npp_rate, color = scenario)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.45) +
  labs(
    title = "Annualized NPP rate through time",
    subtitle = "The Fortran module multiplies this rate by dt_years to obtain carbon input per timestep.",
    x = "Simulation year",
    y = "Annualized NPP rate"
  ) +
  plot_theme

save_plot(plot_npp, "01_npp_rate_time_series.png")

plot_npp_daily <- ggplot(
  daily_plot,
  aes(x = year, y = npp_daily, color = scenario)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.45) +
  labs(
    title = "Carbon input from NPP per timestep",
    x = "Simulation year",
    y = "NPP carbon input per timestep"
  ) +
  plot_theme

save_plot(plot_npp_daily, "02_npp_daily_time_series.png")

# -----------------------------------------------------------------------------
# 8. Individual carbon pools
# -----------------------------------------------------------------------------

carbon_pool_long <- daily_plot %>%
  select(
    year,
    scenario,
    leaf,
    root,
    sapwood,
    heartwood,
    storage
  ) %>%
  pivot_longer(
    cols = c(leaf, root, sapwood, heartwood, storage),
    names_to = "pool",
    values_to = "carbon"
  ) %>%
  mutate(
    pool = recode(
      pool,
      leaf = "Leaf",
      root = "Fine root",
      sapwood = "Sapwood",
      heartwood = "Heartwood",
      storage = "Labile storage"
    )
  )

plot_carbon_pools <- ggplot(
  carbon_pool_long,
  aes(x = year, y = carbon, color = scenario)
) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~pool, scales = "free_y", ncol = 2) +
  labs(
    title = "Carbon pools through time",
    x = "Simulation year",
    y = "Carbon"
  ) +
  plot_theme

save_plot(plot_carbon_pools, "03_carbon_pools_time_series.png", height = 9)

# -----------------------------------------------------------------------------
# 9. Aggregate carbon pools
# -----------------------------------------------------------------------------

aggregate_carbon_long <- daily_plot %>%
  select(
    year,
    scenario,
    living_carbon,
    structural_carbon,
    total_plant_carbon
  ) %>%
  pivot_longer(
    cols = c(living_carbon, structural_carbon, total_plant_carbon),
    names_to = "carbon_group",
    values_to = "carbon"
  ) %>%
  mutate(
    carbon_group = recode(
      carbon_group,
      living_carbon = "Living structural carbon",
      structural_carbon = "Total structural carbon",
      total_plant_carbon = "Total plant carbon including storage"
    )
  )

plot_aggregate_carbon <- ggplot(
  aggregate_carbon_long,
  aes(x = year, y = carbon, color = scenario)
) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~carbon_group, scales = "free_y", ncol = 1) +
  labs(
    title = "Aggregate plant carbon through time",
    x = "Simulation year",
    y = "Carbon"
  ) +
  plot_theme

save_plot(plot_aggregate_carbon, "04_aggregate_carbon_time_series.png", height = 9)

# -----------------------------------------------------------------------------
# 10. Total structural demand and realized allocation
# -----------------------------------------------------------------------------

demand_allocation_long <- daily_plot %>%
  select(
    year,
    scenario,
    total_demand,
    carbon_allocated
  ) %>%
  pivot_longer(
    cols = c(total_demand, carbon_allocated),
    names_to = "variable",
    values_to = "carbon"
  ) %>%
  mutate(
    variable = recode(
      variable,
      total_demand = "Total structural demand",
      carbon_allocated = "Carbon actually allocated"
    )
  )

plot_demand_allocation <- ggplot(
  demand_allocation_long,
  aes(
    x = year,
    y = carbon,
    color = variable,
    linetype = scenario
  )
) +
  geom_line(linewidth = 0.55) +
  facet_wrap(~scenario, scales = "free_y", ncol = 1) +
  labs(
    title = "Structural demand and realized allocation",
    subtitle = "Allocation can be lower than demand because of limited storage or construction capacity.",
    x = "Simulation year",
    y = "Carbon per timestep"
  ) +
  plot_theme +
  theme(legend.title = element_blank())

save_plot(plot_demand_allocation, "05_demand_vs_allocation.png", height = 8)

# -----------------------------------------------------------------------------
# 11. Pool-specific structural demands
# -----------------------------------------------------------------------------

pool_demand_long <- daily_plot %>%
  select(
    year,
    scenario,
    leaf_demand,
    root_demand,
    sapwood_demand
  ) %>%
  pivot_longer(
    cols = c(leaf_demand, root_demand, sapwood_demand),
    names_to = "pool",
    values_to = "demand"
  ) %>%
  mutate(
    pool = recode(
      pool,
      leaf_demand = "Leaf demand",
      root_demand = "Fine-root demand",
      sapwood_demand = "Sapwood demand"
    )
  )

plot_pool_demand <- ggplot(
  pool_demand_long,
  aes(x = year, y = demand, color = scenario)
) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~pool, scales = "free_y", ncol = 1) +
  labs(
    title = "Pool-specific structural demands",
    x = "Simulation year",
    y = "Requested carbon per timestep"
  ) +
  plot_theme

save_plot(plot_pool_demand, "06_pool_specific_demands.png", height = 9)

# -----------------------------------------------------------------------------
# 12. Realized structural increments
# -----------------------------------------------------------------------------

increment_long <- daily_plot %>%
  select(
    year,
    scenario,
    delta_leaf,
    delta_root,
    delta_sapwood
  ) %>%
  pivot_longer(
    cols = c(delta_leaf, delta_root, delta_sapwood),
    names_to = "pool",
    values_to = "increment"
  ) %>%
  mutate(
    pool = recode(
      pool,
      delta_leaf = "Leaf increment",
      delta_root = "Fine-root increment",
      delta_sapwood = "Sapwood increment"
    )
  )

plot_increments <- ggplot(
  increment_long,
  aes(x = year, y = increment, color = scenario)
) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~pool, scales = "free_y", ncol = 1) +
  labs(
    title = "Realized structural increments",
    x = "Simulation year",
    y = "Carbon increment per timestep"
  ) +
  plot_theme

save_plot(plot_increments, "07_realized_structural_increments.png", height = 9)

# -----------------------------------------------------------------------------
# 13. Turnover losses
# -----------------------------------------------------------------------------

turnover_long <- daily_plot %>%
  select(
    year,
    scenario,
    leaf_turnover,
    root_turnover,
    sapwood_turnover,
    storage_turnover,
    heartwood_turnover
  ) %>%
  pivot_longer(
    cols = c(
      leaf_turnover,
      root_turnover,
      sapwood_turnover,
      storage_turnover,
      heartwood_turnover
    ),
    names_to = "pool",
    values_to = "turnover"
  ) %>%
  mutate(
    pool = recode(
      pool,
      leaf_turnover = "Leaf turnover",
      root_turnover = "Fine-root turnover",
      sapwood_turnover = "Sapwood turnover",
      storage_turnover = "Storage turnover",
      heartwood_turnover = "Heartwood turnover"
    )
  )

plot_turnover <- ggplot(
  turnover_long,
  aes(x = year, y = turnover, color = scenario)
) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~pool, scales = "free_y", ncol = 2) +
  labs(
    title = "Turnover fluxes through time",
    x = "Simulation year",
    y = "Carbon turnover per timestep"
  ) +
  plot_theme

save_plot(plot_turnover, "08_turnover_time_series.png", height = 9)

# -----------------------------------------------------------------------------
# 14. Carbon deficits associated with negative NPP
# -----------------------------------------------------------------------------

deficit_long <- daily_plot %>%
  select(
    year,
    scenario,
    unmet_storage_deficit,
    unpaid_carbon_deficit
  ) %>%
  pivot_longer(
    cols = c(unmet_storage_deficit, unpaid_carbon_deficit),
    names_to = "deficit_type",
    values_to = "deficit"
  ) %>%
  mutate(
    deficit_type = recode(
      deficit_type,
      unmet_storage_deficit = "Storage deficit after negative NPP",
      unpaid_carbon_deficit = "Carbon deficit remaining unpaid"
    )
  )

plot_deficits <- ggplot(
  deficit_long,
  aes(x = year, y = deficit, color = scenario)
) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~deficit_type, scales = "free_y", ncol = 1) +
  labs(
    title = "Carbon deficits and starvation-related diagnostics",
    x = "Simulation year",
    y = "Carbon deficit per timestep"
  ) +
  plot_theme

save_plot(plot_deficits, "09_carbon_deficits.png", height = 7)

# -----------------------------------------------------------------------------
# 15. Allometric residuals
# -----------------------------------------------------------------------------

residual_long <- daily_plot %>%
  select(
    year,
    scenario,
    leaf_root_residual,
    pipe_model_residual
  ) %>%
  pivot_longer(
    cols = c(leaf_root_residual, pipe_model_residual),
    names_to = "residual_type",
    values_to = "residual"
  ) %>%
  mutate(
    residual_type = recode(
      residual_type,
      leaf_root_residual = "Leaf-root residual",
      pipe_model_residual = "Pipe-model residual"
    )
  )

plot_residuals <- ggplot(
  residual_long,
  aes(x = year, y = residual, color = scenario)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.55) +
  facet_wrap(~residual_type, scales = "free_y", ncol = 1) +
  labs(
    title = "Allometric residuals through time",
    subtitle = "Zero represents exact agreement with the corresponding target relationship.",
    x = "Simulation year",
    y = "Residual"
  ) +
  plot_theme

save_plot(plot_residuals, "10_allometric_residuals.png", height = 7)

# -----------------------------------------------------------------------------
# 16. Plant height
# -----------------------------------------------------------------------------

plot_height <- ggplot(
  daily_plot,
  aes(x = year, y = height, color = scenario)
) +
  geom_line(linewidth = 0.65) +
  labs(
    title = "Plant height through time",
    x = "Simulation year",
    y = "Height"
  ) +
  plot_theme

save_plot(plot_height, "11_height_time_series.png")

# -----------------------------------------------------------------------------
# 17. Final-state plots
# -----------------------------------------------------------------------------

final_plots <- list()

if (!is.null(final_results)) {
  final_pool_columns <- c(
    "final_leaf",
    "final_root",
    "final_sapwood",
    "final_heartwood",
    "final_storage"
  )

  if (all(final_pool_columns %in% names(final_results))) {
    final_pool_long <- final_results %>%
      select(scenario, all_of(final_pool_columns)) %>%
      pivot_longer(
        cols = all_of(final_pool_columns),
        names_to = "pool",
        values_to = "carbon"
      ) %>%
      mutate(
        pool = recode(
          pool,
          final_leaf = "Leaf",
          final_root = "Fine root",
          final_sapwood = "Sapwood",
          final_heartwood = "Heartwood",
          final_storage = "Labile storage"
        )
      )

    plot_final_pools <- ggplot(
      final_pool_long,
      aes(x = scenario, y = carbon, fill = scenario)
    ) +
      geom_col(show.legend = FALSE) +
      facet_wrap(~pool, scales = "free_y", ncol = 2) +
      labs(
        title = "Final carbon pools by NPP scenario",
        x = NULL,
        y = "Final carbon"
      ) +
      plot_theme +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1)
      )

    save_plot(plot_final_pools, "12_final_carbon_pools.png", height = 9)
    final_plots <- append(final_plots, list(plot_final_pools))
  }

  final_total_columns <- c(
    "final_living_carbon",
    "final_structural_carbon",
    "final_total_plant_carbon"
  )

  if (all(final_total_columns %in% names(final_results))) {
    final_total_long <- final_results %>%
      select(scenario, all_of(final_total_columns)) %>%
      pivot_longer(
        cols = all_of(final_total_columns),
        names_to = "carbon_group",
        values_to = "carbon"
      ) %>%
      mutate(
        carbon_group = recode(
          carbon_group,
          final_living_carbon = "Living structural carbon",
          final_structural_carbon = "Total structural carbon",
          final_total_plant_carbon = "Total plant carbon including storage"
        )
      )

    plot_final_totals <- ggplot(
      final_total_long,
      aes(x = scenario, y = carbon, fill = scenario)
    ) +
      geom_col(show.legend = FALSE) +
      facet_wrap(~carbon_group, scales = "free_y", ncol = 1) +
      labs(
        title = "Final aggregate carbon by NPP scenario",
        x = NULL,
        y = "Final carbon"
      ) +
      plot_theme +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1)
      )

    save_plot(plot_final_totals, "13_final_aggregate_carbon.png", height = 9)
    final_plots <- append(final_plots, list(plot_final_totals))
  }
}

# -----------------------------------------------------------------------------
# 18. Save all figures in one multi-page PDF
# -----------------------------------------------------------------------------

time_series_plots <- list(
  plot_npp,
  plot_npp_daily,
  plot_carbon_pools,
  plot_aggregate_carbon,
  plot_demand_allocation,
  plot_pool_demand,
  plot_increments,
  plot_turnover,
  plot_deficits,
  plot_residuals,
  plot_height
)

all_plots <- c(time_series_plots, final_plots)

pdf(
  file = file.path(output_directory, "all_time_series_plots.pdf"),
  width = 11,
  height = 7,
  onefile = TRUE
)

for (plot_object in all_plots) {
  print(plot_object)
}

dev.off()

# -----------------------------------------------------------------------------
# 19. Report generated files
# -----------------------------------------------------------------------------

message("Plots were written to: ", normalizePath(output_directory))
message(
  "Multi-page PDF: ",
  normalizePath(file.path(output_directory, "all_time_series_plots.pdf"))
)

