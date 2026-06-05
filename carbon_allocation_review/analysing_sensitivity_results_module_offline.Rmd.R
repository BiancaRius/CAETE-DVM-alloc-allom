# Load packages
library(dplyr)
library(readr)
library(ggplot2)

# Define base directory
base_dir <- "~/Desktop/CAETE-DVM-alloc-allom-including_alloc2_Cm2/carbon_allocation_review/"

## __________________________________________________
## Step 1: Check if the text is numerically trustable
## - Before interpreting biology, you need to check
##  if the code ran without any error
## __________________________________________________






# # Load packages
# library(dplyr)
# library(readr)
# library(ggplot2)
# 
# base_dir <- "~/Desktop/CAETE-DVM-alloc-allom-including_alloc2_Cm2/carbon_allocation_review/"
# 
# # Read the scenario-level summary file
# summary_df <- read_csv("~/Desktop/CAETE-DVM-alloc-allom-including_alloc2_Cm2/carbon_allocation_review/storage_allocation_sensitivity_summary_turnover.csv", show_col_types = FALSE)
# 
# # Check numerical status
# summary_df %>%
#   count(status, failure_reason)
# 
# # Create final living carbon if it is missing
# summary_df <- summary_df %>%
#   mutate(
#     final_living_carbon = final_leaf + final_root + final_sapwood
#   )
# 
# # Summarise main outcomes by NPP level
# summary_by_npp <- summary_df %>%
#   group_by(npp_rate) %>%
#   summarise(
#     n_scenarios = n(),
#     n_pass = sum(status == "PASS"),
#     mean_final_living = mean(final_living_carbon, na.rm = TRUE),
#     mean_final_storage = mean(final_storage, na.rm = TRUE),
#     mean_storage_frac_living = mean(final_storage_fraction_living, na.rm = TRUE),
#     max_storage_frac_living = max(max_storage_fraction_living, na.rm = TRUE),
#     mean_turnover_loss = mean(cumulative_turnover_carbon_loss, na.rm = TRUE),
#     mean_starvation_loss = mean(cumulative_starvation_carbon_loss, na.rm = TRUE),
#     mean_unpaid_deficit = mean(cumulative_unpaid_carbon_deficit, na.rm = TRUE),
#     mean_days_with_starvation = mean(days_with_unmet_storage_deficit, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# print(summary_by_npp)
# 
# # Summarise the effect of background demand
# summary_by_npp_background <- summary_df %>%
#   group_by(npp_rate, background_mode) %>%
#   summarise(
#     n_scenarios = n(),
#     mean_final_living = mean(final_living_carbon, na.rm = TRUE),
#     mean_final_storage = mean(final_storage, na.rm = TRUE),
#     mean_storage_frac_living = mean(final_storage_fraction_living, na.rm = TRUE),
#     max_storage_frac_living = max(max_storage_fraction_living, na.rm = TRUE),
#     mean_structural_allocation = mean(cumulative_structural_allocation, na.rm = TRUE),
#     mean_turnover_loss = mean(cumulative_turnover_carbon_loss, na.rm = TRUE),
#     mean_days_storage_above_10pct = mean(days_with_storage_above_10pct_living, na.rm = TRUE),
#     mean_days_storage_above_50pct = mean(days_with_storage_above_50pct_living, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# print(summary_by_npp_background)
# 
# # Classify scenarios into broad outcome classes
# summary_classified <- summary_df %>%
#   mutate(
#     outcome_class = case_when(
#       status != "PASS" ~ "numerical_failure",
#       days_with_unmet_storage_deficit > 0 ~ "carbon_starvation",
#       final_storage_fraction_living > 0.5 ~ "excessive_storage",
#       net_living_carbon_change < 0 ~ "living_biomass_decline",
#       TRUE ~ "plausible_or_positive_growth"
#     )
#   )
# 
# # Count outcome classes by NPP and background mode
# summary_classified %>%
#   count(npp_rate, background_mode, outcome_class)
# 
# # To save this table, ask the name the user wants to put
# file_name <- readline(prompt = "How do you want to name this table (ex: summary_classified_turnover_storageoff.csv): ")
# 
# # create the path
# path_save <- file.path(base_dir, file_name)
# 
# # save dataframe
# write_csv(summary_classified, path_save)
# 
# # Confirm saving
# message("Saved sucssefully: ", path_save)
# 
# # Summarise allocation plausibility by NPP and background mode
# allocation_diagnostics <- summary_df %>%
#   group_by(npp_rate, background_mode) %>%
#   summarise(
#     n_scenarios = n(),
#     
#     # Carbon allocation diagnostics
#     mean_structural_allocation = mean(cumulative_structural_allocation, na.rm = TRUE),
#     mean_positive_npp = mean(cumulative_positive_npp, na.rm = TRUE),
#     mean_structural_fraction = mean(structural_fraction_of_positive_npp, na.rm = TRUE),
#     
#     # Final carbon state
#     mean_final_living = mean(final_living_carbon, na.rm = TRUE),
#     mean_living_change = mean(net_living_carbon_change, na.rm = TRUE),
#     mean_final_storage = mean(final_storage, na.rm = TRUE),
#     mean_storage_frac_living = mean(final_storage_fraction_living, na.rm = TRUE),
#     
#     # Turnover and starvation
#     mean_turnover_loss = mean(cumulative_turnover_carbon_loss, na.rm = TRUE),
#     mean_starvation_loss = mean(cumulative_starvation_carbon_loss, na.rm = TRUE),
#     mean_unpaid_deficit = mean(cumulative_unpaid_carbon_deficit, na.rm = TRUE),
#     
#     # Allometric residuals
#     mean_abs_leaf_root_residual = mean(abs(final_leaf_root_residual), na.rm = TRUE),
#     mean_abs_pipe_residual = mean(abs(final_pipe_residual), na.rm = TRUE),
#     max_abs_leaf_root_residual = max(abs(final_leaf_root_residual), na.rm = TRUE),
#     max_abs_pipe_residual = max(abs(final_pipe_residual), na.rm = TRUE),
#     
#     .groups = "drop"
#   )
# 
# allocation_diagnostics
# 
# # Classify scenarios according to allocation plausibility
# summary_allocation_classes <- summary_df %>%
#   mutate(
#     allocation_class = case_when(
#       status != "PASS" ~ "numerical_failure",
#       
#       abs(final_leaf_root_residual) > 1.0 |
#         abs(final_pipe_residual) > 5.0 ~ "large_allometric_residual",
#       
#       final_storage_fraction_living > 1.0 ~ "excessive_storage",
#       
#       npp_rate > 0 &
#         structural_fraction_of_positive_npp < 0.25 ~ "low_structural_use_of_npp",
#       
#       npp_rate > 0 &
#         net_living_carbon_change < 0 ~ "positive_npp_but_living_decline",
#       
#       cumulative_starvation_carbon_loss > 0 &
#         npp_rate >= 0 ~ "unexpected_starvation",
#       
#       TRUE ~ "apparently_plausible"
#     )
#   )
# 
# summary_allocation_classes %>%
#   count(npp_rate, background_mode, allocation_class)
# 
# # Plot final pipe residual across NPP and background mode
# ggplot(summary_df, aes(x = factor(npp_rate), y = final_pipe_residual)) +
#   geom_boxplot() +
#   facet_wrap(~ background_mode) +
#   theme_bw() +
#   labs(
#     x = "NPP rate",
#     y = "Final pipe residual",
#     title = "Final pipe-model residual across NPP levels"
#   )
# 
# # Plot final leaf-root residual across NPP and background mode
# ggplot(summary_df, aes(x = factor(npp_rate), y = final_leaf_root_residual)) +
#   geom_boxplot() +
#   facet_wrap(~ background_mode) +
#   theme_bw() +
#   labs(
#     x = "NPP rate",
#     y = "Final leaf-root residual",
#     title = "Final leaf-root residual across NPP levels"
#   )
# 
#   # Plot final storage fraction by NPP and background mode
# ggplot(summary_df, aes(x = factor(npp_rate), y = final_storage_fraction_living)) +
#   geom_boxplot() +
#   facet_wrap(~ background_mode) +
#   theme_bw() +
#   labs(
#     x = "NPP rate",
#     y = "Final storage / living carbon",
#     title = "Storage accumulation across NPP levels and background-demand modes"
#   )
# 
# # Plot the fraction of positive NPP allocated to structure
# ggplot(summary_df, aes(x = factor(npp_rate), y = structural_fraction_of_positive_npp)) +
#   geom_boxplot() +
#   facet_wrap(~ background_mode) +
#   theme_bw() +
#   labs(
#     x = "NPP rate",
#     y = "Structural allocation / positive NPP",
#     title = "Fraction of positive NPP allocated to structural growth"
#   )
# # Plot final living carbon by NPP and background mode
# ggplot(summary_df, aes(x = factor(npp_rate), y = final_living_carbon)) +
#   geom_boxplot() +
#   facet_wrap(~ background_mode) +
#   theme_bw() +
#   labs(
#     x = "NPP rate",
#     y = "Final living carbon",
#     title = "Final living carbon across NPP levels and background-demand modes"
#   )
# 
# # Plot turnover loss by NPP
# ggplot(summary_df, aes(x = factor(npp_rate), y = cumulative_turnover_carbon_loss)) +
#   geom_boxplot() +
#   theme_bw() +
#   labs(
#     x = "NPP rate",
#     y = "Cumulative turnover carbon loss",
#     title = "Cumulative tissue turnover across NPP levels"
#   )
# 
