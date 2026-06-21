#!/usr/bin/env Rscript

# Display artifacts for the participant-level ADEMP simulation.
# This script uses the June 19 main simulation outputs. The older June 15
# display script belongs to the earlier maxT-centered simulation branch.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

project_root <- normalizePath(getwd(), mustWork = TRUE)
results_dir <- file.path(project_root, "simulation", "results")
figure_dir <- file.path(project_root, "figures", "simulation")
table_dir <- file.path(project_root, "tables", "simulation")
source_dir <- file.path(figure_dir, "source_data")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

date_tag <- "2026-06-19"
final_table_tag <- "2026-06-20"

read_result <- function(filename) {
  read_csv(file.path(results_dir, filename), show_col_types = FALSE)
}

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

label_scenario <- function(x) {
  dplyr::recode(
    x,
    global_null = "Global null",
    partial_null_AB_minus_A = "Partial null: first contrast",
    partial_null_AB_minus_B = "Partial null: second contrast",
    balanced_positive = "Balanced positive",
    limiting_AB_minus_A = "First contrast harder",
    limiting_AB_minus_B = "Second contrast harder",
    positive_margin = "Positive margin",
    unequal_variance_positive = "Unequal variance",
    heavy_tail_global_null = "Heavy-tailed null",
    overoptimistic_planning = "Overoptimistic planning",
    reversed_limiting_component = "Harder contrast reversed",
    additive_positive_excess = "Additive positive excess",
    .default = x
  )
}

label_method <- function(x) {
  dplyr::recode(
    x,
    equal_iut_full = "Equal allocation, IUT",
    shared_arm_iut_full = "Shared-arm, IUT",
    plugin_iut_full = "Plug-in allocation, IUT",
    conservative_iut_full = "Conservative allocation, IUT",
    marginal_iut_full = "Marginal comparator, IUT",
    equal_maxT_separate = "Equal allocation, maxT separate claims",
    .default = x
  )
}

method_order <- c(
  "Equal allocation, IUT",
  "Shared-arm, IUT",
  "Plug-in allocation, IUT",
  "Conservative allocation, IUT",
  "Marginal comparator, IUT",
  "Equal allocation, maxT separate claims"
)

main_results <- read_result("participant_ademp_main_2026-06-19_operating_characteristics.csv")
scenario_grid <- read_result("participant_ademp_main_2026-06-19_scenario_grid.csv")
method_specs <- read_result("participant_ademp_main_2026-06-19_method_specs.csv")
validation <- read_result("participant_ademp_main_2026-06-19_validation_checks.csv")
error_control <- read_result("participant_ademp_main_2026-06-19_table_error_control.csv")
joint_success <- read_result("participant_ademp_main_2026-06-19_table_joint_success.csv")
fragility <- read_result("participant_ademp_main_2026-06-19_table_allocation_fragility.csv")
robustness <- read_result("participant_ademp_main_2026-06-19_table_robustness.csv")

if (!all(validation$passed)) {
  stop("Main simulation validation checks did not all pass.", call. = FALSE)
}

primary <- main_results %>%
  filter(abs(alpha - 0.025) < 1e-12)

figure_3a_source <- primary %>%
  filter(
    claim_regime == "full_contribution",
    scenario_family %in% c("global_null", "partial_null")
  ) %>%
  transmute(
    scenario_label,
    scenario_display = label_scenario(scenario_label),
    method_id,
    method_display = label_method(method_id),
    empirical_error = empirical_false_full_contribution,
    lower = empirical_false_full_contribution_lcl,
    upper = empirical_false_full_contribution_ucl,
    analytic_error = analytic_false_full_contribution,
    nrep
  ) %>%
  mutate(
    scenario_display = factor(
      scenario_display,
      levels = c(
        "Global null",
        "Partial null: first contrast",
        "Partial null: second contrast"
      )
    ),
    method_display = factor(method_display, levels = method_order)
  ) %>%
  arrange(scenario_display, method_display)

figure_3b_source <- primary %>%
  filter(
    claim_regime == "full_contribution",
    scenario_label %in% c(
      "balanced_positive",
      "limiting_AB_minus_A",
      "limiting_AB_minus_B",
      "positive_margin",
      "additive_positive_excess"
    ),
    method_id %in% c(
      "equal_iut_full",
      "shared_arm_iut_full",
      "plugin_iut_full",
      "conservative_iut_full"
    )
  ) %>%
  transmute(
    scenario_label,
    scenario_display = label_scenario(scenario_label),
    method_id,
    method_display = label_method(method_id),
    empirical_joint_success,
    lower = empirical_joint_success_lcl,
    upper = empirical_joint_success_ucl,
    analytic_joint_success,
    allocation_A,
    allocation_B,
    allocation_AB,
    nrep
  ) %>%
  mutate(
    scenario_display = factor(
      scenario_display,
      levels = c(
        "Balanced positive",
        "First contrast harder",
        "Second contrast harder",
        "Positive margin",
        "Additive positive excess"
      )
    ),
    method_display = factor(method_display, levels = method_order)
  ) %>%
  arrange(scenario_display, method_display)

figure_3c_source <- primary %>%
  filter(
    claim_regime == "full_contribution",
    scenario_label %in% c("overoptimistic_planning", "reversed_limiting_component"),
    method_id %in% c(
      "equal_iut_full",
      "plugin_iut_full",
      "conservative_iut_full"
    )
  ) %>%
  transmute(
    scenario_label,
    scenario_display = label_scenario(scenario_label),
    method_id,
    method_display = label_method(method_id),
    truth_delta_AB_minus_A,
    truth_delta_AB_minus_B,
    plan_delta_AB_minus_A,
    plan_delta_AB_minus_B,
    empirical_joint_success,
    lower = empirical_joint_success_lcl,
    upper = empirical_joint_success_ucl,
    analytic_joint_success,
    allocation_A,
    allocation_B,
    allocation_AB,
    nrep
  ) %>%
  mutate(
    scenario_display = factor(
      scenario_display,
      levels = c("Overoptimistic planning", "Harder contrast reversed")
    ),
    method_display = factor(method_display, levels = method_order)
  ) %>%
  arrange(scenario_display, method_display)

max_iut_error <- primary %>%
  filter(
    claim_regime == "full_contribution",
    scenario_family %in% c("global_null", "partial_null")
  ) %>%
  arrange(desc(empirical_false_full_contribution)) %>%
  slice(1)

max_maxT_error <- primary %>%
  filter(claim_regime == "separate_component_claims", scenario_family == "global_null") %>%
  arrange(desc(empirical_any_true_reject)) %>%
  slice(1)

row_for <- function(label, method) {
  primary %>% filter(scenario_label == label, method_id == method) %>% slice(1)
}

main_table_source <- bind_rows(
  tibble(
    result_type = "Error control",
    comparison = "IUT full contribution, worst primary null row",
    scenario = label_scenario(max_iut_error$scenario_label),
    method = label_method(max_iut_error$method_id),
    estimate = max_iut_error$empirical_false_full_contribution,
    lower = max_iut_error$empirical_false_full_contribution_lcl,
    upper = max_iut_error$empirical_false_full_contribution_ucl,
    reference = max_iut_error$analytic_false_full_contribution
  ),
  tibble(
    result_type = "Error control",
    comparison = "maxT separate claims, global null",
    scenario = label_scenario(max_maxT_error$scenario_label),
    method = label_method(max_maxT_error$method_id),
    estimate = max_maxT_error$empirical_any_true_reject,
    lower = max_maxT_error$empirical_any_true_reject_lcl,
    upper = max_maxT_error$empirical_any_true_reject_ucl,
    reference = max_maxT_error$analytic_any_true_reject
  ),
  bind_rows(
    row_for("balanced_positive", "equal_iut_full"),
    row_for("balanced_positive", "shared_arm_iut_full"),
    row_for("balanced_positive", "plugin_iut_full"),
    row_for("limiting_AB_minus_A", "equal_iut_full"),
    row_for("limiting_AB_minus_A", "shared_arm_iut_full"),
    row_for("limiting_AB_minus_A", "plugin_iut_full"),
    row_for("limiting_AB_minus_B", "equal_iut_full"),
    row_for("limiting_AB_minus_B", "shared_arm_iut_full"),
    row_for("limiting_AB_minus_B", "plugin_iut_full"),
    row_for("reversed_limiting_component", "equal_iut_full"),
    row_for("reversed_limiting_component", "shared_arm_iut_full"),
    row_for("reversed_limiting_component", "plugin_iut_full"),
    row_for("reversed_limiting_component", "conservative_iut_full"),
    row_for("additive_positive_excess", "plugin_iut_full")
  ) %>%
    transmute(
      result_type = "Joint power",
      comparison = case_when(
        scenario_label == "balanced_positive" ~ "Balanced positive",
        scenario_label == "limiting_AB_minus_A" ~ "First contrast harder, correctly planned",
        scenario_label == "limiting_AB_minus_B" ~ "Second contrast harder, correctly planned",
        scenario_label == "reversed_limiting_component" ~ "Harder contrast reversed",
        scenario_label == "additive_positive_excess" ~ "Additive positive excess",
        TRUE ~ label_scenario(scenario_label)
      ),
      scenario = label_scenario(scenario_label),
      method = label_method(method_id),
      estimate = empirical_joint_success,
      lower = empirical_joint_success_lcl,
      upper = empirical_joint_success_ucl,
      reference = analytic_joint_success
    )
)

main_table <- main_table_source %>%
  mutate(
    estimate = format_num(estimate),
    monte_carlo_interval = paste0("(", format_num(lower), ", ", format_num(upper), ")"),
    analytic_reference = format_num(reference)
  ) %>%
  select(
    result_type,
    comparison,
    scenario,
    method,
    estimate,
    monte_carlo_interval,
    analytic_reference
  )

round_for_supp <- function(df) {
  df %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

source_paths <- c(
  figure_3a_source = file.path(source_dir, paste0("figure_3a_error_control_source_", date_tag, ".csv")),
  figure_3b_source = file.path(source_dir, paste0("figure_3b_joint_success_source_", date_tag, ".csv")),
  figure_3c_source = file.path(source_dir, paste0("figure_3c_allocation_fragility_source_", date_tag, ".csv")),
  table_3_source = file.path(source_dir, paste0("table_3_simulation_main_source_", date_tag, ".csv"))
)

write_csv(round_for_supp(figure_3a_source), source_paths[["figure_3a_source"]])
write_csv(round_for_supp(figure_3b_source), source_paths[["figure_3b_source"]])
write_csv(round_for_supp(figure_3c_source), source_paths[["figure_3c_source"]])
write_csv(round_for_supp(main_table_source), source_paths[["table_3_source"]])

main_table_path <- file.path(table_dir, paste0("table_3_simulation_main_", date_tag, ".csv"))
write_csv(main_table, main_table_path)
final_main_table_path <- file.path(table_dir, paste0("table_3_simulation_main_", final_table_tag, ".csv"))
write_csv(main_table, final_main_table_path)

supp_paths <- c(
  methods = file.path(table_dir, paste0("supp_table_simulation_methods_", date_tag, ".csv")),
  scenarios = file.path(table_dir, paste0("supp_table_simulation_scenario_grid_", date_tag, ".csv")),
  error = file.path(table_dir, paste0("supp_table_simulation_error_control_full_", date_tag, ".csv")),
  joint = file.path(table_dir, paste0("supp_table_simulation_joint_success_full_", date_tag, ".csv")),
  fragility = file.path(table_dir, paste0("supp_table_simulation_allocation_fragility_full_", date_tag, ".csv")),
  robustness = file.path(table_dir, paste0("supp_table_simulation_robustness_full_", date_tag, ".csv")),
  validation = file.path(table_dir, paste0("supp_table_simulation_validation_checks_", date_tag, ".csv"))
)

write_csv(method_specs, supp_paths[["methods"]])
write_csv(round_for_supp(scenario_grid), supp_paths[["scenarios"]])
write_csv(round_for_supp(error_control), supp_paths[["error"]])
write_csv(round_for_supp(joint_success), supp_paths[["joint"]])
write_csv(round_for_supp(fragility), supp_paths[["fragility"]])
write_csv(round_for_supp(robustness), supp_paths[["robustness"]])
write_csv(validation, supp_paths[["validation"]])

palette <- c(
  "Equal allocation, IUT" = "#1B5A7A",
  "Shared-arm, IUT" = "#4C4C4C",
  "Plug-in allocation, IUT" = "#7A3E8E",
  "Conservative allocation, IUT" = "#2F6B45",
  "Marginal comparator, IUT" = "#8A5A16",
  "Equal allocation, maxT separate claims" = "#8E2F3F"
)

theme_display <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "#000000"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.title = element_text(color = "#000000"),
      axis.text = element_text(color = "#000000"),
      strip.text = element_text(color = "#000000", face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(color = "#000000"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "#D8D8D8", linewidth = 0.30),
      plot.margin = margin(10, 16, 10, 10)
    )
}

figure_3a <- ggplot(
  figure_3a_source,
  aes(
    x = empirical_error,
    y = method_display,
    xmin = lower,
    xmax = upper,
    color = method_display
  )
) +
  geom_vline(xintercept = 0.025, linewidth = 0.45, linetype = "dashed", color = "#000000") +
  geom_errorbar(orientation = "y", width = 0.16, linewidth = 0.55) +
  geom_point(size = 2.5) +
  facet_wrap(~ scenario_display, ncol = 1) +
  scale_x_continuous(
    limits = c(0, 0.032),
    breaks = seq(0, 0.03, by = 0.01),
    labels = format_num
  ) +
  scale_color_manual(values = palette, drop = TRUE) +
  labs(x = "False full-contribution conclusion", y = NULL) +
  theme_display() +
  theme(legend.position = "none")

figure_3b <- ggplot(
  figure_3b_source,
  aes(
    x = empirical_joint_success,
    y = scenario_display,
    xmin = lower,
    xmax = upper,
    color = method_display
  )
) +
  geom_errorbar(orientation = "y", width = 0.16, linewidth = 0.55,
                position = position_dodge(width = 0.55)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.55)) +
  scale_x_continuous(
    limits = c(0.40, 0.82),
    breaks = seq(0.40, 0.80, by = 0.10),
    labels = format_num
  ) +
  scale_color_manual(values = palette, drop = TRUE) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Joint power", y = NULL) +
  theme_display()

figure_3c <- ggplot(
  figure_3c_source,
  aes(
    x = empirical_joint_success,
    y = method_display,
    xmin = lower,
    xmax = upper,
    color = method_display
  )
) +
  geom_errorbar(orientation = "y", width = 0.16, linewidth = 0.55) +
  geom_point(size = 2.6) +
  facet_wrap(~ scenario_display, ncol = 1) +
  scale_x_continuous(
    limits = c(0.28, 0.40),
    breaks = seq(0.28, 0.40, by = 0.04),
    labels = format_num
  ) +
  scale_color_manual(values = palette, drop = TRUE) +
  labs(x = "Joint power", y = NULL) +
  theme_display() +
  theme(legend.position = "none")

save_plot <- function(plot, stem, width, height) {
  pdf_path <- file.path(figure_dir, paste0(stem, "_", date_tag, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, "_", date_tag, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = grDevices::pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      png_path, plot,
      width = width, height = height, units = "in", dpi = 320,
      device = ragg::agg_png
    )
  } else {
    ggsave(
      png_path, plot,
      width = width, height = height, units = "in", dpi = 320,
      device = grDevices::png
    )
  }
  c(pdf = pdf_path, png = png_path)
}

figure_paths <- c(
  save_plot(figure_3a, "figure_3a_error_control", width = 7.2, height = 5.2),
  save_plot(figure_3b, "figure_3b_joint_success", width = 7.2, height = 5.2),
  save_plot(figure_3c, "figure_3c_allocation_fragility", width = 7.2, height = 3.8)
)

manifest <- tibble(
  artifact_type = c(
    rep("figure", 3),
    rep("source_data", 4),
    "main_table",
    rep("supplement_table", length(supp_paths)),
    "qa"
  ),
  artifact = c(
    "figure_3a_error_control",
    "figure_3b_joint_success",
    "figure_3c_allocation_fragility",
    names(source_paths),
    "table_3_simulation_main",
    names(supp_paths),
    "simulation_display_qa"
  ),
  path = c(
    paste0(file.path("figures", "simulation", paste0("figure_3a_error_control_", date_tag)), ".pdf; .png"),
    paste0(file.path("figures", "simulation", paste0("figure_3b_joint_success_", date_tag)), ".pdf; .png"),
    paste0(file.path("figures", "simulation", paste0("figure_3c_allocation_fragility_", date_tag)), ".pdf; .png"),
    file.path("figures", "simulation", "source_data", basename(unname(source_paths))),
    paste0(file.path("tables", "simulation", basename(main_table_path)), "; ", file.path("tables", "simulation", basename(final_main_table_path))),
    file.path("tables", "simulation", basename(unname(supp_paths))),
    file.path("tables", "simulation", paste0("simulation_display_qa_", date_tag, ".csv"))
  ),
  manuscript_role = c(
    "Main figure candidate: false full-contribution control",
    "Main figure candidate: joint power under allocation choices",
    "Main or supplement figure candidate: allocation fragility",
    "Source data for Figure 3A",
    "Source data for Figure 3B",
    "Source data for Figure 3C",
    "Source data for main simulation table",
    "Main simulation table candidate",
    "Supplement simulation method table",
    "Supplement simulation scenario grid",
    "Supplement full error-control table",
    "Supplement full joint-power table",
    "Supplement full allocation-fragility table",
    "Supplement full robustness table",
    "Supplement validation check table",
    "Display artifact QA"
  )
)

manifest_path <- file.path(table_dir, paste0("simulation_display_manifest_", date_tag, ".csv"))
write_csv(manifest, manifest_path)

qa <- tibble(
  check = c(
    "Main validation checks passed",
    "Figure 3A source rows",
    "Figure 3B source rows",
    "Figure 3C source rows",
    "Main simulation table rows",
    "Figure files exist",
    "Main source rounded to 3 decimals",
    "Main run used as source"
  ),
  value = c(
    paste0(sum(validation$passed), "/", nrow(validation)),
    as.character(nrow(figure_3a_source)),
    as.character(nrow(figure_3b_source)),
    as.character(nrow(figure_3c_source)),
    as.character(nrow(main_table)),
    as.character(length(figure_paths)),
    "yes",
    "participant_ademp_main_2026-06-19"
  ),
  status = c(
    if_else(all(validation$passed), "pass", "fail"),
    if_else(nrow(figure_3a_source) == 15, "pass", "review"),
    if_else(nrow(figure_3b_source) == 20, "pass", "review"),
    if_else(nrow(figure_3c_source) == 6, "pass", "review"),
    if_else(nrow(main_table) == 16, "pass", "review"),
    if_else(all(file.exists(figure_paths)), "pass", "fail"),
    "pass",
    "pass"
  )
)

qa_path <- file.path("tables", "simulation", paste0("simulation_display_qa_", date_tag, ".csv"))
write_csv(qa, qa_path)

cat("Created June 19 simulation display package
")
cat("Figures: ", length(figure_paths), " files\n", sep = "")
cat("Main table: ", main_table_path, "\n", sep = "")
cat("Final main table: ", final_main_table_path, "\n", sep = "")
cat("Manifest: ", manifest_path, "\n", sep = "")
cat("QA: ", qa_path, "\n", sep = "")
print(qa)
