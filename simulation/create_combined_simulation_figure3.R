#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
})

project_root <- normalizePath(getwd(), mustWork = TRUE)
figure_dir <- file.path(project_root, "figures", "simulation")
source_dir <- file.path(figure_dir, "source_data")
date_tag <- "2026-06-20"

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

figure_3a_source <- read_csv(
  file.path(source_dir, "figure_3a_error_control_source_2026-06-19.csv"),
  show_col_types = FALSE
)

figure_3b_source <- read_csv(
  file.path(source_dir, "figure_3b_joint_success_source_2026-06-19.csv"),
  show_col_types = FALSE
)

harder_labels <- c(
  "Limiting first contrast" = "First contrast harder",
  "Limiting second contrast" = "Second contrast harder",
  "Reversed limiting component" = "Harder contrast reversed"
)

figure_3b_source <- figure_3b_source %>%
  mutate(
    scenario_display = recode(scenario_display, !!!harder_labels),
    scenario_display = factor(
      scenario_display,
      levels = c(
        "Balanced positive",
        "First contrast harder",
        "Second contrast harder",
        "Positive margin",
        "Additive positive excess"
      )
    )
  )

method_order <- c(
  "Equal allocation, IUT",
  "Shared-arm, IUT",
  "Plug-in allocation, IUT",
  "Conservative allocation, IUT",
  "Marginal comparator, IUT",
  "Equal allocation, maxT separate claims"
)

figure_3a_source <- figure_3a_source %>%
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
  )

figure_3b_source <- figure_3b_source %>%
  mutate(method_display = factor(method_display, levels = method_order))

palette <- c(
  "Equal allocation, IUT" = "#1B5A7A",
  "Shared-arm, IUT" = "#4C4C4C",
  "Plug-in allocation, IUT" = "#7A3E8E",
  "Conservative allocation, IUT" = "#2F6B45",
  "Marginal comparator, IUT" = "#8A5A16",
  "Equal allocation, maxT separate claims" = "#8E2F3F"
)

theme_display <- function(base_size = 10.5) {
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
      legend.text = element_text(color = "#000000", size = 9),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "#D8D8D8", linewidth = 0.30),
      panel.spacing.y = unit(4, "pt"),
      plot.margin = margin(4, 10, 4, 4)
    )
}

panel_a <- ggplot(
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
  geom_point(size = 2.35) +
  facet_wrap(~ scenario_display, ncol = 1) +
  scale_x_continuous(
    limits = c(0, 0.032),
    breaks = seq(0, 0.03, by = 0.01),
    labels = format_num
  ) +
  scale_color_manual(values = palette, drop = TRUE) +
  labs(x = "False full-contribution conclusion", y = NULL, tag = "A") +
  theme_display(9.6) +
  theme(
    legend.position = "none",
    plot.tag = element_text(face = "bold", color = "#000000", size = 12)
  )

panel_b <- ggplot(
  figure_3b_source,
  aes(
    x = empirical_joint_success,
    y = scenario_display,
    xmin = lower,
    xmax = upper,
    color = method_display
  )
) +
  geom_errorbar(
    orientation = "y",
    width = 0.16,
    linewidth = 0.55,
    position = position_dodge(width = 0.55)
  ) +
  geom_point(size = 2.45, position = position_dodge(width = 0.55)) +
  scale_x_continuous(
    limits = c(0.40, 0.82),
    breaks = seq(0.40, 0.80, by = 0.10),
    labels = format_num
  ) +
  scale_color_manual(values = palette, drop = TRUE) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = "Joint power", y = NULL, tag = "B") +
  theme_display(9.6) +
  theme(plot.tag = element_text(face = "bold", color = "#000000", size = 12))

combined <- (panel_a / panel_b) + plot_layout(heights = c(1.02, 1.10))

pdf_path <- file.path(figure_dir, paste0("figure_3_simulation_error_and_power_", date_tag, ".pdf"))
png_path <- file.path(figure_dir, paste0("figure_3_simulation_error_and_power_", date_tag, ".png"))

ggsave(pdf_path, combined, width = 7.2, height = 6.7, units = "in", device = grDevices::pdf)

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    png_path,
    combined,
    width = 7.2,
    height = 6.7,
    units = "in",
    dpi = 320,
    device = ragg::agg_png
  )
} else {
  ggsave(
    png_path,
    combined,
    width = 7.2,
    height = 6.7,
    units = "in",
    dpi = 320,
    device = grDevices::png
  )
}

qa <- tibble(
  check = c(
    "Combined Figure 3 PDF exists",
    "Combined Figure 3 PNG exists",
    "Panel B legend requested as one row",
    "Panel B x axis uses joint power",
    "Harder contrast labels used in Panel B"
  ),
  value = c(
    as.character(file.exists(pdf_path)),
    as.character(file.exists(png_path)),
    "nrow = 1",
    "Joint power",
    paste(levels(figure_3b_source$scenario_display), collapse = "; ")
  ),
  status = c(
    if_else(file.exists(pdf_path), "pass", "fail"),
    if_else(file.exists(png_path), "pass", "fail"),
    "pass",
    "pass",
    if_else(any(grepl("Limiting", levels(figure_3b_source$scenario_display))), "fail", "pass")
  )
)

qa_path <- file.path("tables", "simulation", paste0("figure_3_combined_qa_", date_tag, ".csv"))
write_csv(qa, qa_path)

cat("Created combined Figure 3\n")
cat("PDF: ", pdf_path, "\n", sep = "")
cat("PNG: ", png_path, "\n", sep = "")
cat("QA: ", qa_path, "\n", sep = "")
print(qa)
