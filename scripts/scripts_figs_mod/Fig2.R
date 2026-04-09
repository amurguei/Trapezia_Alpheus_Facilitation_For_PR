## ============================================================
##  Fig 2 — Unified version
## Colony-level comparisons between Trapezia and Alpheus
## within each host species across the unified control series
## ============================================================

library(tidyverse)
library(ggplot2)
library(ggtext)
library(grid)
library(rstatix)
library(janitor)
library(stringr)

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
repo_dir <- "/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR_rev"
input_file <- file.path(repo_dir, "input", "controls_cleaned_with_months_unified.csv")
out_dir_fig <- file.path(repo_dir, "outputs_unified", "figures")
out_dir_tab <- file.path(repo_dir, "outputs_unified", "output_files")

dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## Colours
## ------------------------------------------------------------
taxon_cols <- c(
  trapezia = "#E69F00",  # Okabe-Ito orange
  alpheus  = "#0072B2"   # Okabe-Ito blue
)

## ------------------------------------------------------------
## 1. Read unified data
## ------------------------------------------------------------
dat_raw <- readr::read_csv(
  input_file,
  col_types = readr::cols(
    collector = readr::col_character(),
    pre_monitoring_date_flag = readr::col_logical(),
    inverts_observed = readr::col_logical(),
    .default = readr::col_guess()
  )
) %>%
  janitor::clean_names()

required_cols <- c(
  "date", "genus", "coral_uid",
  "survival_status", "months_of_monitoring",
  "treatment", "transplantation",
  "collector", "inverts_observed",
  "trapezia", "alpheus"
)

missing_cols <- setdiff(required_cols, names(dat_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

## ------------------------------------------------------------
## 2. Parse dates and exclude problematic survey months
## ------------------------------------------------------------
df_clean <- dat_raw %>%
  mutate(
    collector_clean = stringr::str_to_lower(stringr::str_squish(collector)),
    date = as.Date(date, format = "%Y-%m-%d")
  ) %>%
  filter(!is.na(date)) %>%
  filter(!(format(date, "%Y-%m") %in% c("2009-12", "2010-03")))

## ------------------------------------------------------------
## 3. Build common unified analysis dataset
## ------------------------------------------------------------
dat <- df_clean %>%
  filter(
    genus %in% c("Pocillopora", "Stylophora"),
    treatment == "control" | transplantation == "control_unified",
    survival_status == 0,
    inverts_observed == TRUE
  ) %>%
  mutate(
    trapezia = replace_na(as.numeric(trapezia), 0),
    alpheus  = replace_na(as.numeric(alpheus), 0),
    species = case_match(
      genus,
      "Pocillopora" ~ "Pocillopora favosa",
      "Stylophora"  ~ "Stylophora pistillata"
    ),
    species_md = case_match(
      genus,
      "Pocillopora" ~ "*Pocillopora favosa*",
      "Stylophora"  ~ "*Stylophora pistillata*"
    )
  )

readr::write_csv(
  dat,
  file.path(out_dir_tab, "Fig2_unified_analysis_dataset.csv")
)

## ------------------------------------------------------------
## 4. Per-colony summary metrics for each taxon
## ------------------------------------------------------------
colony_summary_df <- dat %>%
  select(species, species_md, coral_uid, trapezia, alpheus) %>%
  pivot_longer(
    cols = c(trapezia, alpheus),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  group_by(species, species_md, coral_uid, taxon) %>%
  summarise(
    occurrence = mean(count > 0, na.rm = TRUE),
    mean_abundance = mean(count, na.rm = TRUE),
    intensity = if_else(
      any(count > 0, na.rm = TRUE),
      mean(count[count > 0], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )

readr::write_csv(
  colony_summary_df,
  file.path(out_dir_tab, "Fig2_unified_colony_summary.csv")
)

## ------------------------------------------------------------
## 5. Long format table for plotting
## ------------------------------------------------------------
df_colony_long <- colony_summary_df %>%
  pivot_longer(
    cols = c(occurrence, mean_abundance, intensity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = dplyr::recode(
      metric,
      occurrence = "Occurrence",
      mean_abundance = "Mean abundance",
      intensity = "Intensity"
    ),
    host_species = factor(
      species,
      levels = c("Pocillopora favosa", "Stylophora pistillata")
    ),
    host_species_md = factor(
      species_md,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    ),
    taxon = factor(
      taxon,
      levels = c("trapezia", "alpheus")
    )
  )

readr::write_csv(
  df_colony_long,
  file.path(out_dir_tab, "Fig2_unified_df_colony_long.csv")
)

## ------------------------------------------------------------
## 6. Wilcoxon tests within each host species
## ------------------------------------------------------------
test_results_tbl <- df_colony_long %>%
  filter(!is.na(value)) %>%
  group_by(metric, host_species) %>%
  wilcox_test(value ~ taxon, detailed = TRUE) %>%
  ungroup() %>%
  select(metric, host_species, p_value = p)

readr::write_csv(
  test_results_tbl,
  file.path(out_dir_tab, "Fig2_unified_wilcoxon_tests.csv")
)

print(test_results_tbl)

## ------------------------------------------------------------
## 7. Plot data
## ------------------------------------------------------------
df_plot <- df_colony_long %>%
  mutate(
    metric = factor(
      metric,
      levels = c("Occurrence", "Mean abundance", "Intensity")
    ),
    host_species = factor(
      host_species,
      levels = c("Pocillopora favosa", "Stylophora pistillata")
    ),
    host_species_md = factor(
      host_species_md,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    ),
    taxon = factor(
      taxon,
      levels = c("trapezia", "alpheus")
    )
  )

## ------------------------------------------------------------
## 8. Significance table with one bracket per host per metric
## ------------------------------------------------------------
dodge_width <- 0.55

sig_tbl_plot <- test_results_tbl %>%
  mutate(
    metric = factor(
      metric,
      levels = c("Occurrence", "Mean abundance", "Intensity")
    ),
    host_species = factor(
      host_species,
      levels = c("Pocillopora favosa", "Stylophora pistillata")
    ),
    sig_label = case_when(
      p_value < 0.001 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  left_join(
    df_plot %>%
      group_by(metric, host_species) %>%
      summarise(
        y_max = max(value, na.rm = TRUE),
        y_min = min(value, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("metric", "host_species")
  ) %>%
  mutate(
    y_range = y_max - y_min,
    y_bracket = case_when(
      metric == "Occurrence" ~ y_max + 0.18 * if_else(y_range > 0, y_range, 1),
      TRUE ~ y_max + 0.12 * if_else(y_range > 0, y_range, 1)
    ),
    y_text = case_when(
      metric == "Occurrence" ~ y_max + 0.25 * if_else(y_range > 0, y_range, 1),
      TRUE ~ y_max + 0.18 * if_else(y_range > 0, y_range, 1)
    ),
    x_center = case_when(
      host_species == "Pocillopora favosa" ~ 1,
      host_species == "Stylophora pistillata" ~ 2
    ),
    x_left = x_center - dodge_width / 2,
    x_right = x_center + dodge_width / 2,
    x_mid = x_center
  )

## ------------------------------------------------------------
## 9. Figure
## ------------------------------------------------------------
pd <- position_dodge(width = dodge_width)

fig2_unified <- ggplot(
  df_plot,
  aes(x = host_species_md, y = value, colour = taxon)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.40,
    position = pd,
    linewidth = 0.85
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = dodge_width
    ),
    alpha = 0.3,
    size = 1.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    position = pd,
    shape = 23,
    size = 4.2,
    stroke = 1.0,
    fill = "white"
  ) +
  facet_wrap(
    ~ metric,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = taxon_cols,
    labels = c(
      "trapezia" = "*Trapezia*",
      "alpheus"  = "*Alpheus*"
    ),
    name = NULL
  ) +
  geom_segment(
    data = sig_tbl_plot,
    aes(x = x_left, xend = x_right, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.65,
    colour = "black"
  ) +
  geom_segment(
    data = sig_tbl_plot,
    aes(
      x = x_left, xend = x_left,
      y = y_bracket - 0.03 * if_else(y_range > 0, y_range, 1),
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.65,
    colour = "black"
  ) +
  geom_segment(
    data = sig_tbl_plot,
    aes(
      x = x_right, xend = x_right,
      y = y_bracket - 0.03 * if_else(y_range > 0, y_range, 1),
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.65,
    colour = "black"
  ) +
  geom_text(
    data = sig_tbl_plot,
    aes(x = x_mid, y = y_text, label = sig_label),
    inherit.aes = FALSE,
    size = 5.0,
    fontface = "bold",
    colour = "black"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 15) +
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.x = ggtext::element_markdown(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12.5),
    legend.text = ggtext::element_markdown(size = 12),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.8),
    axis.line = element_line(linewidth = 0.6),
    panel.spacing.x = unit(1.2, "lines")
  )

print(fig2_unified)

## ------------------------------------------------------------
## 10. Save figure
## ------------------------------------------------------------
ggsave(
  filename = file.path(out_dir_fig, "Fig2_unified_host_comparison_boxplots.tiff"),
  plot = fig2_unified,
  width = 290,
  height = 110,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig2_unified_host_comparison_boxplots.png"),
  plot = fig2_unified,
  width = 280,
  height = 110,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "Fig2_unified_host_comparison_boxplots.pdf"),
  plot = fig2_unified,
  width = 280,
  height = 110,
  units = "mm"
)

## ------------------------------------------------------------
## 11. Print p-values cleanly
## ------------------------------------------------------------
test_results_tbl %>%
  mutate(p_value = signif(p_value, 3)) %>%
  print()