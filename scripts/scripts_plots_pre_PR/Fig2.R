## ============================================================
##  Fig 2
## ============================================================

library(tidyverse)
library(ggplot2)
library(ggtext)
library(grid)
library(rstatix)

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
setwd("/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR")

input_file <- file.path("input", "controls_cleaned_with_months.csv")
out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## Colours
## ------------------------------------------------------------
taxon_cols <- c(
  trapezia = "#E69F00",  # Okabe-Ito orange
  alpheus  = "#0072B2"   # Okabe-Ito blue
)

## ------------------------------------------------------------
## 1. Read raw data
## ------------------------------------------------------------
dat_raw <- readr::read_csv(
  input_file,
  col_types = readr::cols(
    collector = readr::col_character(),
    .default = readr::col_guess()
  )
)

required_cols <- c(
  "date", "cohort", "genus", "coral_uid",
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
    date = as.Date(date, format = "%d/%m/%Y")
  ) %>%
  filter(!is.na(date)) %>%
  filter(!(format(date, "%Y-%m") %in% c("2009-12", "2010-03")))

## ------------------------------------------------------------
## 3. Build common analysis dataset
## ------------------------------------------------------------
dat <- df_clean %>%
  mutate(
    month_bin = as.integer(floor(months_of_monitoring))
  ) %>%
  filter(
    (treatment == "control") |
      (transplantation %in% c("control_T1", "control_T2")),
    survival_status == 0,
    inverts_observed == TRUE,
    genus %in% c("Pocillopora", "Stylophora")
  ) %>%
  mutate(
    trapezia = replace_na(trapezia, 0),
    alpheus  = replace_na(alpheus, 0),
    species = case_match(
      genus,
      "Pocillopora" ~ "Pocillopora favosa",
      "Stylophora"  ~ "Stylophora pistillata"
    ),
    species_md = case_match(
      genus,
      "Pocillopora" ~ "*Pocillopora favosa*",
      "Stylophora"  ~ "*Stylophora pistillata*"
    ),
    cohort_label = case_match(
      cohort,
      "T1" ~ "N1",
      "T2" ~ "N2",
      .default = as.character(cohort)
    )
  )

## ------------------------------------------------------------
## 4. Per-colony summary metrics for each taxon
## ------------------------------------------------------------
colony_summary_df <- dat %>%
  select(cohort_label, species, species_md, coral_uid, trapezia, alpheus) %>%
  pivot_longer(
    cols = c(trapezia, alpheus),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  group_by(cohort_label, species, species_md, coral_uid, taxon) %>%
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
    ),
    cohort_label = factor(
      cohort_label,
      levels = c("N1", "N2")
    )
  )

## ------------------------------------------------------------
## 6. Wilcoxon tests pooled across cohorts, within host species
## ------------------------------------------------------------
test_results_tbl <- df_colony_long %>%
  filter(!is.na(value)) %>%
  group_by(metric, host_species) %>%
  wilcox_test(value ~ taxon, detailed = TRUE) %>%
  ungroup() %>%
  select(metric, host_species, p_value = p)

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

##hide ns
# sig_tbl_plot <- sig_tbl_plot %>% filter(sig_label != "ns")

## ------------------------------------------------------------
## 9. One-row compact figure
## ------------------------------------------------------------
pd <- position_dodge(width = dodge_width)

fig_host_comparison_stats_v_clean <- ggplot(
  df_plot,
  aes(x = host_species_md, y = value, colour = taxon)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.40,              ## wider boxplots 
    position = pd,
    linewidth = 0.85
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = dodge_width
    ),
    alpha = 0.3,               ## more transparent (was 0.7)
    size = 1.5                 ## lightly smaller (was 2.4)
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

print(fig_host_comparison_stats_v_clean)

## ------------------------------------------------------------
## 10. Save figure
## ------------------------------------------------------------

ggsave(
  filename = file.path(out_dir, "Fig2_host_comparison_boxplots.tiff"),
  plot = fig_host_comparison_stats_v_clean,
  width = 290,     # mm → full-width journal figure
  height = 110,    # nice 1-row proportion
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir, "Fig2_host_comparison_boxplots.png"),
  plot = fig_host_comparison_stats_v_clean,
  width = 290,
  height = 110,
  units = "mm"
)


test_results_tbl %>%
  mutate(p_value = signif(p_value, 3))

test_results_tbl

