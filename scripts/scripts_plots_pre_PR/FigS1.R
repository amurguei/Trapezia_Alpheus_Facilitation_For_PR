## ============================================================
## Fully self-contained code: Supplementary Figure by cohort
## ============================================================

library(tidyverse)
library(ggplot2)
library(ggtext)
library(grid)

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
setwd("/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR")

input_file <- file.path("input", "controls_cleaned_with_months.csv")
out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## Colors (Okabe-Ito palette)
## ------------------------------------------------------------
taxon_cols <- c(
  trapezia = "#E69F00",
  alpheus  = "#0072B2"
)

## ------------------------------------------------------------
## Load and clean data
## ------------------------------------------------------------
dat_raw <- readr::read_csv(
  input_file,
  col_types = readr::cols(
    collector = readr::col_character(),
    .default = readr::col_guess()
  )
)

dat <- dat_raw %>%
  mutate(
    collector_clean = stringr::str_to_lower(stringr::str_squish(collector)),
    date = as.Date(date, format = "%d/%m/%Y"),
    month_bin = as.integer(floor(months_of_monitoring))
  ) %>%
  filter(!is.na(date)) %>%
  filter(!(format(date, "%Y-%m") %in% c("2009-12", "2010-03"))) %>%
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
## Per-colony summary
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
## Long format for plotting
## ------------------------------------------------------------
plot_df_by_cohort <- colony_summary_df %>%
  pivot_longer(
    cols = c(occurrence, mean_abundance, intensity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = case_match(
      metric,
      "occurrence" ~ "Occurrence",
      "mean_abundance" ~ "Mean abundance",
      "intensity" ~ "Intensity"
    ),
    metric = factor(metric, levels = c("Occurrence", "Mean abundance", "Intensity")),
    cohort_label = factor(cohort_label, levels = c("N1", "N2")),
    host_species_md = factor(
      species_md,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    ),
    taxon = factor(taxon, levels = c("trapezia", "alpheus")),
    facet_lab = factor(
      paste(cohort_label, metric, sep = " – "),
      levels = c(
        "N1 – Occurrence",
        "N1 – Mean abundance",
        "N1 – Intensity",
        "N2 – Occurrence",
        "N2 – Mean abundance",
        "N2 – Intensity"
      )
    )
  )

## ------------------------------------------------------------
## Plot
## ------------------------------------------------------------
pd <- position_dodge(width = 0.75)

fig_host_comparison_by_cohort_v_clean <- ggplot(
  plot_df_by_cohort,
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
      dodge.width = 0.75
    ),
    alpha = 0.7,
    size = 2.4
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
    ~ facet_lab,
    nrow = 2,
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
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 15) +
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.x = ggtext::element_markdown(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = ggtext::element_markdown(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.8),
    axis.line = element_line(linewidth = 0.6)
  )

print(fig_host_comparison_by_cohort_v_clean)

## ------------------------------------------------------------
## Export
## ------------------------------------------------------------
ggsave(
  filename = file.path(out_dir, "Fig_S1_clean.png"),
  plot = fig_host_comparison_by_cohort_v_clean,
  width = 11,
  height = 7.0,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_S1_clean.pdf"),
  plot = fig_host_comparison_by_cohort_v_clean,
  width = 11,
  height = 7.0,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, "Fig_S1_clean.tiff"),
  plot = fig_host_comparison_by_cohort_v_clean,
  width = 11,
  height = 7.0,
  units = "in"
)
