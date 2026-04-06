## ============================================================
## Trapezia–Alpheus temporal dynamics in control colonies
##
## Figures:
## 1) incidence  = proportion of colony-months occupied
## 2) abundance  = mean count per colony-month
## 3) intensity  = mean count among occupied colony-months only
##
## Panels:
## - Pocillopora favosa - N1
## - Pocillopora favosa - N2
## - Stylophora pistillata - N1
## - Stylophora pistillata - N2
## ============================================================

library(tidyverse)
library(ggplot2)
library(ggtext)
library(scales)
library(janitor)
library(stringr)

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
setwd("/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR")

input_file <- file.path("input", "controls_cleaned_with_months.csv")
out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## 1. Import raw data from scratch
## ------------------------------------------------------------
dat_raw <- readr::read_csv(
  input_file,
  col_types = readr::cols(
    collector = readr::col_character(),
    pre_monitoring_date_flag = readr::col_logical(),
    fish_observed = readr::col_logical(),
    inverts_observed = readr::col_logical()
  )
) %>%
  janitor::clean_names()

## ------------------------------------------------------------
## 2. Parse date and remove problematic months
## ------------------------------------------------------------
df_clean <- dat_raw %>%
  mutate(
    collector_clean = stringr::str_to_lower(stringr::str_squish(collector)),
    date = as.Date(date, format = "%d/%m/%Y")
  ) %>%
  filter(!is.na(date)) %>%
  filter(!(format(date, "%Y-%m") %in% c("2009-12", "2010-03")))

message("Retained survey months after exclusion:")
print(
  df_clean %>%
    count(format(date, "%Y-%m"), name = "n_rows") %>%
    arrange(`format(date, "%Y-%m")`)
)

## ------------------------------------------------------------
## 3. Build cleaned analysis dataset
## live colonies only, control colonies only, months with
## invertebrate observations only
## ------------------------------------------------------------
focal_hosts <- c("Pocillopora", "Stylophora")

dat <- df_clean %>%
  mutate(
    month_bin = as.integer(floor(months_of_monitoring))
  ) %>%
  filter(
    (treatment == "control") |
      (transplantation %in% c("control_t1", "control_t2", "control_T1", "control_T2")),
    survival_status == 0,
    inverts_observed == TRUE,
    genus %in% focal_hosts
  ) %>%
  mutate(
    cohort = factor(cohort, levels = c("T1", "T2")),
    genus = factor(genus, levels = c("Pocillopora", "Stylophora")),
    coral_uid = factor(coral_uid),
    trapezia = replace_na(as.numeric(trapezia), 0),
    alpheus  = replace_na(as.numeric(alpheus), 0)
  )

readr::write_csv(
  dat,
  file.path(out_dir, "trapezia_alpheus_temporal_dat_v_clean.csv")
)

## ------------------------------------------------------------
## 4. Build long-format dataset
## one row = colony x month x taxon
## ------------------------------------------------------------
dat_long <- dat %>%
  dplyr::select(cohort, genus, coral_uid, month_bin, trapezia, alpheus) %>%
  tidyr::pivot_longer(
    cols = c(trapezia, alpheus),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  dplyr::mutate(
    present = if_else(count > 0, 1L, 0L),
    species = dplyr::case_match(
      as.character(genus),
      "Pocillopora" ~ "*Pocillopora favosa*",
      "Stylophora"  ~ "*Stylophora pistillata*"
    ),
    cohort_label = dplyr::case_match(
      as.character(cohort),
      "T1" ~ "N1",
      "T2" ~ "N2"
    ),
    host_cohort = factor(
      paste(species, cohort_label, sep = " - "),
      levels = c(
        "*Pocillopora favosa* - N1",
        "*Pocillopora favosa* - N2",
        "*Stylophora pistillata* - N1",
        "*Stylophora pistillata* - N2"
      )
    ),
    taxon_label = dplyr::case_match(
      taxon,
      "trapezia" ~ "*Trapezia*",
      "alpheus"  ~ "*Alpheus*"
    )
  )

readr::write_csv(
  dat_long,
  file.path(out_dir, "trapezia_alpheus_temporal_dat_long_v_clean.csv")
)

## ------------------------------------------------------------
## 5. Summaries through time
## ------------------------------------------------------------
time_summary_df <- dat_long %>%
  group_by(host_cohort, species, cohort_label, genus, cohort, month_bin, taxon, taxon_label) %>%
  summarise(
    n_colony_months = n(),
    n_colonies = n_distinct(coral_uid),
    
    ## abundance: all colony-months
    mean_abundance = mean(count, na.rm = TRUE),
    sd_abundance   = sd(count, na.rm = TRUE),
    se_abundance   = sd_abundance / sqrt(n_colony_months),
    
    ## incidence: occupancy frequency
    incidence = mean(present, na.rm = TRUE),
    se_incidence = sqrt((incidence * (1 - incidence)) / n_colony_months),
    
    ## intensity: only occupied colony-months
    n_occupied = sum(present, na.rm = TRUE),
    mean_intensity = if_else(
      n_occupied > 0,
      mean(count[count > 0], na.rm = TRUE),
      NA_real_
    ),
    sd_intensity = if_else(
      n_occupied > 1,
      sd(count[count > 0], na.rm = TRUE),
      NA_real_
    ),
    se_intensity = if_else(
      n_occupied > 1,
      sd_intensity / sqrt(n_occupied),
      NA_real_
    ),
    .groups = "drop"
  )

readr::write_csv(
  time_summary_df,
  file.path(out_dir, "trapezia_alpheus_temporal_summary_v_clean.csv")
)

## ------------------------------------------------------------
## 6. Okabe-Ito colours
## ------------------------------------------------------------
taxon_cols <- c(
  trapezia = "#E69F00",
  alpheus  = "#0072B2"
)

## ------------------------------------------------------------
## 7. Generic plotting function
## larger text, self-contained styling
## ------------------------------------------------------------
plot_temporal_metric <- function(data, y_var, se_var, y_lab, percent_axis = FALSE) {
  
  plot_df <- data %>%
    mutate(
      y = .data[[y_var]],
      se = .data[[se_var]],
      ymin = pmax(y - 1.96 * se, 0),
      ymax = y + 1.96 * se
    ) %>%
    filter(!is.na(y))
  
  p <- ggplot(
    plot_df,
    aes(
      x = month_bin,
      y = y,
      colour = taxon,
      group = taxon
    )
  ) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.8) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0.30,
      linewidth = 0.60
    ) +
    facet_wrap(
      ~ host_cohort,
      nrow = 1,
      scales = "free_x"
    ) +
    scale_colour_manual(
      values = taxon_cols,
      breaks = c("trapezia", "alpheus"),
      labels = c(
        "trapezia" = "*Trapezia*",
        "alpheus"  = "*Alpheus*"
      ),
      name = NULL
    ) +
    labs(
      x = "Monitoring month",
      y = y_lab
    ) +
    theme_bw(base_size = 18) +
    theme(
      strip.text = ggtext::element_markdown(size = 15, face = "plain"),
      legend.text = ggtext::element_markdown(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = ggtext::element_markdown(size = 16),
      axis.text = element_text(size = 13),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  if (percent_axis) {
    p <- p +
      scale_y_continuous(
        labels = scales::label_percent(accuracy = 1),
        expand = expansion(mult = c(0.02, 0.10))
      )
  } else {
    p <- p +
      scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.10))
      )
  }
  
  p
}

## ------------------------------------------------------------
## 8. Build figures
## ------------------------------------------------------------
fig_time_incidence <- plot_temporal_metric(
  data = time_summary_df,
  y_var = "incidence",
  se_var = "se_incidence",
  y_lab = "Proportion of colony-months occupied",
  percent_axis = TRUE
)

fig_time_abundance <- plot_temporal_metric(
  data = time_summary_df,
  y_var = "mean_abundance",
  se_var = "se_abundance",
  y_lab = "Mean individuals per colony-month",
  percent_axis = FALSE
)

fig_time_intensity <- plot_temporal_metric(
  data = time_summary_df,
  y_var = "mean_intensity",
  se_var = "se_intensity",
  y_lab = "Mean individuals per occupied colony",
  percent_axis = FALSE
)

## ------------------------------------------------------------
## 9. Print figures
## ------------------------------------------------------------
print(fig_time_incidence)
print(fig_time_abundance)
print(fig_time_intensity)

## ------------------------------------------------------------
## 10. Save figures
## ------------------------------------------------------------
ggsave(
  filename = file.path(out_dir, "Fig_temporal_incidence_trapezia_alpheus_v_clean.png"),
  plot = fig_time_incidence,
  width = 11,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_temporal_incidence_trapezia_alpheus_v_clean.pdf"),
  plot = fig_time_incidence,
  width = 10.8,
  height = 4.2
)

ggsave(
  filename = file.path(out_dir, "Fig_temporal_abundance_trapezia_alpheus_v_clean.png"),
  plot = fig_time_abundance,
  width = 10.8,
  height = 4.2,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_temporal_abundance_trapezia_alpheus_v_clean.pdf"),
  plot = fig_time_abundance,
  width = 10.8,
  height = 4.2
)

ggsave(
  filename = file.path(out_dir, "Fig_S2.png"),
  plot = fig_time_intensity,
  width = 11,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_S2.pdf"),
  plot = fig_time_intensity,
  width = 11,
  height = 5
)

ggsave(
  filename = file.path(out_dir, "Fig_S2.tiff"),
  plot = fig_time_intensity,
  width = 11,
  height = 5
)
