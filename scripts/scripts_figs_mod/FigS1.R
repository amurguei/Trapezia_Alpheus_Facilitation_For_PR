library(tidyverse)
library(ggplot2)
library(ggtext)
library(scales)
library(janitor)
library(stringr)

repo_dir <- "/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR_rev"
input_file <- file.path(repo_dir, "input", "controls_cleaned_with_months_unified.csv")
out_dir_fig <- file.path(repo_dir, "outputs_unified", "figures")
out_dir_tab <- file.path(repo_dir, "outputs_unified", "output_files")

dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

focal_hosts <- c("Pocillopora", "Stylophora")

taxon_cols <- c(
  trapezia = "#E69F00",
  alpheus  = "#0072B2"
)

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

df_clean <- dat_raw %>%
  mutate(
    collector_clean = stringr::str_to_lower(stringr::str_squish(collector)),
    date = as.Date(date, format = "%Y-%m-%d")
  ) %>%
  filter(!is.na(date)) %>%
  filter(!(format(date, "%Y-%m") %in% c("2009-12", "2010-03")))

message("Retained survey months after exclusion:")
print(
  df_clean %>%
    count(format(date, "%Y-%m"), name = "n_rows") %>%
    arrange(`format(date, "%Y-%m")`)
)

dat <- df_clean %>%
  mutate(
    month_bin = as.integer(floor(months_of_monitoring))
  ) %>%
  filter(
    genus %in% focal_hosts,
    treatment == "control" | transplantation == "control_unified",
    survival_status == 0,
    inverts_observed == TRUE
  ) %>%
  mutate(
    genus = factor(genus, levels = c("Pocillopora", "Stylophora")),
    coral_uid = factor(coral_uid),
    trapezia = replace_na(as.numeric(trapezia), 0),
    alpheus  = replace_na(as.numeric(alpheus), 0),
    host_species = dplyr::case_match(
      as.character(genus),
      "Pocillopora" ~ "Pocillopora favosa",
      "Stylophora"  ~ "Stylophora pistillata"
    ),
    host_species_md = dplyr::case_match(
      as.character(genus),
      "Pocillopora" ~ "*Pocillopora favosa*",
      "Stylophora"  ~ "*Stylophora pistillata*"
    )
  )

readr::write_csv(
  dat,
  file.path(out_dir_tab, "FigS2_unified_analysis_dataset.csv")
)

dat_long <- dat %>%
  dplyr::select(host_species, host_species_md, genus, coral_uid, month_bin, trapezia, alpheus) %>%
  tidyr::pivot_longer(
    cols = c(trapezia, alpheus),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  dplyr::mutate(
    present = if_else(count > 0, 1L, 0L),
    taxon_label = dplyr::case_match(
      taxon,
      "trapezia" ~ "*Trapezia*",
      "alpheus"  ~ "*Alpheus*"
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

readr::write_csv(
  dat_long,
  file.path(out_dir_tab, "FigS2_unified_dat_long.csv")
)

time_summary_df <- dat_long %>%
  group_by(host_species, host_species_md, genus, month_bin, taxon, taxon_label) %>%
  summarise(
    n_colony_observations = n(),
    n_colonies = n_distinct(coral_uid),
    
    mean_abundance = mean(count, na.rm = TRUE),
    sd_abundance   = sd(count, na.rm = TRUE),
    se_abundance   = sd_abundance / sqrt(n_colony_observations),
    
    incidence = mean(present, na.rm = TRUE),
    se_incidence = sqrt((incidence * (1 - incidence)) / n_colony_observations),
    
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
  file.path(out_dir_tab, "FigS2_unified_temporal_summary.csv")
)

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
    geom_line(linewidth = 1.05) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0.25,
      linewidth = 0.55
    ) +
    facet_wrap(
      ~ host_species_md,
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
    scale_x_continuous(
      breaks = seq(
        0,
        max(plot_df$month_bin, na.rm = TRUE),
        by = 10
      )
    ) +
    labs(
      x = "Monitoring month",
      y = y_lab
    ) +
    theme_bw(base_size = 16) +
    theme(
      strip.text = ggtext::element_markdown(size = 15, face = "bold"),
      legend.text = ggtext::element_markdown(size = 15),
      axis.title.x = element_text(size = 14.5),
      axis.title.y = ggtext::element_markdown(size = 14.5),
      axis.text = element_text(size = 12),
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 0.8),
      axis.line = element_line(linewidth = 0.6)
    )
  
  if (percent_axis) {
    p <- p +
      scale_y_continuous(
        labels = scales::label_percent(accuracy = 1),
        expand = expansion(mult = c(0.02, 0.08))
      )
  } else {
    p <- p +
      scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.08))
      )
  }
  
  p
}

fig_time_incidence <- plot_temporal_metric(
  data = time_summary_df,
  y_var = "incidence",
  se_var = "se_incidence",
  y_lab = "Occupancy (%)",
  percent_axis = TRUE
)

fig_time_abundance <- plot_temporal_metric(
  data = time_summary_df,
  y_var = "mean_abundance",
  se_var = "se_abundance",
  y_lab = "Mean abundance",
  percent_axis = FALSE
)

fig_time_intensity <- plot_temporal_metric(
  data = time_summary_df,
  y_var = "mean_intensity",
  se_var = "se_intensity",
  y_lab = "Intensity",
  percent_axis = FALSE
)

print(fig_time_incidence)
print(fig_time_abundance)
print(fig_time_intensity)

ggsave(
  filename = file.path(out_dir_fig, "FigS2_unified_temporal_occupancy.png"),
  plot = fig_time_incidence,
  width = 10.4,
  height = 4.6,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "FigS2_unified_temporal_occupancy.pdf"),
  plot = fig_time_incidence,
  width = 10.4,
  height = 4.6,
  units = "in"
)

ggsave(
  filename = file.path(out_dir_fig, "FigS2_unified_temporal_occupancy.tiff"),
  plot = fig_time_incidence,
  width = 10.4,
  height = 4.6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir_fig, "FigS2_unified_temporal_abundance.png"),
  plot = fig_time_abundance,
  width = 10.4,
  height = 4.6,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "FigS2_unified_temporal_abundance.pdf"),
  plot = fig_time_abundance,
  width = 10.4,
  height = 4.6,
  units = "in"
)

ggsave(
  filename = file.path(out_dir_fig, "FigS2_unified_temporal_abundance.tiff"),
  plot = fig_time_abundance,
  width = 10.4,
  height = 4.6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir_fig, "FigS2_unified_temporal_intensity.png"),
  plot = fig_time_intensity,
  width = 10.4,
  height = 4.6,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "FigS2_unified_temporal_intensity.pdf"),
  plot = fig_time_intensity,
  width = 10.4,
  height = 4.6,
  units = "in"
)

ggsave(
  filename = file.path(out_dir_fig, "FigS2_unified_temporal_intensity.tiff"),
  plot = fig_time_intensity,
  width = 10.4,
  height = 4.6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
