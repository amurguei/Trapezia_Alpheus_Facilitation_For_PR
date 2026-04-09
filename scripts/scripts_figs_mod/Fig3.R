## ============================================================
## Fig 3 — Unified composite figure
##
## Panels:
## a) Conditional probability of Alpheus presence
##    with and without Trapezia
## b) Occupancy through time
## c) Mean abundance through time
##
## Main version uses unified months_of_monitoring
## Backup version uses survey_occasion
## ============================================================

library(tidyverse)
library(ggplot2)
library(ggtext)
library(scales)
library(patchwork)
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
col_trapezia <- "#E69F00"
col_alpheus  <- "#0072B2"
col_probbase <- "#56B4E9"

taxon_cols <- c(
  trapezia = col_trapezia,
  alpheus  = col_alpheus
)

## ------------------------------------------------------------
## Shared theme
## ------------------------------------------------------------
base_theme <- theme_bw(base_size = 16) +
  theme(
    strip.text = ggtext::element_markdown(size = 15, face = "bold"),
    axis.text = element_text(size = 15, colour = "black"),
    axis.text.x = ggtext::element_markdown(size = 15, colour = "black"),
    axis.title.x = element_text(size = 15, colour = "black"),
    axis.title.y = ggtext::element_markdown(size = 15, colour = "black"),
    legend.text = ggtext::element_markdown(size = 15, colour = "black"),
    legend.title = element_blank(),
    plot.title = element_text(size = 17, face = "bold", hjust = 0, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.8),
    axis.line = element_line(linewidth = 0.6),
    legend.position = "bottom"
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
  "survey_occasion",
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
    collector_clean = str_to_lower(str_squish(collector)),
    date = as.Date(date, format = "%Y-%m-%d")
  ) %>%
  filter(!is.na(date)) %>%
  filter(!(format(date, "%Y-%m") %in% c("2009-12", "2010-03")))

message("Date range after exclusion:")
print(range(df_clean$date, na.rm = TRUE))

message("Survey months retained:")
print(
  df_clean %>%
    count(format(date, "%Y-%m"), name = "n_rows") %>%
    arrange(`format(date, "%Y-%m")`)
)

## Optional sanity check on monitoring months
message("Monitoring month range in unified analysis:")

print(range(df_clean$months_of_monitoring, na.rm = TRUE))

message("Number of distinct monitoring months in unified analysis:")
print(n_distinct(df_clean$months_of_monitoring))

message("Monitoring months represented in unified analysis:")
print(sort(unique(df_clean$months_of_monitoring)))

## ------------------------------------------------------------
## 3. Build unified analysis dataset
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
    T_pres = if_else(trapezia > 0, 1L, 0L),
    A_pres = if_else(alpheus  > 0, 1L, 0L),
    months_of_monitoring = as.numeric(months_of_monitoring),
    survey_occasion = as.numeric(survey_occasion),
    host_species = case_match(
      genus,
      "Pocillopora" ~ "Pocillopora favosa",
      "Stylophora"  ~ "Stylophora pistillata"
    ),
    host_species_md = case_match(
      genus,
      "Pocillopora" ~ "*Pocillopora favosa*",
      "Stylophora"  ~ "*Stylophora pistillata*"
    )
  )

readr::write_csv(
  dat,
  file.path(out_dir_tab, "Fig3_unified_analysis_dataset.csv")
)

## ------------------------------------------------------------
## 4. Helper functions
## ------------------------------------------------------------
log_or_ha <- function(A, T) {
  tab <- table(
    factor(A, levels = c(0, 1)),
    factor(T, levels = c(0, 1))
  )
  
  a <- tab[2, 2] + 0.5
  b <- tab[2, 1] + 0.5
  c <- tab[1, 2] + 0.5
  d <- tab[1, 1] + 0.5
  
  log((a * d) / (b * c))
}

or_ha <- function(A, T) {
  tab <- table(
    factor(A, levels = c(0, 1)),
    factor(T, levels = c(0, 1))
  )
  
  a <- tab[2, 2] + 0.5
  b <- tab[2, 1] + 0.5
  c <- tab[1, 2] + 0.5
  d <- tab[1, 1] + 0.5
  
  (a * d) / (b * c)
}

perm_test_group <- function(d, nperm = 10000, seed = 123) {
  set.seed(seed)
  
  stat_obs <- log_or_ha(d$A_pres, d$T_pres)
  
  stat_perm <- replicate(nperm, {
    d_perm <- d %>%
      group_by(coral_uid) %>%
      mutate(A_perm = sample(A_pres, size = n(), replace = FALSE)) %>%
      ungroup()
    
    log_or_ha(d_perm$A_perm, d_perm$T_pres)
  })
  
  tibble(
    logOR_obs = stat_obs,
    p_perm = mean(abs(stat_perm) >= abs(stat_obs), na.rm = TRUE),
    n = nrow(d),
    n_colonies = dplyr::n_distinct(d$coral_uid)
  )
}

plot_temporal_metric <- function(data, x_var, y_var, se_var, x_lab, y_lab,
                                 percent_axis = FALSE, show_legend = TRUE) {
  
  plot_df <- data %>%
    mutate(
      x = .data[[x_var]],
      y = .data[[y_var]],
      se = .data[[se_var]],
      ymin = pmax(y - 1.96 * se, 0),
      ymax = y + 1.96 * se
    ) %>%
    filter(!is.na(y), !is.na(x))
  
  p <- ggplot(
    plot_df,
    aes(x = x, y = y, colour = taxon, group = taxon)
  ) +
    geom_line(linewidth = 1.05) +
    geom_point(size = 2.4) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0.25,
      linewidth = 0.55
    ) +
    facet_wrap(~ host_species_md, nrow = 1, scales = "free_x") +
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
      x = x_lab,
      y = y_lab
    ) +
    base_theme +
    theme(
      legend.position = if (show_legend) "bottom" else "none"
    )
  
  if (x_var == "months_of_monitoring") {
    p <- p +
      scale_x_continuous(
        breaks = seq(0, max(plot_df$x, na.rm = TRUE), by = 10)
      )
  }
  
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

## ------------------------------------------------------------
## 5. Temporal summaries
## ------------------------------------------------------------
dat_long <- dat %>%
  select(
    date, months_of_monitoring, survey_occasion,
    genus, host_species, host_species_md,
    coral_uid, trapezia, alpheus
  ) %>%
  pivot_longer(
    cols = c(trapezia, alpheus),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  mutate(
    present = if_else(count > 0, 1L, 0L),
    taxon_label = case_match(
      taxon,
      "trapezia" ~ "*Trapezia*",
      "alpheus"  ~ "*Alpheus*"
    )
  )

time_summary_month <- dat_long %>%
  group_by(
    host_species, host_species_md, genus,
    months_of_monitoring, taxon, taxon_label
  ) %>%
  summarise(
    n_observations = n(),
    n_colonies = n_distinct(coral_uid),
    
    mean_abundance = mean(count, na.rm = TRUE),
    sd_abundance   = sd(count, na.rm = TRUE),
    se_abundance   = sd_abundance / sqrt(n_observations),
    
    incidence = mean(present, na.rm = TRUE),
    se_incidence = sqrt((incidence * (1 - incidence)) / n_observations),
    
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

time_summary_occ <- dat_long %>%
  group_by(
    host_species, host_species_md, genus,
    survey_occasion, taxon, taxon_label
  ) %>%
  summarise(
    n_observations = n(),
    n_colonies = n_distinct(coral_uid),
    
    mean_abundance = mean(count, na.rm = TRUE),
    sd_abundance   = sd(count, na.rm = TRUE),
    se_abundance   = sd_abundance / sqrt(n_observations),
    
    incidence = mean(present, na.rm = TRUE),
    se_incidence = sqrt((incidence * (1 - incidence)) / n_observations),
    
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
  time_summary_month,
  file.path(out_dir_tab, "Fig3_unified_temporal_summary_by_month.csv")
)

readr::write_csv(
  time_summary_occ,
  file.path(out_dir_tab, "Fig3_unified_temporal_summary_by_occasion.csv")
)

## ------------------------------------------------------------
## 6. Temporal panels
## ------------------------------------------------------------
fig_time_incidence_month <- plot_temporal_metric(
  data = time_summary_month,
  x_var = "months_of_monitoring",
  y_var = "incidence",
  se_var = "se_incidence",
  x_lab = "Monitoring month",
  y_lab = "Occupancy (%)",
  percent_axis = TRUE,
  show_legend = TRUE
)

fig_time_abundance_month <- plot_temporal_metric(
  data = time_summary_month,
  x_var = "months_of_monitoring",
  y_var = "mean_abundance",
  se_var = "se_abundance",
  x_lab = "Monitoring month",
  y_lab = "Mean abundance",
  percent_axis = FALSE,
  show_legend = TRUE
)

fig_time_incidence_occ <- plot_temporal_metric(
  data = time_summary_occ,
  x_var = "survey_occasion",
  y_var = "incidence",
  se_var = "se_incidence",
  x_lab = "Survey occasion",
  y_lab = "Occupancy (%)",
  percent_axis = TRUE,
  show_legend = TRUE
)

fig_time_abundance_occ <- plot_temporal_metric(
  data = time_summary_occ,
  x_var = "survey_occasion",
  y_var = "mean_abundance",
  se_var = "se_abundance",
  x_lab = "Survey occasion",
  y_lab = "Mean abundance",
  percent_axis = FALSE,
  show_legend = TRUE
)

## ------------------------------------------------------------
## 7. Co-occurrence statistics by host species
## ------------------------------------------------------------
desc_by_group <- dat %>%
  group_by(genus, host_species, host_species_md) %>%
  summarise(
    n = n(),
    n_colonies = n_distinct(coral_uid),
    n_T1 = sum(T_pres == 1),
    n_T0 = sum(T_pres == 0),
    pA_given_T = if_else(n_T1 > 0, mean(A_pres[T_pres == 1]), NA_real_),
    pA_given_notT = if_else(n_T0 > 0, mean(A_pres[T_pres == 0]), NA_real_),
    RR = if_else(
      !is.na(pA_given_T) & !is.na(pA_given_notT) & pA_given_notT > 0,
      pA_given_T / pA_given_notT,
      NA_real_
    ),
    OR = or_ha(A_pres, T_pres),
    logOR = log_or_ha(A_pres, T_pres),
    .groups = "drop"
  )

n_perm <- 10000

perm_results <- dat %>%
  group_by(genus, host_species, host_species_md) %>%
  group_modify(~ perm_test_group(.x, nperm = n_perm, seed = 123)) %>%
  ungroup()

summary_table <- desc_by_group %>%
  left_join(
    perm_results %>%
      select(genus, p_perm, logOR_obs),
    by = "genus"
  ) %>%
  mutate(
    p_perm_display = if_else(
      p_perm == 0,
      paste0("<", 1 / n_perm),
      as.character(round(p_perm, 4))
    )
  )

readr::write_csv(
  desc_by_group,
  file.path(out_dir_tab, "Fig3_unified_desc_by_host.csv")
)

readr::write_csv(
  perm_results,
  file.path(out_dir_tab, "Fig3_unified_permutation_results.csv")
)

readr::write_csv(
  summary_table,
  file.path(out_dir_tab, "Fig3_unified_summary_table.csv")
)

## ------------------------------------------------------------
## 8. Co-occurrence panel
## ------------------------------------------------------------
fig_prob_df <- desc_by_group %>%
  select(host_species, host_species_md, pA_given_T, pA_given_notT) %>%
  pivot_longer(
    cols = c(pA_given_T, pA_given_notT),
    names_to = "condition_raw",
    values_to = "prob"
  ) %>%
  mutate(
    condition = case_match(
      condition_raw,
      "pA_given_notT" ~ "No Trapezia",
      "pA_given_T"    ~ "Trapezia present"
    ),
    host_species_md = factor(
      host_species_md,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    ),
    condition = factor(
      condition,
      levels = c("No Trapezia", "Trapezia present")
    )
  )

fig_sig_df <- perm_results %>%
  mutate(
    host_species_md = factor(
      host_species_md,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    ),
    sig_label = case_when(
      p_perm < 0.01 ~ "**",
      p_perm < 0.05 ~ "*",
      TRUE          ~ "ns"
    )
  ) %>%
  select(host_species_md, p_perm, sig_label)

fig_sig_pos_df <- fig_prob_df %>%
  group_by(host_species_md) %>%
  summarise(
    y_max = max(prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(fig_sig_df, by = "host_species_md") %>%
  mutate(
    y_bracket = y_max + 0.04,
    y_text = y_bracket + 0.03,
    x_left = 1,
    x_right = 2,
    x_mid = 1.5
  )

y_upper_barplot <- max(fig_sig_pos_df$y_text, na.rm = TRUE) + 0.03

fig_cooccurrence <- ggplot(
  fig_prob_df,
  aes(x = condition, y = prob, fill = condition)
) +
  geom_col(width = 0.62, colour = "black") +
  facet_wrap(~ host_species_md, nrow = 1) +
  scale_x_discrete(
    labels = c(
      "No Trapezia" = "No *Trapezia*",
      "Trapezia present" = "*Trapezia* present"
    )
  ) +
  scale_y_continuous(
    limits = c(0, y_upper_barplot),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(
      "No Trapezia" = col_probbase,
      "Trapezia present" = col_trapezia
    ),
    guide = "none"
  ) +
  geom_segment(
    data = fig_sig_pos_df,
    aes(x = x_left, xend = x_right, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.65
  ) +
  geom_segment(
    data = fig_sig_pos_df,
    aes(x = x_left, xend = x_left, y = y_bracket - 0.01, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.65
  ) +
  geom_segment(
    data = fig_sig_pos_df,
    aes(x = x_right, xend = x_right, y = y_bracket - 0.01, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.65
  ) +
  geom_text(
    data = fig_sig_pos_df,
    aes(x = x_mid, y = y_text, label = sig_label),
    inherit.aes = FALSE,
    size = 5.0,
    fontface = "bold",
    colour = "black"
  ) +
  labs(
    x = NULL,
    y = "Probability of *Alpheus* presence"
  ) +
  base_theme +
  theme(
    panel.grid.major.x = element_blank()
  )

## ------------------------------------------------------------
## 9. Composite figure — main version
## ------------------------------------------------------------
panel_a_month <- fig_cooccurrence +
  labs(title = "a) Co-occurrence")

panel_b_month <- fig_time_incidence_month +
  labs(title = "b) Occupancy (%)")

panel_c_month <- fig_time_abundance_month +
  labs(title = "c) Mean abundance")

fig3_unified_month <- (
  panel_a_month /
    panel_b_month /
    panel_c_month
) +
  plot_layout(
    heights = c(0.9, 1.2, 1.2),
    guides = "collect"
  ) &
  base_theme

print(fig3_unified_month)

## ------------------------------------------------------------
## 10. Composite figure — backup version
## ------------------------------------------------------------
panel_a_occ <- fig_cooccurrence +
  labs(title = "a) Co-occurrence")

panel_b_occ <- fig_time_incidence_occ +
  labs(title = "b) Occupancy (%)")

panel_c_occ <- fig_time_abundance_occ +
  labs(title = "c) Mean abundance")

fig3_unified_occ <- (
  panel_a_occ /
    panel_b_occ /
    panel_c_occ
) +
  plot_layout(
    heights = c(0.9, 1.2, 1.2),
    guides = "collect"
  ) &
  base_theme

print(fig3_unified_occ)

## ------------------------------------------------------------
## 11. Save outputs
## ------------------------------------------------------------
ggsave(
  filename = file.path(out_dir_fig, "Fig3_unified_composite_monitoring_month.tiff"),
  plot = fig3_unified_month,
  width = 11.5,
  height = 11.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig3_unified_composite_monitoring_month.png"),
  plot = fig3_unified_month,
  width = 11.5,
  height = 11.2,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "Fig3_unified_composite_monitoring_month.pdf"),
  plot = fig3_unified_month,
  width = 11.5,
  height = 11.2,
  units = "in"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig3_unified_composite_survey_occasion.tiff"),
  plot = fig3_unified_occ,
  width = 11.5,
  height = 11.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig3_unified_composite_survey_occasion.png"),
  plot = fig3_unified_occ,
  width = 11.5,
  height = 11.2,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "Fig3_unified_composite_survey_occasion.pdf"),
  plot = fig3_unified_occ,
  width = 11.5,
  height = 11.2,
  units = "in"
)

