## ============================================================
## Trapezia–Alpheus analysis pipeline 
## Re-import raw data, parse dates, exclude problematic survey
## months, recompute statistics, and regenerate figures.
##
## Excluded survey months:
## - 2009-12
## - 2010-03
##
## Outputs are saved with suffix _v_clean so previous versions
## are preserved.
## ============================================================

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
setwd("/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR")

input_file <- file.path("input", "controls_cleaned_with_months.csv")
out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

library(tidyverse)
library(ggplot2)
library(ggtext)
library(scales)
library(patchwork)

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
input_file <- file.path("input", "controls_cleaned_with_months.csv")
out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## Colours
## ------------------------------------------------------------
col_trapezia <- "#E69F00"   # Okabe-Ito orange
col_alpheus  <- "#0072B2"   # Okabe-Ito darker blue
col_probbase <- "#56B4E9"   # Okabe-Ito light blue

taxon_cols <- c(
  trapezia = col_trapezia,
  alpheus  = col_alpheus
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

## Optional sanity checks
message("Date range after exclusion:")
print(range(df_clean$date, na.rm = TRUE))

message("Survey months retained:")
print(
  df_clean %>%
    count(format(date, "%Y-%m"), name = "n_rows") %>%
    arrange(`format(date, "%Y-%m")`)
)

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
    T_pres = if_else(trapezia > 0, 1L, 0L),
    A_pres = if_else(alpheus  > 0, 1L, 0L),
    species = dplyr::case_match(
      genus,
      "Pocillopora" ~ "*Pocillopora favosa*",
      "Stylophora"  ~ "*Stylophora pistillata*"
    ),
    cohort_label = dplyr::case_match(
      cohort,
      "T1" ~ "N1",
      "T2" ~ "N2"
    ),
    host_cohort = factor(
      paste(species, cohort_label, sep = " – "),
      levels = c(
        "*Pocillopora favosa* – N1",
        "*Pocillopora favosa* – N2",
        "*Stylophora pistillata* – N1",
        "*Stylophora pistillata* – N2"
      )
    )
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

plot_temporal_metric <- function(data, y_var, se_var, y_lab, percent_axis = FALSE, show_legend = TRUE) {
  
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
    aes(x = month_bin, y = y, colour = taxon, group = taxon)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0.25,
      linewidth = 0.4
    ) +
    facet_wrap(~ host_cohort, nrow = 1, scales = "free_x", labeller = label_value) +
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
    theme_bw(base_size = 12.5) +
    theme(
      strip.text = ggtext::element_markdown(face = "plain"),
      legend.text = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      panel.grid.minor = element_blank(),
      legend.position = if (show_legend) "bottom" else "none"
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

## ------------------------------------------------------------
## 5. Temporal summaries
## ------------------------------------------------------------
dat_long <- dat %>%
  select(cohort, cohort_label, genus, species, host_cohort, coral_uid, month_bin, trapezia, alpheus) %>%
  pivot_longer(
    cols = c(trapezia, alpheus),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  mutate(
    present = if_else(count > 0, 1, 0),
    taxon_label = dplyr::case_match(
      taxon,
      "trapezia" ~ "*Trapezia*",
      "alpheus"  ~ "*Alpheus*"
    )
  )

time_summary_df <- dat_long %>%
  group_by(host_cohort, species, cohort_label, genus, cohort, month_bin, taxon, taxon_label) %>%
  summarise(
    n_colony_months = n(),
    n_colonies = n_distinct(coral_uid),
    
    mean_abundance = mean(count, na.rm = TRUE),
    sd_abundance   = sd(count, na.rm = TRUE),
    se_abundance   = sd_abundance / sqrt(n_colony_months),
    
    incidence = mean(present, na.rm = TRUE),
    se_incidence = sqrt((incidence * (1 - incidence)) / n_colony_months),
    
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
## 6. Temporal figures
## ------------------------------------------------------------
fig_time_incidence_v_clean <- plot_temporal_metric(
  data = time_summary_df,
  y_var = "incidence",
  se_var = "se_incidence",
  y_lab = "Occupancy (%)",
  percent_axis = TRUE,
  show_legend = TRUE
)

fig_time_abundance_v_clean <- plot_temporal_metric(
  data = time_summary_df,
  y_var = "mean_abundance",
  se_var = "se_abundance",
  y_lab = "Mean abundance (individuals per colony)",
  percent_axis = FALSE,
  show_legend = TRUE
)

fig_time_intensity_v_clean <- plot_temporal_metric(
  data = time_summary_df,
  y_var = "mean_intensity",
  se_var = "se_intensity",
  y_lab = "Mean intensity per occupied colony",
  percent_axis = FALSE,
  show_legend = TRUE
)

print(fig_time_incidence_v_clean)
print(fig_time_abundance_v_clean)
print(fig_time_intensity_v_clean)

ggsave(
  filename = file.path(out_dir, "Fig_temporal_incidence_trapezia_alpheus_v_clean.png"),
  plot = fig_time_incidence_v_clean,
  width = 10.5, height = 3.8, dpi = 600
)
ggsave(
  filename = file.path(out_dir, "Fig_temporal_incidence_trapezia_alpheus_v_clean.pdf"),
  plot = fig_time_incidence_v_clean,
  width = 10.5, height = 3.8
)

ggsave(
  filename = file.path(out_dir, "Fig_temporal_abundance_trapezia_alpheus_v_clean.png"),
  plot = fig_time_abundance_v_clean,
  width = 10.5, height = 3.8, dpi = 600
)
ggsave(
  filename = file.path(out_dir, "Fig_temporal_abundance_trapezia_alpheus_v_clean.pdf"),
  plot = fig_time_abundance_v_clean,
  width = 10.5, height = 3.8
)

ggsave(
  filename = file.path(out_dir, "Fig_temporal_intensity_trapezia_alpheus_v_clean.png"),
  plot = fig_time_intensity_v_clean,
  width = 10.5, height = 3.8, dpi = 600
)
ggsave(
  filename = file.path(out_dir, "Fig_temporal_intensity_trapezia_alpheus_v_clean.pdf"),
  plot = fig_time_intensity_v_clean,
  width = 10.5, height = 3.8
)

## ------------------------------------------------------------
## 7. Co-occurrence statistics
## ------------------------------------------------------------
desc_by_group <- dat %>%
  group_by(cohort, genus) %>%
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
  ) %>%
  mutate(
    species = dplyr::case_match(
      genus,
      "Pocillopora" ~ "*Pocillopora favosa*",
      "Stylophora"  ~ "*Stylophora pistillata*"
    ),
    cohort_label = dplyr::case_match(
      cohort,
      "T1" ~ "N1",
      "T2" ~ "N2"
    )
  )

n_perm <- 10000

perm_results <- dat %>%
  group_by(cohort, genus) %>%
  group_modify(~ perm_test_group(.x, nperm = n_perm, seed = 123)) %>%
  ungroup() %>%
  mutate(
    species = dplyr::case_match(
      genus,
      "Pocillopora" ~ "*Pocillopora favosa*",
      "Stylophora"  ~ "*Stylophora pistillata*"
    ),
    cohort_label = dplyr::case_match(
      cohort,
      "T1" ~ "N1",
      "T2" ~ "N2"
    )
  ) %>%
  arrange(cohort, genus)

summary_table <- desc_by_group %>%
  left_join(
    perm_results %>%
      select(cohort, genus, p_perm, logOR_obs),
    by = c("cohort", "genus")
  ) %>%
  mutate(
    p_perm_display = if_else(
      p_perm == 0,
      paste0("<", 1 / n_perm),
      as.character(round(p_perm, 4))
    )
  )

readr::write_csv(
  summary_table,
  file.path(out_dir, "trapezia_alpheus_cooccurrence_summary_v_clean.csv")
)

print(desc_by_group)
print(perm_results)

## ------------------------------------------------------------
## 8. Co-occurrence figure with significance
## ------------------------------------------------------------
fig_prob_sig_df <- desc_by_group %>%
  select(cohort, cohort_label, genus, species, pA_given_T, pA_given_notT) %>%
  pivot_longer(
    cols = c(pA_given_T, pA_given_notT),
    names_to = "condition_raw",
    values_to = "prob"
  ) %>%
  mutate(
    condition = dplyr::case_match(
      condition_raw,
      "pA_given_notT" ~ "No Trapezia",
      "pA_given_T"    ~ "Trapezia present"
    ),
    host_cohort = factor(
      paste(species, cohort_label, sep = " – "),
      levels = c(
        "*Pocillopora favosa* – N1",
        "*Pocillopora favosa* – N2",
        "*Stylophora pistillata* – N1",
        "*Stylophora pistillata* – N2"
      )
    ),
    condition = factor(
      condition,
      levels = c("No Trapezia", "Trapezia present")
    )
  )

fig_sig_df <- perm_results %>%
  mutate(
    host_cohort = factor(
      paste(species, cohort_label, sep = " – "),
      levels = c(
        "*Pocillopora favosa* – N1",
        "*Pocillopora favosa* – N2",
        "*Stylophora pistillata* – N1",
        "*Stylophora pistillata* – N2"
      )
    ),
    sig_label = case_when(
      p_perm < 0.001 ~ "***",
      p_perm < 0.01  ~ "**",
      p_perm < 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  ) %>%
  select(host_cohort, p_perm, sig_label)

fig_sig_pos_df <- fig_prob_sig_df %>%
  group_by(host_cohort) %>%
  summarise(
    y_max = max(prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(fig_sig_df, by = "host_cohort") %>%
  mutate(
    y_bracket = y_max + 0.04,
    y_text = y_bracket + 0.02,
    x_left = 1,
    x_right = 2,
    x_mid = 1.5
  )

y_upper_barplot <- max(fig_sig_pos_df$y_text, na.rm = TRUE) + 0.03

fig_conditional_sig_v_clean <- ggplot(
  fig_prob_sig_df,
  aes(x = condition, y = prob, fill = condition)
) +
  geom_col(width = 0.62, colour = "black") +
  facet_wrap(~ host_cohort, nrow = 1, labeller = label_value) +
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
    linewidth = 0.5
  ) +
  geom_segment(
    data = fig_sig_pos_df,
    aes(x = x_left, xend = x_left, y = y_bracket - 0.01, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = fig_sig_pos_df,
    aes(x = x_right, xend = x_right, y = y_bracket - 0.01, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = fig_sig_pos_df,
    aes(x = x_mid, y = y_text, label = sig_label),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold"
  ) +
  labs(
    x = NULL,
    y = "Probability of *Alpheus* presence"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text = ggtext::element_markdown(face = "plain"),
    axis.text.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(fig_conditional_sig_v_clean)

ggsave(
  filename = file.path(out_dir, "Fig_Trapezia_Alpheus_conditional_probability_significance_v_clean.png"),
  plot = fig_conditional_sig_v_clean,
  width = 10.5, height = 3.8, dpi = 600
)
ggsave(
  filename = file.path(out_dir, "Fig_Trapezia_Alpheus_conditional_probability_significance_v_clean.pdf"),
  plot = fig_conditional_sig_v_clean,
  width = 10.5, height = 3.8
)

## ------------------------------------------------------------
## 9. Composite figure
## ------------------------------------------------------------
panel_a_v_clean <- fig_conditional_sig_v_clean +
  labs(title = "a)  Co-occurrence") +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 13)
  )

panel_b_v_clean <- fig_time_incidence_v_clean +
  labs(title = "b)  Occupancy (%)") +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 13)
  )

panel_c_v_clean <- fig_time_abundance_v_clean +
  labs(title = "c)  Mean abundance (individuals per colony)") +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 13)
  )

fig_composite_cooccurrence_v_clean <- (
  panel_a_v_clean /
    panel_b_v_clean /
    panel_c_v_clean
) +
  plot_layout(
    heights = c(1, 1.1, 1.1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

print(fig_composite_cooccurrence_v_clean)

ggsave(
  filename = file.path(out_dir, "Fig_composite_cooccurrence_temporal_v_clean.png"),
  plot = fig_composite_cooccurrence_v_clean,
  width = 11, height = 10, dpi = 600
)
ggsave(
  filename = file.path(out_dir, "Fig_composite_cooccurrence_temporal_v_clean.pdf"),
  plot = fig_composite_cooccurrence_v_clean,
  width = 11, height = 10
)

