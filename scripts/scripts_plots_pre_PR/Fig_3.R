## ============================================================
## Trapezia–Alpheus composite figure (Fig 3.)
##
## ============================================================

setwd("/Users/amalia/Documents/GitHub/facilitation_Trapezia_Alpheus")

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
setwd("/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR")

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
## Shared theme for all panels
## ------------------------------------------------------------
base_theme <- theme_bw(base_size = 15) +
  theme(
    strip.text = ggtext::element_markdown(size = 14, face = "bold"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = ggtext::element_markdown(size = 14, colour = "black"),
    axis.title.x = element_text(size = 14, colour = "black"),
    axis.title.y = ggtext::element_markdown(size = 14, colour = "black"),
    legend.text = ggtext::element_markdown(size = 14, colour = "black"),
    legend.title = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.8),
    axis.line = element_line(linewidth = 0.6),
    legend.position = "bottom"
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
) %>%
  janitor::clean_names()

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
      paste(species, cohort_label, sep = " - "),
      levels = c(
        "*Pocillopora favosa* - N1",
        "*Pocillopora favosa* - N2",
        "*Stylophora pistillata* - N1",
        "*Stylophora pistillata* - N2"
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
    geom_line(linewidth = 1.05) +
    geom_point(size = 2.4) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0.25,
      linewidth = 0.55
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
    base_theme +
    theme(
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
      paste(species, cohort_label, sep = " - "),
      levels = c(
        "*Pocillopora favosa* - N1",
        "*Pocillopora favosa* - N2",
        "*Stylophora pistillata* - N1",
        "*Stylophora pistillata* - N2"
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
      paste(species, cohort_label, sep = " - "),
      levels = c(
        "*Pocillopora favosa* - N1",
        "*Pocillopora favosa* - N2",
        "*Stylophora pistillata* - N1",
        "*Stylophora pistillata* - N2"
      )
    ),
    sig_label = case_when(
      p_perm < 0.001 ~ "**",
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
    linewidth = 0.6
  ) +
  geom_segment(
    data = fig_sig_pos_df,
    aes(x = x_left, xend = x_left, y = y_bracket - 0.01, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_segment(
    data = fig_sig_pos_df,
    aes(x = x_right, xend = x_right, y = y_bracket - 0.01, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.6
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
## 9. Composite figure with unified big text
## ------------------------------------------------------------
panel_a_v_clean <- fig_conditional_sig_v_clean +
  labs(title = "a) Co-occurrence")

panel_b_v_clean <- fig_time_incidence_v_clean +
  labs(title = "b) Occupancy (%)")

panel_c_v_clean <- fig_time_abundance_v_clean +
  labs(title = "c) Mean abundance")

fig_composite_cooccurrence_v_clean <- (
  panel_a_v_clean /
    panel_b_v_clean /
    panel_c_v_clean
) +
  plot_layout(
    heights = c(0.9, 1.2, 1.2),
    guides = "collect"
  ) &
  base_theme

print(fig_composite_cooccurrence_v_clean)

## ------------------------------------------------------------
## 10. Save composite figure
## ------------------------------------------------------------
ggsave(
  filename = file.path(out_dir, "Fig_3_mod_composite_cooccurrence_temporal_v_clean_bigtext.png"),
  plot = fig_composite_cooccurrence_v_clean,
  width = 12.7,
  height = 11.5,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_3_mod_composite_cooccurrence_temporal_v_clean_bigtext.pdf"),
  plot = fig_composite_cooccurrence_v_clean,
  width = 12.7,
  height = 11.5,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, "Fig_3_mod_composite_cooccurrence_temporal_v_clean_bigtext.tiff"),
  plot = fig_composite_cooccurrence_v_clean,
  width = 12.7,
  height = 11.5,
  units = "in"
)
