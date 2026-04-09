## ============================================================
## Unified non-random co-occurrence analysis
## Trapezia vs Alpheus in Pocillopora favosa and Stylophora pistillata
## ============================================================

library(tidyverse)
library(ggplot2)
library(ggtext)
library(janitor)
library(stringr)
library(scales)

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
col_probbase <- "#56B4E9"

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
  "survival_status", "treatment", "transplantation",
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

message("Retained survey months after exclusion:")
print(
  df_clean %>%
    count(format(date, "%Y-%m"), name = "n_rows") %>%
    arrange(`format(date, "%Y-%m")`)
)

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
  file.path(out_dir_tab, "cooccurrence_unified_analysis_dataset.csv")
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

## ------------------------------------------------------------
## 5. Unified co-occurrence statistics by host species
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
  file.path(out_dir_tab, "cooccurrence_unified_desc_by_host.csv")
)

readr::write_csv(
  perm_results,
  file.path(out_dir_tab, "cooccurrence_unified_permutation_results.csv")
)

readr::write_csv(
  summary_table,
  file.path(out_dir_tab, "cooccurrence_unified_summary_table.csv")
)

print(desc_by_group)
print(perm_results)
print(summary_table)

## ------------------------------------------------------------
## 6. Conditional probability figure
## ------------------------------------------------------------
fig_prob_df <- desc_by_group %>%
  select(genus, host_species, host_species_md, pA_given_T, pA_given_notT) %>%
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
      TRUE ~ "ns"
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

fig_cooccurrence_unified <- ggplot(
  fig_prob_df,
  aes(x = condition, y = prob, fill = condition)
) +
  geom_col(width = 0.62, colour = "black") +
  facet_wrap(
    ~ host_species_md,
    nrow = 1
  ) +
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
    fontface = "bold"
  ) +
  labs(
    x = NULL,
    y = "Probability of *Alpheus* presence"
  ) +
  theme_bw(base_size = 15) +
  theme(
    strip.text = ggtext::element_markdown(size = 13, face = "bold"),
    axis.text.x = ggtext::element_markdown(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = ggtext::element_markdown(size = 12.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.8),
    axis.line = element_line(linewidth = 0.6)
  )

print(fig_cooccurrence_unified)

## ------------------------------------------------------------
## 7. Save figure
## ------------------------------------------------------------
ggsave(
  filename = file.path(out_dir_fig, "Fig3_unified_cooccurrence.tiff"),
  plot = fig_cooccurrence_unified,
  width = 190,
  height = 95,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig3_unified_cooccurrence.png"),
  plot = fig_cooccurrence_unified,
  width = 190,
  height = 95,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "Fig3_unified_cooccurrence.pdf"),
  plot = fig_cooccurrence_unified,
  width = 190,
  height = 95,
  units = "mm"
)