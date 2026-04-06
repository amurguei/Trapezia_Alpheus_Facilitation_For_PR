library(tidyverse)
library(ggplot2)
library(ggtext)
library(car)
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
## 3. Build cleaned analysis dataset from scratch
## ------------------------------------------------------------
dat <- df_clean %>%
  mutate(
    month_bin = as.integer(floor(months_of_monitoring))
  ) %>%
  filter(
    (treatment == "control") |
      (transplantation %in% c("control_t1", "control_t2", "control_T1", "control_T2")),
    survival_status == 0,
    inverts_observed == TRUE,
    genus %in% c("Pocillopora", "Stylophora")
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
  file.path(out_dir, "trapezia_alpheus_dat_v_clean.csv")
)

## ------------------------------------------------------------
## 4. Build dat_long from dat
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
    present = dplyr::if_else(count > 0, 1L, 0L),
    host_species = dplyr::case_match(
      as.character(genus),
      "Pocillopora" ~ "Pocillopora favosa",
      "Stylophora"  ~ "Stylophora pistillata"
    ),
    host_species_md = dplyr::case_match(
      as.character(genus),
      "Pocillopora" ~ "*Pocillopora favosa*",
      "Stylophora"  ~ "*Stylophora pistillata*"
    ),
    taxon_label = dplyr::case_match(
      taxon,
      "trapezia" ~ "*Trapezia*",
      "alpheus"  ~ "*Alpheus*"
    ),
    cohort_label = dplyr::case_match(
      as.character(cohort),
      "T1" ~ "N1",
      "T2" ~ "N2"
    )
  )

readr::write_csv(
  dat_long,
  file.path(out_dir, "trapezia_alpheus_dat_long_v_clean.csv")
)

print(dat_long)
## ------------------------------------------------------------
## 5. Build colony-level data from dat_long
## This matches the boxplots exactly
## ------------------------------------------------------------
df_colony <- dat_long %>%
  group_by(
    coral_uid,
    host_species,
    host_species_md,
    cohort_label,
    taxon,
    taxon_label
  ) %>%
  summarise(
    n_months = n(),
    n_present = sum(present, na.rm = TRUE),
    total_count = sum(count, na.rm = TRUE),
    
    occurrence = n_present / n_months,
    mean_abundance = total_count / n_months,
    intensity = if_else(
      n_present > 0,
      total_count / n_present,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    taxon = factor(taxon, levels = c("trapezia", "alpheus")),
    host_species = factor(
      host_species,
      levels = c("Pocillopora favosa", "Stylophora pistillata")
    ),
    host_species_md = factor(
      host_species_md,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    )
  )

readr::write_csv(
  df_colony,
  file.path(out_dir, "trapezia_alpheus_df_colony_boxplotmatched_v_clean.csv")
)

print(df_colony)

## ------------------------------------------------------------
## 6. Long version for metric-wise testing and plotting
## ------------------------------------------------------------
df_colony_long <- df_colony %>%
  pivot_longer(
    cols = c(occurrence, mean_abundance, intensity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("occurrence", "mean_abundance", "intensity"),
      labels = c(
        "Occurrence",
        "Mean abundance",
        "Intensity"
      )
    )
  )

readr::write_csv(
  df_colony_long,
  file.path(out_dir, "trapezia_alpheus_df_colony_long_boxplotmatched_v_clean.csv")
)

print(df_colony_long)


## ============================================================
## Colony-level Trapezia vs Alpheus comparisons within each host
## Assumption-based tests matching the boxplots exactly
##
## Input:
##   dat_long
##
## Output:
##   - colony-level summary table
##   - assumption checks
##   - chosen tests + p-values
##   - QQ plots / histograms / boxplots
##
## Decision rule:
##   1) Shapiro-Wilk in each taxon group within host
##   2) Levene test for equal variances
##   3) If both groups normal:
##        - equal variances   -> Student t-test
##        - unequal variances -> Welch t-test
##      Else:
##        - Wilcoxon rank-sum
## ============================================================

library(tidyverse)
library(ggplot2)
library(ggtext)
library(car)

out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## 1. Build colony-level data from dat_long
## This matches the boxplots exactly
## ------------------------------------------------------------
df_colony <- dat_long %>%
  group_by(
    coral_uid,
    host_species,
    host_species_md,
    cohort_label,
    taxon,
    taxon_label
  ) %>%
  summarise(
    n_months = n(),
    n_present = sum(present, na.rm = TRUE),
    total_count = sum(count, na.rm = TRUE),
    
    occurrence = n_present / n_months,
    mean_abundance = total_count / n_months,
    intensity = if_else(
      n_present > 0,
      total_count / n_present,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    taxon = factor(taxon, levels = c("trapezia", "alpheus")),
    host_species = factor(
      host_species,
      levels = c("Pocillopora favosa", "Stylophora pistillata")
    ),
    host_species_md = factor(
      host_species_md,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    )
  )

readr::write_csv(
  df_colony,
  file.path(out_dir, "trapezia_alpheus_df_colony_boxplotmatched_v_clean.csv")
)

print(df_colony)

## ------------------------------------------------------------
## 2. Long version for metric-wise testing and plotting
## ------------------------------------------------------------
df_colony_long <- df_colony %>%
  pivot_longer(
    cols = c(occurrence, mean_abundance, intensity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("occurrence", "mean_abundance", "intensity"),
      labels = c(
        "Occurrence",
        "Mean abundance",
        "Intensity"
      )
    )
  )

readr::write_csv(
  df_colony_long,
  file.path(out_dir, "trapezia_alpheus_df_colony_long_boxplotmatched_v_clean.csv")
)

## ------------------------------------------------------------
## 3. Assumption checks
## - Shapiro within each host x metric x taxon
## - Levene within each host x metric
## ------------------------------------------------------------
safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  if (length(unique(x)) < 3) return(NA_real_)
  shapiro.test(x)$p.value
}

shapiro_tbl <- df_colony_long %>%
  filter(!is.na(value)) %>%
  group_by(host_species, host_species_md, metric, taxon, taxon_label) %>%
  summarise(
    n = n(),
    shapiro_p = safe_shapiro(value),
    .groups = "drop"
  )

levene_tbl <- df_colony_long %>%
  filter(!is.na(value)) %>%
  group_by(host_species, host_species_md, metric) %>%
  group_modify(~{
    lev <- car::leveneTest(value ~ taxon, data = .x)
    tibble(
      levene_p = lev[1, "Pr(>F)"]
    )
  }) %>%
  ungroup()

assumption_tbl <- shapiro_tbl %>%
  left_join(levene_tbl, by = c("host_species", "host_species_md", "metric"))

readr::write_csv(
  assumption_tbl,
  file.path(out_dir, "trapezia_alpheus_assumption_checks_v_clean.csv")
)

print(assumption_tbl)

## ------------------------------------------------------------
## 4. Choose tests automatically
## ------------------------------------------------------------
run_assumption_based_test <- function(d) {
  d <- d %>% filter(!is.na(value))
  
  x <- d %>% filter(taxon == "trapezia") %>% pull(value)
  y <- d %>% filter(taxon == "alpheus") %>% pull(value)
  
  shapiro_x <- safe_shapiro(x)
  shapiro_y <- safe_shapiro(y)
  
  lev <- car::leveneTest(value ~ taxon, data = d)
  lev_p <- lev[1, "Pr(>F)"]
  
  normal_x <- !is.na(shapiro_x) && shapiro_x > 0.05
  normal_y <- !is.na(shapiro_y) && shapiro_y > 0.05
  both_normal <- normal_x && normal_y
  equal_var <- !is.na(lev_p) && lev_p > 0.05
  
  if (both_normal && equal_var) {
    test_obj <- t.test(value ~ taxon, data = d, var.equal = TRUE)
    method <- "Student t-test"
    stat_name <- "t"
    stat_value <- unname(test_obj$statistic)
  } else if (both_normal && !equal_var) {
    test_obj <- t.test(value ~ taxon, data = d, var.equal = FALSE)
    method <- "Welch t-test"
    stat_name <- "t"
    stat_value <- unname(test_obj$statistic)
  } else {
    test_obj <- wilcox.test(value ~ taxon, data = d, exact = FALSE)
    method <- "Wilcoxon rank-sum"
    stat_name <- "W"
    stat_value <- unname(test_obj$statistic)
  }
  
  tibble(
    method = method,
    statistic_name = stat_name,
    statistic_value = stat_value,
    p_value = test_obj$p.value,
    shapiro_trapezia_p = shapiro_x,
    shapiro_alpheus_p = shapiro_y,
    levene_p = lev_p,
    n_trapezia = length(x),
    n_alpheus = length(y),
    mean_trapezia = mean(x, na.rm = TRUE),
    mean_alpheus = mean(y, na.rm = TRUE),
    median_trapezia = median(x, na.rm = TRUE),
    median_alpheus = median(y, na.rm = TRUE),
    sd_trapezia = sd(x, na.rm = TRUE),
    sd_alpheus = sd(y, na.rm = TRUE)
  )
}

test_results_tbl <- df_colony_long %>%
  group_by(host_species, host_species_md, metric) %>%
  group_modify(~ run_assumption_based_test(.x)) %>%
  ungroup()

readr::write_csv(
  test_results_tbl,
  file.path(out_dir, "trapezia_alpheus_boxplotmatched_tests_v_clean.csv")
)

print(test_results_tbl)

## ------------------------------------------------------------
## 5. Descriptive summary table
## ------------------------------------------------------------
descriptive_tbl <- df_colony_long %>%
  filter(!is.na(value)) %>%
  group_by(metric, taxon, taxon_label, host_species, host_species_md) %>%
  summarise(
    n_colonies = n(),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    iqr = IQR(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  descriptive_tbl,
  file.path(out_dir, "trapezia_alpheus_boxplotmatched_descriptives_v_clean.csv")
)

print(descriptive_tbl)

## ------------------------------------------------------------
## 6. Diagnostic plots
## QQ plots + histograms + boxplots
## ------------------------------------------------------------

## 6a. QQ plots
qq_plot <- ggplot(
  df_colony_long %>% filter(!is.na(value)),
  aes(sample = value, colour = taxon)
) +
  stat_qq(size = 1.3, alpha = 0.7) +
  stat_qq_line(linewidth = 0.5) +
  facet_grid(metric ~ host_species_md, scales = "free") +
  scale_colour_manual(
    values = c(
      "trapezia" = "#E69F00",
      "alpheus" = "#0072B2"
    ),
    labels = c(
      "trapezia" = "*Trapezia*",
      "alpheus" = "*Alpheus*"
    ),
    name = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text.x = ggtext::element_markdown(face = "plain"),
    strip.text.y = element_text(face = "bold"),
    legend.text = ggtext::element_markdown(),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Theoretical quantiles",
    y = "Sample quantiles",
    title = "QQ plots for colony-level metrics"
  )

print(qq_plot)

ggsave(
  filename = file.path(out_dir, "Fig_boxplotmatched_QQplots_v_clean.png"),
  plot = qq_plot,
  width = 10.5,
  height = 7.2,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_boxplotmatched_QQplots_v_clean.pdf"),
  plot = qq_plot,
  width = 10.5,
  height = 7.2
)

## 6b. Histograms
hist_plot <- ggplot(
  df_colony_long %>% filter(!is.na(value)),
  aes(x = value, fill = taxon)
) +
  geom_histogram(
    bins = 18,
    alpha = 0.6,
    position = "identity"
  ) +
  facet_grid(metric ~ host_species_md, scales = "free") +
  scale_fill_manual(
    values = c(
      "trapezia" = "#E69F00",
      "alpheus" = "#0072B2"
    ),
    labels = c(
      "trapezia" = "*Trapezia*",
      "alpheus" = "*Alpheus*"
    ),
    name = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text.x = ggtext::element_markdown(face = "plain"),
    strip.text.y = element_text(face = "bold"),
    legend.text = ggtext::element_markdown(),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Value",
    y = "Frequency",
    title = "Histograms for colony-level metrics"
  )

print(hist_plot)

ggsave(
  filename = file.path(out_dir, "Fig_boxplotmatched_histograms_v_clean.png"),
  plot = hist_plot,
  width = 10.5,
  height = 7.2,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_boxplotmatched_histograms_v_clean.pdf"),
  plot = hist_plot,
  width = 10.5,
  height = 7.2
)

## 6c. Boxplots + jitter + mean diamonds
pd <- position_dodge(width = 0.75)

boxplotmatched_fig <- ggplot(
  df_colony_long %>% filter(!is.na(value)),
  aes(
    x = host_species_md,
    y = value,
    colour = taxon
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.22,
    position = pd,
    linewidth = 0.7
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    ),
    alpha = 0.7,
    size = 2.2
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    position = pd,
    shape = 23,
    size = 3.8,
    stroke = 0.9,
    fill = "white"
  ) +
  facet_wrap(
    ~ metric,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = c(
      "trapezia" = "#E69F00",
      "alpheus" = "#0072B2"
    ),
    labels = c(
      "trapezia" = "*Trapezia*",
      "alpheus" = "*Alpheus*"
    ),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(boxplotmatched_fig)

ggsave(
  filename = file.path(out_dir, "Fig_host_comparison_boxplotmatched_tests_v_clean.png"),
  plot = boxplotmatched_fig,
  width = 10.5,
  height = 4.4,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_host_comparison_boxplotmatched_tests_v_clean.pdf"),
  plot = boxplotmatched_fig,
  width = 10.5,
  height = 4.4
)

## ------------------------------------------------------------
## 7.  compact per-metric interpretation table
## ------------------------------------------------------------
interpretation_tbl <- test_results_tbl %>%
  mutate(
    direction = case_when(
      mean_trapezia > mean_alpheus ~ "Trapezia > Alpheus",
      mean_trapezia < mean_alpheus ~ "Alpheus > Trapezia",
      TRUE ~ "Equal means"
    ),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

readr::write_csv(
  interpretation_tbl,
  file.path(out_dir, "trapezia_alpheus_boxplotmatched_interpretation_v_clean.csv")
)

print(interpretation_tbl)


library(tidyverse)
library(ggplot2)
library(ggtext)
library(car)
library(rstatix)
library(ggpubr)

run_auto_test <- function(df) {
  
  x1 <- df$value[df$taxon == "trapezia"]
  x2 <- df$value[df$taxon == "alpheus"]
  
  # Normality
  p1 <- if (length(x1) >= 3) shapiro.test(x1)$p.value else NA
  p2 <- if (length(x2) >= 3) shapiro.test(x2)$p.value else NA
  
  # Homogeneity
  lev <- tryCatch(
    car::leveneTest(value ~ taxon, data = df)$`Pr(>F)`[1],
    error = function(e) NA
  )
  
  if (!is.na(p1) && !is.na(p2) &&
      p1 > 0.05 && p2 > 0.05 &&
      !is.na(lev) && lev > 0.05) {
    
    test <- t.test(value ~ taxon, data = df, var.equal = TRUE)
    method <- "t-test"
    
  } else if (!is.na(p1) && !is.na(p2) &&
             p1 > 0.05 && p2 > 0.05) {
    
    test <- t.test(value ~ taxon, data = df)
    method <- "Welch t-test"
    
  } else {
    
    test <- wilcox.test(value ~ taxon, data = df)
    method <- "Wilcoxon"
  }
  
  tibble(
    method = method,
    statistic = test$statistic,
    p = test$p.value
  )
}

stats_tbl <- stats_tbl %>%
  group_by(metric, host_species) %>%
  mutate(
    y.position = max(df_colony_long$value[
      df_colony_long$metric == unique(metric) &
        df_colony_long$host_species == unique(host_species)
    ], na.rm = TRUE) * 1.08,
    group1 = "trapezia",
    group2 = "alpheus"
  ) %>%
  ungroup()

desc_tbl <- df_colony_long %>%
  group_by(host_species, metric, taxon) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = n(),
    .groups = "drop"
  )

p <- ggplot(df_colony_long,
            aes(x = taxon, y = value, fill = taxon)) +
  
  geom_boxplot(width = 0.65, alpha = 0.7, outlier.shape = NA) +
  
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  
  facet_grid(metric ~ host_species, scales = "free_y") +
  
  scale_fill_manual(values = c("#D55E00", "#0072B2")) +
  
  stat_pvalue_manual(
    stats_tbl,
    label = "p_label",
    tip.length = 0,
    inherit.aes = FALSE
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  
  labs(
    x = NULL,
    y = "Colony-level value",
    title = "Colony-level comparisons of Trapezia and Alpheus across hosts"
  )

print(p)


library(tidyverse)
library(ggplot2)
library(ggtext)

## ------------------------------------------------------------
## Significance labels for annotation
## ------------------------------------------------------------
sig_tbl <- test_results_tbl %>%
  mutate(
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  left_join(
    df_colony_long %>%
      group_by(host_species, metric) %>%
      summarise(
        y_max = max(value, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("host_species", "metric")
  ) %>%
  mutate(
    x_left = 1,
    x_right = 2,
    x_mid = 1.5,
    y_bracket = y_max * 1.06,
    y_text = y_max * 1.11
  )

## ------------------------------------------------------------
## Compact one-line boxplot figure
## ------------------------------------------------------------
fig_host_comparison_stats_v_clean <- ggplot(
  df_colony_long,
  aes(x = taxon, y = value, colour = taxon)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55,
    linewidth = 0.7
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0),
    alpha = 0.55,
    size = 1.8
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3.2,
    stroke = 0.9,
    fill = "white"
  ) +
  facet_grid(
    . ~ metric + host_species_md,
    scales = "free_y",
    space = "free_x"
  ) +
  scale_colour_manual(
    values = c(
      "trapezia" = "#E69F00",
      "alpheus"  = "#0072B2"
    ),
    labels = c(
      "trapezia" = "*Trapezia*",
      "alpheus"  = "*Alpheus*"
    ),
    name = NULL
  ) +
  geom_segment(
    data = sig_tbl,
    aes(x = x_left, xend = x_right, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +
  geom_segment(
    data = sig_tbl,
    aes(x = x_left, xend = x_left, y = y_bracket * 0.985, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +
  geom_segment(
    data = sig_tbl,
    aes(x = x_right, xend = x_right, y = y_bracket * 0.985, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +
  geom_text(
    data = sig_tbl,
    aes(x = x_mid, y = y_text, label = sig_label),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12.5) +
  theme(
    strip.text.x = ggtext::element_markdown(face = "bold"),
    axis.text.x = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.8, "lines")
  )

print(fig_host_comparison_stats_v_clean)

ggsave(
  filename = file.path(out_dir, "Fig_host_comparison_boxplotmatched_stats_v_clean.png"),
  plot = fig_host_comparison_stats_v_clean,
  width = 12.5,
  height = 3.9,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_host_comparison_boxplotmatched_stats_v_clean.pdf"),
  plot = fig_host_comparison_stats_v_clean,
  width = 12.5,
  height = 3.9
)


library(tidyverse)
library(ggplot2)
library(ggtext)

## ------------------------------------------------------------
## Plot data
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
## Significance table with one bracket per host per metric
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
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
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

## optional: hide ns if you prefer
# sig_tbl_plot <- sig_tbl_plot %>% filter(sig_label != "ns")

## ------------------------------------------------------------
## One-row compact figure
## ------------------------------------------------------------
pd <- position_dodge(width = dodge_width)

fig_host_comparison_stats_v_clean <- ggplot(
  df_plot,
  aes(x = host_species_md, y = value, colour = taxon)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.22,
    position = pd,
    linewidth = 0.7
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = dodge_width
    ),
    alpha = 0.7,
    size = 2.1
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    position = pd,
    shape = 23,
    size = 3.8,
    stroke = 0.9,
    fill = "white"
  ) +
  facet_wrap(
    ~ metric,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = c(
      "trapezia" = "#E69F00",
      "alpheus"  = "#0072B2"
    ),
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
    linewidth = 0.5
  ) +
  geom_segment(
    data = sig_tbl_plot,
    aes(
      x = x_left, xend = x_left,
      y = y_bracket - 0.03 * if_else(y_range > 0, y_range, 1),
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = sig_tbl_plot,
    aes(
      x = x_right, xend = x_right,
      y = y_bracket - 0.03 * if_else(y_range > 0, y_range, 1),
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = sig_tbl_plot,
    aes(x = x_mid, y = y_text, label = sig_label),
    inherit.aes = FALSE,
    size = 4.6,
    fontface = "bold"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(1.2, "lines")
  )

print(fig_host_comparison_stats_v_clean)

ggsave(
  filename = file.path(out_dir, "Fig_host_comparison_boxplotmatched_stats_v_clean_v3.png"),
  plot = fig_host_comparison_stats_v_clean,
  width = 9.8,
  height = 4.0,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_host_comparison_boxplotmatched_stats_v_clean_v3.pdf"),
  plot = fig_host_comparison_stats_v_clean,
  width = 9.8,
  height = 4.0
)
