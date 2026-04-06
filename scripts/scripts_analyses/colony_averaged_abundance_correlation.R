## ============================================================
## Trapezia-Alpheus whole-colony correlations
## Outputs:
## - colony means table
## - cohort-specific correlations
## - pooled correlations
## - panelled correlation figure (Combined, N1, N2)
## ============================================================

library(tidyverse)
library(janitor)
library(stringr)
library(ggplot2)
library(ggtext)

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
setwd("/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR")

input_file <- file.path("input", "controls_cleaned_with_months.csv")
out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## 1. Import raw data
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
  file.path(out_dir, "trapezia_alpheus_analysis_dataset_v_clean.csv")
)

## ------------------------------------------------------------
## 4. Whole-colony means
## ------------------------------------------------------------
col_means_v_clean <- dat %>%
  group_by(cohort, genus, coral_uid) %>%
  summarise(
    mean_trap = mean(trapezia, na.rm = TRUE),
    mean_alph = mean(alpheus, na.rm = TRUE),
    n_months = n(),
    .groups = "drop"
  ) %>%
  mutate(
    cohort_label = case_match(
      as.character(cohort),
      "T1" ~ "N1",
      "T2" ~ "N2"
    ),
    genus_label = case_match(
      as.character(genus),
      "Pocillopora" ~ "*Pocillopora favosa*",
      "Stylophora"  ~ "*Stylophora pistillata*"
    )
  )

readr::write_csv(
  col_means_v_clean,
  file.path(out_dir, "trapezia_alpheus_colony_means_v_clean.csv")
)

## ------------------------------------------------------------
## 5. Helper for safe Spearman test extraction
## ------------------------------------------------------------
safe_spearman <- function(x, y) {
  test <- suppressWarnings(
    cor.test(x, y, method = "spearman", exact = FALSE)
  )
  
  ci_low <- NA_real_
  ci_high <- NA_real_
  
  ## cor.test may or may not return conf.int for Spearman depending on context
  if (!is.null(test$conf.int) && length(test$conf.int) == 2) {
    ci_low <- unname(test$conf.int[1])
    ci_high <- unname(test$conf.int[2])
  }
  
  tibble(
    rho = unname(test$estimate),
    p = test$p.value,
    ci_low = ci_low,
    ci_high = ci_high
  )
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ "p = NA",
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", formatC(p, format = "f", digits = 3))
  )
}

## ------------------------------------------------------------
## 6. Cohort-specific correlations
## ------------------------------------------------------------
col_cor_by_cohort_v_clean <- col_means_v_clean %>%
  group_by(cohort_label, genus, genus_label) %>%
  group_modify(~{
    out <- safe_spearman(.x$mean_trap, .x$mean_alph)
    out %>%
      mutate(
        n_colonies = nrow(.x),
        dataset = .y$cohort_label
      )
  }) %>%
  ungroup()

readr::write_csv(
  col_cor_by_cohort_v_clean,
  file.path(out_dir, "trapezia_alpheus_colony_correlations_by_cohort_v_clean.csv")
)

print(col_cor_by_cohort_v_clean)

## ------------------------------------------------------------
## 7. Pooled correlations across cohorts
## ------------------------------------------------------------
col_cor_combined_v_clean <- col_means_v_clean %>%
  group_by(genus, genus_label) %>%
  group_modify(~{
    out <- safe_spearman(.x$mean_trap, .x$mean_alph)
    out %>%
      mutate(
        n_colonies = nrow(.x),
        dataset = "Combined"
      )
  }) %>%
  ungroup()

readr::write_csv(
  col_cor_combined_v_clean,
  file.path(out_dir, "trapezia_alpheus_colony_correlations_combined_v_clean.csv")
)

print(col_cor_combined_v_clean)

## ------------------------------------------------------------
## 8. Merge all correlations
## ------------------------------------------------------------
col_cor_all_v_clean <- bind_rows(
  col_cor_combined_v_clean,
  col_cor_by_cohort_v_clean
) %>%
  mutate(
    dataset = factor(dataset, levels = c("Combined", "N1", "N2")),
    genus = factor(genus, levels = c("Pocillopora", "Stylophora")),
    genus_label = factor(
      genus_label,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    ),
    label_text = paste0(
      "rho = ", sprintf("%.2f", rho), "\n",
      format_p(p), "\n",
      "n = ", n_colonies
    )
  )

readr::write_csv(
  col_cor_all_v_clean,
  file.path(out_dir, "trapezia_alpheus_colony_correlations_all_v_clean.csv")
)

print(col_cor_all_v_clean)


library(tidyverse)
library(ggplot2)
library(ggtext)

out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## 1. Build plotting dataset
## assumes col_means already exists
## ------------------------------------------------------------
df_combined <- col_means_v_clean %>%
  mutate(panel = "Combined")

df_by_cohort <- col_means_v_clean %>%
  mutate(
    panel = case_when(
      as.character(cohort) == "T1" ~ "N1",
      as.character(cohort) == "T2" ~ "N2"
    )
  )

df_plot <- bind_rows(df_combined, df_by_cohort) %>%
  mutate(
    panel = factor(panel, levels = c("Combined", "N1", "N2")),
    host_species = case_when(
      genus == "Pocillopora" ~ "Pocillopora favosa",
      genus == "Stylophora"  ~ "Stylophora pistillata"
    ),
    host_species_md = case_when(
      genus == "Pocillopora" ~ "*Pocillopora favosa*",
      genus == "Stylophora"  ~ "*Stylophora pistillata*"
    ),
    host_species_md = factor(
      host_species_md,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    )
  )

## ------------------------------------------------------------
## 2. Spearman correlations by panel
## ------------------------------------------------------------
safe_spearman <- function(x, y) {
  test <- suppressWarnings(
    cor.test(x, y, method = "spearman", exact = FALSE)
  )
  
  tibble(
    rho = unname(test$estimate),
    p = test$p.value,
    n = length(x)
  )
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ "p = NA",
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", formatC(p, format = "f", digits = 3))
  )
}

cor_tbl <- df_plot %>%
  group_by(host_species_md, panel) %>%
  group_modify(~ safe_spearman(.x$mean_trap, .x$mean_alph)) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      "rho = ", sprintf("%.2f", rho), "\n",
      format_p(p), "\n",
      "n = ", n
    )
  )

## ------------------------------------------------------------
## 3. Stable top-left label placement
## ------------------------------------------------------------
label_pos_tbl <- cor_tbl %>%
  mutate(
    x_pos = -Inf,
    y_pos = Inf
  )

## ------------------------------------------------------------
## 4. Consistent colours
## ------------------------------------------------------------
col_trapezia <- "#E69F00"
col_alpheus  <- "#0072B2"

## ------------------------------------------------------------
## 5. Build BOTH versions
## ------------------------------------------------------------
make_colony_correlation_plot <- function(df_plot,
                                         label_pos_tbl,
                                         smoother = c("lm", "loess")) {
  
  smoother <- match.arg(smoother)
  
  p <- ggplot(
    df_plot,
    aes(x = mean_trap, y = mean_alph)
  ) +
    geom_point(
      colour = col_alpheus,
      alpha = 0.75,
      size = 2.8
    )
  
  if (smoother == "lm") {
    p <- p +
      geom_smooth(
        method = "lm",
        se = TRUE,
        linewidth = 1.2,
        colour = col_trapezia,
        fill = col_trapezia,
        alpha = 0.20
      )
  }
  
  if (smoother == "loess") {
    p <- p +
      geom_smooth(
        method = "loess",
        se = TRUE,
        span = 0.9,
        linewidth = 1.2,
        colour = col_trapezia,
        fill = col_trapezia,
        alpha = 0.20
      )
  }
  
  p +
    geom_text(
      data = label_pos_tbl,
      aes(x = x_pos, y = y_pos, label = label),
      inherit.aes = FALSE,
      hjust = -0.02,
      vjust = 1.25,
      size = 5
    ) +
    facet_grid(
      rows = vars(host_species_md),
      cols = vars(panel),
      scales = "free_y",
      switch = "y"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.18))
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    labs(
      x = "Mean *Trapezia* count per colony",
      y = "Mean *Alpheus* count per colony"
    ) +
    
    ## 🔥 FIGURE 2 STYLE THEME
    theme_classic(base_size = 16) +
    theme(
      ## NO GRID (key change)
      panel.grid = element_blank(),
      
      ## Clean panel borders (subtle, not boxed)
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      
      ## Facet strips
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y.left = ggtext::element_markdown(
        size = 14,
        angle = 90
      ),
      
      ## Axes
      axis.title.x = ggtext::element_markdown(size = 15),
      axis.title.y = ggtext::element_markdown(size = 15),
      axis.text = element_text(size = 13),
      
      ## Spacing (important for clean look)
      panel.spacing = unit(1.2, "lines")
    )
}
fig_colony_correlation_lm_v_clean <- make_colony_correlation_plot(
  df_plot = df_plot,
  label_pos_tbl = label_pos_tbl,
  smoother = "lm"
)

fig_colony_correlation_loess_v_clean <- make_colony_correlation_plot(
  df_plot = df_plot,
  label_pos_tbl = label_pos_tbl,
  smoother = "loess"
)

print(fig_colony_correlation_lm_v_clean)
print(fig_colony_correlation_loess_v_clean)

ggsave(
  filename = file.path(out_dir, "Fig_S3colony_mean_trapezia_alpheus_correlation_lm_v_clean_v2.png"),
  plot = fig_colony_correlation_lm_v_clean,
  width = 10.5,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_colony_mean_trapezia_alpheus_correlation_loess_v_clean_v2.png"),
  plot = fig_colony_correlation_loess_v_clean,
  width = 10.5,
  height = 6.8,
  dpi = 600
)
