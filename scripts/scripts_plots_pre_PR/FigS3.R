library(tidyverse)
library(janitor)
library(stringr)
library(ggplot2)
library(ggtext)

## ------------------------------------------------------------
## Fig S3
##Paths
## ------------------------------------------------------------
setwd("/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR")

input_file <- file.path("input", "controls_cleaned_with_months.csv")
out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 1. Import raw data
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

# 2. Parse date and remove problematic months
df_clean <- dat_raw %>%
  mutate(
    collector_clean = stringr::str_to_lower(stringr::str_squish(collector)),
    date = as.Date(date, format = "%d/%m/%Y")
  ) %>%
  filter(!is.na(date)) %>%
  filter(!(format(date, "%Y-%m") %in% c("2009-12", "2010-03")))

# 3. Build cleaned analysis dataset
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

# 4. Whole-colony means
col_means_v_clean <- dat %>%
  group_by(cohort, genus, coral_uid) %>%
  summarise(
    mean_trap = mean(trapezia, na.rm = TRUE),
    mean_alph = mean(alpheus, na.rm = TRUE),
    n_months = n(),
    .groups = "drop"
  )

# 5. Build plotting dataset
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
    host_species_md = case_when(
      genus == "Pocillopora" ~ "*Pocillopora favosa*",
      genus == "Stylophora"  ~ "*Stylophora pistillata*"
    ),
    host_species_md = factor(
      host_species_md,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    )
  )

# 6. Correlation labels
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

label_pos_tbl <- cor_tbl %>%
  mutate(
    x_pos = -Inf,
    y_pos = Inf
  )

# 7. Plot
col_trapezia <- "#E69F00"
col_alpheus  <- "#0072B2"

make_colony_correlation_plot <- function(df_plot, label_pos_tbl, smoother = "lm") {
  
  p <- ggplot(df_plot, aes(x = mean_trap, y = mean_alph)) +
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
    theme_classic(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y.left = ggtext::element_markdown(size = 14, angle = 90),
      axis.title.x = ggtext::element_markdown(size = 15),
      axis.title.y = ggtext::element_markdown(size = 15),
      axis.text = element_text(size = 13),
      panel.spacing = unit(1.2, "lines")
    )
}

fig_colony_correlation_lm_v_clean <- make_colony_correlation_plot(
  df_plot = df_plot,
  label_pos_tbl = label_pos_tbl,
  smoother = "lm"
)
print(fig_colony_correlation_lm_v_clean)
ggsave(
  filename = file.path(out_dir, "Fig_S3colony_mean_trapezia_alpheus_correlation_lm_v_clean_v2.png"),
  plot = fig_colony_correlation_lm_v_clean,
  width = 10.5,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_S3colony_mean_trapezia_alpheus_correlation_lm_v_clean_v2.tiff"),
  plot = fig_colony_correlation_lm_v_clean,
  width = 10.5,
  height = 6.8,
  dpi = 600
)
