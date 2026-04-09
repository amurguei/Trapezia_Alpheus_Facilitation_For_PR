library(tidyverse)
library(janitor)
library(stringr)
library(ggplot2)
library(ggtext)

repo_dir <- "/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR_rev"
input_file <- file.path(repo_dir, "input", "controls_cleaned_with_months_unified.csv")
out_dir_fig <- file.path(repo_dir, "outputs_unified", "figures")
out_dir_tab <- file.path(repo_dir, "outputs_unified", "output_files")

dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

col_trapezia <- "#E69F00"
col_alpheus  <- "#0072B2"

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

dat <- df_clean %>%
  mutate(
    month_bin = as.integer(floor(months_of_monitoring))
  ) %>%
  filter(
    genus %in% c("Pocillopora", "Stylophora"),
    treatment == "control" | transplantation == "control_unified",
    survival_status == 0,
    inverts_observed == TRUE
  ) %>%
  mutate(
    genus = factor(genus, levels = c("Pocillopora", "Stylophora")),
    coral_uid = factor(coral_uid),
    trapezia = replace_na(as.numeric(trapezia), 0),
    alpheus  = replace_na(as.numeric(alpheus), 0)
  )

readr::write_csv(
  dat,
  file.path(out_dir_tab, "FigS3_unified_analysis_dataset.csv")
)

col_means_unified <- dat %>%
  group_by(genus, coral_uid) %>%
  summarise(
    mean_trap = mean(trapezia, na.rm = TRUE),
    mean_alph = mean(alpheus, na.rm = TRUE),
    n_observations = n(),
    .groups = "drop"
  ) %>%
  mutate(
    host_species = case_match(
      as.character(genus),
      "Pocillopora" ~ "Pocillopora favosa",
      "Stylophora"  ~ "Stylophora pistillata"
    ),
    host_species_md = case_match(
      as.character(genus),
      "Pocillopora" ~ "*Pocillopora favosa*",
      "Stylophora"  ~ "*Stylophora pistillata*"
    ),
    host_species_md = factor(
      host_species_md,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    )
  )

readr::write_csv(
  col_means_unified,
  file.path(out_dir_tab, "FigS3_unified_colony_means.csv")
)

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
    TRUE ~ paste0("p = ", formatC(p, format = "f", digits = 4))
  )
}

cor_tbl_unified <- col_means_unified %>%
  group_by(host_species, host_species_md) %>%
  group_modify(~ safe_spearman(.x$mean_trap, .x$mean_alph)) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      "rho = ", sprintf("%.2f", rho), "\n",
      format_p(p), "\n",
      "n = ", n
    )
  )

readr::write_csv(
  cor_tbl_unified,
  file.path(out_dir_tab, "FigS3_unified_correlation_statistics.csv")
)

print(cor_tbl_unified)

label_pos_tbl <- col_means_unified %>%
  group_by(host_species, host_species_md) %>%
  summarise(
    x_min = min(mean_trap, na.rm = TRUE),
    x_max = max(mean_trap, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    cor_tbl_unified %>%
      select(host_species, host_species_md, label),
    by = c("host_species", "host_species_md")
  ) %>%
  mutate(
    x_pos = x_min + 0.04 * (x_max - x_min),
    y_pos = Inf
  )
make_colony_correlation_plot <- function(df_plot, label_pos_tbl, smoother = c("lm", "loess")) {
  smoother <- match.arg(smoother)
  
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
        linewidth = 1.1,
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
        linewidth = 1.1,
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
      hjust = 0,
      vjust = 1.6,
      size = 4.8
    )+
    facet_wrap(
      ~ host_species_md,
      nrow = 1,
      scales = "free"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.08, 0.18))
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.08, 0.08))
    ) +
    labs(
      x = "Mean *Trapezia* abundance per colony",
      y = "Mean *Alpheus* abundance per colony"
    ) +
    theme_classic(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text = ggtext::element_markdown(size = 14, face = "bold"),
      axis.title.x = ggtext::element_markdown(size = 15),
      axis.title.y = ggtext::element_markdown(size = 15),
      axis.text = element_text(size = 13),
      panel.spacing = unit(1.2, "lines")
    )
}

fig_colony_correlation_unified_lm <- make_colony_correlation_plot(
  df_plot = col_means_unified,
  label_pos_tbl = label_pos_tbl,
  smoother = "lm"
)

fig_colony_correlation_unified_loess <- make_colony_correlation_plot(
  df_plot = col_means_unified,
  label_pos_tbl = label_pos_tbl,
  smoother = "loess"
)

print(fig_colony_correlation_unified_lm)
print(fig_colony_correlation_unified_loess)

ggsave(
  filename = file.path(out_dir_fig, "FigS3_unified_colony_mean_correlation_lm.png"),
  plot = fig_colony_correlation_unified_lm,
  width = 8.2,
  height = 4.8,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "FigS3_unified_colony_mean_correlation_lm.pdf"),
  plot = fig_colony_correlation_unified_lm,
  width = 8.2,
  height = 4.8,
  units = "in"
)

ggsave(
  filename = file.path(out_dir_fig, "FigS3_unified_colony_mean_correlation_lm.tiff"),
  plot = fig_colony_correlation_unified_lm,
  width = 8.2,
  height = 4.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir_fig, "FigS3_unified_colony_mean_correlation_loess.png"),
  plot = fig_colony_correlation_unified_loess,
  width = 8.2,
  height = 4.8,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "FigS3_unified_colony_mean_correlation_loess.pdf"),
  plot = fig_colony_correlation_unified_loess,
  width = 8.2,
  height = 4.8,
  units = "in"
)

ggsave(
  filename = file.path(out_dir_fig, "FigS3_unified_colony_mean_correlation_loess.tiff"),
  plot = fig_colony_correlation_unified_loess,
  width = 8.2,
  height = 4.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
