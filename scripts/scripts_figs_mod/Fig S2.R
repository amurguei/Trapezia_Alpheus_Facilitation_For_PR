## ------------------------------------------------------------
## FIgS2 comparison of symbiont abundance between hosts 
##Paths
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
host_cols <- c(
  "Pocillopora favosa" = "violetred3",
  "Stylophora pistillata" = "chocolate4"
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

## ------------------------------------------------------------
## 2. Parse date and remove problematic survey months
## ------------------------------------------------------------
df_clean <- dat_raw %>%
  mutate(
    collector_clean = stringr::str_to_lower(stringr::str_squish(collector)),
    date = as.Date(date, format = "%Y-%m-%d")
  ) %>%
  filter(!is.na(date)) %>%
  filter(!(format(date, "%Y-%m") %in% c("2009-12", "2010-03")))

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
  file.path(out_dir_tab, "FigS2_host_within_symbiont_analysis_dataset.csv")
)

## ------------------------------------------------------------
## 4. Long-format observation table
## ------------------------------------------------------------
dat_long <- dat %>%
  dplyr::select(
    date,
    genus,
    host_species,
    host_species_md,
    coral_uid,
    trapezia,
    alpheus
  ) %>%
  tidyr::pivot_longer(
    cols = c(trapezia, alpheus),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  mutate(
    present = if_else(count > 0, 1L, 0L),
    taxon = factor(taxon, levels = c("alpheus", "trapezia")),
    taxon_label = dplyr::case_match(
      as.character(taxon),
      "trapezia" ~ "*Trapezia*",
      "alpheus"  ~ "*Alpheus*"
    ),
    host_species = factor(
      host_species,
      levels = c("Pocillopora favosa", "Stylophora pistillata")
    ),
    host_species_md = factor(
      host_species_md,
      levels = c("*Pocillopora favosa*", "*Stylophora pistillata*")
    )
  )

## ------------------------------------------------------------
## 5. Colony-level summaries across the unified series
## ------------------------------------------------------------
df_colony <- dat_long %>%
  group_by(
    coral_uid,
    host_species,
    host_species_md,
    taxon,
    taxon_label
  ) %>%
  summarise(
    n_observations = n(),
    n_present = sum(present, na.rm = TRUE),
    total_count = sum(count, na.rm = TRUE),
    occurrence = n_present / n_observations,
    mean_abundance = total_count / n_observations,
    intensity = if_else(
      n_present > 0,
      total_count / n_present,
      NA_real_
    ),
    .groups = "drop"
  )

readr::write_csv(
  df_colony,
  file.path(out_dir_tab, "FigS2_host_within_symbiont_df_colony.csv")
)

## ------------------------------------------------------------
## 6. Long version for tests and plotting
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
      labels = c("Occurrence", "Mean abundance", "Intensity")
    ),
    taxon = factor(
      taxon,
      levels = c("alpheus", "trapezia")
    ),
    taxon_label = factor(
      taxon_label,
      levels = c("*Alpheus*", "*Trapezia*")
    ),
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
  df_colony_long,
  file.path(out_dir_tab, "FigS2_host_within_symbiont_df_colony_long.csv")
)

## ------------------------------------------------------------
## 7. Host comparisons within each symbiont
## ------------------------------------------------------------
host_test_results_tbl <- df_colony_long %>%
  filter(!is.na(value)) %>%
  group_by(metric, taxon, taxon_label) %>%
  wilcox_test(value ~ host_species, detailed = TRUE) %>%
  ungroup() %>%
  select(
    metric,
    taxon,
    taxon_label,
    p_value = p,
    statistic
  ) %>%
  mutate(
    sig_label = case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    p_label = case_when(
      p_value < 0.001 ~ "p < 0.001",
      TRUE ~ paste0("p = ", signif(p_value, 3))
    )
  )

readr::write_csv(
  host_test_results_tbl,
  file.path(out_dir_tab, "FigS2_host_within_symbiont_wilcoxon_tests.csv")
)

print(host_test_results_tbl)

## ------------------------------------------------------------
## 8. Descriptives
## ------------------------------------------------------------
host_desc_tbl <- df_colony_long %>%
  filter(!is.na(value)) %>%
  group_by(metric, taxon, taxon_label, host_species, host_species_md) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    iqr = IQR(value, na.rm = TRUE),
    n_colonies = n(),
    .groups = "drop"
  )

readr::write_csv(
  host_desc_tbl,
  file.path(out_dir_tab, "FigS2_host_within_symbiont_descriptives.csv")
)

## ------------------------------------------------------------
## 9. Plotting data with facet labels
## ------------------------------------------------------------
host_plot_df <- df_colony_long %>%
  filter(!is.na(value)) %>%
  mutate(
    facet_lab = factor(
      paste(taxon_label, metric, sep = " — "),
      levels = c(
        "*Alpheus* — Occurrence",
        "*Alpheus* — Mean abundance",
        "*Alpheus* — Intensity",
        "*Trapezia* — Occurrence",
        "*Trapezia* — Mean abundance",
        "*Trapezia* — Intensity"
      )
    )
  )

## ------------------------------------------------------------
## 10. Significance positions
## ------------------------------------------------------------
host_sig_tbl_plot <- host_test_results_tbl %>%
  mutate(
    metric = factor(
      metric,
      levels = c("Occurrence", "Mean abundance", "Intensity")
    ),
    taxon = factor(
      taxon,
      levels = c("alpheus", "trapezia")
    ),
    taxon_label = factor(
      taxon_label,
      levels = c("*Alpheus*", "*Trapezia*")
    ),
    facet_lab = factor(
      paste(taxon_label, metric, sep = " — "),
      levels = c(
        "*Alpheus* — Occurrence",
        "*Alpheus* — Mean abundance",
        "*Alpheus* — Intensity",
        "*Trapezia* — Occurrence",
        "*Trapezia* — Mean abundance",
        "*Trapezia* — Intensity"
      )
    )
  ) %>%
  left_join(
    host_plot_df %>%
      group_by(metric, taxon, taxon_label, facet_lab) %>%
      summarise(
        y_max = max(value, na.rm = TRUE),
        y_min = min(value, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("metric", "taxon", "taxon_label", "facet_lab")
  ) %>%
  mutate(
    y_range = y_max - y_min,
    y_bracket = case_when(
      metric == "Occurrence" ~ y_max + 0.22 * if_else(y_range > 0, y_range, 1),
      TRUE ~ y_max + 0.16 * if_else(y_range > 0, y_range, 1)
    ),
    y_text = case_when(
      metric == "Occurrence" ~ y_max + 0.31 * if_else(y_range > 0, y_range, 1),
      TRUE ~ y_max + 0.24 * if_else(y_range > 0, y_range, 1)
    ),
    x_left = 1,
    x_right = 2,
    x_mid = 1.5
  )

## ------------------------------------------------------------
## 11. Plot
## ------------------------------------------------------------
fig_host_within_symbiont <- ggplot(
  host_plot_df,
  aes(x = host_species_md, y = value, colour = host_species, fill = host_species)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.40,
    linewidth = 0.85,
    alpha = 0.18
  ) +
  geom_point(
    position = position_jitter(width = 0.08),
    alpha = 0.30,
    size = 1.5,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 4.2,
    stroke = 1.0,
    fill = "white",
    colour = "black",
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ facet_lab,
    nrow = 2,
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = host_cols,
    labels = c(
      "Pocillopora favosa" = "*Pocillopora favosa*",
      "Stylophora pistillata" = "*Stylophora pistillata*"
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = host_cols,
    labels = c(
      "Pocillopora favosa" = "*Pocillopora favosa*",
      "Stylophora pistillata" = "*Stylophora pistillata*"
    ),
    name = NULL
  ) +
  geom_segment(
    data = host_sig_tbl_plot,
    aes(x = x_left, xend = x_right, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.65,
    colour = "black"
  ) +
  geom_segment(
    data = host_sig_tbl_plot,
    aes(
      x = x_left, xend = x_left,
      y = y_bracket - 0.03 * if_else(y_range > 0, y_range, 1),
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.65,
    colour = "black"
  ) +
  geom_segment(
    data = host_sig_tbl_plot,
    aes(
      x = x_right, xend = x_right,
      y = y_bracket - 0.03 * if_else(y_range > 0, y_range, 1),
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.65,
    colour = "black"
  ) +
  geom_text(
    data = host_sig_tbl_plot,
    aes(x = x_mid, y = y_text, label = sig_label),
    inherit.aes = FALSE,
    size = 5.0,
    fontface = "bold",
    colour = "black"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 15) +
  theme(
    strip.text = ggtext::element_markdown(size = 13, face = "bold"),
    axis.text.x = ggtext::element_markdown(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12.5),
    legend.text = ggtext::element_markdown(size = 12),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.8),
    axis.line = element_line(linewidth = 0.6),
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(1.0, "lines")
  )

print(fig_host_within_symbiont)

## ------------------------------------------------------------
## 12. Save figure
## ------------------------------------------------------------
ggsave(
  filename = file.path(out_dir_fig, "Fig_S_host_comparison_within_symbiont_boxplots.tiff"),
  plot = fig_host_within_symbiont,
  width = 290,
  height = 170,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig_S2_host_comparison_within_symbiont_boxplots.png"),
  plot = fig_host_within_symbiont,
  width = 290,
  height = 170,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "Fig_S2_host_comparison_within_symbiont_boxplots.pdf"),
  plot = fig_host_within_symbiont,
  width = 290,
  height = 170,
  units = "mm"
)
