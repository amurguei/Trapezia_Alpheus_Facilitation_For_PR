# ============================================================
# Trapezia-Alpheus paper heatmap
# matching other figures in the manuscript
# ============================================================

library(tidyverse)
library(janitor)
library(stringr)
library(ggplot2)
library(ggtext)
library(RColorBrewer)
library(patchwork)
library(grid)

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
setwd("/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR")

input_file <- file.path("input", "controls_cleaned_with_months.csv")
out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## Labels and colours
## ------------------------------------------------------------
cohort_map <- c(T1 = "N1", T2 = "N2")

taxon_labels <- c(
  trapezia = "*Trapezia*",
  alpheus  = "*Alpheus*"
)

species_labels <- c(
  Pocillopora = "*Pocillopora favosa*",
  Stylophora  = "*Stylophora pistillata*"
)

focal_genera <- c("Pocillopora", "Stylophora")

## ------------------------------------------------------------
## Shared theme to match other figures
## ------------------------------------------------------------
base_theme_heatmap <- theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks.y = element_blank(),
    
    axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
    
    strip.text.y = ggtext::element_markdown(size = 14.5, face = "bold", colour = "black"),
    
    plot.title = ggtext::element_markdown(size = 17, face = "bold", hjust = 0, colour = "black"),
    
    legend.title = element_text(size = 14, face = "bold", colour = "black"),
    legend.text  = element_text(size = 14, colour = "black"),
    
    legend.key.height = unit(1.3, "cm"),
    legend.key.width  = unit(0.6, "cm"),
    
    plot.margin = margin(6, 6, 6, 6)
  )

## ------------------------------------------------------------
## Read and clean data
## ------------------------------------------------------------
df <- readr::read_csv(
  input_file,
  col_types = readr::cols(
    collector = readr::col_character(),
    pre_monitoring_date_flag = readr::col_logical(),
    fish_observed = readr::col_logical(),
    inverts_observed = readr::col_logical()
  )
) %>%
  janitor::clean_names() %>%
  mutate(
    collector_clean = str_to_lower(str_squish(collector)),
    month_bin = as.integer(floor(months_of_monitoring))
  )

df0 <- df %>%
  filter(
    treatment == "control",
    is.na(collector_clean) | !str_detect(collector_clean, "yariv"),
    !is.na(cohort),
    !is.na(coral_uid),
    !is.na(genus),
    !is.na(month_bin),
    !is.na(survival_status)
  ) %>%
  mutate(
    cohort = factor(cohort, levels = c("T1", "T2")),
    genus = as.character(genus)
  ) %>%
  filter(genus %in% focal_genera) %>%
  mutate(
    genus = factor(genus, levels = focal_genera)
  )

## ------------------------------------------------------------
## Colonies eligible for plotting:
## alive at month 0 within each cohort
## ------------------------------------------------------------
get_eligible_colonies <- function(data, cohort_name) {
  m0_status <- data %>%
    filter(cohort == cohort_name, month_bin == 0) %>%
    group_by(coral_uid) %>%
    summarise(
      has_m0 = TRUE,
      any_nonalive_m0 = any(survival_status != 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  m0_status %>%
    filter(has_m0, !any_nonalive_m0) %>%
    pull(coral_uid)
}

## ------------------------------------------------------------
## Shared colony order within cohort
## rank by combined Trapezia + Alpheus abundance
## across all alive observations
## ------------------------------------------------------------
get_colony_order <- function(data, cohort_name) {
  eligible_uids <- get_eligible_colonies(data, cohort_name)
  
  if (length(eligible_uids) == 0) {
    stop(paste("No eligible colonies found for cohort", cohort_name))
  }
  
  ranking_df <- data %>%
    filter(
      cohort == cohort_name,
      coral_uid %in% eligible_uids,
      survival_status == 0
    ) %>%
    mutate(
      trapezia_val = suppressWarnings(as.numeric(trapezia)),
      alpheus_val  = suppressWarnings(as.numeric(alpheus)),
      combined_symbionts = rowSums(cbind(trapezia_val, alpheus_val), na.rm = TRUE)
    ) %>%
    group_by(coral_uid, genus) %>%
    summarise(
      total_symbionts = sum(combined_symbionts, na.rm = TRUE),
      mean_symbionts = mean(combined_symbionts, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(
      genus,
      desc(total_symbionts),
      desc(mean_symbionts),
      coral_uid
    )
  
  ranking_df$coral_uid
}

## ------------------------------------------------------------
## Prepare one panel dataset
## ------------------------------------------------------------
prepare_heatmap_data <- function(data, taxon_col, cohort_name, colony_order) {
  eligible_uids <- get_eligible_colonies(data, cohort_name)
  
  alive <- data %>%
    filter(
      cohort == cohort_name,
      coral_uid %in% eligible_uids,
      survival_status == 0
    ) %>%
    mutate(
      taxon_value = suppressWarnings(as.numeric(.data[[taxon_col]]))
    ) %>%
    select(cohort, genus, coral_uid, month_bin, taxon_value)
  
  nonalive_keys <- data %>%
    filter(
      cohort == cohort_name,
      coral_uid %in% eligible_uids,
      survival_status != 0
    ) %>%
    distinct(cohort, coral_uid, month_bin) %>%
    mutate(nonalive = TRUE)
  
  sampled_months <- alive %>%
    filter(!is.na(taxon_value)) %>%
    distinct(month_bin) %>%
    pull(month_bin) %>%
    sort()
  
  if (length(sampled_months) == 0) {
    stop(paste("No sampled months found for", taxon_col, "in cohort", cohort_name))
  }
  
  colony_meta <- data %>%
    filter(
      cohort == cohort_name,
      coral_uid %in% eligible_uids
    ) %>%
    distinct(coral_uid, genus) %>%
    mutate(
      coral_uid = factor(coral_uid, levels = colony_order)
    ) %>%
    arrange(coral_uid)
  
  grid_df <- colony_meta %>%
    tidyr::crossing(month_bin = sampled_months) %>%
    mutate(
      cohort = cohort_name,
      genus_lab = species_labels[as.character(genus)]
    ) %>%
    left_join(alive, by = c("cohort", "coral_uid", "month_bin", "genus")) %>%
    left_join(nonalive_keys, by = c("cohort", "coral_uid", "month_bin")) %>%
    mutate(
      nonalive = coalesce(nonalive, FALSE),
      class = case_when(
        !is.na(taxon_value) ~ "value",
        nonalive ~ "nonalive",
        TRUE ~ "missing"
      ),
      month_plot = month_bin + 1,
      taxon_label = taxon_labels[[taxon_col]],
      cohort_label = cohort_map[[cohort_name]],
      coral_uid_ord = factor(as.character(coral_uid), levels = rev(colony_order))
    )
  
  grid_df
}

## ------------------------------------------------------------
## Colony order for each cohort
## ------------------------------------------------------------
colony_order_T1 <- get_colony_order(df0, "T1")
colony_order_T2 <- get_colony_order(df0, "T2")

## ------------------------------------------------------------
## Build the four datasets
## ------------------------------------------------------------
d_A <- prepare_heatmap_data(df0, "trapezia", "T1", colony_order_T1)
d_B <- prepare_heatmap_data(df0, "alpheus",  "T1", colony_order_T1)
d_C <- prepare_heatmap_data(df0, "trapezia", "T2", colony_order_T2)
d_D <- prepare_heatmap_data(df0, "alpheus",  "T2", colony_order_T2)

## ------------------------------------------------------------
## Shared colour scale across all panels
## ------------------------------------------------------------
all_values <- bind_rows(d_A, d_B, d_C, d_D) %>%
  filter(class == "value", !is.na(taxon_value)) %>%
  pull(taxon_value)

fill_limits <- range(all_values, na.rm = TRUE)
fill_palette <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

## ------------------------------------------------------------
## Plot helper
## ------------------------------------------------------------
make_heatmap_panel <- function(dat, panel_letter, show_y_title = TRUE) {
  month_levels <- sort(unique(dat$month_plot))
  
  ggplot(dat, aes(x = factor(month_plot, levels = month_levels), y = coral_uid_ord)) +
    geom_tile(
      data = filter(dat, class == "value"),
      aes(fill = taxon_value),
      width = 0.95,
      height = 0.95
    ) +
    geom_tile(
      data = filter(dat, class == "nonalive"),
      fill = "grey35",
      width = 0.95,
      height = 0.95
    ) +
    geom_tile(
      data = filter(dat, class == "missing"),
      fill = "grey92",
      width = 0.95,
      height = 0.95
    ) +
    facet_grid(genus_lab ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradientn(
      colors = fill_palette,
      limits = fill_limits,
      oob = scales::squish,
      name = "Abundance"
    ) +
    labs(
      title = paste0(panel_letter, ") ", dat$taxon_label[1], " — ", dat$cohort_label[1]),
      x = "Month",
      y = if (show_y_title) "Colony" else NULL
    ) +
    base_theme_heatmap
}

## ------------------------------------------------------------
## Assemble figure
## ------------------------------------------------------------
pA <- make_heatmap_panel(d_A, "a", show_y_title = TRUE)
pB <- make_heatmap_panel(d_B, "b", show_y_title = FALSE)
pC <- make_heatmap_panel(d_C, "c", show_y_title = TRUE)
pD <- make_heatmap_panel(d_D, "d", show_y_title = FALSE)

combined_plot <- ((pA | pB) / (pC | pD)) +
  plot_layout(guides = "collect") &
  base_theme_heatmap &
  theme(legend.position = "right")

print(combined_plot)

## ------------------------------------------------------------
## Save outputs
## ------------------------------------------------------------
ggsave(
  filename = file.path(out_dir, "fig4_heatmap_trapezia_alpheus_N1_N2_pocillopora_stylophora.png"),
  plot = combined_plot,
  width = 12.5,
  height = 10.9,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "heatmap_trapezia_alpheus_N1_N2_pocillopora_stylophora.pdf"),
  plot = combined_plot,
  width = 19,
  height = 15,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, "fig4_heatmap_trapezia_alpheus_N1_N2_pocillopora_stylophora.tiff"),
  plot = combined_plot,
  width = 12.5,
  height = 10.9,
  units = "in",
  dpi = 600
)
