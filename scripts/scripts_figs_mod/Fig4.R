library(tidyverse)
library(janitor)
library(stringr)
library(ggplot2)
library(ggtext)
library(RColorBrewer)
library(patchwork)
library(grid)

repo_dir <- "/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR_rev"
input_file <- file.path(repo_dir, "input", "controls_cleaned_with_months_unified.csv")
out_dir_fig <- file.path(repo_dir, "outputs_unified", "figures")
out_dir_tab <- file.path(repo_dir, "outputs_unified", "output_files")

dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

taxon_labels <- c(
  trapezia = "*Trapezia*",
  alpheus  = "*Alpheus*"
)

species_labels <- c(
  Pocillopora = "*Pocillopora favosa*",
  Stylophora  = "*Stylophora pistillata*"
)

focal_genera <- c("Pocillopora", "Stylophora")

title_size <- 14
legend_title_size <- 12
legend_text_size <- 12

base_theme_heatmap <- theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 13.5, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 13.5, face = "bold", colour = "black"),
    strip.text = ggtext::element_markdown(size = 13.5, face = "bold", colour = "black"),
    plot.title = ggtext::element_markdown(size = title_size, face = "bold", hjust = 0, colour = "black"),
    legend.title = element_text(size = legend_title_size, face = "bold", colour = "black"),
    legend.text = element_text(size = legend_text_size, colour = "black"),
    plot.margin = margin(2, 2, 2, 2)
  )

df <- readr::read_csv(
  input_file,
  col_types = readr::cols(
    collector = readr::col_character(),
    pre_monitoring_date_flag = readr::col_logical(),
    fish_observed = readr::col_logical(),
    inverts_observed = readr::col_logical()
  )
) %>%
  clean_names() %>%
  mutate(
    collector_clean = str_to_lower(str_squish(collector)),
    date = as.Date(date, format = "%Y-%m-%d"),
    month_bin = as.integer(floor(months_of_monitoring))
  )

df0 <- df %>%
  filter(
    genus %in% focal_genera,
    treatment == "control" | transplantation == "control_unified",
    survival_status %in% c(0, 1, 2),
    !is.na(coral_uid),
    !is.na(month_bin),
    !is.na(genus),
    !is.na(inverts_observed),
    inverts_observed == TRUE
  ) %>%
  filter(!(format(date, "%Y-%m") %in% c("2009-12", "2010-03"))) %>%
  mutate(
    genus = factor(genus, levels = focal_genera),
    genus_lab = species_labels[as.character(genus)]
  )

readr::write_csv(
  df0,
  file.path(out_dir_tab, "heatmap_unified_filtered_dataset.csv")
)

get_colony_order <- function(data, genus_name) {
  data %>%
    filter(
      genus == genus_name,
      survival_status == 0
    ) %>%
    mutate(
      trapezia_val = suppressWarnings(as.numeric(trapezia)),
      alpheus_val  = suppressWarnings(as.numeric(alpheus)),
      combined_symbionts = rowSums(cbind(trapezia_val, alpheus_val), na.rm = TRUE)
    ) %>%
    group_by(coral_uid) %>%
    summarise(
      total_symbionts = sum(combined_symbionts, na.rm = TRUE),
      mean_symbionts = mean(combined_symbionts, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(
      desc(total_symbionts),
      desc(mean_symbionts),
      coral_uid
    ) %>%
    pull(coral_uid)
}

colony_order_poc <- get_colony_order(df0, "Pocillopora")
colony_order_sty <- get_colony_order(df0, "Stylophora")

prepare_heatmap_data <- function(data, taxon_col, genus_name, colony_order) {
  alive <- data %>%
    filter(
      genus == genus_name,
      survival_status == 0
    ) %>%
    mutate(
      taxon_value = suppressWarnings(as.numeric(.data[[taxon_col]]))
    ) %>%
    select(genus, coral_uid, month_bin, taxon_value)
  
  nonalive_keys <- data %>%
    filter(
      genus == genus_name,
      survival_status != 0
    ) %>%
    distinct(genus, coral_uid, month_bin) %>%
    mutate(nonalive = TRUE)
  
  sampled_months <- alive %>%
    distinct(month_bin) %>%
    pull(month_bin) %>%
    sort()
  
  eligible_uids <- alive %>%
    filter(month_bin %in% sampled_months) %>%
    distinct(coral_uid) %>%
    pull(coral_uid)
  
  colony_levels <- colony_order[colony_order %in% eligible_uids]
  
  colony_meta <- data %>%
    filter(
      genus == genus_name,
      coral_uid %in% eligible_uids
    ) %>%
    distinct(coral_uid, genus, genus_lab) %>%
    mutate(
      coral_uid = factor(coral_uid, levels = colony_levels)
    ) %>%
    arrange(coral_uid)
  
  grid_df <- colony_meta %>%
    tidyr::crossing(month_bin = sampled_months) %>%
    left_join(alive, by = c("genus", "coral_uid", "month_bin")) %>%
    left_join(nonalive_keys, by = c("genus", "coral_uid", "month_bin")) %>%
    mutate(
      nonalive = coalesce(nonalive, FALSE),
      class = case_when(
        !is.na(taxon_value) ~ "value",
        nonalive ~ "nonalive",
        TRUE ~ "missing"
      ),
      month_display = month_bin + 1L,
      month_f = factor(month_bin, levels = sampled_months),
      taxon_label = taxon_labels[[taxon_col]],
      coral_uid_ord = factor(as.character(coral_uid), levels = rev(colony_levels))
    )
  
  empty_cols <- grid_df %>%
    group_by(coral_uid) %>%
    summarise(any_value = any(class == "value"), .groups = "drop") %>%
    filter(!any_value)
  
  if (nrow(empty_cols) > 0) {
    message("Potentially empty colonies in ", genus_name, " / ", taxon_col, ":")
    print(empty_cols)
  }
  
  attr(grid_df, "sampled_months") <- sampled_months
  attr(grid_df, "sampled_months_display") <- sampled_months + 1L
  grid_df
}

d_A <- prepare_heatmap_data(df0, "trapezia", "Pocillopora", colony_order_poc)
d_B <- prepare_heatmap_data(df0, "alpheus",  "Pocillopora", colony_order_poc)
d_C <- prepare_heatmap_data(df0, "trapezia", "Stylophora",  colony_order_sty)
d_D <- prepare_heatmap_data(df0, "alpheus",  "Stylophora",  colony_order_sty)

all_values <- bind_rows(d_A, d_B, d_C, d_D) %>%
  filter(class == "value", !is.na(taxon_value)) %>%
  pull(taxon_value)

fill_limits <- range(all_values, na.rm = TRUE)
fill_palette <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

make_legend_guide <- function(legend_position = c("right", "bottom")) {
  legend_position <- match.arg(legend_position)
  
  if (legend_position == "right") {
    guide_colorbar(
      title.position = "left",
      title.theme = element_text(
        angle = 90,
        size = legend_title_size,
        face = "bold",
        colour = "black",
        hjust = 0.5,
        vjust = 0.5
      ),
      barheight = unit(8.0, "cm"),
      barwidth = unit(0.65, "cm")
    )
  } else {
    guide_colorbar(
      title.position = "top",
      title.theme = element_text(
        size = legend_title_size,
        face = "bold",
        colour = "black",
        hjust = 0.5
      ),
      barwidth = unit(8.0, "cm"),
      barheight = unit(0.55, "cm")
    )
  }
}

make_heatmap_panel <- function(dat,
                               panel_letter,
                               show_y_title = TRUE,
                               x_axis_mode = c("every2", "all"),
                               legend_position = c("right", "bottom")) {
  x_axis_mode <- match.arg(x_axis_mode)
  legend_position <- match.arg(legend_position)
  
  sampled_months <- attr(dat, "sampled_months")
  sampled_months_display <- attr(dat, "sampled_months_display")
  
  if (x_axis_mode == "every2") {
    label_idx <- seq(1, length(sampled_months), by = 2)
    month_labels <- ifelse(
      seq_along(sampled_months) %in% label_idx,
      as.character(sampled_months_display),
      ""
    )
    axis_text_size <- 10.5
  } else {
    month_labels <- as.character(sampled_months_display)
    axis_text_size <- 9.5
  }
  
  ggplot(dat, aes(x = month_f, y = coral_uid_ord)) +
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
    scale_fill_gradientn(
      colors = fill_palette,
      limits = fill_limits,
      oob = scales::squish,
      name = "Abundance",
      guide = make_legend_guide(legend_position)
    ) +
    scale_x_discrete(
      labels = month_labels,
      drop = FALSE
    ) +
    labs(
      title = paste0(panel_letter, ") ", dat$taxon_label[1], " — ", unique(dat$genus_lab)),
      x = "Monitoring month",
      y = if (show_y_title) "Colony" else NULL
    ) +
    base_theme_heatmap +
    theme(
      axis.text.x = element_text(size = axis_text_size, colour = "black")
    )
}

build_combined_plot <- function(x_axis_mode = c("every2", "all"),
                                legend_position = c("right", "bottom")) {
  x_axis_mode <- match.arg(x_axis_mode)
  legend_position <- match.arg(legend_position)
  
  pA <- make_heatmap_panel(d_A, "a", TRUE,  x_axis_mode, legend_position)
  pB <- make_heatmap_panel(d_B, "b", FALSE, x_axis_mode, legend_position)
  pC <- make_heatmap_panel(d_C, "c", TRUE,  x_axis_mode, legend_position)
  pD <- make_heatmap_panel(d_D, "d", FALSE, x_axis_mode, legend_position)
  
  layout_plot <- ((pA | pB) / (pC | pD)) +
    plot_layout(
      guides = "collect",
      widths = c(1, 1),
      heights = c(1, 1)
    )
  
  if (legend_position == "right") {
    layout_plot &
      base_theme_heatmap &
      theme(
        legend.position = "right",
        legend.box.margin = margin(0, -8, 0, -8),
        legend.margin = margin(0, 0, 0, 0)
      )
  } else {
    layout_plot &
      base_theme_heatmap &
      theme(
        legend.position = "bottom",
        legend.box.margin = margin(-8, 0, -4, 0),
        legend.margin = margin(0, 0, 0, 0)
      )
  }
}

combined_plot_every2_right  <- build_combined_plot("every2", "right")
combined_plot_all_right     <- build_combined_plot("all", "right")
combined_plot_every2_bottom <- build_combined_plot("every2", "bottom")
combined_plot_all_bottom    <- build_combined_plot("all", "bottom")

print(combined_plot_every2_right)
print(combined_plot_all_right)
print(combined_plot_every2_bottom)
print(combined_plot_all_bottom)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_every2_right_legend.png"),
  plot = combined_plot_every2_right,
  width = 11.2,
  height = 8.1,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_every2_right_legend.pdf"),
  plot = combined_plot_every2_right,
  width = 11.2,
  height = 8.1,
  units = "in"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_every2_right_legend.tiff"),
  plot = combined_plot_every2_right,
  width = 11.2,
  height = 8.1,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_all_right_legend.png"),
  plot = combined_plot_all_right,
  width = 12,
  height = 8.1,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_all_right_legend.pdf"),
  plot = combined_plot_all_right,
  width = 12,
  height = 8.1,
  units = "in"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_all_right_legend.tiff"),
  plot = combined_plot_all_right,
  width = 12,
  height = 8.1,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_every2_bottom_legend.png"),
  plot = combined_plot_every2_bottom,
  width = 10.4,
  height = 8.8,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_every2_bottom_legend.pdf"),
  plot = combined_plot_every2_bottom,
  width = 10.4,
  height = 8.8,
  units = "in"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_every2_bottom_legend.tiff"),
  plot = combined_plot_every2_bottom,
  width = 10.4,
  height = 8.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_all_bottom_legend.png"),
  plot = combined_plot_all_bottom,
  width = 11.5,
  height = 8.8,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_all_bottom_legend.pdf"),
  plot = combined_plot_all_bottom,
  width = 11.5,
  height = 8.8,
  units = "in"
)

ggsave(
  filename = file.path(out_dir_fig, "Fig4_unified_heatmap_all_bottom_legend.tiff"),
  plot = combined_plot_all_bottom,
  width = 12,
  height = 8.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
