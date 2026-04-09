library(tidyverse)
library(janitor)
library(stringr)

repo_dir <- "/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR_rev"
input_file <- file.path(repo_dir, "input", "controls_cleaned_with_months_unified.csv")
out_dir_tab <- file.path(repo_dir, "outputs_unified", "output_files")

dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

dat_raw <- readr::read_csv(
  input_file,
  col_types = readr::cols(
    collector = readr::col_character(),
    pre_monitoring_date_flag = readr::col_logical(),
    fish_observed = readr::col_logical(),
    inverts_observed = readr::col_logical()
  )
) %>%
  clean_names()

df_clean <- dat_raw %>%
  mutate(
    collector_clean = str_to_lower(str_squish(collector)),
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
    alpheus  = replace_na(as.numeric(alpheus), 0),
    host_species = case_match(
      as.character(genus),
      "Pocillopora" ~ "Pocillopora favosa",
      "Stylophora"  ~ "Stylophora pistillata"
    )
  )

dat_long <- dat %>%
  select(host_species, genus, coral_uid, trapezia, alpheus) %>%
  pivot_longer(
    cols = c(trapezia, alpheus),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  mutate(
    present = count > 0
  )

colony_intensity <- dat_long %>%
  group_by(host_species, genus, coral_uid, taxon) %>%
  summarise(
    n_present = sum(present, na.rm = TRUE),
    total_count_present = sum(count[present], na.rm = TRUE),
    intensity = if_else(
      n_present > 0,
      total_count_present / n_present,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  select(host_species, genus, coral_uid, taxon, intensity) %>%
  pivot_wider(
    names_from = taxon,
    values_from = intensity,
    names_prefix = "intensity_"
  )

readr::write_csv(
  colony_intensity,
  file.path(out_dir_tab, "colony_intensity_by_taxon_unified.csv")
)

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  
  if (sum(keep) < 3) {
    return(tibble(
      rho = NA_real_,
      p = NA_real_,
      n_colonies = sum(keep)
    ))
  }
  
  test <- suppressWarnings(
    cor.test(x[keep], y[keep], method = "spearman", exact = FALSE)
  )
  
  tibble(
    rho = unname(test$estimate),
    p = test$p.value,
    n_colonies = sum(keep)
  )
}

intensity_cor_tbl <- colony_intensity %>%
  group_by(host_species, genus) %>%
  group_modify(~ safe_spearman(.x$intensity_trapezia, .x$intensity_alpheus)) %>%
  ungroup()

print(intensity_cor_tbl)

readr::write_csv(
  intensity_cor_tbl,
  file.path(out_dir_tab, "colony_intensity_correlations_unified.csv")
)
