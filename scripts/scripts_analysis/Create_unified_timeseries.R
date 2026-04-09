## ============================================================
## Unify T1 + T2 control series for Pocillopora and Stylophora
## Creates a canonical input file for re-analysis
##
## Main idea:
## - Keep focal control colonies only (Pocillopora, Stylophora)
## - Use coral_id_raw as the biological colony identifier
## - Remove duplicated overlap rows (same colony + same date)
## - Recompute months_of_monitoring from a single global start
## - Replace split coral_uid values with unified IDs
## - Preserve provenance columns for traceability
## ============================================================

library(tidyverse)
library(janitor)
library(stringr)
library(lubridate)

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
repo_dir <- "/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR_rev"
input_dir <- file.path(repo_dir, "input")
output_dir <- file.path(repo_dir, "outputs", "output_files")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(input_dir, "controls_cleaned_with_months.csv")

## canonical unified output
output_file_unified <- file.path(input_dir, "controls_cleaned_with_months_unified.csv")

## audit outputs
output_dup_full <- file.path(output_dir, "audit_duplicate_rows_full.csv")
output_dup_consistency <- file.path(output_dir, "audit_duplicate_rows_consistency.csv")
output_removed <- file.path(output_dir, "audit_removed_duplicate_rows.csv")
output_summary <- file.path(output_dir, "audit_unified_controls_summary.csv")

## ------------------------------------------------------------
## 1. Read raw data
## ------------------------------------------------------------
raw <- readr::read_csv(
  input_file,
  col_types = readr::cols(
    collector = readr::col_character(),
    pre_monitoring_date_flag = readr::col_logical(),
    inverts_observed = readr::col_logical(),
    .default = readr::col_guess()
  )
) %>%
  janitor::clean_names() %>%
  mutate(
    date = dmy(date),
    monitoring_start = dmy(monitoring_start),
    collector_clean = str_to_lower(str_squish(collector))
  )

## ------------------------------------------------------------
## 2. Split focal controls vs everything else
## ------------------------------------------------------------
focal_controls <- raw %>%
  filter(
    genus %in% c("Pocillopora", "Stylophora"),
    treatment == "control" | transplantation %in% c("control_T1", "control_T2", "control_t1", "control_t2")
  )

non_focal_or_other <- raw %>%
  anti_join(focal_controls, by = names(raw))

## ------------------------------------------------------------
## 3. Build biological colony identifier
##    Prefer coral_id_raw; fall back to coral_id if needed
## ------------------------------------------------------------
focal_controls <- focal_controls %>%
  mutate(
    bio_colony_id = case_when(
      !is.na(coral_id_raw) & coral_id_raw != "" ~ as.character(coral_id_raw),
      !is.na(coral_id) & coral_id != "" ~ as.character(coral_id),
      TRUE ~ as.character(coral_uid)
    ),
    source_coral_uid = coral_uid,
    source_cohort = cohort,
    source_monitoring_start = monitoring_start,
    source_transplantation = transplantation
  )

## ------------------------------------------------------------
## 4. Inspect duplicated overlap rows
##    Duplicates are defined as same biological colony + same date
## ------------------------------------------------------------
dup_full <- focal_controls %>%
  group_by(bio_colony_id, genus, date) %>%
  filter(n() > 1) %>%
  arrange(bio_colony_id, date, source_cohort, source_monitoring_start) %>%
  ungroup()

readr::write_csv(dup_full, output_dup_full)

dup_consistency <- dup_full %>%
  group_by(bio_colony_id, genus, date) %>%
  summarise(
    n_rows = n(),
    n_survival_status = n_distinct(survival_status),
    n_trapezia = n_distinct(trapezia),
    n_alpheus = n_distinct(alpheus),
    n_spirobranchus = n_distinct(spirobranchus_giganteus),
    n_inverts_observed = n_distinct(inverts_observed),
    identical_biological_content =
      n_survival_status == 1 &
      n_trapezia == 1 &
      n_alpheus == 1 &
      n_spirobranchus == 1 &
      n_inverts_observed == 1,
    .groups = "drop"
  )

readr::write_csv(dup_consistency, output_dup_consistency)

## stop if any duplicated rows disagree in key fields
if (any(!dup_consistency$identical_biological_content)) {
  stop("Some duplicated overlap rows are not identical in biological content. Check audit_duplicate_rows_consistency.csv before proceeding.")
}

## ------------------------------------------------------------
## 5. Deduplicate focal controls
##    Keep earliest monitoring_start when duplicates exist
##    (this will prefer T1 rows in the overlap)
## ------------------------------------------------------------
focal_controls_unified <- focal_controls %>%
  group_by(bio_colony_id, genus, date) %>%
  arrange(source_monitoring_start, source_cohort, .by_group = TRUE) %>%
  mutate(row_rank_within_duplicate = row_number()) %>%
  ungroup()

removed_duplicate_rows <- focal_controls_unified %>%
  filter(row_rank_within_duplicate > 1)

readr::write_csv(removed_duplicate_rows, output_removed)

focal_controls_unified <- focal_controls_unified %>%
  filter(row_rank_within_duplicate == 1) %>%
  select(-row_rank_within_duplicate)

## ------------------------------------------------------------
## 6. Set unified series metadata
## ------------------------------------------------------------
global_start <- as.Date("2005-12-01")

focal_controls_unified <- focal_controls_unified %>%
  mutate(
    cohort = "Unified",
    transplantation = "control_unified",
    treatment = "control",
    monitoring_start = global_start,
    months_of_monitoring = interval(global_start, date) %/% months(1),
    pre_monitoring_date_flag = date < global_start,
    coral_uid = paste0(bio_colony_id, "_control_unified")
  )

## ------------------------------------------------------------
## 7. Add a survey occasion index based on actual dates
##    Useful for plots/tests that should use ordered observations
## ------------------------------------------------------------
survey_calendar <- focal_controls_unified %>%
  distinct(date) %>%
  arrange(date) %>%
  mutate(survey_occasion = row_number())

focal_controls_unified <- focal_controls_unified %>%
  left_join(survey_calendar, by = "date")

## ------------------------------------------------------------
## 8. Optional species relabeling check
## ------------------------------------------------------------
focal_controls_unified <- focal_controls_unified %>%
  mutate(
    species = case_when(
      genus == "Pocillopora" ~ species,
      genus == "Stylophora" ~ species,
      TRUE ~ species
    )
  )

## ------------------------------------------------------------
## 9. Recombine with everything else
## ------------------------------------------------------------
final_unified <- bind_rows(
  focal_controls_unified,
  non_focal_or_other
)

## ------------------------------------------------------------
## 10. Save canonical unified file
## ------------------------------------------------------------
readr::write_csv(final_unified, output_file_unified)

## ------------------------------------------------------------
## 11. Save audit summary
## ------------------------------------------------------------
audit_summary <- bind_rows(
  focal_controls %>%
    summarise(
      dataset = "focal_controls_before_dedup",
      n_rows = n(),
      n_colonies = n_distinct(bio_colony_id),
      min_date = min(date, na.rm = TRUE),
      max_date = max(date, na.rm = TRUE)
    ),
  focal_controls_unified %>%
    summarise(
      dataset = "focal_controls_after_dedup",
      n_rows = n(),
      n_colonies = n_distinct(bio_colony_id),
      min_date = min(date, na.rm = TRUE),
      max_date = max(date, na.rm = TRUE)
    ),
  removed_duplicate_rows %>%
    summarise(
      dataset = "removed_duplicate_rows",
      n_rows = n(),
      n_colonies = n_distinct(bio_colony_id),
      min_date = min(date, na.rm = TRUE),
      max_date = max(date, na.rm = TRUE)
    )
)

readr::write_csv(audit_summary, output_summary)

## ------------------------------------------------------------
## 12. Console checks
## ------------------------------------------------------------
message("Unified file written to:")
message(output_file_unified)

message("\nAudit summary:")
print(audit_summary)

message("\nUnique survey dates in unified focal controls:")
print(
  focal_controls_unified %>%
    distinct(date, survey_occasion, months_of_monitoring) %>%
    arrange(date)
)

message("\nCounts by genus after unification:")
print(
  focal_controls_unified %>%
    count(genus, name = "n_rows") %>%
    arrange(genus)
)

message("\nUnique biological colonies by genus:")
print(
  focal_controls_unified %>%
    count(genus, bio_colony_id, name = "n_obs") %>%
    count(genus, name = "n_colonies")
)
