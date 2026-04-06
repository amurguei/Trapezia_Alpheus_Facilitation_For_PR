## ============================================================
## Trapezia–Alpheus abundance coupling
##
## Goal:
## Test whether Alpheus abundance scales with Trapezia abundance
## separately for each host genus and cohort, 
## Steps:
## 1) import raw file from scratch
## 2) parse date correctly
## 3) exclude problematic months (2009-12, 2010-03)
## 4) restrict to control colonies, alive months, valid invertebrate surveys
## 5) fit NB-GLMMs separately by cohort x genus
## 6) export results
## ============================================================

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
setwd("/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR")

input_file <- file.path("input", "controls_cleaned_with_months.csv")
out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

library(tidyverse)
library(janitor)
library(stringr)
library(glmmTMB)
library(broom.mixed)
library(DHARMa)
library(performance)

out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path("input", "controls_cleaned_with_months.csv")

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

## Optional sanity checks
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

## Save a cleaned analysis snapshot
readr::write_csv(
  dat,
  file.path(out_dir, "trapezia_alpheus_analysis_dataset_v_clean.csv")
)

## ------------------------------------------------------------
## 4. Create dat2 explicitly for modeling
## ------------------------------------------------------------
dat2 <- dat %>%
  mutate(
    cohort = factor(cohort, levels = c("T1", "T2")),
    genus = factor(genus, levels = c("Pocillopora", "Stylophora")),
    coral_uid = factor(coral_uid),
    month_num = as.numeric(months_of_monitoring),
    trapezia = as.numeric(trapezia),
    alpheus = as.numeric(alpheus)
  )

## ------------------------------------------------------------
## 5. Fit one NB-GLMM per host x cohort
## Model: alpheus ~ trapezia + (1 | coral_uid)
## ------------------------------------------------------------
fit_nb_group <- function(d) {
  if (nrow(d) < 20) return(NULL)
  if (dplyr::n_distinct(d$coral_uid) < 2) return(NULL)
  
  tryCatch(
    glmmTMB(
      alpheus ~ trapezia + (1 | coral_uid),
      family = nbinom2(),
      data = d
    ),
    error = function(e) NULL
  )
}

results_tbl <- dat2 %>%
  group_by(cohort, genus) %>%
  group_modify(~{
    m <- fit_nb_group(.x)
    
    if (is.null(m)) {
      return(tibble(
        beta = NA_real_,
        se = NA_real_,
        z = NA_real_,
        p = NA_real_,
        irr = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        n_obs = nrow(.x),
        n_colonies = dplyr::n_distinct(.x$coral_uid),
        converged = FALSE
      ))
    }
    
    s <- summary(m)$coefficients$cond
    beta <- s["trapezia", "Estimate"]
    se   <- s["trapezia", "Std. Error"]
    z    <- s["trapezia", "z value"]
    p    <- s["trapezia", "Pr(>|z|)"]
    
    tibble(
      beta = beta,
      se = se,
      z = z,
      p = p,
      irr = exp(beta),
      conf_low = exp(beta - 1.96 * se),
      conf_high = exp(beta + 1.96 * se),
      n_obs = nrow(.x),
      n_colonies = dplyr::n_distinct(.x$coral_uid),
      converged = TRUE
    )
  }) %>%
  ungroup()

print(results_tbl)

readr::write_csv(
  results_tbl,
  file.path(out_dir, "trapezia_alpheus_abundance_coupling_results_v_clean.csv")
)

## ------------------------------------------------------------
## 6. Save model objects for diagnostics / summaries
## ------------------------------------------------------------
model_list <- dat2 %>%
  group_by(cohort, genus) %>%
  group_split() %>%
  setNames(
    dat2 %>%
      group_by(cohort, genus) %>%
      group_keys() %>%
      tidyr::unite("label", cohort, genus, sep = "_") %>%
      pull(label)
  ) %>%
  purrr::map(fit_nb_group)

## ------------------------------------------------------------
## 7. Model diagnostics
## ------------------------------------------------------------
check_model_group <- function(model) {
  if (is.null(model)) return(NULL)
  
  res <- DHARMa::simulateResiduals(model, n = 1000)
  
  list(
    overdispersion = performance::check_overdispersion(model),
    zeroinflation = performance::check_zeroinflation(model),
    singularity = performance::check_singularity(model),
    uniformity_test = DHARMa::testUniformity(res),
    dispersion_test = DHARMa::testDispersion(res),
    zero_test = DHARMa::testZeroInflation(res),
    residual_object = res
  )
}

diagnostics_list <- purrr::map(model_list, check_model_group)

## Optional printed summaries
message("Model summaries:")
for (nm in names(model_list)) {
  message("---- ", nm, " ----")
  if (is.null(model_list[[nm]])) {
    message("Model not fit")
  } else {
    print(summary(model_list[[nm]]))
  }
}

## ------------------------------------------------------------
## 8. Colony-level mean correlations 
## ------------------------------------------------------------
col_means <- dat2 %>%
  group_by(cohort, genus, coral_uid) %>%
  summarise(
    mean_trap = mean(trapezia, na.rm = TRUE),
    mean_alph = mean(alpheus, na.rm = TRUE),
    .groups = "drop"
  )

col_cor <- col_means %>%
  group_by(cohort, genus) %>%
  summarise(
    rho = suppressWarnings(cor(mean_trap, mean_alph, method = "spearman")),
    p = tryCatch(
      cor.test(mean_trap, mean_alph, method = "spearman")$p.value,
      error = function(e) NA_real_
    ),
    n_colonies = n(),
    .groups = "drop"
  )

print(col_cor)

readr::write_csv(
  col_means,
  file.path(out_dir, "trapezia_alpheus_colony_means_v_clean.csv")
)

readr::write_csv(
  col_cor,
  file.path(out_dir, "trapezia_alpheus_colony_correlations_v_clean.csv")
)

## ------------------------------------------------------------
## 9. prediction plots by host x cohort
## ------------------------------------------------------------
plot_prediction_group <- function(model, data_used, cohort_name, genus_name) {
  if (is.null(model)) return(NULL)
  
  newdat <- tibble(
    trapezia = seq(
      from = min(data_used$trapezia, na.rm = TRUE),
      to   = max(data_used$trapezia, na.rm = TRUE),
      length.out = 100
    ),
    coral_uid = data_used$coral_uid[1]
  )
  
  pred <- predict(model, newdata = newdat, type = "response", se.fit = TRUE)
  
  pred_df <- newdat %>%
    mutate(
      predicted = pred$fit,
      se = pred$se.fit,
      conf_low = pmax(predicted - 1.96 * se, 0),
      conf_high = predicted + 1.96 * se,
      cohort = cohort_name,
      genus = genus_name
    )
  
  ggplot(pred_df, aes(x = trapezia, y = predicted)) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.2) +
    geom_line(linewidth = 1) +
    theme_bw(base_size = 12) +
    labs(
      x = "Trapezia count",
      y = "Predicted Alpheus count",
      title = paste0(genus_name, " — ", cohort_name)
    )
}

## Build prediction plots if desired
prediction_plots <- list()
group_keys_tbl <- dat2 %>% group_by(cohort, genus) %>% group_keys()

for (i in seq_len(nrow(group_keys_tbl))) {
  coh <- as.character(group_keys_tbl$cohort[i])
  gen <- as.character(group_keys_tbl$genus[i])
  label <- paste(coh, gen, sep = "_")
  
  d_sub <- dat2 %>%
    filter(cohort == coh, genus == gen)
  
  prediction_plots[[label]] <- plot_prediction_group(
    model = model_list[[label]],
    data_used = d_sub,
    cohort_name = coh,
    genus_name = gen
  )
}

## Example print
print(prediction_plots$T1_Pocillopora)

## ------------------------------------------------------------
## 10. Save a compact diagnostics summary table
## ------------------------------------------------------------
diagnostics_summary <- purrr::imap_dfr(diagnostics_list, ~{
  if (is.null(.x)) {
    return(tibble(
      model_id = .y,
      overdispersion_ratio = NA_real_,
      overdispersion_p = NA_real_,
      zeroinflation_ratio = NA_real_,
      zeroinflation_p = NA_real_,
      singular = NA
    ))
  }
  
  tibble(
    model_id = .y,
    overdispersion_ratio = tryCatch(.x$overdispersion$dispersion_ratio, error = function(e) NA_real_),
    overdispersion_p = tryCatch(.x$overdispersion$p.value, error = function(e) NA_real_),
    zeroinflation_ratio = tryCatch(.x$zeroinflation$ratio, error = function(e) NA_real_),
    zeroinflation_p = tryCatch(.x$zero_test$p.value, error = function(e) NA_real_),
    singular = tryCatch(.x$singularity$singular, error = function(e) NA)
  )
})

print(diagnostics_summary)

readr::write_csv(
  diagnostics_summary,
  file.path(out_dir, "trapezia_alpheus_abundance_coupling_diagnostics_v_clean.csv")
)
