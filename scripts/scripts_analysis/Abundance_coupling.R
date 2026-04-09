## ============================================================
## Unified Trapezia–Alpheus abundance coupling
##
## Goal:
## Test whether Alpheus abundance scales with Trapezia abundance
## separately for each host species using the unified control series
##
## Outputs:
## - model coefficients
## - DHARMa / performance diagnostics
## - colony-level mean correlations
## - optional prediction plots
## ============================================================

library(tidyverse)
library(janitor)
library(stringr)
library(glmmTMB)
library(broom.mixed)
library(DHARMa)
library(performance)
library(ggplot2)
library(ggtext)

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
repo_dir <- "/Users/amalia/Documents/GitHub/Trapezia_Alpheus_Facilitation_For_PR_rev"
input_file <- file.path(repo_dir, "input", "controls_cleaned_with_months_unified.csv")
out_dir_fig <- file.path(repo_dir, "outputs_unified", "figures")
out_dir_tab <- file.path(repo_dir, "outputs_unified", "output_files")

dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tab, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## 1. Import unified data
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
    date = as.Date(date, format = "%Y-%m-%d")
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
## 3. Build cleaned unified analysis dataset
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
    coral_uid = factor(coral_uid),
    month_num = as.numeric(months_of_monitoring),
    survey_occasion = as.numeric(survey_occasion),
    trapezia = replace_na(as.numeric(trapezia), 0),
    alpheus  = replace_na(as.numeric(alpheus), 0)
  )

readr::write_csv(
  dat,
  file.path(out_dir_tab, "abundance_coupling_unified_analysis_dataset.csv")
)

message("Counts by host species:")
print(
  dat %>%
    count(host_species, name = "n_obs") %>%
    left_join(
      dat %>%
        group_by(host_species) %>%
        summarise(n_colonies = n_distinct(coral_uid), .groups = "drop"),
      by = "host_species"
    )
)

## ------------------------------------------------------------
## 4. Fit one NB-GLMM per host species
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
    error = function(e) {
      message("Model fit failed: ", conditionMessage(e))
      NULL
    }
  )
}

results_tbl <- dat %>%
  group_by(genus, host_species, host_species_md) %>%
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
  file.path(out_dir_tab, "abundance_coupling_unified_results.csv")
)

## ------------------------------------------------------------
## 5. Save model objects for diagnostics / summaries
## ------------------------------------------------------------
group_keys_tbl <- dat %>%
  group_by(genus, host_species) %>%
  group_keys() %>%
  mutate(label = paste(genus, sep = "_"))

model_list <- dat %>%
  group_by(genus, host_species) %>%
  group_split() %>%
  setNames(group_keys_tbl$label) %>%
  purrr::map(fit_nb_group)

## ------------------------------------------------------------
## 6. Model diagnostics
## ------------------------------------------------------------
extract_diagnostics <- function(model, label) {
  if (is.null(model)) {
    return(tibble(
      label = label,
      overdispersion_ratio = NA_real_,
      overdispersion_p = NA_real_,
      zero_inflation_ratio = NA_real_,
      zero_inflation_p = NA_real_,
      singular = NA,
      dharma_uniformity_p = NA_real_,
      dharma_dispersion_p = NA_real_,
      dharma_zero_p = NA_real_
    ))
  }
  
  res <- DHARMa::simulateResiduals(model, n = 1000)
  
  od <- tryCatch(performance::check_overdispersion(model), error = function(e) NULL)
  zi <- tryCatch(performance::check_zeroinflation(model), error = function(e) NULL)
  sg <- tryCatch(performance::check_singularity(model), error = function(e) NULL)
  
  tu <- tryCatch(DHARMa::testUniformity(res), error = function(e) NULL)
  td <- tryCatch(DHARMa::testDispersion(res), error = function(e) NULL)
  tz <- tryCatch(DHARMa::testZeroInflation(res), error = function(e) NULL)
  
  tibble(
    label = label,
    overdispersion_ratio = if (!is.null(od) && "dispersion_ratio" %in% names(od)) od$dispersion_ratio else NA_real_,
    overdispersion_p = if (!is.null(od) && "p_value" %in% names(od)) od$p_value else NA_real_,
    zero_inflation_ratio = if (!is.null(zi) && "ratio" %in% names(zi)) zi$ratio else NA_real_,
    zero_inflation_p = if (!is.null(zi) && "p_value" %in% names(zi)) zi$p_value else NA_real_,
    singular = if (!is.null(sg) && "singular" %in% names(sg)) sg$singular else NA,
    dharma_uniformity_p = if (!is.null(tu)) tu$p.value else NA_real_,
    dharma_dispersion_p = if (!is.null(td)) td$p.value else NA_real_,
    dharma_zero_p = if (!is.null(tz)) tz$p.value else NA_real_
  )
}

diagnostics_tbl <- purrr::imap_dfr(model_list, extract_diagnostics) %>%
  separate(label, into = c("genus"), sep = "_", remove = FALSE) %>%
  mutate(
    host_species = dplyr::case_match(
      genus,
      "Pocillopora" ~ "Pocillopora favosa",
      "Stylophora"  ~ "Stylophora pistillata"
    )
  ) %>%
  select(host_species, everything(), -label, -genus)

print(diagnostics_tbl)

readr::write_csv(
  diagnostics_tbl,
  file.path(out_dir_tab, "abundance_coupling_unified_diagnostics.csv")
)

## ------------------------------------------------------------
## 7. Print full model summaries and diagnostic tests
## ------------------------------------------------------------
message("Model summaries:")
for (nm in names(model_list)) {
  message("---- ", nm, " ----")
  if (is.null(model_list[[nm]])) {
    message("Model not fit")
  } else {
    print(summary(model_list[[nm]]))
  }
}

message("Diagnostics summary:")
print(diagnostics_tbl)

## ------------------------------------------------------------
## 8. Save DHARMa diagnostic plots
## ------------------------------------------------------------
for (nm in names(model_list)) {
  model <- model_list[[nm]]
  if (is.null(model)) next
  
  png(
    filename = file.path(out_dir_fig, paste0("DHARMa_", nm, ".png")),
    width = 1800,
    height = 1800,
    res = 220
  )
  res <- DHARMa::simulateResiduals(model, n = 1000)
  plot(res)
  dev.off()
}

## ------------------------------------------------------------
## 9. Colony-level mean correlations
## ------------------------------------------------------------
col_means <- dat %>%
  group_by(genus, host_species, host_species_md, coral_uid) %>%
  summarise(
    mean_trap = mean(trapezia, na.rm = TRUE),
    mean_alph = mean(alpheus, na.rm = TRUE),
    n_observations = n(),
    .groups = "drop"
  )

safe_spearman <- function(x, y) {
  test <- suppressWarnings(
    cor.test(x, y, method = "spearman", exact = FALSE)
  )
  
  tibble(
    rho = unname(test$estimate),
    p = test$p.value
  )
}

col_cor <- col_means %>%
  group_by(genus, host_species, host_species_md) %>%
  group_modify(~{
    out <- safe_spearman(.x$mean_trap, .x$mean_alph)
    out %>%
      mutate(
        n_colonies = nrow(.x)
      )
  }) %>%
  ungroup()

print(col_cor)

readr::write_csv(
  col_means,
  file.path(out_dir_tab, "abundance_coupling_unified_colony_means.csv")
)

readr::write_csv(
  col_cor,
  file.path(out_dir_tab, "abundance_coupling_unified_colony_correlations.csv")
)

## ------------------------------------------------------------
## 10. Prediction plots by host species
## ------------------------------------------------------------
plot_prediction_group <- function(model, data_used, host_name) {
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
      host_species = host_name
    )
  
  ggplot(pred_df, aes(x = trapezia, y = predicted)) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.2) +
    geom_line(linewidth = 1) +
    theme_bw(base_size = 14) +
    labs(
      x = "Trapezia count",
      y = "Predicted Alpheus count",
      title = host_name
    )
}

prediction_plots <- list()

for (host_name in unique(dat$host_species)) {
  genus_name <- ifelse(host_name == "Pocillopora favosa", "Pocillopora", "Stylophora")
  label <- genus_name
  
  d_sub <- dat %>%
    filter(host_species == host_name)
  
  prediction_plots[[label]] <- plot_prediction_group(
    model = model_list[[label]],
    data_used = d_sub,
    host_name = host_name
  )
}

if (!is.null(prediction_plots$Pocillopora)) {
  print(prediction_plots$Pocillopora)
}
if (!is.null(prediction_plots$Stylophora)) {
  print(prediction_plots$Stylophora)
}

if (!is.null(prediction_plots$Pocillopora)) {
  ggsave(
    filename = file.path(out_dir_fig, "abundance_coupling_prediction_Pocillopora.png"),
    plot = prediction_plots$Pocillopora,
    width = 6,
    height = 4.5,
    dpi = 600
  )
}

if (!is.null(prediction_plots$Stylophora)) {
  ggsave(
    filename = file.path(out_dir_fig, "abundance_coupling_prediction_Stylophora.png"),
    plot = prediction_plots$Stylophora,
    width = 6,
    height = 4.5,
    dpi = 600
  )
}
