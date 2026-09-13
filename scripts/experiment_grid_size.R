# =============================================================================
# Appendix experiment: sensitivity of CostNP+ and CostNP to the cost grid step size
# Includes code for generating corresponding latex table
# =============================================================================

library(costnp)
source("scripts/paper_utils.R")

out_dir   <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
file_name <- paste0("experiment_grid_size_", format(Sys.time(), "%Y%m%d_%H%M"), ".RData")
out_path  <- file.path(out_dir, file_name)

library(dplyr)
library(tidyr)
library(furrr)
library(future)
plan(multisession, workers = parallel::detectCores() - 1)

# =============================================================================
# EXPERIMENT CONFIGURATION
# =============================================================================

SAMPLE_SIZE     <- 1000
DATA_MODELS     <- c("gaussian", "planet_2")
N_REPLICATIONS  <- 100
DIMENSIONS      <- 20
CLASS_BALANCE   <- 0.5
POPULATION_SIZE <- 1000000
TARGET_ALPHA    <- 0.1
TARGET_DELTA    <- 0.1

COST_START <- 0.99
COST_END   <- 0.01
STEP_SIZES <- c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2)

# =============================================================================
# RUN SIMULATION
# =============================================================================

scenario_groups <- expand.grid(
  cost_by = STEP_SIZES,
  model   = DATA_MODELS,
  stringsAsFactors = FALSE
)

costnp_plus_results_list <- list()
costnp_results_list      <- list()

for (g in seq_len(nrow(scenario_groups))) {
  g_by    <- scenario_groups$cost_by[g]
  g_model <- scenario_groups$model[g]

  cat(sprintf("\n--- Running: model=%s, cost_by=%.3f (%d/%d) ---\n",
              g_model, g_by, g, nrow(scenario_groups)))

  rep_results <- future_map(seq_len(N_REPLICATIONS), function(rep) {
    data <- gen_data(model = g_model, n = SAMPLE_SIZE, d = DIMENSIONS, pi = CLASS_BALANCE)
    pop  <- gen_data(model = g_model, n = POPULATION_SIZE, d = DIMENSIONS, pi = 0.5)

    # costnp_plus (interpolation)
    res_costnp_plus <- costnp_plus(
      data, pop, alpha = TARGET_ALPHA, delta = TARGET_DELTA, method = "LR",
      cost_start = COST_START, cost_end = COST_END, cost_by = g_by
    )
    costnp_plus_row <- data.frame(cost   = res_costnp_plus$cost,
                                  pop_t1 = res_costnp_plus$population_type1,
                                  pop_t2 = res_costnp_plus$population_type2)

    # costnp (conservative, no interpolation)
    res_costnp <- costnp(
      data, pop, alpha = TARGET_ALPHA, delta = TARGET_DELTA, method = "LR",
      cost_start = COST_START, cost_end = COST_END, cost_by = g_by
    )
    costnp_row <- data.frame(cost   = res_costnp$cost,
                             pop_t1 = res_costnp$population_type1,
                             pop_t2 = res_costnp$population_type2)

    list(costnp_plus = costnp_plus_row, costnp = costnp_row)
  }, .options = furrr_options(seed = TRUE), .progress = TRUE)

  add_meta <- function(df, model, by) { df$model <- model; df$cost_by <- by; df }

  costnp_plus_results_list[[g]] <- add_meta(bind_rows(lapply(rep_results, `[[`, "costnp_plus")), g_model, g_by)
  costnp_results_list[[g]]      <- add_meta(bind_rows(lapply(rep_results, `[[`, "costnp")),      g_model, g_by)

  cat(sprintf("  CostNP+  — violation rate: %.1f%%, avg type2: %.3f\n",
              mean(costnp_plus_results_list[[g]]$pop_t1 > TARGET_ALPHA) * 100,
              mean(costnp_plus_results_list[[g]]$pop_t2)))
  cat(sprintf("  CostNP   — violation rate: %.1f%%, avg type2: %.3f\n",
              mean(costnp_results_list[[g]]$pop_t1      > TARGET_ALPHA) * 100,
              mean(costnp_results_list[[g]]$pop_t2)))
}

costnp_plus_results <- bind_rows(costnp_plus_results_list)
costnp_results      <- bind_rows(costnp_results_list)

plan(sequential)

save(costnp_plus_results, costnp_results, file = out_path)
message("Success! Results saved to: ", out_path)

# =============================================================================
# LaTeX TABLE
# =============================================================================
# Move/rename the chosen run from output/ to data/experiment_grid_size.RData first.

results_env <- new.env()
load("data/experiment_grid_size.RData", envir = results_env)

MODEL_LABELS <- c(gaussian = "Tri-diagonal", planet_2 = "Elliptical")

summarise_method <- function(df, method_name) {
  df %>%
    group_by(model, cost_by) %>%
    summarise(
      viol    = mean(pop_t1 > TARGET_ALPHA) * 100,
      mean_t2 = mean(pop_t2) * 100,
      .groups = "drop"
    ) %>%
    mutate(method = method_name)
}

summary_wide <- bind_rows(
  summarise_method(results_env$costnp_plus_results, "costnp_plus"),
  summarise_method(results_env$costnp_results,      "costnp")
) %>%
  pivot_wider(id_cols = c(model, cost_by), names_from = method,
              values_from = c(viol, mean_t2)) %>%
  arrange(cost_by)

model_blocks <- lapply(intersect(names(MODEL_LABELS), unique(summary_wide$model)), function(m) {
  rows <- summary_wide %>% filter(model == m)
  c(sprintf("\\multicolumn{6}{l}{\\textbf{%s}} \\\\", MODEL_LABELS[[m]]),
    "\\midrule",
    sprintf("& %.3f & %.1f & %.1f & %.1f & %.1f \\\\",
            rows$cost_by, rows$viol_costnp_plus, rows$viol_costnp,
            rows$mean_t2_costnp_plus, rows$mean_t2_costnp))
})
body_rows <- Reduce(function(a, b) c(a, "\\midrule", b), model_blocks)

latex_lines <- c(
  "% Requires \\usepackage{booktabs} in preamble.",
  "\\begin{table}[ht]",
  "\\centering",
  "\\begin{tabular}{ll cc cc}",
  "\\toprule",
  "Model & Step & \\multicolumn{2}{c}{Type I Violation Rate (\\%)} & \\multicolumn{2}{c}{Avg.\\ Type II Error (\\%)} \\\\",
  "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6}",
  " & & CostNP+ & CostNP & CostNP+ & CostNP \\\\",
  "\\midrule",
  body_rows,
  "\\bottomrule",
  "\\end{tabular}",
  paste0("\\caption{Type I violation rates and average type II error rates for CostNP+ and CostNP ",
         "across data models and cost grid step sizes. All entries are percentages.}"),
  "\\label{Tab::ExperimentAdditional}",
  "\\end{table}"
)

cat("\n\n=== LaTeX Table ===\n\n")
cat(paste(latex_lines, collapse = "\n"), "\n")
