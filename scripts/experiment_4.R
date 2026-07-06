# =============================================================================
# Experiment 4: Comparing CostNP and NP under equal class balance
# =============================================================================

library(costnp)
library(dplyr)
library(furrr)
library(future)
library(nproc)
source("scripts/paper_utils.R")

out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

plan(multisession, workers = parallel::detectCores() - 1)

# =============================================================================
# SHARED EXPERIMENT CONFIGURATION
# =============================================================================

SAMPLE_SIZES    <- c(1000, 5000)
DATA_MODELS     <- c("gaussian", "t", "mixture")
N_REPLICATIONS  <- 500
DIMENSIONS      <- 20
CLASS_BALANCE   <- 0.5
POPULATION_SIZE <- 1000000
TARGET_ALPHA    <- 0.1
TARGET_DELTA    <- 0.1

# Which classifiers to run (subset to re-run a single one).
METHODS_TO_RUN  <- c("LDA", "LR", "NB")

# --- Cost grids -------------------------------------------------------------
# GB and NB share the "default" grids; LR uses its own.
COST_GRID_DEFAULT <- list(
  gaussian = list(
    "1000" = list(start = 0.99, end = 0.2, by = 0.02),
    "5000" = list(start = 0.95, end = 0.2, by = 0.01)
  ),
  t = list(
    "1000" = list(start = 0.99, end = 0.2, by = 0.02),
    "5000" = list(start = 0.95, end = 0.2, by = 0.01)
  ),
  mixture = list(
    "1000" = list(start = 0.99, end = 0.2, by = 0.02),
    "5000" = list(start = 0.95, end = 0.2, by = 0.01)
  )
)

NAIVE_COST_GRID_DEFAULT <- list(
  gaussian = list(
    "1000" = list(start = 0.45, end = 0.995, by = 0.02),
    "5000" = list(start = 0.45, end = 0.9,   by = 0.01)
  ),
  t = list(
    "1000" = list(start = 0.15, end = 0.995, by = 0.02),
    "5000" = list(start = 0.25, end = 0.9,   by = 0.01)
  ),
  mixture = list(
    "1000" = list(start = 0.5,  end = 0.995, by = 0.02),
    "5000" = list(start = 0.5,  end = 0.9,   by = 0.01)
  )
)

# --- Per-classifier configuration -------------------------------------------
# np_strategy: "manual" => npc_manual with this classifier; otherwise the string
# is passed as nproc::npc(method = ...).
METHOD_CONFIGS <- list(
  LR = list(method = "LR", np_strategy = "logistic",   threshold_fun = function(cost) 0.5,
            cost_grid = COST_GRID_DEFAULT, naive_cost_grid = NAIVE_COST_GRID_DEFAULT),
  LDA = list(method = "LDA", np_strategy = "lda", threshold_fun = function(cost) 0.5,
            cost_grid = COST_GRID_DEFAULT,      naive_cost_grid = NAIVE_COST_GRID_DEFAULT),
  NB = list(method = "NB", np_strategy = "nb",       threshold_fun = function(cost) 0.5,
            cost_grid = COST_GRID_DEFAULT, naive_cost_grid = NAIVE_COST_GRID_DEFAULT)
)

# =============================================================================
# RUN SIMULATION
# =============================================================================

# Compute the NPROC comparison row for one replication.
compute_np_row <- function(cfg, data, pop) {
  if (cfg$np_strategy == "manual") {
    res <- npc_manual(data = data, alpha = TARGET_ALPHA, delta = TARGET_DELTA,
                      population_data = pop, method = cfg$method)
    data.frame(cost = 0.5, pop_t1 = res$population_type1, pop_t2 = res$population_type2)
  } else {
    fit <- npc(x = data$x, y = data$y, alpha = TARGET_ALPHA,
               delta = TARGET_DELTA, method = cfg$np_strategy)
    data.frame(
      cost   = 0.5,
      pop_t1 = sum(predict(fit, newx = pop$x[pop$y == 0, ])$pred.label == 1) / sum(pop$y == 0),
      pop_t2 = sum(predict(fit, newx = pop$x[pop$y == 1, ])$pred.label == 0) / sum(pop$y == 1)
    )
  }
}

run_experiment <- function(method_name, cfg) {
  cat(sprintf("\n================================================================\n"))
  cat(sprintf("  Classifier: %s\n", method_name))
  cat(sprintf("================================================================\n"))
  
  scenario_groups <- expand.grid(
    n     = SAMPLE_SIZES,
    model = DATA_MODELS,
    stringsAsFactors = FALSE
  )
  
  costnp_plus_results_list <- list()
  costnp_results_list      <- list()
  naive_results_list       <- list()
  np_results_list          <- list()
  
  for (g in seq_len(nrow(scenario_groups))) {
    g_n     <- scenario_groups$n[g]
    g_model <- scenario_groups$model[g]
    cost_settings       <- cfg$cost_grid[[g_model]][[as.character(g_n)]]
    naive_cost_settings <- cfg$naive_cost_grid[[g_model]][[as.character(g_n)]]
    
    cat(sprintf("\n--- Running: model=%s, n=%d (%d/%d) ---\n",
                g_model, g_n, g, nrow(scenario_groups)))
    
    rep_results <- future_map(seq_len(N_REPLICATIONS), function(rep) {
      data <- gen_data(model = g_model, n = g_n, d = DIMENSIONS, pi = CLASS_BALANCE)
      pop  <- gen_data(model = g_model, n = POPULATION_SIZE, d = DIMENSIONS, pi = 0.5)
      
      # costnp_plus (interpolation)
      res_costnp_plus <- costnp_plus(
        data, pop, alpha = TARGET_ALPHA, method = cfg$method,
        cost_start = cost_settings$start, cost_end = cost_settings$end,
        cost_by = cost_settings$by, threshold_fun = cfg$threshold_fun
      )
      costnp_plus_row <- data.frame(cost   = res_costnp_plus$cost,
                                    pop_t1 = res_costnp_plus$population_type1,
                                    pop_t2 = res_costnp_plus$population_type2)
      
      # costnp (conservative, no interpolation)
      res_costnp <- costnp(
        data, pop, alpha = TARGET_ALPHA, method = cfg$method,
        cost_start = cost_settings$start, cost_end = cost_settings$end,
        cost_by = cost_settings$by, threshold_fun = cfg$threshold_fun
      )
      costnp_row <- data.frame(cost   = res_costnp$cost,
                               pop_t1 = res_costnp$population_type1,
                               pop_t2 = res_costnp$population_type2)
      
      # Naive
      res_naive <- get_cost_naive(
        data, alpha = TARGET_ALPHA, population_data = pop,
        cost_start = naive_cost_settings$start, cost_end = naive_cost_settings$end,
        cost_by = naive_cost_settings$by, method = cfg$method
      )
      naive_row <- data.frame(cost   = res_naive$cost,
                              pop_t1 = res_naive$population_type1,
                              pop_t2 = res_naive$population_type2)
      
      # NPROC
      np_row <- compute_np_row(cfg, data, pop)
      
      list(costnp_plus = costnp_plus_row, costnp = costnp_row, naive = naive_row, np = np_row)
    }, .options = furrr_options(seed = TRUE), .progress = TRUE)
    
    add_meta <- function(df, model, n) { df$model <- model; df$n <- n; df }
    
    costnp_plus_results_list[[g]] <- add_meta(bind_rows(lapply(rep_results, `[[`, "costnp_plus")), g_model, g_n)
    costnp_results_list[[g]]      <- add_meta(bind_rows(lapply(rep_results, `[[`, "costnp")),      g_model, g_n)
    naive_results_list[[g]]       <- add_meta(bind_rows(lapply(rep_results, `[[`, "naive")),       g_model, g_n)
    np_results_list[[g]]          <- add_meta(bind_rows(lapply(rep_results, `[[`, "np")),          g_model, g_n)
    
    cat(sprintf("  CostNP+  — violation rate: %.1f%%, avg type2: %.3f\n",
                mean(costnp_plus_results_list[[g]]$pop_t1 > TARGET_ALPHA) * 100,
                mean(costnp_plus_results_list[[g]]$pop_t2)))
    cat(sprintf("  CostNP   — violation rate: %.1f%%, avg type2: %.3f\n",
                mean(costnp_results_list[[g]]$pop_t1      > TARGET_ALPHA) * 100,
                mean(costnp_results_list[[g]]$pop_t2)))
    cat(sprintf("  Naive    — violation rate: %.1f%%, avg type2: %.3f\n",
                mean(naive_results_list[[g]]$pop_t1  > TARGET_ALPHA) * 100,
                mean(naive_results_list[[g]]$pop_t2)))
    cat(sprintf("  NPROC    — violation rate: %.1f%%, avg type2: %.3f\n",
                mean(np_results_list[[g]]$pop_t1     > TARGET_ALPHA) * 100,
                mean(np_results_list[[g]]$pop_t2)))
  }
  
  costnp_plus_results <- bind_rows(costnp_plus_results_list)
  costnp_results      <- bind_rows(costnp_results_list)
  naive_results       <- bind_rows(naive_results_list)
  np_results          <- bind_rows(np_results_list)
  
  file_name <- paste0("experiment4_", tolower(method_name), "_",
                      format(Sys.time(), "%Y%m%d_%H%M"), ".RData")
  out_path  <- file.path(out_dir, file_name)
  save(costnp_plus_results, costnp_results, naive_results, np_results, file = out_path)
  message("Success! ", method_name, " results saved to: ", out_path)
}

for (method_name in METHODS_TO_RUN) {
  run_experiment(method_name, METHOD_CONFIGS[[method_name]])
}

plan(sequential)
