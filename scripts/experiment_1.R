# =============================================================================
# Experiment 1: Comparing CostNP and NP under varying class imbalance
# =============================================================================

library(costnp)
source("scripts/paper_utils.R")

out_dir   <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
file_name <- paste0("experiment1_", format(Sys.time(), "%Y%m%d_%H%M"), ".RData")
out_path  <- file.path(out_dir, file_name)

library(dplyr)
library(furrr)
library(future)
library(nproc)
plan(multisession, workers = parallel::detectCores() - 1)

# =============================================================================
# EXPERIMENT CONFIGURATION
# =============================================================================

SAMPLE_SIZES    <- c(1000)
DATA_MODELS     <- c("planet")
N_REPLICATIONS  <- 500
CLASS_BALANCES  <- seq(0.1, 0.9, by = 0.1)
POPULATION_SIZE <- 1000000
TARGET_ALPHA    <- 0.1
TARGET_DELTA    <- 0.1

COST_GRID <- list(
  planet = list(
    "1000"  = list(start = 0.99, end = 0.01,  by = 0.02)
  )
)

# =============================================================================
# RUN SIMULATION
# =============================================================================

set.seed(1357)

scenario_groups <- expand.grid(
  n     = SAMPLE_SIZES,
  model = DATA_MODELS,
  pi    = CLASS_BALANCES,
  stringsAsFactors = FALSE
)

costnp_plus_results_list <- list()
costnp_results_list       <- list()
naive_results_list     <- list()
np_results_list        <- list()
np_manual_results_list <- list()

for (g in seq_len(nrow(scenario_groups))) {
  g_n     <- scenario_groups$n[g]
  g_model <- scenario_groups$model[g]
  g_pi    <- scenario_groups$pi[g]
  cost_settings <- COST_GRID[[g_model]][[as.character(g_n)]]

  cat(sprintf("\n--- Running: model=%s, n=%d, pi=%.2f (%d/%d) ---\n",
              g_model, g_n, g_pi, g, nrow(scenario_groups)))

  rep_results <- future_map(seq_len(N_REPLICATIONS), function(rep) {
    data <- gen_data(model = g_model, n = g_n, d = DIMENSIONS, pi = g_pi)
    pop  <- gen_data(model = g_model, n = POPULATION_SIZE, d = DIMENSIONS, pi = 0.5)

    # costnp_plus (interpolation)
    res_costnp_plus <- costnp_plus(
      data, pop, alpha = TARGET_ALPHA, method = "LR",
      cost_start = cost_settings$start, cost_end = cost_settings$end,
      cost_by = cost_settings$by
    )
    costnp_plus_row <- data.frame(cost   = res_costnp_plus$cost,
                                  pop_t1 = res_costnp_plus$population_type1,
                                  pop_t2 = res_costnp_plus$population_type2)

    # costnp (conservative, no interpolation)
    res_costnp <- costnp(
      data, pop, alpha = TARGET_ALPHA, method = "LR",
      cost_start = cost_settings$start, cost_end = cost_settings$end,
      cost_by = cost_settings$by
    )
    costnp_row <- data.frame(cost   = res_costnp$cost,
                            pop_t1 = res_costnp$population_type1,
                            pop_t2 = res_costnp$population_type2)

    # Naive
    res_naive <- get_cost_naive(
      data, alpha = TARGET_ALPHA, population_data = pop,
      cost_start = cost_settings$end, cost_end = cost_settings$start,
      cost_by = cost_settings$by
    )
    naive_row <- data.frame(cost   = res_naive$cost,
                            pop_t1 = res_naive$population_type1,
                            pop_t2 = res_naive$population_type2)

    # NPROC_logistic
    nproc_fit <- npc(x = data$x, y = data$y, alpha = TARGET_ALPHA,
                     delta = TARGET_DELTA, method = "logistic")
    np_row <- data.frame(
      cost   = 0.5,
      pop_t1 = sum(predict(nproc_fit, newx = pop$x[pop$y == 0, ])$pred.label == 1) /
               sum(pop$y == 0),
      pop_t2 = sum(predict(nproc_fit, newx = pop$x[pop$y == 1, ])$pred.label == 0) /
               sum(pop$y == 1)
    )

    list(costnp_plus = costnp_plus_row, costnp = costnp_row, naive = naive_row,
         np = np_row)
  }, .options = furrr_options(seed = TRUE), .progress = TRUE)

  add_meta <- function(df, model, n, pi) { df$model <- model; df$n <- n; df$pi <- pi; df }

  costnp_plus_results_list[[g]] <- add_meta(bind_rows(lapply(rep_results, `[[`, "costnp_plus")), g_model, g_n, g_pi)
  costnp_results_list[[g]]       <- add_meta(bind_rows(lapply(rep_results, `[[`, "costnp")),       g_model, g_n, g_pi)
  naive_results_list[[g]]     <- add_meta(bind_rows(lapply(rep_results, `[[`, "naive")),     g_model, g_n, g_pi)
  np_results_list[[g]]        <- add_meta(bind_rows(lapply(rep_results, `[[`, "np")),        g_model, g_n, g_pi)

  cat(sprintf("  costnp_plus — violation rate: %.1f%%, avg type2: %.3f\n",
              mean(costnp_plus_results_list[[g]]$pop_t1 > TARGET_ALPHA) * 100,
              mean(costnp_plus_results_list[[g]]$pop_t2)))
  cat(sprintf("  costnp       — violation rate: %.1f%%, avg type2: %.3f\n",
              mean(costnp_results_list[[g]]$pop_t1       > TARGET_ALPHA) * 100,
              mean(costnp_results_list[[g]]$pop_t2)))
  cat(sprintf("  Naive     — violation rate: %.1f%%, avg type2: %.3f\n",
              mean(naive_results_list[[g]]$pop_t1     > TARGET_ALPHA) * 100,
              mean(naive_results_list[[g]]$pop_t2)))
  cat(sprintf("  NPROC     — violation rate: %.1f%%, avg type2: %.3f\n",
              mean(np_results_list[[g]]$pop_t1        > TARGET_ALPHA) * 100,
              mean(np_results_list[[g]]$pop_t2)))
}

costnp_plus_results <- bind_rows(costnp_plus_results_list)
costnp_results       <- bind_rows(costnp_results_list)
naive_results     <- bind_rows(naive_results_list)
np_results        <- bind_rows(np_results_list)

plan(sequential)

save(costnp_plus_results, costnp_results, naive_results, np_results,
     file = out_path)
message("Success! Results saved to: ", out_path)
