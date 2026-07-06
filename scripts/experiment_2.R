# =============================================================================
# Experiment 2: Comparing CostNP, NP and eLDA under class imbalance
# Includes code for generating corresponding latex table
# =============================================================================

library(costnp)
source("scripts/paper_utils.R")

out_dir   <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
file_name <- paste0("experiment2_", format(Sys.time(), "%Y%m%d_%H%M"), ".RData")
out_path  <- file.path(out_dir, file_name)

library(dplyr)
library(furrr)
library(future)
library(nproc)
plan(multisession, workers = parallel::detectCores() - 1)

# =============================================================================
# EXPERIMENT CONFIGURATION
# =============================================================================

SAMPLE_SIZES    <- c(1000, 5000)
DATA_MODELS     <- c("homgaussian", "gaussian", "planet_2")
N_REPLICATIONS  <- 500
DIMENSIONS      <- 20
CLASS_BALANCE   <- 0.2
POPULATION_SIZE <- 1000000
TARGET_ALPHA    <- 0.1
TARGET_DELTA    <- 0.1

COST_GRID <- list(
  gaussian = list(
    "1000" = list(start = 0.99, end = 0.5, by = 0.02),
    "5000" = list(start = 0.99, end = 0.5, by = 0.01)
  ),
  homgaussian = list(
    "1000" = list(start = 0.99, end = 0.5, by = 0.02),
    "5000" = list(start = 0.99, end = 0.5, by = 0.01)
  ),
  planet_2 = list(
    "1000" = list(start = 0.99, end = 0.5, by = 0.02),
    "5000" = list(start = 0.99, end = 0.5, by = 0.01)
  )
)

# =============================================================================
# RUN SIMULATION
# =============================================================================

scenario_groups <- expand.grid(
  n     = SAMPLE_SIZES,
  model = DATA_MODELS,
  stringsAsFactors = FALSE
)

costnp_plus_results_list <- list()
costnp_results_list      <- list()
naive_results_list  <- list()
np_results_list     <- list()
elda_results_list   <- list()

for (g in seq_len(nrow(scenario_groups))) {
  g_n     <- scenario_groups$n[g]
  g_model <- scenario_groups$model[g]
  cost_settings <- COST_GRID[[g_model]][[as.character(g_n)]]

  cat(sprintf("\n--- Running: model=%s, n=%d (%d/%d) ---\n",
              g_model, g_n, g, nrow(scenario_groups)))

  rep_results <- future_map(seq_len(N_REPLICATIONS), function(rep) {
    data <- gen_data(model = g_model, n = g_n, d = DIMENSIONS, pi = CLASS_BALANCE)
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

    # NPROC
    nproc_fit  <- npc(x = data$x, y = data$y, alpha = TARGET_ALPHA,
                      delta = TARGET_DELTA, method = "logistic")
    np_row <- data.frame(
      cost   = 0.5,
      pop_t1 = sum(predict(nproc_fit, newx = pop$x[pop$y == 0, ])$pred.label == 1) /
               sum(pop$y == 0),
      pop_t2 = sum(predict(nproc_fit, newx = pop$x[pop$y == 1, ])$pred.label == 0) /
               sum(pop$y == 1)
    )

    # eLDA
    res_elda <- elda_npc(data, alpha = TARGET_ALPHA, population_data = pop,
                         delta = TARGET_DELTA, variant = "eLDA")
    elda_row <- data.frame(cost      = NA_real_,
                           pop_t1    = res_elda$population_type1,
                           pop_t2    = res_elda$population_type2,
                           threshold = res_elda$threshold)

    list(costnp_plus = costnp_plus_row, costnp = costnp_row, naive = naive_row, np = np_row, elda = elda_row)
  }, .options = furrr_options(seed = TRUE), .progress = TRUE)

  add_meta <- function(df, model, n) { df$model <- model; df$n <- n; df }

  costnp_plus_results_list[[g]] <- add_meta(bind_rows(lapply(rep_results, `[[`, "costnp_plus")), g_model, g_n)
  costnp_results_list[[g]]      <- add_meta(bind_rows(lapply(rep_results, `[[`, "costnp")),      g_model, g_n)
  naive_results_list[[g]]  <- add_meta(bind_rows(lapply(rep_results, `[[`, "naive")),  g_model, g_n)
  np_results_list[[g]]     <- add_meta(bind_rows(lapply(rep_results, `[[`, "np")),     g_model, g_n)
  elda_results_list[[g]]   <- add_meta(bind_rows(lapply(rep_results, `[[`, "elda")),   g_model, g_n)

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
  cat(sprintf("  eLDA     — violation rate: %.1f%%, avg type2: %.3f\n",
              mean(elda_results_list[[g]]$pop_t1   > TARGET_ALPHA) * 100,
              mean(elda_results_list[[g]]$pop_t2)))
}

costnp_plus_results <- bind_rows(costnp_plus_results_list)
costnp_results      <- bind_rows(costnp_results_list)
naive_results  <- bind_rows(naive_results_list)
np_results     <- bind_rows(np_results_list)
elda_results   <- bind_rows(elda_results_list)

plan(sequential)

save(costnp_plus_results, costnp_results, naive_results, np_results, elda_results, file = out_path)
message("Success! Results saved to: ", out_path)

# =============================================================================
# LaTeX TABLE
# =============================================================================

build_summary <- function(df, method_name) {
  df %>%
    group_by(model, n) %>%
    summarise(
      viol_rate = mean(pop_t1 > TARGET_ALPHA) * 100,
      mean_t2   = mean(pop_t2) * 100,
      sd_t2     = sd(pop_t2) * 100,
      .groups   = "drop"
    ) %>%
    mutate(method = method_name)
}

summary_all <- bind_rows(list(
  build_summary(naive_results,  "Naive"),
  build_summary(np_results,     "NP"),
  build_summary(costnp_plus_results, "CostNP+"),
  build_summary(costnp_results,      "CostNP"),
  build_summary(elda_results,   "eLDA")
))

model_labels <- c(homgaussian = "Shared", gaussian = "Tri-diagonal", planet_2 = "Elliptical")
method_order <- c("Naive", "NP", "CostNP+", "CostNP", "eLDA")
model_order  <- c("homgaussian", "gaussian", "planet_2")
n_order      <- c(1000, 5000)
VIOLATION_SLACK <- 2*sqrt(TARGET_ALPHA * (1 - TARGET_ALPHA) / N_REPLICATIONS)
viol_threshold  <- (TARGET_ALPHA + VIOLATION_SLACK) * 100

header <- paste0(
  "\\begin{table}[ht]\n",
  "\\centering\n",
  "\\setlength{\\tabcolsep}{5pt}\n",
  "\\begin{tabular}{llccccc}\n",
  "\\hline\n",
  "Model & $n$ & Naive & NP & CostNP+ & CostNP & eLDA \\\\\n",
  "\\hline\n"
)

body <- "\\multicolumn{7}{l}{\\textit{Type I violation rate (\\%)}} \\\\\n"

for (mod in model_order) {
  for (i_n in seq_along(n_order)) {
    n_val     <- n_order[i_n]
    mod_label <- if (i_n == 1) model_labels[mod] else ""
    row_data  <- summary_all %>% filter(model == mod, n == n_val)

    t1_vals <- sapply(method_order, function(m) {
      v   <- row_data$viol_rate[row_data$method == m]
      val <- sprintf("%.1f", v)
      if (v > viol_threshold) val <- paste0(val, "$^*$")
      val
    })

    body <- paste0(body, sprintf(
      "%s & %d & %s \\\\\n",
      mod_label, n_val, paste(t1_vals, collapse = " & ")
    ))
  }
}

body <- paste0(body, "\\hline\n",
               "\\multicolumn{7}{l}{\\textit{Average type II error (\\%)}} \\\\\n")

for (mod in model_order) {
  for (i_n in seq_along(n_order)) {
    n_val     <- n_order[i_n]
    mod_label <- if (i_n == 1) model_labels[mod] else ""
    row_data  <- summary_all %>% filter(model == mod, n == n_val)

    means   <- sapply(method_order, function(m) row_data$mean_t2[row_data$method == m])
    sds     <- sapply(method_order, function(m) row_data$sd_t2[row_data$method == m])
    viols   <- sapply(method_order, function(m) row_data$viol_rate[row_data$method == m])
    starred <- viols > viol_threshold

    # Lowest average type II error among non-asterisked (valid) entries
    best_idx <- if (any(!starred)) which(!starred)[which.min(means[!starred])] else NA_integer_

    t2_vals <- sapply(seq_along(method_order), function(j) {
      val <- sprintf("%.1f (%.1f)", means[j], sds[j])
      if (starred[j]) val <- paste0(val, "$^*$")
      if (!is.na(best_idx) && j == best_idx) val <- paste0("\\textbf{", val, "}")
      val
    })

    body <- paste0(body, sprintf(
      "%s & %d & %s \\\\\n",
      mod_label, n_val, paste(t2_vals, collapse = " & ")
    ))
  }
}

footer <- paste0(
  "\\hline\n",
  "\\end{tabular}\n",
  "\\caption{Average type I violation rates and mean (sd) Type II error rates (in percentage) across classifiers, ",
  "models, and sample sizes. An asterisk indicates that the type I violation rate exceeds ",
  "the target violation level of $10\\%$ (within Monte Carlo tolerance). The lowest ",
  "average type II error in each row, among non-asterisked entries, is shown in bold.}\n",
  "\\label{Tab::Experiment2Table}\n",
  "\\end{table}"
)

cat(paste0(header, body, footer))
