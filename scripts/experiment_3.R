# =============================================================================
# Experiment 3: Comparing CostNP and NP under many learning methods
# Includes code for generating corresponding latex table
# =============================================================================

library(costnp)
source("scripts/paper_utils.R")

out_dir   <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
file_name <- paste0("experiment3_test_", format(Sys.time(), "%Y%m%d_%H%M"), ".RData")
out_path  <- file.path(out_dir, file_name)

library(dplyr)
library(furrr)
library(future)
library(nproc)
plan(multisession, workers = parallel::detectCores() - 1)

# =============================================================================
# EXPERIMENT CONFIGURATION
# =============================================================================
set.seed(1)
SAMPLE_SIZES    <- c(1000)
DATA_MODELS     <- c("gaussian")
N_REPLICATIONS  <- 500
DIMENSIONS      <- 20
CLASS_BALANCE   <- 0.2
POPULATION_SIZE <- 100000
TARGET_ALPHA    <- 0.1
TARGET_DELTA    <- 0.1

METHODS <- c("QDA", "LDA", "SVM", "LR", "NB", "penLR")

# Maps our method names to nproc's method argument; NULL means use npc_manual
NPROC_METHOD_MAP <- list(
  QDA   = NULL,
  LDA   = "lda",
  SVM   = "svm",
  LR    = "logistic",
  NB    = "nb",
  penLR = "penlog"
)

COST_GRID <- list(
  gaussian = list(
    "1000" = list(start = 0.99, end = 0.5, by = 0.01)
  ),
  planet = list(
    "1000" = list(start = 0.99, end = 0.5, by = 0.01)
  ),
  homgaussian = list(
    "1000" = list(start = 0.99, end = 0.5, by = 0.01)
  )
)

NAIVE_COST_GRID <- list(
  gaussian = list(
    "1000" = list(start = 0.5, end = 0.995, by = 0.01)
  ),
  planet = list(
    "1000" = list(start = 0.5, end = 0.995, by = 0.01)
  ),
  homgaussian = list(
    "1000" = list(start = 0.5, end = 0.995, by = 0.01)
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

all_results <- list()

for (method in METHODS) {
  cat(sprintf("\n============================================================\n"))
  cat(sprintf("  METHOD: %s\n", method))
  cat(sprintf("============================================================\n"))

  nproc_method <- NPROC_METHOD_MAP[[method]]
  cfun <- if (method == "SVM") classify_fun_svm_weighted else classify_fun_stratified
  tfun <- function(cost) cost

  costnp_plus_results_list <- list()
  costnp_results_list      <- list()
  naive_results_list  <- list()
  np_results_list     <- list()

  for (g in seq_len(nrow(scenario_groups))) {
    g_n     <- scenario_groups$n[g]
    g_model <- scenario_groups$model[g]
    cost_settings       <- COST_GRID[[g_model]][[as.character(g_n)]]
    naive_cost_settings <- NAIVE_COST_GRID[[g_model]][[as.character(g_n)]]

    cat(sprintf("\n--- Running: model=%s, n=%d (%d/%d) ---\n",
                g_model, g_n, g, nrow(scenario_groups)))

    rep_results <- future_map(seq_len(N_REPLICATIONS), function(rep) {
      data <- gen_data(model = g_model, n = g_n, d = DIMENSIONS, pi = CLASS_BALANCE)
      pop  <- gen_data(model = g_model, n = POPULATION_SIZE, d = DIMENSIONS, pi = 0.5)

      # costnp_plus (interpolation)
      res_costnp_plus <- costnp_plus(
        data, pop, alpha = TARGET_ALPHA, method = method,
        cost_start = cost_settings$start, cost_end = cost_settings$end,
        cost_by = cost_settings$by,
        classify_fun = cfun, threshold_fun = tfun
      )
      costnp_plus_row <- data.frame(cost   = res_costnp_plus$cost,
                                    pop_t1 = res_costnp_plus$population_type1,
                                    pop_t2 = res_costnp_plus$population_type2)

      # costnp (conservative, no interpolation)
      res_costnp <- costnp(
        data, pop, alpha = TARGET_ALPHA, method = method,
        cost_start = cost_settings$start, cost_end = cost_settings$end,
        cost_by = cost_settings$by,
        classify_fun = cfun, threshold_fun = tfun
      )
      costnp_row <- data.frame(cost   = res_costnp$cost,
                               pop_t1 = res_costnp$population_type1,
                               pop_t2 = res_costnp$population_type2)

      # Naive
      res_naive <- get_cost_naive(
        data, alpha = TARGET_ALPHA, population_data = pop,
        cost_start = naive_cost_settings$start, cost_end = naive_cost_settings$end,
        cost_by = naive_cost_settings$by, method = method,
        classify_fun = cfun
      )
      naive_row <- data.frame(cost   = res_naive$cost,
                              pop_t1 = res_naive$population_type1,
                              pop_t2 = res_naive$population_type2)

      # NPROC — use npc() for supported methods, npc_manual() for GB
      if (!is.null(nproc_method)) {
        nproc_fit <- npc(x = data$x, y = data$y, alpha = TARGET_ALPHA,
                         delta = TARGET_DELTA, method = nproc_method)
        np_row <- data.frame(
          cost   = 0.5,
          pop_t1 = sum(predict(nproc_fit, newx = pop$x[pop$y == 0, ])$pred.label == 1) /
                   sum(pop$y == 0),
          pop_t2 = sum(predict(nproc_fit, newx = pop$x[pop$y == 1, ])$pred.label == 0) /
                   sum(pop$y == 1)
        )
      } else {
        res_np <- npc_manual(
          data = data, alpha = TARGET_ALPHA, delta = TARGET_DELTA,
          population_data = pop, method = method
        )
        np_row <- data.frame(
          cost   = 0.5,
          pop_t1 = res_np$population_type1,
          pop_t2 = res_np$population_type2
        )
      }

      list(costnp_plus = costnp_plus_row, costnp = costnp_row, naive = naive_row, np = np_row)
    }, .options = furrr_options(seed = TRUE), .progress = TRUE)

    add_meta <- function(df, model, n) { df$model <- model; df$n <- n; df }

    costnp_plus_results_list[[g]] <- add_meta(bind_rows(lapply(rep_results, `[[`, "costnp_plus")), g_model, g_n)
    costnp_results_list[[g]]      <- add_meta(bind_rows(lapply(rep_results, `[[`, "costnp")),      g_model, g_n)
    naive_results_list[[g]]  <- add_meta(bind_rows(lapply(rep_results, `[[`, "naive")),  g_model, g_n)
    np_results_list[[g]]     <- add_meta(bind_rows(lapply(rep_results, `[[`, "np")),     g_model, g_n)

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

  all_results[[method]] <- list(
    costnp_plus = bind_rows(costnp_plus_results_list),
    costnp      = bind_rows(costnp_results_list),
    naive       = bind_rows(naive_results_list),
    np          = bind_rows(np_results_list)
  )
}

plan(sequential)

save(all_results, file = out_path)
message("Success! Results saved to: ", out_path)

# =============================================================================
# LaTeX TABLE
# =============================================================================
load_results <- function(path) {
  e <- new.env()
  load(path, envir = e)
  e$all_results
}

MODEL_RESULTS <- list(

  planet   = load_results("data/experiment_3_planet.RData"),
  gaussian = load_results("data/experiment_3_gaussian.RData")
)

MODEL_LABELS <- c(
  planet   = "Elliptical",
  gaussian = "Tri-diagonal"
)


TABLE_ALPHA <- 0.1
TABLE_DELTA <- 0.1
ALGOS       <- c("naive", "np", "costnp_plus", "costnp")
T2_ALGOS    <- c("np", "costnp_plus", "costnp")
ALGO_LABELS <- c(naive = "Naive", np = "NP", costnp_plus = "CostNP+", costnp = "CostNP")


METHOD_LABELS <- c(
  QDA   = "QDA",
  LDA   = "LDA",
  SVM   = "SVM",
  LR    = "LR",
  NB    = "NB",
  GB    = "XGB",
  penLR = "penLR"
)

fmt_t1 <- function(v, viol_star_thr) {
  if (is.na(v)) return("---")
  star <- if (v > viol_star_thr) "$^*$" else ""
  sprintf("%.1f%s", v * 100, star)
}

fmt_t2_plain <- function(v, s) {
  if (is.na(v)) return("---")
  sprintf("%.1f (%.1f)", v * 100, s * 100)
}


n_t1   <- length(ALGOS)
n_t2   <- length(T2_ALGOS)
n_cols <- 1 + n_t1 + n_t2

build_model_block <- function(all_results, model_label) {
  res_methods <- names(all_results)
  res_n_reps  <- nrow(all_results[[1]][[1]])

  mc_se         <- sqrt(TABLE_DELTA * (1 - TABLE_DELTA) / res_n_reps)
  viol_star_thr <- TABLE_DELTA + 2 * mc_se

  viol_mat <- matrix(NA_real_, nrow = length(res_methods), ncol = length(ALGOS),
                     dimnames = list(res_methods, ALGOS))
  avg_t2   <- matrix(NA_real_, nrow = length(res_methods), ncol = length(ALGOS),
                     dimnames = list(res_methods, ALGOS))
  sd_t2    <- matrix(NA_real_, nrow = length(res_methods), ncol = length(ALGOS),
                     dimnames = list(res_methods, ALGOS))

  for (meth in res_methods) {
    for (alg in ALGOS) {
      df <- all_results[[meth]][[alg]]
      if (is.null(df) || nrow(df) == 0) next
      viol_mat[meth, alg] <- mean(df$pop_t1 > TABLE_ALPHA)
      avg_t2  [meth, alg] <- mean(df$pop_t2)
      sd_t2   [meth, alg] <- sd(df$pop_t2)
    }
  }

  cat(sprintf("\n=== [%s] TYPE I VIOLATION RATES (%%) ===\n", model_label))
  print(round(viol_mat * 100, 1))
  cat(sprintf("\n=== [%s] AVERAGE TYPE II ERROR (%%) ===\n", model_label))
  print(round(avg_t2 * 100, 1))

  row_labels <- ifelse(res_methods %in% names(METHOD_LABELS),
                       METHOD_LABELS[res_methods], res_methods)

  starred_t2 <- matrix(FALSE, nrow = length(res_methods), ncol = length(T2_ALGOS),
                       dimnames = list(res_methods, T2_ALGOS))
  for (meth in res_methods)
    for (alg in T2_ALGOS)
      starred_t2[meth, alg] <- !is.na(viol_mat[meth, alg]) && viol_mat[meth, alg] > viol_star_thr

  t2_sub <- matrix(NA_real_, nrow = length(res_methods), ncol = length(T2_ALGOS),
                   dimnames = list(res_methods, T2_ALGOS))
  for (i in seq_along(res_methods))
    for (alg in T2_ALGOS)
      t2_sub[i, alg] <- if (!starred_t2[res_methods[i], alg]) avg_t2[res_methods[i], alg] else NA_real_

  row_min_j <- apply(t2_sub, 1, function(x) if (all(is.na(x))) NA_integer_ else which.min(x))
  col_min_i <- apply(t2_sub, 2, function(x) if (all(is.na(x))) NA_integer_ else which.min(x))

  data_rows <- vapply(seq_along(res_methods), function(i) {
    meth <- res_methods[i]
    t1_cells <- vapply(ALGOS, function(alg) fmt_t1(viol_mat[meth, alg], viol_star_thr), character(1))
    t2_cells <- vapply(seq_along(T2_ALGOS), function(j) {
      alg     <- T2_ALGOS[j]
      is_star <- starred_t2[meth, alg]
      cell    <- fmt_t2_plain(avg_t2[meth, alg], sd_t2[meth, alg])
      if (is_star) cell <- paste0(cell, "$^*$")
      is_row_min <- !is_star && !is.na(row_min_j[meth]) && j == row_min_j[meth]
      is_col_min <- !is_star && !is.na(col_min_i[alg])  && i == col_min_i[alg]
      if (is_col_min && is_row_min) {
        sprintf("\\textcolor{green!60!black}{\\textbf{%s}}", cell)
      } else if (is_col_min) {
        sprintf("\\textcolor{green!60!black}{%s}", cell)
      } else if (is_row_min) {
        sprintf("\\textbf{%s}", cell)
      } else {
        cell
      }
    }, character(1))
    sprintf("%s & %s & %s \\\\",
            row_labels[i],
            paste(t1_cells, collapse = " & "),
            paste(t2_cells, collapse = " & "))
  }, character(1))

  model_header <- sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", n_cols, model_label)
  c(model_header, "\\hline", data_rows)
}

# ---- LaTeX table ------------------------------------------------------------

col_spec <- paste0("l|", paste(rep("c", n_t1), collapse = ""),
                   "|", paste(rep("c", n_t2), collapse = ""))

header_row1 <- sprintf(
  " & \\multicolumn{%d}{c|}{Type I Violation Rate (\\%%)} & \\multicolumn{%d}{c}{Avg.\\ Type II Error (\\%%, SD)} \\\\",
  n_t1, n_t2
)
header_row2 <- sprintf(
  "Clf. & %s & %s \\\\",
  paste(ALGO_LABELS[ALGOS], collapse = " & "),
  paste(ALGO_LABELS[T2_ALGOS], collapse = " & ")
)

model_blocks <- lapply(names(MODEL_RESULTS), function(m) {
  build_model_block(MODEL_RESULTS[[m]], MODEL_LABELS[m])
})
body_rows <- Reduce(function(a, b) c(a, "\\hline", b), model_blocks)

latex_lines <- c(
  "% Requires \\usepackage{xcolor} in preamble.",
  "% $^*$ = type I violation rate exceeds target $\\delta$.",
  "\\begin{table}[ht]",
  "\\centering",
  sprintf("\\begin{tabular}{%s}", col_spec),
  "\\hline",
  header_row1,
  header_row2,
  "\\hline",
  body_rows,
  "\\hline",
  "\\end{tabular}",
  sprintf(paste0(
    "\\caption{Type I violation rates and mean (sd) Type II error rates (in percentage) across classifiers, ",
    "models, and sample sizes. An asterisk indicates that the type I violation rate exceeds ",
    "the target violation level of $%.0f\\%%$ (within Monte Carlo tolerance). The lowest ",
    "average type II error in each row, among non-asterisked entries, is shown in bold.}"
  ), TABLE_ALPHA * 100),
  "\\label{Tab::Experiment3Results}",
  "\\end{table}"
)

cat("\n\n=== LaTeX Table ===\n\n")
cat(paste(latex_lines, collapse = "\n"), "\n")
