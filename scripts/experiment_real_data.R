# =============================================================================
# Real Data Analysis
# =============================================================================

library(costnp)
library(kernlab)
library(dplyr)
library(furrr)
library(future)
library(patchwork)
library(pbapply)
library(nproc)
source("scripts/paper_utils.R")

out_dir <- "output"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

DATASETS <- c("mammography", "phoneme", "wine_quality", "magic",
              "diabetes", "transfusion", "heart")
plan(multisession, workers = parallel::detectCores() - 1)

# =============================================================================
# DATA & CONFIGURATION
# =============================================================================
TARGET_ALPHA <- 0.1
TARGET_DELTA <- 0.1
REPS <- 500
split_ratio <- 0.7

ALGORITHMS <- c("costnp_plus", "costnp", "nproc", "naive", "np_man", "elda")

split_data <- function(x, y) {
  ind0 <- which(y == 0)
  ind1 <- which(y == 1)
  train_ind0 <- sample(ind0, floor(split_ratio * length(ind0)))
  train_ind1 <- sample(ind1, floor(split_ratio * length(ind1)))
  train_ind <- seq_along(y) %in% c(train_ind0, train_ind1)
  return(list(
    train_data = list(x = x[train_ind, ], y = y[train_ind]),
    test_data  = list(x = x[!train_ind, ], y = y[!train_ind])
  ))
}

all_results <- list()

for (DATASET in DATASETS) {
  cat(sprintf("\n================================================================\n"))
  cat(sprintf("  Dataset: %s\n", DATASET))
  cat(sprintf("================================================================\n"))

  ld              <- load_dataset(DATASET)
  x_all           <- ld$x_all
  y_all           <- ld$y_all
  subsample_data  <- ld$subsample_data
  population_data <- ld$population_data

  set.seed(1)
  uses_subsampling <- DATASET %in% c("mammography", "phoneme", "wine_quality", "magic")
  if (!uses_subsampling) {
    cat(sprintf("Pre-generating %d shared train/test splits...\n", REPS))
    data_splits <- lapply(1:REPS, function(i) split_data(x_all, y_all))
  }

  rep_results <- future_map(seq_len(REPS), function(i) {
    dl  <- if (uses_subsampling) subsample_data() else data_splits[[i]]
    pop <- if (uses_subsampling) dl$pop_data else dl$test_data
    train_data <- dl$train_data

    # costnp_plus (linear interpolation)
    costnp_plus_row <- if ("costnp_plus" %in% ALGORITHMS) {
      res <- costnp_plus(
        train_data, pop,
        alpha = TARGET_ALPHA, delta = TARGET_DELTA,
        cost_start = 0.99, cost_end = 0.2, cost_by = 0.02,
        method = "LR"
      )
      data.frame(type1 = res$population_type1, cost = res$cost, type2 = res$population_type2)
    } else NULL

    # costnp (no interpolation)
    costnp_row <- if ("costnp" %in% ALGORITHMS) {
      res <- costnp(
        train_data, pop,
        alpha = TARGET_ALPHA, delta = TARGET_DELTA,
        cost_start = 0.99, cost_end = 0.2, cost_by = 0.02,
        method = "LR"
      )
      data.frame(type1 = res$population_type1, cost = res$cost, type2 = res$population_type2)
    } else NULL

    # Naive
    naive_row <- if ("naive" %in% ALGORITHMS) {
      res <- get_cost_naive(
        data = train_data,
        alpha = TARGET_ALPHA,
        population_data = pop,
        verbose = FALSE,
        cost_start = 0.4, cost_end = 0.99, cost_by = 0.01, data_pct = 0.5,
        method = "LR"
      )
      data.frame(type1 = res$population_type1, cost = res$cost, type2 = res$population_type2)
    } else NULL

    # NPROC
    nproc_row <- if ("nproc" %in% ALGORITHMS) {
      fit <- npc(x = train_data$x, y = train_data$y,
                 alpha = TARGET_ALPHA, delta = TARGET_DELTA,
                 split = 1, method = "logistic", randSeed = i)
      data.frame(
        type1 = sum(predict(fit, newx = pop$x[pop$y == 0, ])$pred.label == 1) / sum(pop$y == 0),
        type2 = sum(predict(fit, newx = pop$x[pop$y == 1, ])$pred.label == 0) / sum(pop$y == 1)
      )
    } else NULL

    # eLDA
    elda_row <- if ("elda" %in% ALGORITHMS) {
      res <- elda_npc(dl$train_data, alpha = TARGET_ALPHA, population_data = pop,
                      delta = TARGET_DELTA, variant = "eLDA")
      data.frame(type1 = res$population_type1, type2 = res$population_type2)
    } else NULL

    # npMan
    np_man_row <- if ("np_man" %in% ALGORITHMS) {
      res <- npc_manual(
        data = train_data,
        alpha = TARGET_ALPHA,
        delta = TARGET_DELTA,
        population_data = pop, method = "LR"
      )
      data.frame(type1 = res$population_type1, cost = 0.5, type2 = res$population_type2)
    } else NULL

    list(costnp_plus = costnp_plus_row, costnp = costnp_row, naive = naive_row,
         nproc = nproc_row, elda = elda_row, np_man = np_man_row)
  }, .options = furrr_options(seed = TRUE), .progress = TRUE)

  if ("costnp_plus" %in% ALGORITHMS) {
    costnp_plus_results   <- bind_rows(lapply(rep_results, `[[`, "costnp_plus"))
    costnp_plus_type1errs <- costnp_plus_results$type1; costnp_plus_costs <- costnp_plus_results$cost; costnp_plus_type2errs <- costnp_plus_results$type2
    cat(sprintf("  CostNP+  — violation rate: %.1f%%, avg type2: %.3f\n",
                mean(costnp_plus_type1errs > TARGET_ALPHA) * 100, mean(costnp_plus_type2errs)))
  }
  if ("costnp" %in% ALGORITHMS) {
    costnp_results   <- bind_rows(lapply(rep_results, `[[`, "costnp"))
    costnp_type1errs <- costnp_results$type1; costnp_costs <- costnp_results$cost; costnp_type2errs <- costnp_results$type2
    cat(sprintf("  CostNP   — violation rate: %.1f%%, avg type2: %.3f\n",
                mean(costnp_type1errs > TARGET_ALPHA) * 100, mean(costnp_type2errs)))
  }
  if ("naive" %in% ALGORITHMS) {
    naive_results   <- bind_rows(lapply(rep_results, `[[`, "naive"))
    naive_type1errs <- naive_results$type1; naive_costs <- naive_results$cost; naive_type2errs <- naive_results$type2
    cat(sprintf("  Naive    — violation rate: %.1f%%, avg type2: %.3f\n",
                mean(naive_type1errs > TARGET_ALPHA) * 100, mean(naive_type2errs)))
  }
  if ("nproc" %in% ALGORITHMS) {
    nproc_results   <- bind_rows(lapply(rep_results, `[[`, "nproc"))
    nproc_type1errs <- nproc_results$type1; nproc_type2errs <- nproc_results$type2
    cat(sprintf("  NPROC    — violation rate: %.1f%%, avg type2: %.3f\n",
                mean(nproc_type1errs > TARGET_ALPHA) * 100, mean(nproc_type2errs)))
  }
  if ("elda" %in% ALGORITHMS) {
    elda_results   <- bind_rows(lapply(rep_results, `[[`, "elda"))
    elda_type1errs <- elda_results$type1; elda_type2errs <- elda_results$type2
    cat(sprintf("  eLDA     — violation rate: %.1f%%, avg type2: %.3f\n",
                mean(elda_type1errs > TARGET_ALPHA) * 100, mean(elda_type2errs)))
  }
  if ("np_man" %in% ALGORITHMS) {
    np_man_results   <- bind_rows(lapply(rep_results, `[[`, "np_man"))
    np_man_type1errs <- np_man_results$type1; np_man_type2errs <- np_man_results$type2
    cat(sprintf("  npMan    — violation rate: %.1f%%, avg type2: %.3f\n",
                mean(np_man_type1errs > TARGET_ALPHA) * 100, mean(np_man_type2errs)))
  }

  all_results[[DATASET]] <- rep_results

  out_path <- file.path(out_dir, paste0("real_data_", DATASET, "_", format(Sys.time(), "%Y%m%d_%H%M"), ".RData"))
  save(rep_results, file = out_path)
  message("Saved: ", out_path)
}

plan(sequential)
