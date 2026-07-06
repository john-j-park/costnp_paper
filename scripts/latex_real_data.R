# =============================================================================
# Real Data Latex Table
# =============================================================================

DATA_DIR     <- "data"
TARGET_ALPHA <- 0.1
TARGET_DELTA <- 0.1
REPS         <- 500

VIOLATION_SLACK     <- 2 * sqrt(TARGET_DELTA * (1 - TARGET_DELTA) / REPS)
VIOLATION_THRESHOLD <- TARGET_DELTA + VIOLATION_SLACK

fixed_datasets <- c("heart", "diabetes", "transfusion")
fixed_labels   <- c(
  heart       = "Cardiotocography",
  diabetes    = "Diabetes",
  transfusion = "Transfusion"
)

sub_datasets <- c("mammography", "phoneme", "wine_quality", "magic")
sub_labels   <- c(
  mammography = "Mammography",
  phoneme     = "Phoneme",
  wine_quality = "Wine",
  magic       = "Telescope"
)

t1_methods <- c("naive", "nproc", "costnp_plus", "costnp")
t2_methods <- c("nproc", "costnp_plus", "costnp")

t1_labels  <- c(naive = "Naive", nproc = "NP", costnp_plus = "CostNP+", costnp = "CostNP")
t2_labels  <- c(nproc = "NP",   costnp_plus = "CostNP+",   costnp = "CostNP")

# =============================================================================
# Helpers
# =============================================================================

# Results are read from data/<key>.RData. The experiment writes timestamped
# files to output/; move/rename the chosen run into data/ before building the
# table. Most datasets use their own name; a few committed caches use a short key.
FILE_KEY <- c(wine_quality = "wine")
file_key <- function(dataset) if (dataset %in% names(FILE_KEY)) FILE_KEY[[dataset]] else dataset

extract_summary <- function(rep_results, alpha) {
  algos <- names(rep_results[[1]])
  out   <- list()
  for (alg in algos) {
    rows <- Filter(Negate(is.null), lapply(rep_results, `[[`, alg))
    if (length(rows) == 0) next
    t1 <- sapply(rows, `[[`, "type1")
    t2 <- sapply(rows, `[[`, "type2")
    out[[alg]] <- list(
      viol_rate = mean(t1 > alpha),
      avg_type2 = mean(t2),
      sd_type2  = sd(t2)
    )
  }
  out
}

load_results <- function(datasets, dir, alpha) {
  out <- list()
  for (dname in datasets) {
    fpath <- file.path(dir, paste0(file_key(dname), ".RData"))
    if (!file.exists(fpath)) { warning("No RData found for: ", dname); next }
    message("Loading ", basename(fpath))
    env <- new.env(parent = emptyenv())
    load(fpath, envir = env)
    if (!exists("rep_results", envir = env)) {
      warning("rep_results not found in: ", fpath); next
    }
    out[[dname]] <- extract_summary(env$rep_results, alpha)
  }
  out
}

build_matrices <- function(results, datasets, methods) {
  viol  <- matrix(NA_real_, nrow = length(datasets), ncol = length(methods),
                  dimnames = list(datasets, methods))
  type2 <- matrix(NA_real_, nrow = length(datasets), ncol = length(methods),
                  dimnames = list(datasets, methods))
  sd2   <- matrix(NA_real_, nrow = length(datasets), ncol = length(methods),
                  dimnames = list(datasets, methods))
  for (dname in datasets) {
    if (is.null(results[[dname]])) next
    for (m in methods) {
      if (!is.null(results[[dname]][[m]])) {
        viol[dname,  m] <- results[[dname]][[m]]$viol_rate
        type2[dname, m] <- results[[dname]][[m]]$avg_type2
        sd2[dname,   m] <- results[[dname]][[m]]$sd_type2
      }
    }
  }
  list(viol = viol, type2 = type2, sd2 = sd2)
}

fmt_t1 <- function(v, star) {
  if (is.na(v)) return("---")
  sprintf("%.1f%s", v * 100, if (star) "$^*$" else "")
}

fmt_t2 <- function(v, sd, bold) {
  if (is.na(v)) return("---")
  s <- sprintf("%.1f (%.1f)", v * 100, sd * 100)
  if (bold) sprintf("\\textbf{%s}", s) else s
}

build_latex_table <- function(viol, type2, sd_type2, dataset_labels,
                               t1_methods, t2_methods, t1_labels, t2_labels,
                               viol_threshold,
                               bold_min_t2 = FALSE, caption, label) {
  n_t1    <- length(t1_methods)
  n_t2    <- length(t2_methods)
  dnames  <- rownames(viol)

  col_spec <- paste0("l|", paste(rep("c", n_t1), collapse = ""), "|",
                     paste(rep("c", n_t2), collapse = ""))

  lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    "\\resizebox{\\textwidth}{!}{%",
    sprintf("\\begin{tabular}{%s}", col_spec),
    "\\hline",
    sprintf(" & \\multicolumn{%d}{c|}{Type I Violation Rate (\\%%)} & \\multicolumn{%d}{c}{Average Type II Error (\\%%)} \\\\",
            n_t1, n_t2),
    sprintf("Dataset & %s & %s \\\\",
            paste(t1_labels[t1_methods], collapse = " & "),
            paste(t2_labels[t2_methods], collapse = " & ")),
    "\\hline"
  )

  for (dname in dnames) {
    violating <- sapply(t1_methods, function(m)
      !is.na(viol[dname, m]) && viol[dname, m] > viol_threshold)

    t1_cells <- mapply(fmt_t1, viol[dname, t1_methods], violating)

    # For bold: among t2_methods, find the one with lowest avg_type2 where the
    # corresponding t1 method is not asterisked (only applied when bold_min_t2=TRUE)
    t2_viol <- sapply(t2_methods, function(m)
      !is.na(viol[dname, m]) && viol[dname, m] > viol_threshold)
    t2_vals <- type2[dname, t2_methods]
    if (bold_min_t2 && any(!t2_viol & !is.na(t2_vals))) {
      eligible <- which(!t2_viol & !is.na(t2_vals))
      best     <- eligible[which.min(t2_vals[eligible])]
      is_bold  <- seq_along(t2_methods) == best
    } else {
      is_bold <- rep(FALSE, length(t2_methods))
    }

    t2_cells <- mapply(fmt_t2, type2[dname, t2_methods], sd_type2[dname, t2_methods], is_bold)

    lines <- c(lines,
      sprintf("%s & %s & %s \\\\",
              dataset_labels[dname],
              paste(t1_cells, collapse = " & "),
              paste(t2_cells, collapse = " & ")))
  }

  lines <- c(lines,
    "\\hline",
    "\\end{tabular}",
    "}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\end{table}"
  )
  paste(lines, collapse = "\n")
}

# =============================================================================
# Load all results
# =============================================================================

all_methods <- union(t1_methods, t2_methods)
fixed_results <- load_results(fixed_datasets, DATA_DIR, TARGET_ALPHA)
sub_results   <- load_results(sub_datasets,   DATA_DIR, TARGET_ALPHA)

fixed_mats <- build_matrices(fixed_results, fixed_datasets, all_methods)
sub_mats   <- build_matrices(sub_results,   sub_datasets,   all_methods)

# =============================================================================
# Table 1: fixed-split datasets (heart, diabetes, transfusion) — no bold
# =============================================================================

latex1 <- build_latex_table(
  viol           = fixed_mats$viol,
  type2          = fixed_mats$type2,
  sd_type2       = fixed_mats$sd2,
  dataset_labels = fixed_labels,
  t1_methods     = t1_methods,
  t2_methods     = t2_methods,
  t1_labels      = t1_labels,
  t2_labels      = t2_labels,
  viol_threshold = VIOLATION_THRESHOLD,
  bold_min_t2    = FALSE,
  caption = paste0(
    "Average type I violation rates and mean (sd) Type II error rates (in percentage) ",
    "across classifiers, models, and sample sizes. ",
    "An asterisk indicates that the type I violation rate exceeds ",
    "the target violation level of $10\\%$ (within Monte Carlo tolerance)."
  ),
  label = "Tab::RealDataFixed"
)

cat("\n\n=== LaTeX Table 1 (fixed-split datasets) ===\n\n")
cat(latex1, "\n")

# =============================================================================
# Table 2: subsampled datasets (mammography, phoneme, wine_quality, magic) — bold min t2
# =============================================================================

latex2 <- build_latex_table(
  viol           = sub_mats$viol,
  type2          = sub_mats$type2,
  sd_type2       = sub_mats$sd2,
  dataset_labels = sub_labels,
  t1_methods     = t1_methods,
  t2_methods     = t2_methods,
  t1_labels      = t1_labels,
  t2_labels      = t2_labels,
  viol_threshold = VIOLATION_THRESHOLD,
  bold_min_t2    = TRUE,
  caption = paste0(
    "Average type I violation rates and mean (sd) Type II error rates (in percentage) ",
    "across classifiers, models, and sample sizes. ",
    "An asterisk indicates that the type I violation rate exceeds ",
    "the target violation level of $10\\%$ (within Monte Carlo tolerance). The lowest ",
    "average type II error in each row, among non-asterisked entries, is shown in bold."
  ),
  label = "Tab::RealDataSubs"
)

cat("\n\n=== LaTeX Table 2 (subsampled datasets) ===\n\n")
cat(latex2, "\n")
