# =============================================================================
# Real Data Plots
# =============================================================================

library(ggplot2)
library(dplyr)

TARGET_ALPHA <- 0.15
TARGET_DELTA <- 0.1
REPS         <- 500

VIOLATION_SLACK     <- 2 * sqrt(TARGET_DELTA * (1 - TARGET_DELTA) / REPS)
VIOLATION_THRESHOLD <- TARGET_DELTA + VIOLATION_SLACK

METHOD_COLORS <- c(
  Naive     = "#D9534F",
  NP        = "#A67C00",
  `CostNP+` = "#0072B2",
  CostNP    = "#6BAED6"
)

ALGO_LABELS <- c(
  naive = "Naive", nproc = "NP", np = "NP",
  costnp_plus = "CostNP+", calcs = "CostNP+",
  costnp = "CostNP",       strict = "CostNP"
)
T1_ALGOS <- c("naive", "nproc", "np", "costnp_plus", "calcs", "costnp", "strict")
T2_ALGOS <- c("nproc", "np", "costnp_plus", "calcs", "costnp", "strict")

# Results are read from data/<key>.RData. The experiment writes timestamped
# files to output/; move/rename the chosen run into data/ before plotting.
DATA_DIR <- "data"

extract <- function(rep_results, algo, field) {
  rows <- Filter(Negate(is.null), lapply(rep_results, `[[`, algo))
  if (length(rows) == 0) return(numeric(0))
  as.numeric(sapply(rows, `[[`, field))
}

make_plots <- function(rep_results, dataset_name, target_alpha, viol_threshold) {
  t1_data <- lapply(setNames(T1_ALGOS, T1_ALGOS), extract,
                    rep_results = rep_results, field = "type1")
  t2_data <- lapply(setNames(T2_ALGOS, T2_ALGOS), extract,
                    rep_results = rep_results, field = "type2")

  present  <- function(algos, data) algos[sapply(data[algos], length) > 0]
  dedup    <- function(algos) algos[!duplicated(ALGO_LABELS[algos])]
  t1_avail <- dedup(present(T1_ALGOS, t1_data))
  t2_avail <- dedup(present(T2_ALGOS, t2_data))

  if (length(t1_avail) == 0) { warning("No type1 data for: ", dataset_name); return(NULL) }

  t1_labels <- ALGO_LABELS[t1_avail]
  t2_labels <- ALGO_LABELS[t2_avail]

  df_type1 <- data.frame(
    pop_t1 = unlist(t1_data[t1_avail], use.names = FALSE),
    method  = factor(
      rep(t1_labels, times = sapply(t1_data[t1_avail], length)),
      levels = t1_labels
    )
  )

  summary_type1 <- data.frame(
    method         = factor(t1_labels, levels = t1_labels),
    violation_rate = sapply(t1_data[t1_avail], function(v) mean(v > target_alpha))
  )

  p1 <- ggplot(df_type1, aes(x = method, y = pop_t1, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    geom_hline(yintercept = target_alpha, linetype = "dashed", color = "red", linewidth = 1) +
    geom_text(data = summary_type1,
              aes(x = method, y = Inf,
                  label = sprintf("%.1f%%", violation_rate * 100)),
              vjust = 1.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = METHOD_COLORS) +
    labs(title = dataset_name, x = "Method", y = "Type I Error") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.margin = ggplot2::margin(t = 25, r = 10, b = 10, l = 10)) +
    coord_cartesian(clip = "off")

  if (length(t2_avail) == 0) {
    warning("No type2 data for: ", dataset_name)
    return(list(p1 = p1, p2 = NULL))
  }

  df_type2 <- data.frame(
    pop_t2 = unlist(t2_data[t2_avail], use.names = FALSE),
    method  = factor(
      rep(t2_labels, times = sapply(t2_data[t2_avail], length)),
      levels = t2_labels
    )
  )

  summary_type2 <- data.frame(
    method     = factor(t2_labels, levels = t2_labels),
    mean_type2 = sapply(t2_data[t2_avail], mean)
  )

  non_trivial <- df_type2$pop_t2[df_type2$pop_t2 < 1]
  y_max <- quantile(non_trivial, 0.99, na.rm = TRUE) * 1.3

  p2 <- ggplot(df_type2, aes(x = method, y = pop_t2, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    geom_text(data = summary_type2,
              aes(x = method, y = Inf,
                  label = sprintf("Mean: %.1f%%", mean_type2 * 100)),
              vjust = 1.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = METHOD_COLORS) +
    labs(title = dataset_name, x = "Method", y = "Type II Error") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.margin = ggplot2::margin(t = 25, r = 10, b = 10, l = 10)) +
    coord_cartesian(ylim = c(0, y_max), clip = "off")

  list(p1 = p1, p2 = p2)
}

# =============================================================================
# Loop over datasets
# =============================================================================

datasets <- list(
  list(key = "heart",        name = "Heart"),
  list(key = "diabetes",     name = "Diabetes"),
  list(key = "transfusion",  name = "Transfusion"),
  list(key = "mammography",  name = "Mammography"),
  list(key = "phoneme",      name = "Phoneme"),
  list(key = "wine", name = "Wine"),
  list(key = "magic",        name = "Telescope")
)

for (ds in datasets) {
  fpath <- file.path(DATA_DIR, paste0(ds$key, ".RData"))
  if (!file.exists(fpath)) { warning("No RData found for: ", ds$key); next }
  message("Loading ", basename(fpath))
  env <- new.env(parent = emptyenv())
  load(fpath, envir = env)
  plots <- make_plots(env$rep_results, ds$name, TARGET_ALPHA, VIOLATION_THRESHOLD)
  if (is.null(plots)) next
  ggsave(paste0("images/", ds$name, "_T1.pdf"), plot = plots$p1, width = 8, height = 4)
  if (!is.null(plots$p2))
    ggsave(paste0("images/", ds$name, "_T2.pdf"), plot = plots$p2, width = 8, height = 4)
}

# =============================================================================
# LaTeX subfigure code
# =============================================================================
figures <- paste(sapply(datasets, function(ds) {
  name    <- ds$name
  t1_file <- sprintf("Images/%s_T1.pdf", name)
  t2_file <- sprintf("Images/%s_T2.pdf", name)
  sprintf(
    paste0(
      "\\begin{figure}[ht]\n  \\centering\n",
      "  \\begin{subfigure}[t]{0.48\\textwidth}\n    \\centering\n",
      "    \\includegraphics[width=\\linewidth]{%s}\n",
      "    \\caption{Type I Error}\n    \\label{fig:%s_t1}\n  \\end{subfigure}\n",
      "  \\hfill\n",
      "  \\begin{subfigure}[t]{0.48\\textwidth}\n    \\centering\n",
      "    \\includegraphics[width=\\linewidth]{%s}\n",
      "    \\caption{Type II Error}\n    \\label{fig:%s_t2}\n  \\end{subfigure}\n",
      "  \\caption{Type I and Type II error distributions for the %s dataset.}\n",
      "  \\label{fig:%s_errors}\n\\end{figure}"
    ),
    t1_file, tolower(gsub(" ", "_", name)),
    t2_file, tolower(gsub(" ", "_", name)),
    name,    tolower(gsub(" ", "_", name))
  )
}), collapse = "\n\n")

cat(sprintf("\n\n%% Add to preamble: \\usepackage{subcaption}\n\n%s\n", figures))

