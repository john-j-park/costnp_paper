# =============================================================================
# Experiment 4 Plots and latex tables
# =============================================================================

library(dplyr)
library(ggplot2)

METHOD_COLORS <- c(
  Naive     = "#D55E00",
  NP        = "#A67C00",
  eLDA      = "#66A61E",
  `CostNP+` = "#0072B2",
  CostNP    = "#6BAED6"
)

load_classifier <- function(file, classifier_label) {
  env <- new.env()
  load(file, envir = env)
  costnp_plus <- if (exists("costnp_plus_results", envir = env)) env$costnp_plus_results else env$calcs_results
  costnp      <- if (exists("costnp_results",      envir = env)) env$costnp_results      else env$strict_results
  costnp_plus$method <- "CostNP+"
  costnp$method      <- "CostNP"
  env$np_results$method    <- "NP"
  env$naive_results$method <- "Naive"
  combined <- bind_rows(costnp_plus, costnp, env$np_results, env$naive_results)
  combined$classifier <- classifier_label
  combined
}

all_results <- bind_rows(
  load_classifier("data/experiment_4_nb.RData", "NB"),
  load_classifier("data/experiment_4_lr.RData", "LR"),
  load_classifier("data/experiment_4_lda.RData", "LDA")
)

preferred_models      <- c("homgaussian", "gaussian", "mixture", "t", "planet_2")
preferred_methods     <- c("Naive", "NP", "CostNP+", "CostNP")
preferred_classifiers <- c("LR", "LDA", "NB")
all_results$model      <- factor(all_results$model,
                                 levels = intersect(preferred_models,      unique(all_results$model)))
all_results$method     <- factor(all_results$method,
                                 levels = intersect(preferred_methods,     unique(all_results$method)))
all_results$classifier <- factor(all_results$classifier,
                                 levels = intersect(preferred_classifiers, unique(all_results$classifier)))
TARGET_ALPHA <- 0.1

summary_all <- all_results %>%
  group_by(method, model, n, classifier) %>%
  summarise(
    violation_rate = mean(pop_t1 > TARGET_ALPHA),
    .groups = 'drop'
  )


LABEL_HEADROOM <- 0.10
make_plot <- function(clf, n_val, x_max, model_val = NULL) {
  sub_data    <- all_results %>% filter(classifier == clf, n == n_val)
  sub_summary <- summary_all %>% filter(classifier == clf, n == n_val)
  if (!is.null(model_val)) {
    sub_data    <- sub_data    %>% filter(model == model_val)
    sub_summary <- sub_summary %>% filter(model == model_val)
  }

  y_lo  <- 0
  y_hi  <- x_max
  y_top <- y_hi + (y_hi - y_lo) * LABEL_HEADROOM
  y_lab <- y_hi + (y_hi - y_lo) * (LABEL_HEADROOM / 2)

  p <- ggplot(sub_data, aes(x = method, y = pop_t1, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5, width = 0.6) +
    geom_hline(yintercept = TARGET_ALPHA, linetype = "dashed",
               color = "red", linewidth = 0.8) +
    geom_text(data = sub_summary,
              aes(x = method, y = y_lab,
                  label = sprintf("%.1f%%", round(violation_rate * 100, 1))),
              hjust = 0.5, vjust = 0.5, size = 2.6, fontface = "bold",
              inherit.aes = FALSE) +
    scale_fill_manual(values = METHOD_COLORS) +
    labs(x = NULL, y = "Type I Error") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(ylim = c(y_lo, y_top))
  if (is.null(model_val)) {
    p <- p + facet_wrap(~ model, ncol = 1)
  } else {
    p <- p + ggtitle(tools::toTitleCase(model_val))
  }
  p
}

make_plot_t2 <- function(clf, n_val, x_max, model_val = NULL) {
  sub_data    <- all_results %>% filter(classifier == clf, n == n_val)
  sub_summary <- sub_data %>%
    group_by(method, model) %>%
    summarise(mean_t2 = mean(pop_t2, na.rm = TRUE), .groups = "drop")
  if (!is.null(model_val)) {
    sub_data    <- sub_data    %>% filter(model == model_val)
    sub_summary <- sub_summary %>% filter(model == model_val)
  }

  y_lo  <- max(0, quantile(sub_data$pop_t2, 0.005, na.rm = TRUE) - 0.02)
  y_hi  <- x_max
  y_top <- y_hi + (y_hi - y_lo) * LABEL_HEADROOM
  y_lab <- y_hi + (y_hi - y_lo) * (LABEL_HEADROOM / 2)

  p <- ggplot(sub_data, aes(x = method, y = pop_t2, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5, width = 0.6) +
    geom_text(data = sub_summary,
              aes(x = method, y = y_lab,
                  label = sprintf("%.3f", round(mean_t2, 3))),
              hjust = 0.5, vjust = 0.5, size = 2.6, fontface = "bold",
              inherit.aes = FALSE) +
    scale_fill_manual(values = METHOD_COLORS) +
    labs(x = NULL, y = "Type II Error") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(ylim = c(y_lo, y_top))
  if (is.null(model_val)) {
    p <- p + facet_wrap(~ model, ncol = 1)
  } else {
    p <- p + ggtitle(tools::toTitleCase(model_val))
  }
  p
}

make_grid <- function(n_val, x_max, metric = c("t1", "t2")) {
  metric  <- match.arg(metric)
  err_col <- if (metric == "t1") "pop_t1" else "pop_t2"

  sub_data <- all_results %>% filter(n == n_val)

  y_lo <- if (metric == "t1") 0 else
          max(0, quantile(sub_data$pop_t2, 0.005, na.rm = TRUE) - 0.02)
  y_hi  <- x_max
  y_top <- y_hi + (y_hi - y_lo) * LABEL_HEADROOM
  y_lab <- y_hi + (y_hi - y_lo) * (LABEL_HEADROOM / 2)

  if (metric == "t1") {
    sub_summary <- summary_all %>% filter(n == n_val)
    label_aes   <- aes(x = method, y = y_lab,
                       label = sprintf("%.1f%%", round(violation_rate * 100, 1)))
    y_title <- "Type I Error"
  } else {
    sub_summary <- sub_data %>%
      group_by(method, model, classifier) %>%
      summarise(mean_t2 = mean(pop_t2, na.rm = TRUE), .groups = "drop")
    label_aes <- aes(x = method, y = y_lab,
                     label = sprintf("%.3f", round(mean_t2, 3)))
    y_title <- "Type II Error"
  }

  p <- ggplot(sub_data, aes(x = method, y = .data[[err_col]], fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5, width = 0.6)
  if (metric == "t1") {
    p <- p + geom_hline(yintercept = TARGET_ALPHA, linetype = "dashed",
                        color = "red", linewidth = 0.8)
  }
  p +
    geom_text(data = sub_summary, label_aes,
              hjust = 0.5, vjust = 0.5, size = 2.4, fontface = "bold",
              inherit.aes = FALSE) +
    scale_fill_manual(values = METHOD_COLORS) +
    facet_grid(model ~ classifier, labeller = labeller(model = tools::toTitleCase)) +
    labs(x = NULL, y = y_title) +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(ylim = c(y_lo, y_top))
}


p_LR_1000 <- make_plot("LR", 1000, 0.18)
print(p_LR_1000)
ggsave("Images/Synthetic_1000_LR.pdf", plot = p_LR_1000, width = 5, height = 5)

p_LR_5000 <- make_plot("LR", 5000, 0.15)
print(p_LR_5000)
ggsave("Images/Synthetic_5000_LR.pdf", plot = p_LR_5000, width = 5, height = 5)

p_LDA_1000 <- make_plot("LDA", 1000, 0.18)
print(p_LDA_1000)
ggsave("Images/Synthetic_1000_LDA.pdf", plot = p_LDA_1000, width = 5, height = 5)

p_LDA_5000 <- make_plot("LDA", 5000, 0.15)
print(p_LDA_5000)
ggsave("Images/Synthetic_5000_LDA.pdf", plot = p_LDA_5000, width = 5, height = 5)

p_NB_1000 <- make_plot("NB", 1000, 0.18)
print(p_NB_1000)
ggsave("Images/Synthetic_1000_NB.pdf", plot = p_NB_1000, width = 5, height = 5)

p_NB_5000 <- make_plot("NB", 5000, 0.15)
print(p_NB_5000)
ggsave("Images/Synthetic_5000_NB.pdf", plot = p_NB_5000, width = 5, height = 5)

p_t2_LR_1000 <- make_plot_t2("LR", 1000, 0.6)
print(p_t2_LR_1000)
ggsave("Images/Synthetic_1000_LR_t2.pdf", plot = p_t2_LR_1000, width = 5, height = 5)

p_t2_LR_5000 <- make_plot_t2("LR", 5000, 0.6)
print(p_t2_LR_5000)
ggsave("Images/Synthetic_5000_LR_t2.pdf", plot = p_t2_LR_5000, width = 5, height = 5)

p_t2_LDA_1000 <- make_plot_t2("LDA", 1000, 0.9)
print(p_t2_LDA_1000)
ggsave("Images/Synthetic_1000_LDA_t2.pdf", plot = p_t2_LDA_1000, width = 5, height = 5)

p_t2_LDA_5000 <- make_plot_t2("LDA", 5000, 0.8)
print(p_t2_LDA_5000)
ggsave("Images/Synthetic_5000_LDA_t2.pdf", plot = p_t2_LDA_5000, width = 5, height = 5)

p_t2_NB_1000 <- make_plot_t2("NB", 1000, 0.6)
print(p_t2_NB_1000)
ggsave("Images/Synthetic_1000_NB_t2.pdf", plot = p_t2_NB_1000, width = 5, height = 5)

p_t2_NB_5000 <- make_plot_t2("NB", 5000, 0.6)
print(p_t2_NB_5000)
ggsave("Images/Synthetic_5000_NB_t2.pdf", plot = p_t2_NB_5000, width = 5, height = 5)

p_grid_t1_1000 <- make_grid(1000, 0.18, "t1")
print(p_grid_t1_1000)
ggsave("Images/Experiment_4_grid_1000_t1.pdf", plot = p_grid_t1_1000, width = 10, height = 5)

p_grid_t1_5000 <- make_grid(5000, 0.15, "t1")
print(p_grid_t1_5000)
ggsave("Images/Experiment_4_grid_5000_t1.pdf", plot = p_grid_t1_5000, width = 10, height = 5)

p_grid_t2_1000 <- make_grid(1000, 0.9, "t2")
print(p_grid_t2_1000)
ggsave("Images/Experiment_4_grid_1000_t2.pdf", plot = p_grid_t2_1000, width = 10, height = 5)

p_grid_t2_5000 <- make_grid(5000, 0.8, "t2")
print(p_grid_t2_5000)
ggsave("Images/Experiment_4_grid_5000_t2.pdf", plot = p_grid_t2_5000, width = 10, height = 5)

# =============================================================================
# LaTeX TABLE: Violation rates (%) and mean Type II error by method
# =============================================================================

t1_methods <- levels(all_results$method)

t1_wide <- summary_all %>%
  mutate(vr = sprintf("%.1f", violation_rate * 100)) %>%
  tidyr::pivot_wider(
    id_cols = c(classifier, model, n),
    names_from = method,
    values_from = vr
  )

vr_wide <- summary_all %>%
  tidyr::pivot_wider(
    id_cols = c(classifier, model, n),
    names_from = method,
    values_from = violation_rate
  ) %>%
  rename_with(~ paste0(.x, "_vr"), all_of(t1_methods))

t2_summary <- all_results %>%
  group_by(method, model, n, classifier) %>%
  summarise(
    mean_t2 = mean(pop_t2, na.rm = TRUE),
    sd_t2   = sd(pop_t2,   na.rm = TRUE),
    .groups = "drop"
  )

t2_wide <- t2_summary %>%
  mutate(mt2 = sprintf("%.1f (%.1f)", mean_t2 * 100, sd_t2 * 100)) %>%
  tidyr::pivot_wider(
    id_cols = c(classifier, model, n),
    names_from = method,
    values_from = mt2
  ) %>%
  rename_with(~ paste0(.x, "_t2"), all_of(t1_methods))

t2_raw_wide <- t2_summary %>%
  tidyr::pivot_wider(
    id_cols = c(classifier, model, n),
    names_from = method,
    values_from = mean_t2
  ) %>%
  rename_with(~ paste0(.x, "_t2_raw"), all_of(t1_methods))

latex_data <- left_join(t1_wide, t2_wide,     by = c("classifier", "model", "n")) %>%
              left_join(t2_raw_wide,            by = c("classifier", "model", "n")) %>%
              left_join(vr_wide,                by = c("classifier", "model", "n"))

VIOLATION_SLACK <- 0.028


t2_methods   <- t1_methods[t1_methods != "Naive"]
n_t1         <- length(t1_methods)
n_t2         <- length(t2_methods)
multi_clf    <- nlevels(all_results$classifier) > 1
clf_col_spec <- if (multi_clf) "lll|" else "ll|"
clf_col_hdr  <- if (multi_clf) "Clf. & Model & $n$" else "Model & $n$"

# Build LaTeX table string
if (n_t1 == 1) {
  header <- paste0(
    "\\begin{table}[ht]\n",
    "\\centering\n",
    "\\begin{tabular}{", clf_col_spec, "c|c}\n",
    "\\hline\n",
    clf_col_hdr, " & Type I Viol.\\ (\\%) & Avg.\\ $R_1$ (\\%) \\\\\n",
    "\\hline\n"
  )
} else {
  header <- paste0(
    "\\begin{table}[ht]\n",
    "\\centering\n",
    "\\small\n",
    "\\begin{tabular}{", clf_col_spec,
    paste(rep("c", n_t1), collapse = ""),
    "|", paste(rep("c", n_t2), collapse = ""), "}\n",
    "\\hline\n",
    " & & & \\multicolumn{", n_t1, "}{c|}{Type I Violation Rate (\\%)} &",
    " \\multicolumn{", n_t2, "}{c}{Average type II error (\\%, mean (sd))} \\\\\n",
    clf_col_hdr, " & ",
    paste(t1_methods, collapse = " & "), " & ",
    paste(t2_methods, collapse = " & "), " \\\\\n",
    "\\hline\n"
  )
}

body <- ""
for (clf in levels(all_results$classifier)) {
  clf_rows <- latex_data %>% filter(classifier == clf)
  for (mod in levels(all_results$model)) {
    mod_rows <- clf_rows %>% filter(model == mod)
    for (i in seq_len(nrow(mod_rows))) {
      clf_label <- if (multi_clf) {
        if (mod == levels(all_results$model)[1] && i == 1) paste0(clf, " & ") else " & "
      } else ""
      mod_label <- if (i == 1) mod else ""
      t1_vals <- paste(sapply(t1_methods, function(m) mod_rows[[m]][i]), collapse = " & ")
      t2_vals <- paste(sapply(t2_methods, function(m) {
        t2_val <- mod_rows[[paste0(m, "_t2")]][i]
        vr_val <- mod_rows[[paste0(m, "_vr")]][i]
        flagged <- !is.na(vr_val) && vr_val > TARGET_ALPHA + VIOLATION_SLACK
        if (flagged) paste0(t2_val, "$^*$") else t2_val
      }), collapse = " & ")
      body <- paste0(body, sprintf(
        "%s%s & %s & %s & %s \\\\\n",
        clf_label, mod_label, mod_rows$n[i], t1_vals, t2_vals
      ))
    }
  }
  body <- paste0(body, "\\hline\n")
}

footer <- paste0(
  "\\end{tabular}\n",
  "\\caption{Type I violation rates and mean (sd) Type II error rates (in percentage) across classifiers, ",
  "models, and sample sizes.}\n",
  "\\label{Tab::Experiment4Results}\n",
  "\\end{table}"
)
latex_table <- paste0(header, body, footer)
cat(latex_table)
