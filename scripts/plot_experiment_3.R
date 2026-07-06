# =============================================================================
# Experiment 3 — boxplots (varying classifier)
# Boxplot analog of the Experiment 3 table (experiment_3.R). The x axis is the
# cost-selection *algorithm* (Naive / NP / CostNP+ / CostNP) and panels are
# faceted by the base *classifier* (QDA, LDA, SVM, LR, NB, XGB, penLR).
#
# Loads, per data model, an .RData file containing `all_results`: a named list
# keyed by classifier, each element a named list (costnp_plus, costnp, naive,
# np) of data frames with columns cost, pop_t1, pop_t2, model, n.
# =============================================================================

library(ggplot2)
library(dplyr)

TARGET_ALPHA <- 0.1

# Files keyed by data model (mirrors MODEL_RESULTS in experiment_3.R)
MODEL_FILES <- list(
  planet   = "data/experiment_3_planet.RData",
  gaussian = "data/experiment_3_gaussian.RData"
)
MODEL_LABELS <- c(planet = "Elliptical", gaussian = "Tri-diagonal")

# Algorithm (x axis) — labels, colours, left-to-right order
ALGO_LABELS <- c(naive = "Naive", np = "NP", costnp_plus = "CostNP+", costnp = "CostNP")
ALGO_LEVELS <- c("Naive", "NP", "CostNP+", "CostNP")
ALGO_COLORS <- c(
  Naive      = "#D55E00",  # vermillion / warning color
  NP         = "#A67C00",  # strong green
  `CostNP+`  = "#0072B2",  # strong blue
  CostNP     = "#6BAED6"   # light blue, visually related to CostNP+
)

# Classifier (facet) — labels and order
CLF_LABELS <- c(QDA = "QDA", LDA = "LDA", SVM = "SVM", LR = "LR",
                NB = "NB", penLR = "penLR")
CLF_LEVELS <- c("QDA", "LDA", "SVM", "LR", "NB", "penLR")

# --- Assemble long-format results from a single model's all_results list ---
make_all_results <- function(all_results) {
  clfs <- names(all_results)
  bind_rows(lapply(clfs, function(clf) {
    bind_rows(lapply(names(ALGO_LABELS), function(alg) {
      df <- all_results[[clf]][[alg]]
      if (is.null(df) || nrow(df) == 0) return(NULL)
      data.frame(clf    = clf,
                 algo   = ALGO_LABELS[[alg]],
                 pop_t1 = df$pop_t1,
                 pop_t2 = df$pop_t2)
    }))
  })) %>%
    mutate(
      algo = factor(algo, levels = ALGO_LEVELS),
      clf  = factor(ifelse(clf %in% names(CLF_LABELS), CLF_LABELS[clf], clf),
                    levels = CLF_LEVELS)
    )
}

# Fraction of the visible y-window reserved above the data for the annotation row
LABEL_HEADROOM <- 0.10

# --- Type I error boxplots, faceted by classifier ---
make_box_t1 <- function(all_results, violation_summary, y_lo, y_hi) {
  y_top <- y_hi + (y_hi - y_lo) * LABEL_HEADROOM
  y_lab <- y_hi + (y_hi - y_lo) * (LABEL_HEADROOM / 2)
  ggplot(all_results, aes(x = algo, y = pop_t1, fill = algo)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5, width = 0.6) +
    geom_hline(yintercept = TARGET_ALPHA, linetype = "dashed",
               color = "red", linewidth = 0.8) +
    geom_text(data = violation_summary,
              aes(x = algo, y = y_lab,
                  label = sprintf("%.1f%%", violation_rate * 100)),
              hjust = 0.5, vjust = 0.5, size = 2.2, fontface = "bold",
              inherit.aes = FALSE) +
    scale_fill_manual(values = ALGO_COLORS) +
    scale_x_discrete(limits = ALGO_LEVELS) +
    facet_wrap(~ clf) +
    coord_cartesian(ylim = c(y_lo, y_top)) +
    labs(x = NULL, y = "Type I Error Rate") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(hjust = 0.5))
}

# --- Type II error boxplots, faceted by classifier ---
make_box_t2 <- function(all_results, type2_summary, y_lo, y_hi) {
  y_top <- y_hi + (y_hi - y_lo) * LABEL_HEADROOM
  y_lab <- y_hi + (y_hi - y_lo) * (LABEL_HEADROOM / 2)
  ggplot(all_results, aes(x = algo, y = pop_t2, fill = algo)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5, width = 0.6) +
    geom_text(data = type2_summary,
              aes(x = algo, y = y_lab,
                  label = sprintf("%.3f", mean_t2)),
              hjust = 0.5, vjust = 0.5, size = 2.2, fontface = "bold",
              inherit.aes = FALSE) +
    scale_fill_manual(values = ALGO_COLORS) +
    scale_x_discrete(limits = ALGO_LEVELS) +
    facet_wrap(~ clf) +
    coord_cartesian(ylim = c(y_lo, y_top)) +
    labs(x = NULL, y = "Type II Error Rate") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(hjust = 0.5))
}

# Load `all_results` from a file into an isolated environment (each file defines
# the same object name, so loading into new.env() avoids clobbering).
load_results <- function(path) {
  e <- new.env()
  load(path, envir = e)
  e$all_results
}

# =============================================================================
# BUILD — one figure pair per data model
# =============================================================================

for (model in names(MODEL_FILES)) {
  all_results <- make_all_results(load_results(MODEL_FILES[[model]]))

  violation_summary <- all_results %>%
    group_by(algo, clf) %>%
    summarise(violation_rate = mean(pop_t1 > TARGET_ALPHA, na.rm = TRUE),
              .groups = "drop")

  type2_summary <- all_results %>%
    group_by(algo, clf) %>%
    summarise(mean_t2 = mean(pop_t2, na.rm = TRUE), .groups = "drop")

  # Data windows (trimmed to the bulk so outliers don't drive the range).
  y_lo_t1 <- 0
  y_hi_t1 <- quantile(all_results$pop_t1, 0.995, na.rm = TRUE)

  t2_pad  <- 0.02
  y_lo_t2 <- max(0, quantile(all_results$pop_t2, 0.005, na.rm = TRUE) - t2_pad)
  y_hi_t2 <- quantile(all_results$pop_t2, 0.995, na.rm = TRUE)

  p_t1 <- make_box_t1(all_results, violation_summary, y_lo_t1, y_hi_t1)
  p_t2 <- make_box_t2(all_results, type2_summary,     y_lo_t2, y_hi_t2)

  print(p_t1)
  print(p_t2)

  ggsave(sprintf("images/experiment_3_test_%s_type1_box.pdf", model),
         plot = p_t1, width = 12, height = 6, dpi = 120)
  ggsave(sprintf("images/experiment_3_test_%s_type2_box.pdf", model),
         plot = p_t2, width = 12, height = 6, dpi = 120)
}

