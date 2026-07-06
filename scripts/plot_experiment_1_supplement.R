# =============================================================================
# Supplemental Plots for Experiment 1
# =============================================================================

library(ggplot2)
library(dplyr)

TARGET_ALPHA  <- 0.1
bl_type2_theo <- 0.1943861   # theoretical best-linear Type II baseline

METHOD_COLORS <- c(
  Naive      = "#D55E00",  # vermillion / warning color
  NP         = "#A67C00",  # strong green
  `CostNP+`  = "#0072B2",  # strong blue
  CostNP     = "#6BAED6"   # light blue, visually related to CostNP+
)
# Left-to-right order along the x axis
METHOD_LEVELS <- c("Naive", "NP", "CostNP+", "CostNP")

# --- Assemble long-format results from the four method frames ---
make_all_results <- function(costnp_plus_results, costnp_results,
                             naive_results, np_results) {
  bind_rows(
    costnp_plus_results %>% mutate(method = "CostNP+"),
    costnp_results      %>% mutate(method = "CostNP"),
    naive_results  %>% mutate(method = "Naive"),
    np_results     %>% mutate(method = "NP")
  ) %>%
    mutate(method = factor(method, levels = METHOD_LEVELS))
}

# Fraction of the visible y-window reserved above the data for the annotation row
LABEL_HEADROOM <- 0.10

# --- Type I error boxplots, faceted by pi (vertical) ---
# y_lo/y_hi define the data window; labels sit in headroom added above y_hi.
make_box_t1 <- function(all_results, violation_summary, y_lo, y_hi) {
  y_top <- y_hi + (y_hi - y_lo) * LABEL_HEADROOM
  y_lab <- y_hi + (y_hi - y_lo) * (LABEL_HEADROOM / 2)
  ggplot(all_results, aes(x = method, y = pop_t1, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5, width = 0.6) +
    geom_hline(yintercept = TARGET_ALPHA, linetype = "dashed",
               color = "red", linewidth = 0.8) +
    geom_text(data = violation_summary,
              aes(x = method, y = y_lab,
                  label = sprintf("%.1f%%", violation_rate * 100)),
              hjust = 0.5, vjust = 0.5, size = 2.6, fontface = "bold",
              inherit.aes = FALSE) +
    scale_fill_manual(values = METHOD_COLORS) +
    scale_x_discrete(limits = METHOD_LEVELS) +
    facet_wrap(~ pi, labeller = as_labeller(function(x) paste0("pi == ", x),
                                            label_parsed)) +
    coord_cartesian(ylim = c(y_lo, y_top)) +
    labs(x = NULL, y = "Type I Error Rate") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(hjust = 0.5))
}

# --- Type II error boxplots, faceted by pi (vertical) ---
make_box_t2 <- function(all_results, type2_summary, y_lo, y_hi) {
  y_top <- y_hi + (y_hi - y_lo) * LABEL_HEADROOM
  y_lab <- y_hi + (y_hi - y_lo) * (LABEL_HEADROOM / 2)
  ggplot(all_results, aes(x = method, y = pop_t2, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5, width = 0.6) +
    geom_hline(yintercept = bl_type2_theo, linetype = "dashed",
               color = "black", linewidth = 0.8) +
    geom_text(data = type2_summary,
              aes(x = method, y = y_lab,
                  label = sprintf("%.3f", mean_t2)),
              hjust = 0.5, vjust = 0.5, size = 2.6, fontface = "bold",
              inherit.aes = FALSE) +
    scale_fill_manual(values = METHOD_COLORS) +
    scale_x_discrete(limits = METHOD_LEVELS) +
    facet_wrap(~ pi, labeller = as_labeller(function(x) paste0("pi == ", x),
                                            label_parsed)) +
    coord_cartesian(ylim = c(y_lo, y_top)) +
    labs(x = NULL, y = "Type II Error Rate") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(hjust = 0.5))
}

# =============================================================================
# BUILD
# =============================================================================

local({
  load("data/experiment_1.RData")

  all_results <- make_all_results(costnp_plus_results, costnp_results,
                                  naive_results, np_results)

  violation_summary <- all_results %>%
    group_by(method, pi) %>%
    summarise(violation_rate = mean(pop_t1 > TARGET_ALPHA, na.rm = TRUE),
              .groups = "drop")

  type2_summary <- all_results %>%
    group_by(method, pi) %>%
    summarise(mean_t2 = mean(pop_t2, na.rm = TRUE), .groups = "drop")

  # Data windows (trimmed to the bulk so outliers don't drive the range).
  # Type I keeps 0 as a meaningful floor; Type II zooms to a padded data band.
  y_lo_t1 <- 0
  y_hi_t1 <- quantile(all_results$pop_t1, 0.995, na.rm = TRUE)

  t2_pad  <- 0.02
  y_lo_t2 <- max(0, quantile(all_results$pop_t2, 0.005, na.rm = TRUE) - t2_pad)
  y_hi_t2 <- quantile(all_results$pop_t2, 0.995, na.rm = TRUE)

  p_t1 <- make_box_t1(all_results, violation_summary, y_lo_t1, y_hi_t1)
  p_t2 <- make_box_t2(all_results, type2_summary,     y_lo_t2, y_hi_t2)

  print(p_t1)
  print(p_t2)

  ggsave("images/experiment_1_type1_box.pdf", plot = p_t1, width = 10, height = 5, dpi = 120)
  ggsave("images/experiment_1_type2_box.pdf", plot = p_t2, width = 10, height = 5, dpi = 120)
})
