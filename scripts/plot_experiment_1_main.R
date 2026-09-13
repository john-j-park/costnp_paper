# =============================================================================
# Main Plots for Experiment 1
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

TARGET_ALPHA <- 0.1
bl_type2_theo <- 0.1943861

METHOD_COLORS <- c(
  Naive      = "#D9534F",
  NP         = "#A67C00",
  `CostNP+`  = "#0072B2",
  CostNP     = "#6BAED6"
)

make_all_results <- function(sim_results, conservative_sim_results,
                             naive_sim_results, np_results) {
  bind_rows(
    sim_results              %>% mutate(method = "CostNP+"),
    conservative_sim_results %>% mutate(method = "CostNP"),
    naive_sim_results        %>% mutate(method = "Naive"),
    np_results               %>% mutate(method = "NP")
  ) %>%
    mutate(method = factor(method, levels = c("Naive", "NP", "CostNP+", "CostNP")))
}


plot_violation_tile <- function(violation_summary, split_var) {
  df <- violation_summary %>%
    mutate(
      label     = sprintf("%.1f%%", violation_rate * 100),
      exceeds   = violation_rate > TARGET_ALPHA,
      col_value = pmin(pmax(violation_rate, 0), 0.3)
    )

  x_labeller <- if (split_var == "pi") {
    function(x) parse(text = paste0("pi == ", x))
  } else {
    function(x) paste0("n = ", x)
  }

  ggplot(df, aes(x = factor(.data[[split_var]]), y = method)) +
    geom_tile(fill = "white", color = "grey60", linewidth = 0.6) +
    geom_text(aes(label = label), color = "black", size = 3.7, fontface = "bold") +
    scale_x_discrete(labels = x_labeller) +
    scale_y_discrete(limits = rev(levels(df$method))) +
    labs(x = NULL, y = NULL, title = "Violation Rate") +
    theme_minimal() +
    theme(
      panel.grid    = element_blank(),
      plot.title    = element_text(hjust = 0.5, size = 11),
      axis.text.x   = element_text(size = 11),
      axis.text.y   = element_text(size = 11),
      panel.border  = element_blank()
    )
}

plot_line_summary <- function(all_results, x_var, x_label, title_suffix) {
  line_summary <- all_results %>%
    group_by(method, .data[[x_var]]) %>%
    summarise(
      mean_t1 = mean(pop_t1, na.rm = TRUE),
      mean_t2 = mean(pop_t2, na.rm = TRUE),
      .groups = "drop"
    )

  n_methods <- length(METHOD_COLORS)

  # Shared color scale for both panels so their legends are identical and
  # patchwork can collect them into a single legend. It lists BOTH dashed
  # reference lines even though each panel only draws its own.
  legend_values <- c(METHOD_COLORS,
                     "Target Type I Error Rate" = "red",
                     "Best linear"              = "black")
  legend_breaks <- c(names(METHOD_COLORS), "Target Type I Error Rate", "Best linear")
  method_overrides <- list(
    linetype  = c(rep("solid", n_methods), "2323", "2323"),
    shape     = c(rep(16, n_methods), NA, NA),
    linewidth = c(rep(0.9, n_methods), 0.8, 0.8)
  )
  shared_color_scale <- scale_color_manual(
    values = legend_values,
    limits = legend_breaks,  # force ALL entries into each panel's legend (not just those present in its data)
    guide  = guide_legend(title = NULL, nrow = 1, override.aes = method_overrides)
  )

  p_t1 <- ggplot(line_summary, aes(x = .data[[x_var]], y = mean_t1, color = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_hline(aes(yintercept = TARGET_ALPHA, color = "Target Type I Error Rate"),
               linetype = "dashed", linewidth = 0.8) +
    shared_color_scale +
    scale_x_continuous(breaks = sort(unique(line_summary[[x_var]]))) +
    labs(x = x_label, y = "Average Type I Error Rate") +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), plot.title = element_text(hjust = 0.5), legend.position = "bottom")

  p_t2 <- ggplot(line_summary, aes(x = .data[[x_var]], y = mean_t2, color = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_hline(aes(yintercept = bl_type2_theo, color = "Best linear"),
               linetype = "dashed", linewidth = 0.8) +
    shared_color_scale +
    scale_x_continuous(breaks = sort(unique(line_summary[[x_var]]))) +
    labs(x = x_label, y = "Average Type II Error Rate") +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), plot.title = element_text(hjust = 0.5), legend.position = "bottom")

  list(t1 = p_t1, t2 = p_t2)
}

# =============================================================================
# BY-PI DATASET
# =============================================================================

local({
  load("data/experiment_1.RData")

  all_results <- make_all_results(costnp_plus_results, costnp_results,
                                  naive_results, np_results)

  violation_summary <- all_results %>%
    group_by(method, pi) %>%
    summarise(violation_rate = mean(pop_t1 > TARGET_ALPHA, na.rm = TRUE),
              .groups = "drop")

  violation_tile <- plot_violation_tile(violation_summary, "pi")

  line_plots     <- plot_line_summary(all_results, "pi", expression("Class 0 Proportion (" * pi * ")"), "Class Balance")
  combined_lines <- (line_plots$t1 | line_plots$t2) / guide_area() / violation_tile +
    plot_layout(heights = c(3, 0.3, 1), guides = "collect") &
    theme(legend.position   = "bottom",
          legend.direction  = "horizontal",
          legend.text       = element_text(size = 12),
          legend.key.width  = unit(1.0, "cm"))

  print(combined_lines)

  ggsave("images/experiment_1_line.pdf", plot = combined_lines, width = 10, height = 6, dpi = 100)

  # Slide variant: line panels only. The paper figure keeps the violation-rate
  # table, but at projection size its numbers are too small to read, so the deck
  # uses this version and the violation rates are given verbally instead.
  # Height drops from 6 to 4.6 = 6 * 3.3/4.3, the share of the layout that the
  # two line panels and the legend strip occupied.
  slide_lines <- (line_plots$t1 | line_plots$t2) / guide_area() +
    plot_layout(heights = c(3, 0.3), guides = "collect") &
    theme(legend.position   = "bottom",
          legend.direction  = "horizontal",
          legend.text       = element_text(size = 12),
          legend.key.width  = unit(1.0, "cm"))

  ggsave("images/experiment_1_line_slides.pdf", plot = slide_lines,
         width = 10, height = 4.6, dpi = 100)
})
