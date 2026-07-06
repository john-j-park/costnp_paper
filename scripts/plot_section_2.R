# =============================================================================
# Section 2 Figure Generation
# =============================================================================

library(ggplot2)
library(data.table)

# --- Beta Distribution Analysis ---
n_simulations <- 50000
mu_0 <- -1.5
sigma_0 <- 1.0
target_error <- 0.1

n_emp_values <- c(50, 500)

for (n_emp in n_emp_values) {

  population_errors <- numeric(n_simulations)

  for (i in 1:n_simulations) {
    x_emp_0 <- rnorm(n_emp, mean = mu_0, sd = sigma_0)
    k <- ceiling(n_emp * (1 - target_error))
    sorted_emp <- sort(x_emp_0)
    chosen_c <- sorted_emp[k]
    true_error <- pnorm(chosen_c, mean = mu_0, sd = sigma_0, lower.tail = FALSE)
    population_errors[i] <- true_error
  }

  results_df <- data.frame(True_Error = population_errors)

  k <- ceiling(n_emp * (1 - target_error))
  print(k)
  alpha_beta <- n_emp - k + 1
  beta_beta <- k
  median <- qbeta(0.5, alpha_beta, beta_beta)

  pct_over <- mean(population_errors > target_error) * 100
  pct_label <- paste0(sprintf("%.1f", pct_over), "%")

  x_pos <- max(population_errors) * 0.95
  y_max <- max(dbeta(seq(0.01, 0.99, length.out = 500), alpha_beta, beta_beta))
  y_pos <- y_max * 0.05

  p <- ggplot(results_df, aes(x = True_Error)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "steelblue", color = "black", alpha = 0.7) +
    stat_function(fun = dbeta, args = list(shape1 = alpha_beta, shape2 = beta_beta),
                  color = "darkorange", linewidth = 1.2) +
    geom_vline(xintercept = target_error, color = "darkred", linetype = "dashed", linewidth = 1.2) +
    geom_vline(xintercept = alpha_beta / (alpha_beta + beta_beta), color = "darkorange", linetype = "dashed", linewidth = 1.2) +
    annotate("text", x = x_pos, y = y_pos, label = pct_label, size = 4, fontface = "bold", hjust = 1) +
    labs(
      x = "Population Type I Error",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  print(p)

  ggsave(paste0("images/beta_distribution_plot_", n_emp, ".pdf"), plot = p, width = 5, height = 5)
}

# --- Theoretical Mean Heatmap ---
mu_0 <- -0.5
mu_1 <- 0.5
sigma <- 1
midpoint <- (mu_0 + mu_1) / 2
delta <- mu_1 - mu_0
n <- 50

alpha_grid <- seq(0.001, 0.25, length = 1000)
step_grid  <- seq(0.001, 0.15, length = 1000)
dt <- as.data.table(expand.grid(alpha_star = alpha_grid, step_size = step_grid))

dt[, expected_R0 := {
  s <- step_size[1]
  cost_grid <- seq(s, 1 - s, by = s)
  t_grid <- midpoint + sigma^2 * log(cost_grid / (1 - cost_grid)) / delta
  p_grid <- pnorm(t_grid, mean = mu_0, sd = sigma, lower.tail = FALSE)

  ks <- floor(alpha_star * n)

  cdf_matrix <- outer(ks, p_grid, function(k, p) pbinom(k, size = n, prob = p))

  pmf_matrix <- cdf_matrix
  if (ncol(cdf_matrix) > 1) {
    pmf_matrix[, 2:ncol(pmf_matrix)] <- cdf_matrix[, 2:ncol(cdf_matrix)] - cdf_matrix[, 1:(ncol(cdf_matrix)-1)]
  }

  as.vector(pmf_matrix %*% p_grid)

}, by = step_size]

dt[, bias := expected_R0 - alpha_star]
dt$over <- ifelse(dt$bias >= 0, "Aggressive", "Conservative")

p1 <- ggplot(dt, aes(x = step_size, y = alpha_star, fill = over)) +
  geom_raster() +
  scale_fill_manual(values = c("Aggressive" = "darkred", "Conservative" = "steelblue")) +
  labs(
    x = "Grid step size",
    y = expression(alpha ~ "(target Type I error)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA)
  )
p1

ggsave("images/TheoreticalMeanHeatmap50.pdf", plot = p1, width = 5, height = 4)
