# =============================================================================
# Intro Figure Generation
#
# Produces two figures:
#   images/np_planet_data.pdf    (plot_npc_boundary)
#   images/calcs_planet_data.pdf (plot_cost_grid)
# =============================================================================

library(costnp)
library(nproc)
library(ggplot2)
source("scripts/paper_utils.R")

TARGET_ALPHA <- 0.1
TARGET_DELTA <- 0.1
DIMENSIONS   <- 2

set.seed(7)
lr_boundary_df <- function(fit_obj, data, margin = 1) {
  if (inherits(fit_obj, "npc")) {
    lr_model     <- fit_obj$fits[[1]]$fit
    lr_threshold <- fit_obj$fits[[1]]$cutoff
  } else {
    lr_model     <- fit_obj$classifier
    lr_threshold <- 0.5
  }
  beta    <- coef(lr_model)
  logit_t <- log(lr_threshold / (1 - lr_threshold))
  x1      <- seq(min(data$x[, 1]) - margin, max(data$x[, 1]) + margin, length.out = 400)
  data.frame(x1 = x1, x2 = (logit_t - beta[1] - beta[2] * x1) / beta[3])
}

best_linear_params <- function(mu0, Sigma0, mu1, Sigma1, alpha, n_grid = 3600) {
  z_alpha    <- qnorm(1 - alpha)
  delta      <- mu1 - mu0
  theta_grid <- seq(0, 2 * pi, length.out = n_grid + 1)[-(n_grid + 1)]
  scores     <- sapply(theta_grid, function(theta) {
    w   <- c(cos(theta), sin(theta))
    sd0 <- sqrt(as.numeric(t(w) %*% Sigma0 %*% w))
    sd1 <- sqrt(as.numeric(t(w) %*% Sigma1 %*% w))
    (sum(w * delta) - z_alpha * sd0) / sd1
  })
  theta_opt <- theta_grid[which.max(scores)]
  w_opt     <- c(cos(theta_opt), sin(theta_opt))
  sd0_opt   <- sqrt(as.numeric(t(w_opt) %*% Sigma0 %*% w_opt))
  c_opt     <- sum(w_opt * mu0) + z_alpha * sd0_opt
  list(w = w_opt, c = c_opt, power = pnorm(max(scores)))
}

best_linear_df <- function(mu0, Sigma0, mu1, Sigma1, alpha, data, n_grid = 3600, margin = 1) {
  params <- best_linear_params(mu0, Sigma0, mu1, Sigma1, alpha, n_grid)
  x1     <- seq(min(data$x[, 1]) - margin, max(data$x[, 1]) + margin, length.out = 400)
  data.frame(x1 = x1, x2 = (params$c - params$w[1] * x1) / params$w[2])
}

plot_npc_boundary <- function(nproc_fit, data, mu0, Sigma0, mu1, Sigma1, save_path = NULL) {
  npc_bnd  <- lr_boundary_df(nproc_fit, data, margin = 20)
  bl_bnd   <- best_linear_df(mu0, Sigma0, mu1, Sigma1, TARGET_ALPHA, data, margin = 20)

  beta_npc <- coef(nproc_fit$fits[[1]]$fit)
  x1_seq   <- seq(min(data$x[, 1]) - 1, max(data$x[, 1]) + 1, length.out = 400)
  x2_range <- range(data$x[, 2])
  corners  <- expand.grid(x1 = range(x1_seq), x2 = x2_range)
  c_scores <- beta_npc[1] + beta_npc[2] * corners$x1 + beta_npc[3] * corners$x2
  c_vals   <- seq(min(c_scores), max(c_scores), length.out = 30)
  pal      <- colorRampPalette(c("#A67C00", "lightgray"))(length(c_vals))
  cand_col <- pal[ceiling(length(pal) / 2)]  # mid-gradient swatch for the legend
  grad_lines <- do.call(rbind, lapply(seq_along(c_vals), function(i) {
    data.frame(x1  = x1_seq,
               x2  = (c_vals[i] - beta_npc[1] - beta_npc[2] * x1_seq) / beta_npc[3],
               grp = i,
               col = pal[i])
  }))

  logit_t <- log(nproc_fit$fits[[1]]$cutoff / (1 - nproc_fit$fits[[1]]$cutoff))
  slope   <- unname(beta_npc[2:3])
  int     <- unname(beta_npc[1])
  sd0     <- sqrt(as.numeric(t(slope) %*% Sigma0 %*% slope))
  sd1     <- sqrt(as.numeric(t(slope) %*% Sigma1 %*% slope))
  print(1 - pnorm((logit_t - int - sum(slope * mu0)) / sd0))  # npc_type1_theo
  print(pnorm((logit_t - int - sum(slope * mu1)) / sd1))       # npc_type2_theo

  scatter   <- data.frame(x1  = data$x[, 1], x2  = data$x[, 2],
                          col = ifelse(data$y == 0, "#D9534F", "#00BFC4"))
  bnd_lines <- rbind(data.frame(npc_bnd, grp = "NPC",         col = "#A67C00", lty = "solid"),
                     data.frame(bl_bnd,  grp = "best linear", col = "black",   lty = "dashed"))

  p <- ggplot() +
    geom_line(data = grad_lines,
              aes(x = x1, y = x2, group = grp, color = col), linewidth = 0.5) +
    geom_point(data = scatter,
               aes(x = x1, y = x2, color = col), shape = 19, size = 0.8) +
    geom_line(data = bnd_lines,
              aes(x = x1, y = x2, group = grp, color = col, linetype = lty), linewidth = 1.2) +
    scale_linetype_identity(guide = "none") +
    scale_color_identity(
      guide  = guide_legend(
        title = NULL,
        override.aes = list(shape     = c(19,      19,      NA,      NA,      NA),
                            linetype  = c("blank", "blank", "solid", "2122",  "solid"),
                            linewidth = c(NA,      NA,      0.5,     1.2,     1.2))
      ),
      breaks = c("#D9534F", "#00BFC4", cand_col, "black", "#A67C00"),
      labels = c("Class 0 Points", "Class 1 Points", "Candidate classifiers", "Best linear boundary", "NP boundary")
    ) +
    coord_cartesian(xlim = range(data$x[, 1]), ylim = range(data$x[, 2])) +
    theme_minimal() +
    theme(panel.border         = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid           = element_blank(),
          axis.title           = element_blank(),
          legend.position      = c(1, 0),
          legend.justification = c(1, 0),
          legend.background    = element_rect(fill = "white", colour = "black", linewidth = 0.5),
          legend.key           = element_rect(fill = "white", colour = NA))
  if (!is.null(save_path)) ggsave(save_path, plot = p, width = 5, height = 4, dpi = 100)
  p
}

plot_cost_grid <- function(data, res, mu0, Sigma0, mu1, Sigma1, classify_fn, save_path = NULL) {
  costvec     <- seq(0.99, 0.01, by = -abs(0.05))
  idx_0       <- which(data$y == 0)
  train_idx_0 <- sample(idx_0, floor(0.5 * length(idx_0)))
  train_idx   <- (data$y == 1) | (seq_along(data$y) %in% train_idx_0)
  xtrain_vis  <- data$x[train_idx, ]
  ytrain_vis  <- data$y[train_idx]
  x1_seq      <- seq(min(data$x[, 1]) - 2, max(data$x[, 1]) + 2, length.out = 400)
  pal         <- colorRampPalette(c("lightgray", "#4292C6"))(length(costvec))
  cand_col    <- pal[ceiling(length(pal) / 2)]  # mid-gradient swatch for the legend

  cost_lines <- do.call(rbind, lapply(seq_along(costvec), function(i) {
    b <- coef(classify_fn(xtrain_vis, ytrain_vis, cost = costvec[i],
                          xnew = xtrain_vis[1, , drop = FALSE], method = "LR")$obj)
    data.frame(x1 = x1_seq, x2 = (-b[1] - b[2] * x1_seq) / b[3], grp = i, col = pal[i])
  }))

  chosen_line     <- lr_boundary_df(res, data, margin = 2)
  chosen_line$grp <- "chosen"
  chosen_line$col <- "#0072B2"

  bl_bnd     <- best_linear_df(mu0, Sigma0, mu1, Sigma1, TARGET_ALPHA, data)
  bl_bnd$col <- "black"
  scatter    <- data.frame(x1  = data$x[, 1], x2  = data$x[, 2],
                           col = ifelse(data$y == 0, "#D9534F", "#00BFC4"))

  p <- ggplot() +
    geom_line(data = cost_lines,
              aes(x = x1, y = x2, group = grp, color = col), linewidth = 0.5) +
    geom_point(data = scatter,
               aes(x = x1, y = x2, color = col), shape = 19, size = 0.8) +
    geom_line(data = chosen_line,
              aes(x = x1, y = x2, group = grp, color = col), linewidth = 1.2) +
    geom_line(data = bl_bnd,
              aes(x = x1, y = x2, color = col), linetype = "dashed", linewidth = 1.2) +
    coord_cartesian(xlim = range(data$x[, 1]), ylim = range(data$x[, 2])) +
    scale_color_identity(
      guide  = guide_legend(
        title = NULL,
        override.aes = list(shape     = c(19,      19,      NA,      NA,      NA),
                            linetype  = c("blank", "blank", "solid", "2122",  "solid"),
                            linewidth = c(NA,      NA,      0.5,     1.2,     1.2))
      ),
      breaks = c("#D9534F", "#00BFC4", cand_col, "black", "#0072B2"),
      labels = c("Class 0 Points", "Class 1 Points", "Candidate classifiers", "Best linear boundary", "CS boundary")
    ) +
    theme_minimal() +
    theme(panel.border         = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid           = element_blank(),
          axis.title           = element_blank(),
          legend.position      = c(1, 0),
          legend.justification = c(1, 0),
          legend.background    = element_rect(fill = "white", colour = "black", linewidth = 0.5),
          legend.key           = element_rect(fill = "white", colour = NA))
  if (!is.null(save_path)) ggsave(save_path, plot = p, width = 5, height = 4, dpi = 100)
  p
}

# ── Elliptical model ──────────────────────────────────────────────────────────────

PLANET_MU0    <- c(2.5, 0)
PLANET_SIGMA0 <- matrix(c(3, 1.5, 1.5, 1), nrow = 2)
PLANET_MU1    <- c(0, 0)
PLANET_SIGMA1 <- diag(2)

# NOTE: this draw is retained (though its figure is no longer produced) so the
# RNG stream feeding data/pop/nproc_fit/res and plot_cost_grid stays identical.
plot_data <- gen_data(model = "planet", n = 3000, d = DIMENSIONS, pi = 0.5)

data <- gen_data(model = "planet", n = 1000, d = DIMENSIONS, pi = 0.1)
pop  <- gen_data(model = "planet", n = 500000, d = DIMENSIONS, pi = 0.5)
nproc_fit <- npc(x = data$x, y = data$y, alpha = 0.1, delta = 0.1, model = "logistic")

plot_npc_boundary(nproc_fit, data, PLANET_MU0, PLANET_SIGMA0, PLANET_MU1, PLANET_SIGMA1, save_path = "images/np_planet_data.pdf")

res <- costnp_plus(
  data, population_data = pop, alpha = TARGET_ALPHA, method = "LR",
  cost_start = 0.99, cost_end = 0.5, cost_by = 0.02
)

plot_cost_grid(data, res, PLANET_MU0, PLANET_SIGMA0, PLANET_MU1, PLANET_SIGMA1, classify_fn = classify_fun_stratified, save_path = "images/calcs_planet_data.pdf")
