# =============================================================================
# Helper functions required for the other scripts
# =============================================================================
library(MASS)
library(mvtnorm)

# =============================================================================
# DATA GENERATION
# =============================================================================

gen_data <- function(model, seed = NULL, n, d, pi = 0.5, nleave = 0, ribbon_angle = 3.14/6) {
  if (!is.null(seed)) set.seed(seed)

  modellist <- c("gaussian", "t", "mixture","planet","planet_2","homgaussian")
  if (!(model %in% modellist)) stop("Unrecognized data model!")

  if (model == "gaussian") {
    s <- 5
    mu0 <- c(rep(0,s), rep(0, d - s))
    mu1 <- c(rep(1,s), rep(0, d - s))
    cov <- diag(d)
    #cov[abs(row(cov) - col(cov)) == 1] <- 0.5
    cov[abs(row(cov) - col(cov)) == 1] <- 0.5
    sigma1 <- cov
    sigma0 <- diag(d)

    y <- sample(0:1, size = n, replace = TRUE, prob = c(pi, 1 - pi))
    n1 <- sum(y)
    n0 <- n - n1 + nleave

    if (d != length(mu0) | d != nrow(sigma0)) stop("Dimension mismatch!")

    S0 <- list(x = MASS::mvrnorm(n = n0, mu0, sigma0), y = rep(0, n0))
    S1 <- list(x = MASS::mvrnorm(n = n1, mu1, sigma1), y = rep(1, n1))
  }

  if (model == "homgaussian") {
    p0 <- min(5, d)
    rho <- 0.5
    signal <- 0.6

    Sigma <- outer(
      1:d, 1:d,
      FUN = function(i, j) rho^abs(i - j)
    )

    beta_bayes <- signal * c(rep(1, p0), rep(0, d - p0))
    mu0 <- rep(0, d)
    mu1 <- as.vector(Sigma %*% beta_bayes)
    sigma0 <- Sigma
    sigma1 <- Sigma

    y <- sample(0:1, size = n, replace = TRUE, prob = c(pi, 1 - pi))
    n1 <- sum(y)
    n0 <- n - n1 + nleave

    if (d != length(mu0) | d != nrow(sigma0)) stop("Dimension mismatch!")

    S0 <- list(x = MASS::mvrnorm(n = n0, mu0, sigma0), y = rep(0, n0))
    S1 <- list(x = MASS::mvrnorm(n = n1, mu1, sigma1), y = rep(1, n1))
  }

  if (model == "t") {
    if (d < 3) stop("d must be >= 3 for t distribution!")

    s <- 2
    mu0 <- rep(0,s)
    mu1 <- rep(2,s)

    y <- sample(0:1, size = n, replace = TRUE, prob = c(pi, 1 - pi))
    n1 <- sum(y)
    n0 <- n - n1 + nleave

    x01 <- mvtnorm::rmvt(n0, sigma = diag(s), df = 3, delta = mu0, type = "shifted")
    x02 <- matrix(rnorm(n0 * (d - s)), nrow = n0)
    x0 <- cbind(x01, x02)

    x11 <- mvtnorm::rmvt(n1, sigma = diag(s), df = 3, delta = mu1, type = "shifted")
    x12 <- matrix(rnorm(n1 * (d - s)), nrow = n1)
    x1 <- cbind(x11, x12)

    S0 <- list(x = x0, y = rep(0, n0))
    S1 <- list(x = x1, y = rep(1, n1))
  }

  if (model == "mixture") {
    a <- 2 / sqrt(d)
    mu01 <- rep(a, d)
    mu02 <- rep(-a, d)
    mu1 <- rep_len(c(a, -a), length.out = d)

    nn <- rmultinom(1, n, c(pi, 1 - pi))
    n0 <- rmultinom(1, nn[1, 1] + nleave, c(1/2, 1/2))
    n1 <- nn[2, 1]

    X1 <- rbind(t(matrix(mu01, d, n0[1])), t(matrix(mu02, d, n0[2]))) +
          matrix(rnorm(sum(n0) * d), sum(n0), d)
    X2 <- t(matrix(mu1, d, n1)) + matrix(rnorm(n1 * d), n1, d)

    S0 <- list(x = X1, y = rep(0, sum(n0)))
    S1 <- list(x = X2, y = rep(1, n1))
  }

  if(model == "planet"){
    y <- sample(0:1, size = n, replace = TRUE, prob = c(pi, 1 - pi))
    n1 <- sum(y)
    n1 <- n1 + n1 %% 2
    n0 <- n - n1 + nleave

    class1_center = c(0.0, 0.0)
    class1_cov = matrix(c(1.0, 0.0, 0.0, 1.0), nrow = 2)
    ribbon_center = c(2.5, 0.0)
    ribbon_angle = 3.14/ 6
    ribbon_length = 4.0
    ribbon_width = 0.15

    R <- matrix(c(cos(ribbon_angle), -sin(ribbon_angle),
                  sin(ribbon_angle),  cos(ribbon_angle)), nrow = 2, byrow = TRUE)
    D <- diag(c(ribbon_length, ribbon_width))
    ribbon_cov <- R %*% D %*% t(R)

    ribbon_cov <- matrix(c(3,1.5,1.5,1),nrow=2)

    X1 <- MASS::mvrnorm(n1, mu = class1_center,  Sigma = class1_cov)
    X0 <- MASS::mvrnorm(n0, mu = ribbon_center,  Sigma = ribbon_cov)

    X <- rbind(X0, X1)
    y <- c(rep(0L, n0), rep(1L, n1))

    perm <- sample.int(length(y))
    return(list(x = X[perm, , drop = FALSE], y = y[perm]))
  }

  # Elliptical model: 2 signal dimensions plus d - 2 independent N(0, 1) noise dimensions
  if (model == "planet_2") {
    if (d < 2) stop("d must be >= 2 for planet_2!")

    nn <- rmultinom(1, n, c(pi, 1 - pi))
    n0 <- nn[1, 1] + nleave
    n1 <- nn[2, 1]

    class1_center <- c(0.0, 0.0)
    class1_cov    <- diag(2)

    ribbon_center <- c(2,2)
    ribbon_cov    <- matrix(c(3, 1.5, 1.5, 1), nrow = 2)

    X1_signal <- MASS::mvrnorm(n1, mu = class1_center, Sigma = class1_cov)
    X0_signal <- MASS::mvrnorm(n0, mu = ribbon_center, Sigma = ribbon_cov)

    if (d > 2) {
      X1 <- cbind(X1_signal, matrix(rnorm(n1 * (d - 2)), nrow = n1))
      X0 <- cbind(X0_signal, matrix(rnorm(n0 * (d - 2)), nrow = n0))
    } else {
      X1 <- X1_signal
      X0 <- X0_signal
    }

    S0 <- list(x = X0, y = rep(0L, n0))
    S1 <- list(x = X1, y = rep(1L, n1))
  }

  x <- rbind(S0$x, S1$x)
  y <- c(S0$y, S1$y)
  perm <- sample.int(length(y))
  return(list(x = x[perm, , drop = FALSE], y = y[perm]))
}

# =============================================================================
# NPC (Neyman-Pearson Classifier) - comparison method (for using NP with XGB)
# =============================================================================

find_np_order_statistic <- function(n, alpha, delta) {
  k_seq <- 1:n
  quantiles <- qbeta(1 - delta, shape1 = n - k_seq + 1, shape2 = k_seq)
  valid_k <- which(quantiles <= alpha)
  if (length(valid_k) == 0) {
    min_n <- ceiling(log(delta) / log(1 - alpha))
    stop(sprintf("Calibration set too small. For alpha=%g and delta=%g, you need at least n=%d.",
                 alpha, delta, min_n))
  }
  return(min(valid_k))
}

npc_manual <- function(data, alpha = 0.1, population_data, delta = 0.1,
                      data_pct = 0.5, method = "LR", strat_pct = 0.5,
                      seeded = FALSE) {
  if (seeded) set.seed(1)
  idx_0 <- which(data$y == 0)
  train_idx_0 <- sample(idx_0, floor(data_pct * length(idx_0)))
  train_idx <- (data$y == 1) | (seq_along(data$y) %in% train_idx_0)
  xtrain <- data$x[train_idx, , drop = FALSE]
  ytrain <- data$y[train_idx]
  xcalib <- data$x[!train_idx, , drop = FALSE]
  n_calib <- nrow(xcalib)

  base_result <- costnp:::classify_base(xtrain, ytrain, xnew = xcalib, method = method)

  k <- find_np_order_statistic(n_calib, alpha, delta)
  thresh <- sort(base_result$prob)[k]

  pop_is_0 <- population_data$y == 0
  pop_score_on_0 <- costnp:::score_base(base_result$obj,
                               population_data$x[pop_is_0, , drop = FALSE], method)
  pop_type1 <- sum(pop_score_on_0 > thresh) / sum(pop_is_0)

  pop_is_1 <- population_data$y == 1
  pop_score_on_1 <- costnp:::score_base(base_result$obj,
                               population_data$x[pop_is_1, , drop = FALSE], method)
  pop_type2 <- sum(pop_score_on_1 <= thresh) / sum(pop_is_1)

  return(list(
    population_type1 = pop_type1,
    population_type2 = pop_type2,
    n_calib = n_calib,
    classifier = base_result$obj,
    threshold = thresh,
    method = method
  ))
}

# =============================================================================
# eLDA / feLDA NP CLASSIFIERS - comparison methods
# =============================================================================

estimate_lda_params <- function(x, y) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (length(unique(y)) != 2 || !all(sort(unique(y)) == c(0, 1)))
    stop("y must be a 0/1 vector with both classes present.")
  x0 <- x[y == 0, , drop = FALSE]
  x1 <- x[y == 1, , drop = FALSE]
  n0 <- nrow(x0); n1 <- nrow(x1); n <- n0 + n1; p <- ncol(x)
  if (p >= n - 1) stop(sprintf("eLDA requires p < n-1 (got p=%d, n=%d).", p, n))
  mu0 <- colMeans(x0); mu1 <- colMeans(x1); mu_d <- mu1 - mu0
  Sigma_hat <- ((n0 - 1) * cov(x0) + (n1 - 1) * cov(x1)) / (n - 2)
  A <- tryCatch(
    solve(Sigma_hat, mu_d),
    error = function(e) {
      # Singular covariance (e.g. zero-variance features in a split); add ridge
      eps <- 1e-6 * mean(diag(Sigma_hat))
      solve(Sigma_hat + eps * diag(p), mu_d)
    }
  )
  list(mu0 = mu0, mu1 = mu1, mu_d = mu_d, Sigma_hat = Sigma_hat, A = A,
       n0 = n0, n1 = n1, n = n, p = p, r = p / n,
       Ahat_Sigma_Ahat = as.numeric(t(A) %*% Sigma_hat %*% A),
       Ahat_mu0 = as.numeric(sum(A * mu0)),
       mahalanobis_d = as.numeric(sum(A * mu_d)))
}

elda_score <- function(x, params) {
  if (!is.matrix(x)) x <- as.matrix(x)
  as.numeric(x %*% params$A)
}

felda_threshold <- function(params, alpha = 0.05, delta = 0.1) {
  with(params, {
    Phi_a <- qnorm(1 - alpha); Phi_d <- qnorm(1 - delta)
    V_tilde <- Phi_a^2 / 2 + n / n0
    F_tilde <- sqrt(Ahat_Sigma_Ahat) * Phi_a + Ahat_mu0
    C_tilde <- F_tilde + sqrt(Ahat_Sigma_Ahat * V_tilde / n) * Phi_d
    list(threshold = C_tilde, F_hat = F_tilde, V_hat = V_tilde,
         alpha = alpha, delta = delta, variant = "feLDA")
  })
}

elda_threshold <- function(params, alpha = 0.05, delta = 0.1) {
  with(params, {
    if (r >= 1) stop("eLDA requires r = p/n < 1.")
    Phi_a <- qnorm(1 - alpha); Phi_d <- qnorm(1 - delta)
    v1n <- n / n0 + n / n1
    q <- (1 - r) * Ahat_Sigma_Ahat - r * v1n
    if (q <= 0) {
      warning(sprintf("eLDA: (1-r) A'Sigma_hat A - r ||v1||^2 = %.4g is non-positive; falling back to feLDA threshold.", q))
      return(felda_threshold(params, alpha = alpha, delta = delta))
    }
    C <- (1 - r) / (2 * sqrt(mahalanobis_d))
    F_hat <- sqrt(Ahat_Sigma_Ahat / (1 - r)) * Phi_a + Ahat_mu0 + (n / n0) * r / (1 - r)
    V1 <- q * C^2 * Phi_a^2 * 2 * (1 + r) / (1 - r)^7
    V2 <- C^2 * Phi_a^2 * v1n * 4 * r * (1 + r) / (1 - r)^7 +
          n / (n0 * (1 - r)^3) +
          2 * C * Phi_a * sqrt(v1n) * sqrt(n1 / n0) * 2 * r / (1 - r)^5
    V3 <- (v1n / q) * (C^2 * Phi_a^2 * v1n * 2 * r^2 * (1 + r) / (1 - r)^7 +
           (n + n1) * r / (n0 * (1 - r)^3) +
           2 * C * Phi_a * sqrt(v1n) * sqrt(n1 / n0) * 2 * r^2 / (1 - r)^5)
    V_hat <- V1 + V2 + V3
    C_alpha <- F_hat + sqrt(q * V_hat / n) * Phi_d
    list(threshold = C_alpha, F_hat = F_hat, V_hat = V_hat, q = q,
         v1_norm_sq = v1n, alpha = alpha, delta = delta, variant = "eLDA")
  })
}

elda_npc <- function(data, alpha = 0.1, population_data, delta = 0.1,
                     variant = c("eLDA", "feLDA"), seeded = FALSE) {
  if (seeded) set.seed(1)
  variant <- match.arg(variant)
  params <- estimate_lda_params(data$x, data$y)
  thr_info <- if (variant == "eLDA") elda_threshold(params, alpha = alpha, delta = delta) else
                                     felda_threshold(params, alpha = alpha, delta = delta)
  thresh <- thr_info$threshold
  pop_is_0 <- population_data$y == 0
  pop_type1 <- sum(elda_score(population_data$x[pop_is_0, , drop = FALSE], params) > thresh) / sum(pop_is_0)
  pop_is_1 <- population_data$y == 1
  pop_type2 <- sum(elda_score(population_data$x[pop_is_1, , drop = FALSE], params) <= thresh) / sum(pop_is_1)
  list(population_type1 = pop_type1, population_type2 = pop_type2,
       n_train = params$n, threshold = thresh, method = variant,
       params = params, threshold_info = thr_info)
}

felda_npc <- function(data, alpha = 0.1, population_data, delta = 0.1, seeded = FALSE) {
  elda_npc(data, alpha = alpha, population_data = population_data,
           delta = delta, variant = "feLDA", seeded = seeded)
}

get_cost_score_matrix <- function(data, cost_by = 0.01, data_pct = 0.5,
                                  method = "LR", classify_fun = classify_fun_stratified,
                                  seeded = FALSE, cost_start = 0.5, cost_end = 0.99) {
  if (seeded) set.seed(1)

  costvec <- seq(cost_start, cost_end, by = cost_by)

  idx_0       <- which(data$y == 0)
  train_idx_0 <- sample(idx_0, floor(data_pct * length(idx_0)))
  train_idx   <- (data$y == 1) | (seq_along(data$y) %in% train_idx_0)
  xtrain      <- data$x[train_idx, , drop = FALSE]
  ytrain      <- data$y[train_idx]

  train_is_0 <- ytrain == 0
  xtrain_0   <- xtrain[train_is_0, , drop = FALSE]
  ytrain_0   <- ytrain[train_is_0]
  n_train_0  <- sum(train_is_0)

  score_mat <- matrix(NA_real_, nrow = n_train_0, ncol = length(costvec))
  m_vec     <- numeric(length(costvec))

  for (i in seq_along(costvec)) {
    result         <- classify_fun(xtrain, ytrain, cost = costvec[i], xnew = xtrain_0, method = method)
    score_mat[, i] <- result$prob
    m_vec[i]       <- sum(result$prob > 0.5)
  }

  list(scores = score_mat, costs = costvec, m = m_vec, y = ytrain_0)
}

# =============================================================================
# NAIVE COST SELECTION - comparison method
# =============================================================================

.npcs_build_cost_grid <- function(cost_start, cost_end, cost_by) {
  step <- if (cost_start > cost_end) -abs(cost_by) else abs(cost_by)
  seq(cost_start, cost_end, by = step)
}

.npcs_eval_population <- function(calib, population_data, chosen_cost, get_threshold, obj = NULL) {
  pop_obj    <- if (is.null(obj)) calib$fit_pop(chosen_cost) else obj
  threshold  <- get_threshold(chosen_cost)
  pop_is_0   <- population_data$y == 0
  pop_prob_0 <- calib$predict_pop(pop_obj, population_data$x[pop_is_0, , drop = FALSE])
  pop_is_1   <- population_data$y == 1
  pop_prob_1 <- calib$predict_pop(pop_obj, population_data$x[pop_is_1, , drop = FALSE])
  list(obj = pop_obj, type1 = mean(pop_prob_0 > threshold), type2 = mean(pop_prob_1 <= threshold))
}

.npcs_setup_calibration <- function(data, method = "LR", data_pct = 0.5,
                                    classify_fun = costnp::classify_fun_stratified,
                                    threshold_fun = function(cost) 0.5) {
  idx_0       <- which(data$y == 0)
  train_idx_0 <- sample(idx_0, floor(data_pct * length(idx_0)))
  train_idx   <- (data$y == 1) | (seq_along(data$y) %in% train_idx_0)
  xtrain      <- data$x[train_idx, , drop = FALSE]
  ytrain      <- data$y[train_idx]
  xcalib      <- data$x[!train_idx, , drop = FALSE]
  n_calib     <- nrow(xcalib)
  get_threshold <- threshold_fun
  fit_pop     <- function(cost) classify_fun(xtrain, ytrain, cost = cost,
                                             xnew = xtrain[1, , drop = FALSE],
                                             method = method)$obj
  get_result  <- function(cost) classify_fun(xtrain, ytrain, cost = cost,
                                             xnew = xcalib, method = method)
  predict_pop <- function(obj, xnew) costnp:::score_base(obj, xnew, method)
  list(get_result = get_result, get_threshold = get_threshold,
       fit_pop = fit_pop, predict_pop = predict_pop, n_calib = n_calib)
}

get_cost_naive <- function(data, alpha = 0.1, population_data, delta = 0.1,
                           cost_start = 0.99, cost_end = 0.5, cost_by = 0.01,
                           data_pct = 0.5, method = "LR",
                           verbose = FALSE,
                           classify_fun = costnp::classify_fun_stratified,
                           threshold_fun = function(cost) 0.5) {
  calib <- .npcs_setup_calibration(data, method = method, data_pct = data_pct,
                                   classify_fun = classify_fun, threshold_fun = threshold_fun)
  get_result    <- calib$get_result
  get_threshold <- calib$get_threshold
  n_calib       <- calib$n_calib
  costvec       <- .npcs_build_cost_grid(cost_start, cost_end, cost_by)
  chosen_cost   <- costvec[1]
  m             <- NA
  alpha_observed <- NA

  for (i in seq_along(costvec)) {
    current_cost   <- costvec[i]
    result         <- get_result(current_cost)
    yhat           <- as.numeric(result$prob > get_threshold(current_cost))
    m              <- sum(yhat == 1)
    alpha_observed <- m / n_calib
    chosen_cost    <- current_cost
    if (verbose) cat("Cost:", chosen_cost, "\n", "Alpha Observed:", alpha_observed, "\n")
    if (alpha_observed <= alpha) break
  }

  pop <- .npcs_eval_population(calib, population_data, chosen_cost, get_threshold)

  if (verbose) {
    cat("Alpha Observed:", alpha_observed, "\n")
    cat("Population Type 1:", pop$type1, "\n")
    cat("Cost:", chosen_cost, "\n")
  }

  list(
    cost              = chosen_cost,
    population_type1  = pop$type1,
    class0_probs      = calib$predict_pop(pop$obj,
                                          population_data$x[population_data$y == 0, , drop = FALSE]),
    n_calib           = n_calib,
    m_errors_in_calib = m,
    population_type2  = pop$type2,
    classifier        = pop$obj
  )
}

# =============================================================================
# REAL DATA: dataset loading
# =============================================================================
load_dataset <- function(DATASET) {
  x_all <- NULL
  y_all <- NULL
  subsample_data <- NULL
  population_data <- NULL

  if (DATASET == "diabetes") {
    diabetes_data   <- read.csv("data/diabetes.csv")
    x_all           <- as.matrix(diabetes_data[, 1:8])
    y_all           <- 1 - diabetes_data[, 9]

  } else if (DATASET == "transfusion") {
    transfusion_data <- read.csv("data/transfusion.data", header = TRUE)
    x_all            <- as.matrix(transfusion_data[, 1:4])
    y_all            <- 1 - transfusion_data[, 5]

  } else if (DATASET == "heart") {
    heart_data      <- read.csv("data/CTG.csv", header = TRUE)
    x_all           <- as.matrix(heart_data[, -22])
    y_all           <- ifelse(heart_data$NSP == 1, 1L, 0L)
    keep            <- complete.cases(x_all) & !is.na(y_all)
    x_all           <- scale(x_all[keep, , drop = FALSE])
    y_all           <- y_all[keep]

  } else if (DATASET == "mammography") {
    MAMMO_N  <- 1000
    MAMMO_PI <- 0.2

    mammo_raw    <- read.csv("data/mammography.csv")
    x_mammo_full <- as.matrix(mammo_raw[, -ncol(mammo_raw)])
    y_raw        <- mammo_raw[[ncol(mammo_raw)]]
    minority_label <- names(which.min(table(as.character(y_raw))))
    y_mammo_full <- ifelse(as.character(y_raw) == minority_label, 0, 1)
    keep         <- complete.cases(x_mammo_full) & !is.na(y_mammo_full)
    x_mammo_full <- scale(x_mammo_full[keep, , drop = FALSE])
    y_mammo_full <- y_mammo_full[keep]

    subsample_data <- function() {
      idx0      <- which(y_mammo_full == 0)
      idx1      <- which(y_mammo_full == 1)
      n0        <- round(MAMMO_N * MAMMO_PI)
      n1        <- MAMMO_N - n0
      idx_train <- c(sample(idx0, min(n0, length(idx0))), sample(idx1, min(n1, length(idx1))))
      idx_pop   <- setdiff(seq_along(y_mammo_full), idx_train)
      list(
        train_data = list(x = x_mammo_full[idx_train, , drop = FALSE], y = y_mammo_full[idx_train]),
        pop_data   = list(x = x_mammo_full[idx_pop,   , drop = FALSE], y = y_mammo_full[idx_pop])
      )
    }

  } else if (DATASET == "phoneme") {
    PHONEME_N  <- 1000
    PHONEME_PI <- 0.2

    phoneme_raw    <- read.csv("data/phoneme.csv")
    label_col      <- "Class"
    y_raw          <- phoneme_raw[[label_col]]
    x_raw          <- phoneme_raw[, setdiff(names(phoneme_raw), label_col), drop = FALSE]
    x_phoneme_full <- model.matrix(~ . - 1, data = x_raw)
    minority_label <- names(which.min(table(as.character(y_raw))))
    y_phoneme_full <- ifelse(as.character(y_raw) == minority_label, 0, 1)
    keep           <- complete.cases(x_phoneme_full) & !is.na(y_phoneme_full)
    x_phoneme_full <- scale(x_phoneme_full[keep, , drop = FALSE])
    y_phoneme_full <- y_phoneme_full[keep]

    subsample_data <- function() {
      idx0      <- which(y_phoneme_full == 0)
      idx1      <- which(y_phoneme_full == 1)
      n0        <- round(PHONEME_N * PHONEME_PI)
      n1        <- PHONEME_N - n0
      idx_train <- c(sample(idx0, min(n0, length(idx0))), sample(idx1, min(n1, length(idx1))))
      idx_pop   <- setdiff(seq_along(y_phoneme_full), idx_train)
      list(
        train_data = list(x = x_phoneme_full[idx_train, , drop = FALSE], y = y_phoneme_full[idx_train]),
        pop_data   = list(x = x_phoneme_full[idx_pop,   , drop = FALSE], y = y_phoneme_full[idx_pop])
      )
    }

  } else if (DATASET == "wine_quality") {
    WINE_N              <- 1000
    WINE_PI             <- 0.2
    WINE_HIGH_THRESHOLD <- 7

    wine_raw <- read.csv("data/winequality-white.csv", sep = ";", stringsAsFactors = FALSE)

    y_wine_full <- ifelse(wine_raw$quality >= WINE_HIGH_THRESHOLD, 0, 1)
    x_raw       <- wine_raw[, setdiff(names(wine_raw), "quality"), drop = FALSE]
    x_wine_full <- model.matrix(~ . - 1, data = x_raw)
    x_wine_full <- x_wine_full[, apply(x_wine_full, 2, sd) > 0, drop = FALSE]

    subsample_data <- function() {
      idx0 <- which(y_wine_full == 0)
      idx1 <- which(y_wine_full == 1)
      n0   <- round(WINE_N * WINE_PI)
      n1   <- WINE_N - n0
      if (n0 > length(idx0)) { warning("Capping n0."); n0 <- length(idx0); n1 <- min(WINE_N - n0, length(idx1)) }
      if (n1 > length(idx1)) { warning("Capping n1."); n1 <- length(idx1); n0 <- min(WINE_N - n1, length(idx0)) }
      idx_train <- c(sample(idx0, n0), sample(idx1, n1))
      idx_pop   <- setdiff(seq_along(y_wine_full), idx_train)
      list(
        train_data = list(x = x_wine_full[idx_train, , drop = FALSE], y = y_wine_full[idx_train]),
        pop_data   = list(x = x_wine_full[idx_pop,   , drop = FALSE], y = y_wine_full[idx_pop])
      )
    }

  } else if (DATASET == "magic") {
    MAGIC_N  <- 1000
    MAGIC_PI <- 0.2

    magic_names <- c("fLength", "fWidth", "fSize", "fConc", "fConc1",
                     "fAsym", "fM3Long", "fM3Trans", "fAlpha", "fDist", "class")
    magic_raw    <- read.csv("data/magic04.data", header = FALSE, col.names = magic_names)
    x_magic_full <- as.matrix(magic_raw[, magic_names[-length(magic_names)]])
    y_magic_full <- ifelse(as.character(magic_raw$class) == "g", 0L, 1L)
    keep         <- complete.cases(x_magic_full) & !is.na(y_magic_full)
    x_magic_full <- scale(x_magic_full[keep, , drop = FALSE])
    y_magic_full <- y_magic_full[keep]

    subsample_data <- function() {
      idx0      <- which(y_magic_full == 0)
      idx1      <- which(y_magic_full == 1)
      n0        <- round(MAGIC_N * MAGIC_PI)
      n1        <- MAGIC_N - n0
      idx_train <- c(sample(idx0, min(n0, length(idx0))), sample(idx1, min(n1, length(idx1))))
      idx_pop   <- setdiff(seq_along(y_magic_full), idx_train)
      list(
        train_data = list(x = x_magic_full[idx_train, , drop = FALSE], y = y_magic_full[idx_train]),
        pop_data   = list(x = x_magic_full[idx_pop,   , drop = FALSE], y = y_magic_full[idx_pop])
      )
    }

  } else {
    stop("Unknown DATASET: ", DATASET)
  }

  list(x_all = x_all, y_all = y_all,
       subsample_data = subsample_data, population_data = population_data)
}
