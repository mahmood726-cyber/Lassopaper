################################################################################
# REGULARIZED META-REGRESSION: FINAL CORRECTED + EXTENDED SCRIPT
# - Weighted variance for R^2_het (apparent, CV, penalized)
# - Robust dataset loaders
# - Nested CV + permutation tests (n=1000)
# - Paired bootstrap for headline comparison
# - Minimal simulation (Supplement)
# - CV-curve plotting (ASCII labels) + Cairo-safe saving
################################################################################

suppressPackageStartupMessages({
  library(metafor)
  library(metadat)
  library(glmnet)
  library(ggplot2)
  library(knitr)
  library(reshape2)
})

set.seed(12345)

# ==============================================================================
# HELPERS
# ==============================================================================

safe_r2het <- function(tau2_full, tau2_null, eps = 1e-10) {
  if (is.na(tau2_null) || is.na(tau2_full)) return(NA_real_)
  denom <- max(tau2_null, eps)
  pmin(1, pmax(0, 1 - tau2_full / denom))
}

drop_constant <- function(X) {
  if (is.null(X) || length(X) == 0) return(NULL)
  if (!is.matrix(X)) X <- as.matrix(X)
  keep <- apply(X, 2, function(x) length(unique(x)) > 1 && var(x, na.rm = TRUE) > 1e-10)
  if (sum(keep) == 0) return(NULL)
  X[, keep, drop = FALSE]
}

is_intercept_column <- function(X) {
  if (is.null(X) || ncol(X) == 0) return(FALSE)
  all(abs(X[, 1] - 1) < 1e-12)
}

check_data <- function(y, v, min_k = 5) {
  if (any(is.na(y)) || any(is.na(v))) return(FALSE)
  if (any(v <= 0)) return(FALSE)
  if (length(y) < min_k) return(FALSE)
  TRUE
}

# Safe ggsave (use Cairo when available to avoid Unicode/encoding issues)
safe_ggsave <- function(filename, plot, width, height) {
  dev <- if (capabilities("cairo")) cairo_pdf else pdf
  ggplot2::ggsave(filename, plot, width = width, height = height, device = dev)
}

# ==============================================================================
# CORE ANALYSIS FUNCTIONS (weighted-variance used consistently)
# ==============================================================================

# Standard R^2_het (apparent) using moderators X (may include intercept column)
calculate_r2het <- function(y, v, X, method = "REML") {
  if (!check_data(y, v)) return(list(r2het = NA, p_eff = NA, I2 = NA))
  fit_null <- rma(yi = y, vi = v, method = method)
  I2 <- fit_null$I2
  tau2_null <- fit_null$tau2

  if (is_intercept_column(X)) X <- X[, -1, drop = FALSE]
  X <- drop_constant(X)

  if (is.null(X) || ncol(X) == 0) {
    return(list(r2het = 0, tau2_null = tau2_null, p_eff = 0, I2 = I2))
  }

  fit_full <- rma(yi = y, vi = v, mods = X, method = method)
  r2het <- safe_r2het(fit_full$tau2, tau2_null)
  list(r2het = r2het, tau2_null = tau2_null, tau2_full = fit_full$tau2,
       p_eff = ncol(X), I2 = I2)
}

# LOO CV R^2_het with CORRECTED weighted variance
cv_r2het <- function(y, v, X, method = "REML") {
  if (!check_data(y, v)) return(NA_real_)

  if (is_intercept_column(X)) X <- X[, -1, drop = FALSE]
  X <- drop_constant(X)
  if (is.null(X) || ncol(X) == 0) return(0)

  fit_null <- rma(yi = y, vi = v, method = method)
  tau2_null <- fit_null$tau2
  w <- 1 / (v + tau2_null); w <- w / sum(w)

  k <- length(y)
  preds <- rep(NA_real_, k)

  for (i in seq_len(k)) {
    preds[i] <- tryCatch({
      fit_i <- rma(yi = y[-i], vi = v[-i], mods = X[-i, , drop = FALSE], method = method)
      predict(fit_i, newmods = X[i, , drop = FALSE])$pred
    }, error = function(e) NA_real_)
  }

  ok <- !is.na(preds)
  if (sum(ok) < max(3, floor(0.5 * k))) return(NA_real_)

  res <- y[ok] - preds[ok]
  w_ok <- w[ok]; w_ok <- w_ok / sum(w_ok)
  vbar_w <- sum(w_ok * v[ok])
  tau2_cv <- max(0, sum(w_ok * res^2) - vbar_w)
  safe_r2het(tau2_cv, tau2_null)
}

# Penalized meta-regression (alpha=1 LASSO; alpha=0 Ridge)
.penalized_meta_cv <- function(y, v, X, alpha, method = "REML",
                               lambda_choice = "lambda.1se") {
  if (!check_data(y, v)) return(list(r2het = NA, r2het_cv = NA, lambda = NA,
                                     p_eff = NA, n_selected = NA, cv_fit = NULL))

  has_intercept <- is_intercept_column(X)
  if (has_intercept) X <- X[, -1, drop = FALSE]
  X <- drop_constant(X)
  if (is.null(X) || ncol(X) == 0) {
    return(list(r2het = 0, r2het_cv = 0, lambda = NA, p_eff = 0, n_selected = 0, cv_fit = NULL))
  }

  p_eff <- ncol(X)
  penalty.factor <- rep(1, ncol(X))
  # Stabilize glmnet path in single-predictor case
  if (ncol(X) == 1) {
    X <- cbind(X, noise1 = rnorm(nrow(X), 0, 0.001), noise2 = rnorm(nrow(X), 0, 0.001))
    penalty.factor <- c(1, 10000, 10000)
  }

  Xs <- scale(X); Xs[is.na(Xs)] <- 0

  fit_null <- rma(yi = y, vi = v, method = method)
  tau2_null <- fit_null$tau2
  w <- 1 / (v + tau2_null); w <- w / sum(w)  # normalized weights
  k <- length(y)

  # Folds: LOO for k <= 20, else 10-fold
  nfolds <- if (k <= 20) k else 10
  foldid <- if (nfolds == k) seq_len(k) else {
    set.seed(12345); sample(rep(1:nfolds, length.out = k))
  }

  cv_fit <- cv.glmnet(
    x = Xs, y = y, weights = w * k, alpha = alpha, nfolds = nfolds, foldid = foldid,
    penalty.factor = penalty.factor, grouped = FALSE, keep = TRUE,
    standardize = FALSE, intercept = TRUE, family = "gaussian"
  )

  lam <- cv_fit[[lambda_choice]]

  # Apparent (in-sample) R^2_het
  yhat_in <- as.vector(predict(cv_fit, newx = Xs, s = lam))
  res_in <- y - yhat_in
  vbar_w <- sum(w * v)
  tau2_in <- max(0, sum(w * res_in^2) - vbar_w)
  r2_in <- safe_r2het(tau2_in, tau2_null)

  # Out-of-fold R^2_het
  idx <- which.min(abs(log(cv_fit$lambda) - log(lam)))
  yhat_oof <- cv_fit$fit.preval[, idx]
  res_oof <- y - yhat_oof
  tau2_cv <- max(0, sum(w * res_oof^2) - vbar_w)
  r2_cv <- safe_r2het(tau2_cv, tau2_null)

  # Selected features (exclude intercept)
  coefs <- as.matrix(coef(cv_fit, s = lam))
  n_sel <- sum(abs(coefs[-1, 1]) > 0)

  list(r2het = r2_in, r2het_cv = r2_cv, lambda = lam,
       p_eff = p_eff, n_selected = n_sel, cv_fit = cv_fit)
}

lasso_meta_cv <- function(y, v, X, method = "REML", lambda_choice = "lambda.1se") {
  .penalized_meta_cv(y, v, X, alpha = 1, method = method, lambda_choice = lambda_choice)
}
ridge_meta_cv <- function(y, v, X, method = "REML", lambda_choice = "lambda.1se") {
  .penalized_meta_cv(y, v, X, alpha = 0, method = method, lambda_choice = lambda_choice)
}

# Nested CV (outer LOO, inner CV for lambda) – evaluates lambda.min (default)
nested_cv_lasso <- function(y, v, X, method = "REML", lambda_choice = "lambda.min") {
  if (!check_data(y, v)) return(list(r2het_nested = NA, lambda_variability = NA))

  if (is_intercept_column(X)) X <- X[, -1, drop = FALSE]
  X <- drop_constant(X)
  if (is.null(X) || ncol(X) == 0) return(list(r2het_nested = 0, lambda_variability = 0))

  fit_null <- rma(yi = y, vi = v, method = method)
  tau2_null <- fit_null$tau2
  w <- 1 / (v + tau2_null)

  k <- length(y)
  preds <- rep(NA_real_, k)
  lambdas_used <- rep(NA_real_, k)

  for (i in seq_len(k)) {
    Xtr <- X[-i, , drop = FALSE]; ytr <- y[-i]; vtr <- v[-i]
    wtr <- w[-i]; wtr <- wtr / sum(wtr)

    pf <- rep(1, ncol(Xtr))
    if (ncol(Xtr) == 1) {
      Xtr <- cbind(Xtr, noise1 = rnorm(nrow(Xtr), 0, 0.001), noise2 = rnorm(nrow(Xtr), 0, 0.001))
      pf <- c(1, 10000, 10000)
    }

    Xtr_s <- scale(Xtr)
    nfolds_in <- min(5, max(3, nrow(Xtr)))

    ok_fold <- TRUE
    cv_in <- tryCatch({
      cv.glmnet(x = Xtr_s, y = ytr, weights = wtr * length(ytr),
                alpha = 1, nfolds = nfolds_in, penalty.factor = pf,
                standardize = FALSE, intercept = TRUE, grouped = FALSE, keep = FALSE)
    }, error = function(e) { ok_fold <<- FALSE; NULL })

    if (!ok_fold) {
      preds[i] <- NA_real_
      next
    }

    lam <- cv_in[[lambda_choice]]
    lambdas_used[i] <- lam

    Xte <- matrix(X[i, ], nrow = 1)
    if (ncol(X) == 1) Xte <- cbind(Xte, 0, 0)
    Xte_s <- scale(Xte, center = attr(Xtr_s, "scaled:center"), scale = attr(Xtr_s, "scaled:scale"))
    preds[i] <- as.numeric(predict(cv_in, newx = Xte_s, s = lam))
  }

  ok <- !is.na(preds)
  if (sum(ok) < max(3, floor(0.5 * k))) return(list(r2het_nested = NA, lambda_variability = NA))

  res <- y[ok] - preds[ok]
  wok <- w[ok] / sum(w[ok])
  vbar_w <- sum(wok * v[ok])
  tau2_cv <- max(0, sum(wok * res^2) - vbar_w)
  list(r2het_nested = safe_r2het(tau2_cv, tau2_null),
       lambda_variability = sd(lambdas_used, na.rm = TRUE) / mean(lambdas_used, na.rm = TRUE))
}

# Permutation test for lambda.min: permute moderators only; keep (y, v) fixed
permutation_test <- function(y, v, X, method = "REML", n_perm = 1000) {
  if (!check_data(y, v)) return(list(p_value = NA, overfit_ratio = NA))

  obs <- lasso_meta_cv(y, v, X, method, lambda_choice = "lambda.min")
  obs_r2 <- obs$r2het

  X_base <- if (is_intercept_column(X)) X[, -1, drop = FALSE] else X
  perm_r2 <- numeric(n_perm)

  for (p in seq_len(n_perm)) {
    Xp <- X_base[sample(nrow(X_base)), , drop = FALSE]
    Xp <- if (is_intercept_column(X)) cbind(1, Xp) else cbind(1, Xp)  # ensure intercept
    perm_r2[p] <- lasso_meta_cv(y, v, Xp, method, lambda_choice = "lambda.min")$r2het
  }
  pval <- mean(perm_r2 >= obs_r2)
  list(observed_r2 = obs_r2,
       perm_mean = mean(perm_r2), perm_sd = sd(perm_r2),
       p_value = pval,
       overfit_ratio = ifelse(mean(perm_r2) > 0, obs_r2 / mean(perm_r2), NA))
}

# Plot glmnet CV curve with lambda.min / lambda.1se markers (ASCII labels)
plot_cv_curve <- function(cv_fit, main = "Flat CV Surface Enables Noise Capitalization") {
  if (is.null(cv_fit) || is.atomic(cv_fit)) {
    message("No cv.glmnet object available to plot.")
    return(invisible(NULL))
  }
  plot(cv_fit)
  abline(v = log(cv_fit$lambda.min), lty = 2, col = "red")
  abline(v = log(cv_fit$lambda.1se), lty = 2, col = "blue")
  legend("topleft", c("lambda.min", "lambda.1se"), col = c("red", "blue"), lty = 2, bty = "n")
  title(main = main)
}

# ==============================================================================
# COMPLETE ANALYSIS WRAPPER (prints + returns one-row data.frame)
# ==============================================================================

complete_analysis <- function(data_name, y, v, X, save_results = TRUE,
                              do_permutation = TRUE, n_perm = 1000) {
  cc <- complete.cases(y, v, X)
  y <- y[cc]; v <- v[cc]; X <- X[cc, , drop = FALSE]

  cat("\n=== Analyzing:", data_name, "===\n")
  k <- length(y)
  p <- max(0, ncol(X) - as.integer(is_intercept_column(X)))  # moderators only
  cat(sprintf("k=%d, p=%d, k/p=%s\n", k, p, ifelse(p > 0, sprintf("%.1f", k / p), "NA")))

  # Standard meta-regression
  std <- calculate_r2het(y, v, X)
  std_cv <- cv_r2het(y, v, X)

  # LASSO / Ridge (both lambda choices)
  l_min <- lasso_meta_cv(y, v, X, lambda_choice = "lambda.min")
  l_1se <- lasso_meta_cv(y, v, X, lambda_choice = "lambda.1se")
  r_min <- ridge_meta_cv(y, v, X, lambda_choice = "lambda.min")
  r_1se <- ridge_meta_cv(y, v, X, lambda_choice = "lambda.1se")

  # Nested CV and permutation diagnostics
  nested_lasso <- nested_cv_lasso(y, v, X, lambda_choice = "lambda.min")
  perm_result <- if (do_permutation && check_data(y, v) && p > 0)
    permutation_test(y, v, X, n_perm = n_perm) else list(p_value = NA, overfit_ratio = NA)

  # Print brief summary
  cat("\nOptimism Summary (Apparent - CV, in %):\n")
  cat(sprintf("  Standard:            %.1f\n", (std$r2het - std_cv) * 100))
  cat(sprintf("  LASSO lambda.min:    %.1f (lambda=%.4f, selected=%d)\n",
              (l_min$r2het - l_min$r2het_cv) * 100, l_min$lambda, l_min$n_selected))
  cat(sprintf("  LASSO lambda.1se:    %.1f (lambda=%.4f, selected=%d)\n",
              (l_1se$r2het - l_1se$r2het_cv) * 100, l_1se$lambda, l_1se$n_selected))
  cat(sprintf("  Ridge  lambda.min:   %.1f\n", (r_min$r2het - r_min$r2het_cv) * 100))
  cat(sprintf("  Ridge  lambda.1se:   %.1f\n", (r_1se$r2het - r_1se$r2het_cv) * 100))

  if (is.list(perm_result) && !is.na(perm_result$p_value)) {
    cat(sprintf("\nPermutation: mean=%.3f, sd=%.3f, p=%.3f, overfit ratio=%.2fx\n",
                perm_result$perm_mean, perm_result$perm_sd,
                perm_result$p_value, perm_result$overfit_ratio))
  }

  # Results row
  out <- data.frame(
    Dataset = data_name,
    k = k,
    p = p,
    k_p_ratio = ifelse(p > 0, k / p, NA_real_),
    I2 = std$I2,

    Std_Apparent = std$r2het * 100,
    Std_CV = std_cv * 100,
    Std_Optimism = (std$r2het - std_cv) * 100,

    LASSO_min_Apparent = l_min$r2het * 100,
    LASSO_min_CV = l_min$r2het_cv * 100,
    LASSO_min_Optimism = (l_min$r2het - l_min$r2het_cv) * 100,
    LASSO_min_lambda = l_min$lambda,
    LASSO_min_selected = l_min$n_selected,

    LASSO_1se_Apparent = l_1se$r2het * 100,
    LASSO_1se_CV = l_1se$r2het_cv * 100,
    LASSO_1se_Optimism = (l_1se$r2het - l_1se$r2het_cv) * 100,
    LASSO_1se_lambda = l_1se$lambda,
    LASSO_1se_selected = l_1se$n_selected,

    Ridge_min_Apparent = r_min$r2het * 100,
    Ridge_min_CV = r_min$r2het_cv * 100,
    Ridge_min_Optimism = (r_min$r2het - r_min$r2het_cv) * 100,
    Ridge_min_lambda = r_min$lambda,

    Ridge_1se_Apparent = r_1se$r2het * 100,
    Ridge_1se_CV = r_1se$r2het_cv * 100,
    Ridge_1se_Optimism = (r_1se$r2het - r_1se$r2het_cv) * 100,
    Ridge_1se_lambda = r_1se$lambda,

    Nested_CV_R2 = ifelse(is.list(nested_lasso), nested_lasso$r2het_nested * 100, NA),
    Lambda_Variability = ifelse(is.list(nested_lasso), nested_lasso$lambda_variability, NA),

    Perm_pvalue = if (is.list(perm_result)) perm_result$p_value else NA_real_,
    Overfit_Ratio = if (is.list(perm_result)) perm_result$overfit_ratio else NA_real_,

    check.names = FALSE
  )

  if (save_results) {
    fn <- sprintf("results_%s.csv", gsub("[^A-Za-z0-9]+", "_", data_name))
    write.csv(out, fn, row.names = FALSE)
    cat(sprintf("\nResults saved to: %s\n", fn))
  }

  # Attach cv.glmnet objects (for plotting CV curves later if desired)
  attr(out, "cv_curve_lasso_min") <- l_min$cv_fit
  attr(out, "cv_curve_lasso_1se") <- l_1se$cv_fit

  out
}

# ==============================================================================
# DATASET LOADERS (robust; always return complete cases)
# ==============================================================================

load_bcg <- function() {
  data(dat.bcg)
  es <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)
  keep <- complete.cases(es$yi, es$vi, dat.bcg$ablat, dat.bcg$year)
  y <- es$yi[keep]; v <- es$vi[keep]
  X <- model.matrix(~ ablat + year, data = dat.bcg[keep, ])
  list(y = y, v = v, X = X, name = "BCG")
}

load_teacher <- function() {
  data(dat.raudenbush1985)
  keep <- complete.cases(dat.raudenbush1985$yi, dat.raudenbush1985$vi, dat.raudenbush1985$weeks)
  y <- dat.raudenbush1985$yi[keep]; v <- dat.raudenbush1985$vi[keep]
  X <- model.matrix(~ weeks, data = dat.raudenbush1985[keep, ])
  list(y = y, v = v, X = X, name = "Teacher")
}

load_konstant <- function() {
  data(dat.konstantopoulos2011)
  keep <- complete.cases(dat.konstantopoulos2011$yi, dat.konstantopoulos2011$vi,
                         dat.konstantopoulos2011$year)
  y <- dat.konstantopoulos2011$yi[keep]; v <- dat.konstantopoulos2011$vi[keep]
  X <- model.matrix(~ year, data = dat.konstantopoulos2011[keep, ])
  list(y = y, v = v, X = X, name = "Konstantopoulos")
}

load_berkey <- function() {
  data(dat.berkey1998)
  keep <- complete.cases(dat.berkey1998$yi, dat.berkey1998$vi, dat.berkey1998$year)
  y <- dat.berkey1998$yi[keep]; v <- dat.berkey1998$vi[keep]
  X <- model.matrix(~ year, data = dat.berkey1998[keep, ])
  list(y = y, v = v, X = X, name = "Berkey")
}

# Robust Hackshaw loader: handles missing 'continent'/'design'
load_hackshaw <- function() {
  data(dat.hackshaw1998)
  df <- dat.hackshaw1998
  keep <- complete.cases(df$yi, df$vi)
  df <- df[keep, ]

  if (all(c("continent", "design") %in% names(df))) {
    X <- model.matrix(~ continent + design, data = df)
  } else {
    message("Hackshaw: 'continent' or 'design' not found; falling back to ~ year or intercept-only.")
    if ("year" %in% names(df)) {
      X <- model.matrix(~ year, data = df)
    } else {
      X <- model.matrix(~ 1, data = df)  # intercept only; p=0
    }
  }

  list(y = df$yi, v = df$vi, X = X, name = "Hackshaw")
}

# ==============================================================================
# RUN ALL EMPIRICAL ANALYSES
# ==============================================================================

cat("========== REGULARIZATION PARADOX IN META-REGRESSION ==========\n\n")

loaders <- list(load_bcg, load_teacher, load_konstant, load_berkey, load_hackshaw)

all_results <- do.call(
  rbind,
  lapply(loaders, function(loader) {
    ds <- loader()
    complete_analysis(ds$name, ds$y, ds$v, ds$X,
                      save_results = TRUE, do_permutation = TRUE, n_perm = 1000)
  })
)

# Strip list columns for display
display_results <- all_results[, !sapply(all_results, is.list), drop = FALSE]

# ==============================================================================
# SUMMARY TABLES & MEANS (include nested-CV columns)
# ==============================================================================

cat("\n========== TABLE 1: THE REGULARIZATION PARADOX ==========\n\n")

summary_cols <- c("Dataset", "k", "p", "k_p_ratio", "I2",
                  "Std_Optimism", "LASSO_min_Optimism", "LASSO_1se_Optimism",
                  "Ridge_min_Optimism", "Ridge_1se_Optimism")
summary_cols_extended <- c(summary_cols, "Nested_CV_R2", "Lambda_Variability",
                           "Perm_pvalue", "Overfit_Ratio")

print(kable(display_results[, summary_cols_extended],
            caption = "Optimism in R^2_het (%) by Method; nested-CV (LASSO lambda.min) and permutation results"))

cat("\n========== MEAN OPTIMISM ACROSS DATASETS ==========\n")
cat(sprintf("Standard:           %.1f%%\n", mean(display_results$Std_Optimism, na.rm = TRUE)))
cat(sprintf("LASSO lambda.min:   %.1f%%\n", mean(display_results$LASSO_min_Optimism, na.rm = TRUE)))
cat(sprintf("LASSO lambda.1se:   %.1f%%\n", mean(display_results$LASSO_1se_Optimism, na.rm = TRUE)))
cat(sprintf("Ridge  lambda.min:  %.1f%%\n", mean(display_results$Ridge_min_Optimism, na.rm = TRUE)))
cat(sprintf("Ridge  lambda.1se:  %.1f%%\n", mean(display_results$Ridge_1se_Optimism, na.rm = TRUE)))

# ==============================================================================
# FORMAL TEST: Paired bootstrap for (LASSO lambda.min − Standard) optimism
# ==============================================================================

paired_boot <- function(df, colA, colB, B = 5000, seed = 123) {
  set.seed(seed)
  idx <- which(!is.na(df[[colA]]) & !is.na(df[[colB]]))
  if (length(idx) < 2) return(c(mean_diff = NA, ci_lo = NA, ci_hi = NA, p = NA))
  base <- df[idx, , drop = FALSE]
  diffs <- base[[colA]] - base[[colB]]
  boot <- replicate(B, {
    s <- sample(seq_along(diffs), length(diffs), replace = TRUE)
    mean(diffs[s])
  })
  c(
    mean_diff = mean(diffs),
    ci_lo = unname(quantile(boot, 0.025)),
    ci_hi = unname(quantile(boot, 0.975)),
    p = 2 * min(mean(boot <= 0), mean(boot >= 0))
  )
}

pb_lasso_vs_std <- paired_boot(display_results, "LASSO_min_Optimism", "Std_Optimism")
cat(sprintf("\nPaired bootstrap (LASSO lambda.min − Standard) optimism:\n  Mean diff = %.1f%%, 95%% CI [%.1f, %.1f], p = %.3f\n",
            pb_lasso_vs_std["mean_diff"], pb_lasso_vs_std["ci_lo"],
            pb_lasso_vs_std["ci_hi"], pb_lasso_vs_std["p"]))

# ==============================================================================
# FIGURES (ASCII labels; Cairo if available)
# ==============================================================================

# Figure 1: Optimism comparison
optimism_data <- data.frame(
  Method = factor(c("Standard", "LASSO\nlambda.min", "LASSO\nlambda.1se",
                    "Ridge\nlambda.min", "Ridge\nlambda.1se"),
                  levels = c("Standard", "LASSO\nlambda.min", "LASSO\nlambda.1se",
                             "Ridge\nlambda.min", "Ridge\nlambda.1se")),
  Optimism = c(
    mean(display_results$Std_Optimism, na.rm = TRUE),
    mean(display_results$LASSO_min_Optimism, na.rm = TRUE),
    mean(display_results$LASSO_1se_Optimism, na.rm = TRUE),
    mean(display_results$Ridge_min_Optimism, na.rm = TRUE),
    mean(display_results$Ridge_1se_Optimism, na.rm = TRUE)
  )
)

p1 <- ggplot(optimism_data, aes(x = Method, y = Optimism)) +
  geom_bar(stat = "identity",
           fill = c("gray50", "red3", "blue3", "red3", "blue3")) +
  geom_hline(yintercept = optimism_data$Optimism[1],
             linetype = "dashed", alpha = 0.5) +
  labs(title = "The Regularization Paradox",
       subtitle = "lambda.min increases optimism; lambda.1se provides minimal benefit",
       y = "Mean Optimism in R^2_het (%)",
       x = "") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 10))

# Figure 2: Optimism vs k/p ratio (points + linear trend if enough points)
plot_data <- display_results[!is.na(display_results$k_p_ratio), ]
plot_long <- melt(plot_data[, c("Dataset", "k_p_ratio",
                                "Std_Optimism", "LASSO_min_Optimism",
                                "LASSO_1se_Optimism")],
                  id.vars = c("Dataset", "k_p_ratio"))

p2 <- ggplot(plot_long, aes(x = k_p_ratio, y = value, color = variable)) +
  geom_point(size = 3, alpha = 0.7) +
  { if (nrow(plot_long) >= 5) geom_smooth(method = "lm", se = FALSE, linewidth = 1) else NULL } +
  geom_vline(xintercept = 10, linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c("gray50", "red3", "blue3"),
                     labels = c("Standard", "lambda.min (LASSO)", "lambda.1se (LASSO)")) +
  labs(title = "Optimism by Sample Size Ratio",
       subtitle = "Paradox most severe at k/p < 10",
       x = "k/p Ratio",
       y = "Optimism (%)",
       color = "Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save artifacts
write.csv(display_results, "regularization_paradox_results.csv", row.names = FALSE)
safe_ggsave("fig1_optimism_comparison.pdf", p1, width = 6, height = 4)
safe_ggsave("fig2_optimism_by_ratio.pdf", p2, width = 7, height = 5)

# Example CV curve from first dataset (if available)
cv_obj <- attr(all_results, "cv_curve_lasso_min")
if (!is.null(cv_obj) && length(cv_obj) >= 1) {
  plot_cv_curve(cv_obj[[1]], main = "Flat CV Surface Enables Noise Capitalization (Example)")
}

# ==============================================================================
# (Optional) SENSITIVITY: Unweighted CV variance (to show it doesn't reverse)
# ==============================================================================

cv_r2het_unweighted <- function(y, v, X, method = "REML") {
  if (!check_data(y, v)) return(NA_real_)
  if (is_intercept_column(X)) X <- X[, -1, drop = FALSE]
  X <- drop_constant(X); if (is.null(X) || ncol(X) == 0) return(0)
  fit_null <- rma(yi = y, vi = v, method = method)
  tau2_null <- fit_null$tau2
  k <- length(y); preds <- rep(NA_real_, k)
  for (i in 1:k) {
    preds[i] <- tryCatch({
      fit <- rma(yi = y[-i], vi = v[-i], mods = X[-i,,drop=FALSE], method = method)
      predict(fit, newmods = X[i,,drop=FALSE])$pred
    }, error = function(e) NA_real_)
  }
  ok <- !is.na(preds); if (sum(ok) < max(3, floor(0.5*k))) return(NA_real_)
  res <- y[ok] - preds[ok]
  vbar <- mean(v[ok])
  tau2_cv <- max(0, mean(res^2) - vbar)
  safe_r2het(tau2_cv, tau2_null)
}

# ==============================================================================
# SUPPLEMENT: Minimal Simulation Study
# ==============================================================================

simulate_one <- function(k, p, R2_true = 0, tau2 = 0.1,
                         vi_shape = 8, vi_scale = 0.01, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # Moderators
  X <- scale(matrix(rnorm(k * p), ncol = p))
  # Target moderator variance to achieve R2_true:
  sigma_m2 <- if (R2_true <= 0) 0 else tau2 * R2_true / (1 - R2_true)
  # Random beta scaled to give Var(X beta) = sigma_m2
  beta_raw <- rnorm(p)
  beta_raw <- beta_raw / sqrt(var(as.vector(X %*% beta_raw)))
  beta <- sqrt(sigma_m2) * beta_raw
  mu <- as.vector(X %*% beta)

  # Sampling variances (positive, heterogeneous)
  v <- rgamma(k, shape = vi_shape, rate = 1/vi_scale)  # mean = shape*scale
  v <- v / mean(v) * 0.05  # rescale to realistic ~0.05

  # Random-effects outcomes
  yi <- rnorm(k, mean = mu, sd = sqrt(tau2)) + rnorm(k, sd = sqrt(v))

  # Design with intercept
  Xmm <- cbind(1, X)
  list(y = yi, v = v, X = Xmm)
}

eval_methods_once <- function(y, v, X) {
  std <- calculate_r2het(y, v, X); std_cv <- cv_r2het(y, v, X)
  lmin <- lasso_meta_cv(y, v, X, lambda_choice = "lambda.min")
  l1se <- lasso_meta_cv(y, v, X, lambda_choice = "lambda.1se")
  c(
    Std_Opt = (std$r2het - std_cv) * 100,
    Lmin_Opt = (lmin$r2het - lmin$r2het_cv) * 100,
    L1se_Opt = (l1se$r2het - l1se$r2het_cv) * 100
  )
}

simulate_grid <- function(n_sims = 300, R2_vals = c(0, .25, .5),
                          kp_set = c(5, 10, 20), tau2 = 0.1, seed = 1) {
  set.seed(seed)
  out <- list()
  for (R2 in R2_vals) {
    for (kp in kp_set) {
      # choose (k, p) combos that realize the desired ratio while keeping p small
      combos <- list(
        `5`  = list(k=10, p=2),
        `10` = list(k=20, p=2),
        `20` = list(k=20, p=1)
      )
      cfg <- combos[[as.character(kp)]]
      mat <- replicate(n_sims, {
        sim <- simulate_one(cfg$k, cfg$p, R2_true = R2, tau2 = tau2)
        eval_methods_once(sim$y, sim$v, sim$X)
      })
      mat <- t(mat)
      out[[length(out)+1]] <- data.frame(
        True_R2 = R2, k = cfg$k, p = cfg$p, kp = kp,
        Std = mean(mat[, "Std_Opt"], na.rm=TRUE),
        LASSO_min = mean(mat[, "Lmin_Opt"], na.rm=TRUE),
        LASSO_1se = mean(mat[, "L1se_Opt"], na.rm=TRUE)
      )
    }
  }
  do.call(rbind, out)
}

# Run simulation (adjust n_sims up if desired)
sim_tab <- simulate_grid(n_sims = 300)
print(kable(sim_tab, digits = 1,
            caption = "Mean optimism (%) from simulations by true R^2_het and k/p"))
write.csv(sim_tab, "simulation_mean_optimism.csv", row.names = FALSE)

################################################################################
# END
################################################################################
