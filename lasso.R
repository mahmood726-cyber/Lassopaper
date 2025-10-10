################################################################################
# REGULARIZED META-REGRESSION: FINAL FIXED & PAPER-READY SCRIPT (v2)
# - Correct weighted variance for R^2_het (apparent, CV, penalized)
# - Robust dataset loaders (fixes `continent` error for Hackshaw)
# - Nested CV and permutation diagnostics
# - CV-curve plotting with lambda.min / lambda.1se markers (ASCII-safe)
# - Figures save cleanly on Windows; smoothing avoids NaNs with tiny N
# - FIX: store cv.glmnet objects in separate lists (no attribute loss on rbind)
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
  if (!any(keep)) return(NULL)
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

# Safe ggsave (use cairo when available to avoid encoding issues)
safe_ggsave <- function(filename, plot, width, height) {
  dev <- if (capabilities("cairo")) cairo_pdf else pdf
  ggplot2::ggsave(filename, plot, width = width, height = height, device = dev)
}

# ==============================================================================
# CORE ANALYSIS FUNCTIONS (Corrected weighted-variance everywhere)
# ==============================================================================

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
  list(r2het = r2het, tau2_null = tau2_null, tau2_full = fit_full$tau2, p_eff = ncol(X), I2 = I2)
}

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

.penalized_meta_cv <- function(y, v, X, alpha, method = "REML", lambda_choice = "lambda.1se") {
  if (!check_data(y, v)) return(list(r2het = NA, r2het_cv = NA, lambda = NA,
                                     p_eff = NA, n_selected = NA, cv_fit = NULL))

  if (is_intercept_column(X)) X <- X[, -1, drop = FALSE]
  X <- drop_constant(X)
  if (is.null(X) || ncol(X) == 0) {
    return(list(r2het = 0, r2het_cv = 0, lambda = NA, p_eff = 0, n_selected = 0, cv_fit = NULL))
  }

  p_eff <- ncol(X)
  penalty.factor <- rep(1, ncol(X))
  if (ncol(X) == 1) {
    X <- cbind(X, noise1 = rnorm(nrow(X), 0, 0.001), noise2 = rnorm(nrow(X), 0, 0.001))
    penalty.factor <- c(1, 10000, 10000)
  }

  Xs <- scale(X); Xs[is.na(Xs)] <- 0

  fit_null <- rma(yi = y, vi = v, method = method)
  tau2_null <- fit_null$tau2
  w <- 1 / (v + tau2_null); w <- w / sum(w)
  k <- length(y)

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

  yhat_in <- as.vector(predict(cv_fit, newx = Xs, s = lam))
  res_in <- y - yhat_in
  vbar_w <- sum(w * v)
  tau2_in <- max(0, sum(w * res_in^2) - vbar_w)
  r2_in <- safe_r2het(tau2_in, tau2_null)

  idx <- which.min(abs(log(cv_fit$lambda) - log(lam)))
  yhat_oof <- cv_fit$fit.preval[, idx]
  res_oof <- y - yhat_oof
  tau2_cv <- max(0, sum(w * res_oof^2) - vbar_w)
  r2_cv <- safe_r2het(tau2_cv, tau2_null)

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

    if (!ok_fold) { preds[i] <- NA_real_; next }

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

permutation_test <- function(y, v, X, method = "REML", n_perm = 100) {
  if (!check_data(y, v)) return(list(p_value = NA, overfit_ratio = NA))

  obs <- lasso_meta_cv(y, v, X, method, lambda_choice = "lambda.min")
  obs_r2 <- obs$r2het

  Xperm_base <- if (is_intercept_column(X)) X[, -1, drop = FALSE] else X
  perm_r2 <- numeric(n_perm)

  for (p in seq_len(n_perm)) {
    Xp <- Xperm_base[sample(nrow(Xperm_base)), , drop = FALSE]
    Xp <- if (is_intercept_column(X)) cbind(1, Xp) else Xp
    perm_r2[p] <- lasso_meta_cv(y, v, Xp, method, lambda_choice = "lambda.min")$r2het
  }
  pval <- mean(perm_r2 >= obs_r2)
  list(observed_r2 = obs_r2,
       perm_mean = mean(perm_r2), perm_sd = sd(perm_r2),
       p_value = pval,
       overfit_ratio = ifelse(mean(perm_r2) > 0, obs_r2 / mean(perm_r2), NA))
}

plot_cv_curve <- function(cv_fit, main = "Flat CV Surface Enables Noise Capitalization") {
  plot(cv_fit)
  abline(v = log(cv_fit$lambda.min), lty = 2, col = "red")
  abline(v = log(cv_fit$lambda.1se), lty = 2, col = "blue")
  legend("topleft", c("lambda.min", "lambda.1se"), col = c("red", "blue"), lty = 2, bty = "n")
  title(main = main)
}

# ==============================================================================
# COMPLETE ANALYSIS WRAPPER
# ==============================================================================

complete_analysis <- function(data_name, y, v, X, save_results = TRUE) {
  cc <- complete.cases(y, v, X)
  y <- y[cc]; v <- v[cc]; X <- X[cc, , drop = FALSE]

  cat("\n=== Analyzing:", data_name, "===\n")
  k <- length(y)
  p <- max(0, ncol(X) - as.integer(is_intercept_column(X)))
  cat(sprintf("k=%d, p=%d, k/p=%s\n", k, p, ifelse(p > 0, sprintf("%.1f", k / p), "NA")))

  std <- calculate_r2het(y, v, X)
  std_cv <- cv_r2het(y, v, X)

  l_min <- lasso_meta_cv(y, v, X, lambda_choice = "lambda.min")
  l_1se <- lasso_meta_cv(y, v, X, lambda_choice = "lambda.1se")
  r_min <- ridge_meta_cv(y, v, X, lambda_choice = "lambda.min")
  r_1se <- ridge_meta_cv(y, v, X, lambda_choice = "lambda.1se")

  nested_lasso <- nested_cv_lasso(y, v, X, lambda_choice = "lambda.min")
  perm_result <- if (check_data(y, v) && p > 0) permutation_test(y, v, X, n_perm = 100) else list(p_value = NA, overfit_ratio = NA)

  cat("\nOptimism Summary (Apparent - CV, in %):\n")
  cat(sprintf("  Standard:            %.1f\n", (std$r2het - std_cv) * 100))
  cat(sprintf("  LASSO lambda.min:    %.1f (lambda=%.4f, selected=%d)\n",
              (l_min$r2het - l_min$r2het_cv) * 100, l_min$lambda, l_min$n_selected))
  cat(sprintf("  LASSO lambda.1se:    %.1f (lambda=%.4f, selected=%d)\n",
              (l_1se$r2het - l_1se$r2het_cv) * 100, l_1se$lambda, l_1se$n_selected))
  cat(sprintf("  Ridge  lambda.min:   %.1f\n", (r_min$r2het - r_min$r2het_cv) * 100))
  cat(sprintf("  Ridge  lambda.1se:   %.1f\n", (r_1se$r2het - r_1se$r2het_cv) * 100))

  if (is.list(perm_result) && !is.na(perm_result$p_value)) {
    cat(sprintf("\nPermutation test p-value: %.3f\n", perm_result$p_value))
    cat(sprintf("Overfit ratio: %.2fx null expectation\n", perm_result$overfit_ratio))
  }

  out <- data.frame(
    Dataset = data_name, k = k, p = p, k_p_ratio = ifelse(p > 0, k / p, NA_real_), I2 = std$I2,
    Std_Apparent = std$r2het * 100, Std_CV = std_cv * 100, Std_Optimism = (std$r2het - std_cv) * 100,
    LASSO_min_Apparent = l_min$r2het * 100, LASSO_min_CV = l_min$r2het_cv * 100,
    LASSO_min_Optimism = (l_min$r2het - l_min$r2het_cv) * 100, LASSO_min_lambda = l_min$lambda,
    LASSO_min_selected = l_min$n_selected,
    LASSO_1se_Apparent = l_1se$r2het * 100, LASSO_1se_CV = l_1se$r2het_cv * 100,
    LASSO_1se_Optimism = (l_1se$r2het - l_1se$r2het_cv) * 100, LASSO_1se_lambda = l_1se$lambda,
    LASSO_1se_selected = l_1se$n_selected,
    Ridge_min_Apparent = r_min$r2het * 100, Ridge_min_CV = r_min$r2het_cv * 100,
    Ridge_min_Optimism = (r_min$r2het - r_min$r2het_cv) * 100, Ridge_min_lambda = r_min$lambda,
    Ridge_1se_Apparent = r_1se$r2het * 100, Ridge_1se_CV = r_1se$r2het_cv * 100,
    Ridge_1se_Optimism = (r_1se$r2het - r_1se$r2het_cv) * 100, Ridge_1se_lambda = r_1se$lambda,
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

  # Return results plus cv objects for external storage
  attr(out, "cv_curve_lasso_min") <- l_min$cv_fit
  attr(out, "cv_curve_lasso_1se") <- l_1se$cv_fit
  out
}

# ==============================================================================
# DATASET LOADERS
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
  keep <- complete.cases(dat.konstantopoulos2011$yi, dat.konstantopoulos2011$vi, dat.konstantopoulos2011$year)
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

load_hackshaw <- function() {
  data(dat.hackshaw1998)
  df <- dat.hackshaw1998
  keep <- complete.cases(df$yi, df$vi)
  df <- df[keep, ]

  if (all(c("continent", "design") %in% names(df))) {
    X <- model.matrix(~ continent + design, data = df)
  } else {
    message("Hackshaw: 'continent' or 'design' not found; falling back to ~ year or intercept-only.")
    if ("year" %in% names(df)) X <- model.matrix(~ year, data = df) else X <- model.matrix(~ 1, data = df)
  }
  list(y = df$yi, v = df$vi, X = X, name = "Hackshaw")
}

# ==============================================================================
# RUN ALL ANALYSES (store CV curves in side lists)
# ==============================================================================

cat("========== REGULARIZATION PARADOX IN META-REGRESSION ==========\n\n")

loaders <- list(load_bcg, load_teacher, load_konstant, load_berkey, load_hackshaw)

res_list <- list()
cv_curves_min  <- list()
cv_curves_1se  <- list()

for (loader in loaders) {
  ds <- loader()
  res <- complete_analysis(ds$name, ds$y, ds$v, ds$X, save_results = TRUE)
  res_list[[ds$name]] <- res
  cv_curves_min[[ds$name]] <- attr(res, "cv_curve_lasso_min")
  cv_curves_1se[[ds$name]] <- attr(res, "cv_curve_lasso_1se")
}

all_results <- do.call(rbind, res_list)

# Strip list columns for display
display_results <- all_results[, !sapply(all_results, is.list), drop = FALSE]

# ==============================================================================
# SUMMARY TABLE & MEANS
# ==============================================================================

cat("\n========== TABLE 1: THE REGULARIZATION PARADOX ==========\n\n")
summary_cols <- c("Dataset", "k", "p", "k_p_ratio", "I2",
                  "Std_Optimism", "LASSO_min_Optimism", "LASSO_1se_Optimism",
                  "Ridge_min_Optimism", "Ridge_1se_Optimism")
print(kable(display_results[, summary_cols], caption = "Optimism in R^2_het (%) by Method"))

cat("\n========== MEAN OPTIMISM ACROSS DATASETS ==========\n")
cat(sprintf("Standard:          %.1f%%\n", mean(display_results$Std_Optimism, na.rm = TRUE)))
cat(sprintf("LASSO lambda.min:  %.1f%%\n", mean(display_results$LASSO_min_Optimism, na.rm = TRUE)))
cat(sprintf("LASSO lambda.1se:  %.1f%%\n", mean(display_results$LASSO_1se_Optimism, na.rm = TRUE)))
cat(sprintf("Ridge  lambda.min: %.1f%%\n", mean(display_results$Ridge_min_Optimism, na.rm = TRUE)))
cat(sprintf("Ridge  lambda.1se: %.1f%%\n", mean(display_results$Ridge_1se_Optimism, na.rm = TRUE)))

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
  geom_bar(stat = "identity", fill = c("gray50", "red3", "blue3", "red3", "blue3")) +
  geom_hline(yintercept = optimism_data$Optimism[1], linetype = "dashed", alpha = 0.5) +
  labs(title = "The Regularization Paradox",
       subtitle = "lambda.min increases optimism; lambda.1se provides minimal benefit",
       y = "Mean Optimism in R^2_het (%)",
       x = "") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 10))

print(p1)

# Figure 2: Optimism vs k/p ratio (linear smoothing only if >=5 points)
plot_data <- display_results[!is.na(display_results$k_p_ratio), ]
plot_long <- melt(plot_data[, c("Dataset", "k_p_ratio", "Std_Optimism",
                                "LASSO_min_Optimism", "LASSO_1se_Optimism")],
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

print(p2)

# Save artifacts
write.csv(display_results, "regularization_paradox_results.csv", row.names = FALSE)
safe_ggsave("fig1_optimism_comparison.pdf", p1, width = 6, height = 4)
safe_ggsave("fig2_optimism_by_ratio.pdf", p2, width = 7, height = 5)

# ==============================================================================
# CV CURVE PLOT (pull from stored list, not from rbind-ed data)
# ==============================================================================

# Grab the first non-null lambda.min CV object and plot it
cv_min_first <- NULL
for (nm in names(cv_curves_min)) {
  if (!is.null(cv_curves_min[[nm]])) { cv_min_first <- cv_curves_min[[nm]]; break }
}
if (!is.null(cv_min_first)) {
  plot_cv_curve(cv_min_first, main = "Flat CV Surface Enables Noise Capitalization (Example)")
} else {
  message("No available cv.glmnet object to plot.")
}

cat("\n========== DONE ==========\n")
