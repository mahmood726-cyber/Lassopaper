# Validation helpers that do not depend on metafor/glmnet.
# Out-of-fold skill estimates must be allowed to go negative.
oof_r2het <- function(tau2_oof, tau2_null, eps = 1e-10) {
  if (length(tau2_oof) != 1L || length(tau2_null) != 1L ||
      !is.finite(tau2_oof) || !is.finite(tau2_null) || tau2_null <= eps ||
      tau2_oof < 0) {
    return(NA_real_)
  }
  1 - tau2_oof / tau2_null
}

# The usual permutation plus-one correction avoids impossible p=0.
permutation_pvalue <- function(observed, permuted) {
  if (length(observed) != 1L || !is.finite(observed)) return(NA_real_)
  good <- permuted[is.finite(permuted)]
  if (!length(good)) return(NA_real_)
  (1 + sum(good >= observed)) / (length(good) + 1)
}
