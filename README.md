# Regularized meta-regression (“Lassopaper”) — exploratory methods project

**Status: not publication-ready.** The R script compares conventional,
LASSO and ridge meta-regression using a small set of `metadat` example
datasets and simulations. Its working title/conclusions are **hypotheses**:
the claimed “regularization paradox” has not been independently validated.

## Reproducible checks

```bash
# Basic syntax and pure-math regression tests (no CRAN packages required):
Rscript -e "parse(file='lasso.R')"
Rscript tests/test_validation_metrics.R

# Full analyses require compatible installed R packages listed in lasso.R:
Rscript lasso.R
```

The script currently auto-runs its analyses when invoked and may require
substantial CPU/time. It is **not** currently a package with isolated,
independently exercised modelling functions. Dependencies are not pinned.

## Specific audit findings and changes

- Previously, the out-of-fold explained-heterogeneity estimate used
  `safe_r2het`, clipping negative out-of-fold values to zero. This
  obscured performance *worse than the null*. An `oof_r2het` function now
  allows negative skill and returns NA when the null heterogeneity is zero.
  Apparent `R²_het` remains a separate, bounded descriptive quantity.
- The nested LASSO outer loop previously computed training weights from a
  null-model heterogeneity estimate fitted to **all** outcomes, including the
  held-out outcome. Training-fold weights now use a null model fitted only
  to each training split.
- The empirical permutation tail calculation previously could return
  `p=0` with a finite number of shuffles. It now uses the finite-sample
  `(1+b)/(B+1)` correction. This is not sufficient by itself to validate
  an exchangeable null hypothesis or deal with model-selection issues.

## Remaining publication blockers

1. **Independent statistical validation:** compare every metric, intercept
   convention, design matrix and fold prediction with independent R
   implementations, including null `tau²=0`, p near k and rank deficiency.
2. **Unresolved selection/weighting optimism:** the ordinary penalized
   `cv.glmnet` calculation still uses full-data heterogeneity for weights and
   derives tuning and apparent/out-of-fold summaries from the same process.
   Treat these as exploratory, *not* as fully external validation. The
   nested-CV estimator also needs larger-sample benchmarking and
   preprocessing-leakage review.
3. **Estimand clarity:** the residual-moment `tau²_cv` approximation and
   bounded apparent `R²_het` are not interchangeable estimands. Provide
   a formal definition and uncertainty, demonstrate its behavior when
   true heterogeneity is small or zero, and report negative validation skill.
4. **Small study-count evidence:** empirical comparisons cover only five
   datasets and can be unstable. Independent datasets, multi-seed simulation
   sweeps and sensitivity to missing moderators are needed.
5. **Reproducibility and reporting:** pin R/package versions, lock datasets,
   predefine evaluation design and primary outcomes, produce figures/tables
   from a clean checkout, and include a versioned protocol and test archive.
6. **No inferential conclusions yet:** multiple testing, simulation design,
   null exchangeability, uncorrected selection effects and robustness of any
   claimed method superiority need independent assessment.

The new GitHub Actions workflow **only checks R syntax and pure-math helper
functions**. It does *not* run full computational analyses; do not cite a green
check as proof of a reproducible published result.
