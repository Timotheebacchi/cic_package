  # cic: Changes-in-Changes Estimator

  `cic` is an R package for the Changes-in-Changes estimator and asymptotic inference for empirical quantile-based estimators. It computes the plug-in estimate, provides several confidence interval methods, and includes a diagnostic helper to check whether the input data look compatible with the model assumptions.

  ## Features

  - Point estimation of the CiC parameter from `Y`, `X`, and `Z`
  - Confidence intervals with five methods:
    - `"no-split"`
    - `"split"`
    - `"kde"`
    - `"bse"`
    - `"bpc"`
  - Diagnostic checks via `check_cic_assumptions()`
  - Rcpp-backed computation for the core routines, with a pure R fallback where available
  - Simulation helpers `sim_dgp()`, `qY_dgp()`, and `theta_true()`

  ## Installation

  ```r
  devtools::install_github("Timotheebacchi/cic_package")
  ```

  ## Quick Start

  ```r
  library(cic)

  set.seed(42)
  d <- sim_dgp(500)

  diag <- check_cic_assumptions(d$Y, d$X, d$Z)
  diag$pass_all

  fit <- cic(d$Y, d$X, d$Z, method = c("no-split", "split", "kde", "bse", "bpc"))

  fit
  fit$ci
  summary(fit)
  ```

  ## Diagnostics

  `check_cic_assumptions()` returns a list with:

  - `pass_all`: overall logical result
  - `metrics`: sample ratios, tail indices, boundary estimates, and related checks
  - `messages`: warnings or success messages

  It is designed as a quick pre-check before running `cic()` on empirical data.

  ## Mathematical Reminder for the Methods

  Let $F_Z$ denote the empirical distribution function of $Z$ and let $Q_Y(u)=F_Y^{-1}(u)$ be the left-continuous quantile of $Y$. The common point estimator used by all methods is

  $$
  \hat\theta = \frac{1}{n}\sum_{i=1}^n Q_Y(F_Z(X_i)).
  $$

  The confidence interval is then built around $\hat\theta$ using either variance estimators or bootstrap. In compact notation, the asymptotic methods use a form of

  $$
  \hat\theta \pm z_{1-\alpha}\,\widehat{\mathrm{se}},
  $$

  where $\alpha = (1-\text{level})/2$ and $z_{1-\alpha}$ is the standard normal quantile. The code also separates the variance contributions via

  $$
  \hat\sigma^2 = \lambda_{1,3}\,\hat\eta + \lambda_2\,\hat\varepsilon,
  $$

  with $\hat\varepsilon = n^{-1}\sum_i (\hat\theta - Q_Y(F_Z(X_i)))^2$, $\lambda_{1,3} = N(n_1+n_3)/(n_1n_3)$, $\lambda_2 = N/n_2$, and $N = \min(n_1,n_2,n_3)$.

  | Method | What it does | Main difference |
  |---|---|---|
  | "no-split" | Estimates the nonparametric variance component on the full sample, with smoothing controlled by `h`. | It is the most direct asymptotic method and uses all observations at once. |
  | "split" | Splits the samples into two halves, re-estimates the density components on each half, and then builds the interval with a sample-splitting variance. | It reduces data reuse bias, but typically loses information and increases variance. |
  | "kde" | Approximates the density using an Epanechnikov kernel applied to the scores $Q_Y(F_Z(X_i))$. | It replaces local counting with kernel density smoothing and provides a smooth alternative. |
  | "bse" | Resamples the triplets $(Y,X,Z)$, computes $\hat\theta^*$ at each replication, and uses the bootstrap standard deviation as the standard error. | The interval is centered on $\hat\theta$ and reflects the empirical bootstrap dispersion. |
  | "bpc" | Resamples as in "bse", but uses the empirical quantiles of $\hat\theta^*$ directly to build the interval. | The interval is not forced to be symmetric around $\hat\theta$ and better reflects finite-sample asymmetry. |

  In practice, `cic()` can receive one method or several methods at once, and it returns the intervals in the requested order. As a quick reading guide: `no-split`, `split`, and `kde` are the three main asymptotic variants, while `bse` and `bpc` are bootstrap comparisons that help assess robustness in small samples.

  ## Package Contents

  - `cic()` for estimation and confidence intervals
  - `check_cic_assumptions()` for diagnostics
  - `sim_dgp()` for a reproducible data-generating process
  - `qY_dgp()` and `theta_true()` for simulation support

  ## Reference

  Chhor, J., D'Haultfoeuille, X., L'Hour, J., & Mugnier, M. (2026). Asymptotic Properties of Empirical Quantile-Based Estimators. Manuscript.

  ## License

  MIT + file LICENSE

  ## Author

  Timothée Bacchi
