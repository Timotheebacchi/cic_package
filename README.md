  # cic.newassumptions.newvarianceestimator: Changes-in-Changes Estimator

  `cic.newassumptions.newvarianceestimator` is an R package for the Changes-in-Changes estimator and asymptotic inference for empirical quantile-based estimators. It computes the plug-in estimate, provides several confidence interval methods. This package is based on the inference methods proposed in https://arxiv.org/abs/2607.00219

  ## Features

  - Point estimation of the CiC parameter from `Y`, `X`, and `Z`
  - Confidence intervals with five methods:
    - `"no-split"`
    - `"split"`
    - `"kde"`
    - `"bse"`
    - `"bpc"`
  - Optional `panel_data = TRUE` workflow for paired `(Y, Z)` samples
  - Rcpp-backed computation for the core routines, with a pure R fallback where available
  - Simulation helpers `sim_dgp()`, `qY_dgp()`, and `theta_true()`

  ## Installation

  ```r
  devtools::install_github("Timotheebacchi/cic_package")
  ```

  ## Quick Start

  ```r
  library(cic.newassumptions.newvarianceestimator)

  set.seed(2026) #To match the code of the manuscript
  d <- sim_dgp(2000)
  
  fit <- cic(d$Y, d$X, d$Z, method = c("no-split", "split", "kde", "bse", "bpc"))
  summary(fit)
  ```

  ## Validation and warnings

  The package performs lightweight input validation and emits clear warnings for common
  issues (e.g., non-numeric inputs, mismatched lengths, or invalid bandwidth
  values). Use `sim_dgp()` and the estimation `cic()` for simulation and
  inference; inspect warnings to help debug input problems.

  ## Assumptions Reminder

  The package is built for the CiC setup described in [the manuscript](https://arxiv.org/abs/2607.00219). Before
  interpreting the output, the input data should be checked against the main
  assumptions used by the estimator:

  - the observations should look approximately i.i.d. and continuous enough for
    the rank transformation to make sense;
  - the empirical quantile function of $Y$ should stay below a boundary envelope
    of the form $C_Y t^{-d_1}(1-t)^{-d_2}$ on the interior of $(0,1)$;
  - the transformed covariate ranks $U = F_Z(X)$ should admit a smooth density
    that can be screened against $C_U u^{-b_1}(1-u)^{-b_2}$;
  - the outcome distribution should have tails that are not too heavy;
  - the combined tail and boundary behavior should stay within the rate
    conditions required by the asymptotic theory;
  - the smoothing bandwidth should be reasonable for the sample size.

  The package provides
  only lightweight validation warnings; inspect warnings produced by `cic()`
  when running your data.
 Use `summary(fit)` to view
  the estimation and inference table.

  In practice, `cic()` can receive one method of inference or several methods at once, and it
  returns the intervals in the requested order. The asymptotic methods are
  `no-split`, `split`, and `kde`, while `bse` and `bpc` are bootstrap
  comparisons that help assess robustness in small samples.

  
  ## Package Contents

  - `cic()` for estimation and confidence intervals
  - `sim_dgp()` for a reproducible data-generating process
  - `qY_dgp()` and `theta_true()` for simulation support
  

  ## Reference

  Chhor, J., D'Haultfoeuille, X., L'Hour, J., & Mugnier, M. (2026). Asymptotic Properties of Empirical Quantile-Based Estimators. Manuscript.

  ## License

  MIT + file LICENSE

  ## Author

  Timothée Bacchi
