# cic: Changes-in-Changes Estimator

An R package implementing the **Changes-in-Changes (CiC) estimator** from Athey & Imbens (2006) with asymptotic confidence intervals and bootstrap alternatives.

## Features

- **Core estimator**: Average of quantile-transformed ranks
- **Confidence intervals via three methods**:
  - Split-sample (nonparametric, fast)
  - Bootstrap standard-error (bse)
  - Bootstrap percentile (bpc)
- **Performance**: C++ backend via Rcpp
- **Validated**: Comprehensive Monte Carlo tests with known true parameter

## Installation

```r
# Install from GitHub
devtools::install_github("timbak/cic_package")
```

## Quick Start

```r
library(cic)

# Simulate from the DGP (Athey & Imbens 2006)
set.seed(42)
n <- 500
W <- runif(n)
Y <- -W^(-0.2) + (1-W)^(-0.05)  # Quantile function
X <- qnorm(rbeta(n, 1, 1.05))   # Covariate
Z <- rnorm(n)                     # Instrument

# Estimate with split-sample CI
fit <- cic(Y, X, Z, method = "split")

# View results
print(fit)

# Get confidence interval
fit$ci
```

## Methods

The `cic()` function supports:

| Method | Speed | Notes |
|--------|-------|-------|
| `"split"` | Fast | Nonparametric, requires even n |
| `"bse"` | Medium | Bootstrap standard-error |
| `"bpc"` | Medium | Bootstrap percentile |

## References

- Athey, S., & Imbens, G. W. (2006). Identification and inference in nonlinear difference-in-differences models. Econometrica, 74(2), 431–497.

## License

MIT + file LICENSE

## Author

Timothée Bacchi
