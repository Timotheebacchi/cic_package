#' @useDynLib cic, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import stats
NULL

#' Changes-in-Changes estimator
#'
#' Computes the Changes-in-Changes (CiC) estimator
#' \deqn{\hat\theta = \frac{1}{n_2}\sum_{j=1}^{n_2} \hat F_Y^{-1}(\hat F_Z(X_j))}
#' and returns asymptotic confidence intervals based on the chosen variance
#' estimator(s).
#'
#' @param Y Numeric vector. Outcome sample of length \eqn{n_1}.
#' @param X Numeric vector. Covariate sample of length \eqn{n_2}.
#' @param Z Numeric vector. Instrument sample of length \eqn{n_3}.
#' @param method Character vector. One or more of `"no-split"`, `"bse"`, `"bpc"`.
#'   Defaults to all three.
#'   * `"no-split"` : asymptotic CI using the full-sample (no-split) variance estimator.
#'   * `"bse"`      : bootstrap CI using the standard-error bootstrap.
#'   * `"bpc"`      : bootstrap CI using the percentile bootstrap.
#' @param B Integer. Number of bootstrap replications. Ignored if neither
#'   `"bse"` nor `"bpc"` is in `method`. Must be >= 200. Default: 999.
#' @param h Numeric or `NULL`. Bandwidth for the Epanechnikov KDE used in the
#'   `"no-split"` (full-sample) variance estimator. If `NULL` (default), Silverman's
#'   rule-of-thumb is applied automatically.
#' @param level Numeric in (0, 1). Confidence level. Default: 0.95.
#'
#' @return An object of class `"cic"`, a list containing:
#'   * `theta_hat`  : point estimate \eqn{\hat\theta}.
#'   * `ci`         : data frame with columns `method`, `lower`, `upper`,
#'     `length` for each requested method.
#'   * `level`      : confidence level used.
#'   * `n`          : named integer vector `c(n1, n2, n3)`.
#'   * `method`     : methods requested.
#'   * `h`          : bandwidth used (NA if only bootstrap methods requested).
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' Y <- rnorm(n)
#' Z <- rnorm(n)
#' X <- qnorm(pbeta(pnorm(Z), 0.95, 0.95))
#' fit <- cic(Y, X, Z)
#' print(fit)
#' summary(fit)
#'
#' @export
cic <- function(Y, X, Z,
                method = c("no-split", "bse", "bpc"),
                B      = 999L,
                h      = NULL,
                level  = 0.95) {

  # ── Input checks ─────────────────────────────────────────────────────────────
  method <- match.arg(method, several.ok = TRUE)
  stopifnot(
    is.numeric(Y), is.numeric(X), is.numeric(Z),
    length(Y) >= 4, length(X) >= 4, length(Z) >= 4,
    is.null(h) || (is.numeric(h) && length(h) == 1 && h > 0),
    is.numeric(level), level > 0, level < 1
  )
  if (any(c("bse", "bpc") %in% method)) {
    B <- as.integer(B)
    if (B < 200L) {
      warning("B < 200: the bootstrap standard error may be unstable. ",
              "Using B = 200.")
      B <- 200L
    }
  }

  n1 <- length(Y); n2 <- length(X); n3 <- length(Z)
  alpha <- (1 - level) / 2
  z_a   <- stats::qnorm(1 - alpha)

  # ── Point estimate ────────────────────────────────────────────────────────────
  FZ             <- stats::ecdf(Z)
  Uhat           <- FZ(X)
  qcdf_transform <- .prepare_left_quantile(Y)(Uhat)
  theta_hat      <- mean(qcdf_transform)
  eps_hat        <- mean((theta_hat - qcdf_transform)^2)

  # ── Bandwidth ─────────────────────────────────────────────────────────────────
  if (is.null(h)) h <- .default_bandwidth(Y)

  # ── Variance estimation ───────────────────────────────────────────────────────
  ci_rows <- list()
  N        <- min(n1, n2, n3)
  lbda1_3  <- N * (n1 + n3) / (n1 * n3)
  lbda2    <- N / n2

  if ("no-split" %in% method) {
    # Use the user-provided or default Silverman bandwidth in the no-split
    # adaptive density estimator. This makes fit$h the actual smoothing
    # parameter used by the full-sample variance estimator.
    eps_bw <- h

    # No-split: use full samples, not divided into halves
    Ysortdiff <- diff(sort(Y))
    FYhat     <- seq_len(n1 - 1) / n1
    
    est_full   <- .make_density_estimator(sort(Uhat), FYhat)
    fUhat      <- est_full$estimate(eps_bw, pointwise = 1)
    fUhat_unif <- est_full$estimate(eps_bw, pointwise = 0)
    fUhat_t2   <- est_full$estimate(2 * eps_bw, pointwise = 1)
    fUhat_d2   <- est_full$estimate(eps_bw / 2, pointwise = 1)

    eta_hat    <- .fast_eta(Ysortdiff, fUhat,      Ysortdiff, fUhat,      FYhat)
    eta_unif   <- .fast_eta(Ysortdiff, fUhat_unif, Ysortdiff, fUhat_unif, FYhat)
    eta_t2     <- .fast_eta(Ysortdiff, fUhat_t2,   Ysortdiff, fUhat_t2,   FYhat)
    eta_d2     <- .fast_eta(Ysortdiff, fUhat_d2,   Ysortdiff, fUhat_d2,   FYhat)

    eta_nosplit <- (eta_hat + eta_unif + eta_t2 + eta_d2) / 4
    sigma_sq    <- lbda1_3 * eta_nosplit + lbda2 * eps_hat
    se          <- sqrt(sigma_sq / N)

    ci_rows[["no-split"]] <- data.frame(
      method = "no-split",
      lower  = theta_hat - z_a * se,
      upper  = theta_hat + z_a * se,
      length = 2 * z_a * se
    )
    rm(Ysortdiff, FYhat, fUhat, fUhat_unif, fUhat_t2, fUhat_d2, est_full)
  }

  if (any(c("bse", "bpc") %in% method)) {
    boot_vals <- boot_core(sort(Y), X, sort(Z), B = B)

    if ("bpc" %in% method) {
      q_lo <- stats::quantile(boot_vals, probs = alpha,     names = FALSE)
      q_hi <- stats::quantile(boot_vals, probs = 1 - alpha, names = FALSE)
      ci_rows[["bpc"]] <- data.frame(
        method = "bpc",
        lower  = q_lo,
        upper  = q_hi,
        length = q_hi - q_lo
      )
    }

    if ("bse" %in% method) {
      se_boot <- stats::sd(boot_vals)
      ci_rows[["bse"]] <- data.frame(
        method = "bse",
        lower  = theta_hat - z_a * se_boot,
        upper  = theta_hat + z_a * se_boot,
        length = 2 * z_a * se_boot
      )
    }
  }

  ci <- do.call(rbind, ci_rows[method])  # preserve requested order
  rownames(ci) <- NULL

  # ── Output ───────────────────────────────────────────────────────────────────
  structure(
    list(
      theta_hat = theta_hat,
      ci        = ci,
      level     = level,
      n         = c(n1 = n1, n2 = n2, n3 = n3),
        method    = method,
      h         = if ("no-split" %in% method) h else NA_real_
    ),
    class = "cic"
  )
}
