# ── print.cic ──────────────────────────────────────────────────────────────────

#' @export
print.cic <- function(x, digits = 4, ...) {
  cat("Changes-in-Changes Estimator\n")
  cat(rep("-", 44), "\n", sep = "")
  cat(sprintf("Estimate (theta_hat) : %.4f\n", x$theta_hat))
  cat(sprintf("Confidence level     : %.0f%%\n", x$level * 100))
  cat(sprintf("Sample sizes         : n1 = %d, n2 = %d, n3 = %d\n",
              x$n["n1"], x$n["n2"], x$n["n3"]))
  if (!is.na(x$h))
    cat(sprintf("Bandwidth (h)        : %.4f\n", x$h))
  cat(rep("-", 44), "\n", sep = "")
  cat("Confidence intervals:\n\n")

  # Format the CI table
  ci <- x$ci
  ci$lower  <- round(ci$lower,  digits)
  ci$upper  <- round(ci$upper,  digits)
  ci$length <- round(ci$length, digits)

  label_map <- c(`no-split` = "No-split", bpc = "Bootstrap (pct.)",
                 bse   = "Bootstrap (SE)")
  ci$method <- label_map[ci$method]

  print(ci, row.names = FALSE, right = FALSE)
  invisible(x)
}

# ── summary.cic ────────────────────────────────────────────────────────────────

#' @export
summary.cic <- function(object, digits = 4, ...) {
  cat("Changes-in-Changes - Summary\n")
  cat(rep("=", 44), "\n", sep = "")
  cat(sprintf("Point estimate : %.4f\n\n", object$theta_hat))

  cat("Sample sizes:\n")
  cat(sprintf("  Y sample (n1) : %d\n", object$n["n1"]))
  cat(sprintf("  X sample (n2) : %d\n", object$n["n2"]))
  cat(sprintf("  Z sample (n3) : %d\n", object$n["n3"]))

  cat(sprintf("\nConfidence level : %.0f%%\n", object$level * 100))

  if (!is.na(object$h))
    cat(sprintf("Bandwidth (h)    : %.4f  [Silverman rule-of-thumb]\n",
                object$h))

  cat(rep("-", 44), "\n", sep = "")
  cat("Confidence intervals:\n\n")

  ci <- object$ci
  label_map <- c(`no-split` = "No-split        ",
                 bpc   = "Bootstrap (pct.)",
                 bse   = "Bootstrap (SE)  ")
  for (i in seq_len(nrow(ci))) {
    cat(sprintf("  %s  [%.*f, %.*f]  (length %.4f)\n",
                label_map[ci$method[i]],
                digits, ci$lower[i],
                digits, ci$upper[i],
                ci$length[i]))
  }
  cat(rep("=", 44), "\n", sep = "")
  invisible(object)
}

#' @export
coef.cic <- function(object, ...) {
  object$theta_hat
}

#' @export
confint.cic <- function(object, parm = NULL, level = object$level, ...) {
  if (!is.null(level) && !identical(level, object$level)) {
    warning(
      "`level` differs from the confidence level used to fit the object. ",
      "Returning stored intervals from `object$level`."
    )
  }
  ci <- object$ci[, c("lower", "upper")]
  rownames(ci) <- object$ci$method
  as.matrix(ci)
}

#' @export
plot.cic <- function(x, y = NULL, ..., col = c("black", "blue"), pch = 19, lwd = 2) {
  ci <- x$ci
  n <- nrow(ci)
  y_pos <- seq_len(n)
  x_min <- min(ci$lower, x$theta_hat)
  x_max <- max(ci$upper, x$theta_hat)
  x_range <- if (x_max > x_min) x_max - x_min else 1
  x_lim <- c(x_min - 0.05 * x_range, x_max + 0.05 * x_range)

  plot(
    NA, NA,
    xlim = x_lim,
    ylim = c(0.5, n + 0.5),
    yaxt = "n",
    ylab = "",
    xlab = "Estimate",
    main = "CiC confidence intervals",
    ...
  )
  axis(2, at = y_pos, labels = ci$method)
  segments(ci$lower, y_pos, ci$upper, y_pos, lwd = lwd, col = col[1])
  points(rep(x$theta_hat, n), y_pos, pch = pch, col = col[2])
  abline(v = x$theta_hat, lty = 2, col = "gray70")
  invisible(x)
}
