## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(cic.newassumptions.newvarianceestimator)

## ----simulate-data------------------------------------------------------------
set.seed(42)
d <- sim_dgp(500)

str(d)
head(d$Y)
head(d$X)
head(d$Z)

## ----diagnostics--------------------------------------------------------------
diag <- check_cic_assumptions(d$Y, d$X, d$Z)
diag$pass_all
diag$messages

## ----estimation---------------------------------------------------------------
fit <- cic(
  d$Y, d$X, d$Z,
  method = c("no-split", "split", "kde", "bse", "bpc"),
  B = 500
)

fit
fit$theta_hat
fit$ci

## ----compare-methods----------------------------------------------------------
fit$ci[, c("method", "lower", "upper", "length")]

## ----helpers------------------------------------------------------------------
coef(fit)
confint(fit)
summary(fit)

## ----plot-fit, fig.width=7, fig.height=3.5------------------------------------
plot(fit)

