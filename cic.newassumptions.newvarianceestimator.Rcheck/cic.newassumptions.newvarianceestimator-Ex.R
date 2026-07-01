pkgname <- "cic.newassumptions.newvarianceestimator"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('cic.newassumptions.newvarianceestimator')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cic")
### * cic

flush(stderr()); flush(stdout())

### Name: cic
### Title: Changes-in-Changes Estimator
### Aliases: cic

### ** Examples

set.seed(42)
d <- sim_dgp(500)
fit <- cic(d$Y, d$X, d$Z, method = "no-split")
summary(fit)



cleanEx()
nameEx("qY_dgp")
### * qY_dgp

flush(stderr()); flush(stdout())

### Name: qY_dgp
### Title: Quantile Function for the CiC Simulation DGP
### Aliases: qY_dgp

### ** Examples

qY_dgp(0.5, d1 = 0, d2 = 0.05)



cleanEx()
nameEx("sim_dgp")
### * sim_dgp

flush(stderr()); flush(stdout())

### Name: sim_dgp
### Title: Simulate Data from the CiC Simulation DGP
### Aliases: sim_dgp

### ** Examples

set.seed(42)
d <- sim_dgp(500)
head(d)

set.seed(42)
d_panel <- sim_dgp(500, panel_data = TRUE)
diag_panel <- check_cic_assumptions(d_panel$Y, d_panel$X, d_panel$Z, panel_data = TRUE)
diag_panel$pass_all

set.seed(42)
d_fail <- sim_dgp(500)
check_cic_assumptions(d_fail$Y, d_fail$X, d_fail$Z)$pass_all



cleanEx()
nameEx("theta_true")
### * theta_true

flush(stderr()); flush(stdout())

### Name: theta_true
### Title: True CiC Parameter for the Simulation DGP
### Aliases: theta_true

### ** Examples

theta_true(b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
