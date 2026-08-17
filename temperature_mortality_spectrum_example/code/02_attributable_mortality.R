# 02_attributable_mortality.R
# Example backward attributable-risk calculation using the model from script 01.
# Uses synthetic data only; outputs are not study estimates.

rm(list = ls())
set.seed(20260817)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(tibble)
  library(dlnm)
  library(splines)
})

### ---------- 1. Paths and settings ----------

DATA_FILE <- file.path("data", "example_daily_data.csv")
MODEL_FILE <- file.path("output", "example_model_package.rds")
OUT_DIR <- "output"
NSIM <- 1000
CI_PROBS <- c(0.025, 0.975)

if (!file.exists(MODEL_FILE)) {
  stop("Missing model package. Run code/01_temperature_mortality_model.R first.")
}
if (!file.exists(DATA_FILE)) stop("Missing: ", DATA_FILE)

pkg <- readRDS(MODEL_FILE)
TEMP_VAR <- pkg$temp_var
OUTCOME <- pkg$outcome
MAIN_LAG <- pkg$main_lag
argvar_main <- pkg$argvar_main
arglag_main <- pkg$arglag_main

### ---------- 2. Rebuild the fitted analysis rows ----------

dat <- read_csv(DATA_FILE, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(date),
    pref = factor(pref),
    holiday = factor(holiday),
    year = year(date),
    month = month(date),
    dow = factor(wday(date, label = TRUE, week_start = 1)),
    stratum = interaction(pref, year, month, dow, drop = TRUE),
    y = as.numeric(.data[[OUTCOME]])
  ) %>%
  arrange(pref, date) %>%
  group_by(stratum) %>%
  mutate(ind = sum(y, na.rm = TRUE)) %>%
  ungroup()

make_lag_matrix_by_group <- function(x, group, lag) {
  out <- matrix(NA_real_, nrow = length(x), ncol = lag + 1L)
  colnames(out) <- paste0("lag", 0:lag)
  for (idx in split(seq_along(x), group, drop = TRUE)) {
    x_i <- x[idx]
    for (lag_i in 0:lag) out[idx, lag_i + 1L] <- dplyr::lag(x_i, n = lag_i)
  }
  out
}

cb <- crossbasis(
  x = dat[[TEMP_VAR]],
  lag = MAIN_LAG,
  group = dat$pref,
  argvar = argvar_main,
  arglag = arglag_main
)

x_lag_all <- make_lag_matrix_by_group(dat[[TEMP_VAR]], dat$pref, MAIN_LAG)
fit_rows <- dat$ind > 0 & complete.cases(cb) & complete.cases(x_lag_all) & !is.na(dat$y)

if (sum(fit_rows) != pkg$n_fit_rows) {
  warning("Rebuilt fitted rows differ from script 01: ", sum(fit_rows), " vs ", pkg$n_fit_rows)
}

x_lag_fit <- x_lag_all[fit_rows, , drop = FALSE]
cases_fit <- dat$y[fit_rows]

### ---------- 3. Attributable-risk helpers ----------

rmvnorm_base <- function(n, mean, sigma) {
  p <- length(mean)
  sigma <- (sigma + t(sigma)) / 2
  eg <- eigen(sigma, symmetric = TRUE)
  values <- pmax(eg$values, 0)
  A <- eg$vectors %*% diag(sqrt(values), nrow = p, ncol = p)
  z <- matrix(rnorm(p * n), nrow = p, ncol = n)
  sims <- t(A %*% z)
  sweep(sims, 2, mean, "+")
}

calc_af_an_from_basis <- function(X, cases, beta) {
  eta <- drop(X %*% beta)
  eta <- pmin(pmax(eta, -50), 50)
  af_i <- 1 - exp(-eta)
  an <- sum(cases * af_i, na.rm = TRUE)
  af <- an / sum(cases, na.rm = TRUE)
  c(af = af, an = an)
}

attr_one <- function(x_lag, cases, cb_coef, cb_vcov, cen, range, component) {
  x_use <- x_lag
  if (!is.null(range)) x_use[x_use < range[1] | x_use > range[2]] <- cen
  x_cen <- matrix(cen, nrow = nrow(x_use), ncol = ncol(x_use))

  b_use <- crossbasis(x_use, lag = c(0, MAIN_LAG), argvar = argvar_main, arglag = arglag_main)
  b_cen <- crossbasis(x_cen, lag = c(0, MAIN_LAG), argvar = argvar_main, arglag = arglag_main)
  X <- as.matrix(b_use) - as.matrix(b_cen)

  if (ncol(X) != length(cb_coef)) stop("Cross-basis coefficient length mismatch.")

  point <- calc_af_an_from_basis(X, cases, cb_coef)
  beta_sim <- rmvnorm_base(NSIM, cb_coef, cb_vcov)
  sim_mat <- t(apply(beta_sim, 1, function(beta) calc_af_an_from_basis(X, cases, beta)))
  if (is.null(dim(sim_mat))) sim_mat <- matrix(sim_mat, nrow = NSIM, byrow = TRUE)
  colnames(sim_mat) <- c("af", "an")

  af_ci <- quantile(sim_mat[, "af"], CI_PROBS, na.rm = TRUE, names = FALSE)
  an_ci <- quantile(sim_mat[, "an"], CI_PROBS, na.rm = TRUE, names = FALSE)

  tibble(
    component = component,
    af = as.numeric(point["af"]),
    af_low = af_ci[1],
    af_high = af_ci[2],
    af_percent = 100 * as.numeric(point["af"]),
    af_percent_low = 100 * af_ci[1],
    af_percent_high = 100 * af_ci[2],
    an = as.numeric(point["an"]),
    an_low = an_ci[1],
    an_high = an_ci[2]
  )
}

### ---------- 4. Cold, heat, and total burden ----------

mmt <- pkg$mmt
temp_min <- min(dat[[TEMP_VAR]][fit_rows], na.rm = TRUE)
temp_max <- max(dat[[TEMP_VAR]][fit_rows], na.rm = TRUE)

result <- bind_rows(
  attr_one(x_lag_fit, cases_fit, pkg$cb_coef, pkg$cb_vcov, mmt, c(temp_min, mmt), "cold"),
  attr_one(x_lag_fit, cases_fit, pkg$cb_coef, pkg$cb_vcov, mmt, c(mmt, temp_max), "heat"),
  attr_one(x_lag_fit, cases_fit, pkg$cb_coef, pkg$cb_vcov, mmt, NULL, "total")
) %>%
  mutate(
    total_deaths_fitted = sum(cases_fit),
    mmt = mmt,
    n_simulations = NSIM
  )

write_csv(result, file.path(OUT_DIR, "example_attributable_burden.csv"))
message("Finished. Attributable-burden output written to ", OUT_DIR, "/example_attributable_burden.csv")
