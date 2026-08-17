# 01_temperature_mortality_model.R
# Example one-stage space-time-stratified case-crossover DLNM.
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
  library(gnm)
})

### ---------- 1. Paths and settings ----------

DATA_FILE <- file.path("data", "example_daily_data.csv")
OUT_DIR <- "output"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(DATA_FILE)) {
  stop("Run this script from the repository root. Missing: ", DATA_FILE)
}

TEMP_VAR <- "temperature"
OUTCOME <- "deaths"
MAIN_LAG <- 21

### ---------- 2. Read and prepare synthetic data ----------

dat <- read_csv(DATA_FILE, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(date),
    pref = factor(pref),
    holiday = factor(holiday),
    year = year(date),
    month = month(date),
    dow = factor(wday(date, label = TRUE, week_start = 1)),
    stratum = interaction(pref, year, month, dow, drop = TRUE)
  ) %>%
  arrange(pref, date)

make_analysis_data <- function(data, outcome) {
  data %>%
    filter(!is.na(.data[[TEMP_VAR]]), !is.na(.data[[outcome]])) %>%
    mutate(
      y = as.numeric(.data[[outcome]]),
      stratum = droplevels(factor(stratum)),
      holiday = droplevels(factor(holiday)),
      pref = droplevels(factor(pref))
    ) %>%
    group_by(stratum) %>%
    mutate(ind = sum(y, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(pref, date)
}

make_lag_matrix_by_group <- function(x, group, lag) {
  out <- matrix(NA_real_, nrow = length(x), ncol = lag + 1L)
  colnames(out) <- paste0("lag", 0:lag)
  
  for (idx in split(seq_along(x), group, drop = TRUE)) {
    x_i <- x[idx]
    for (lag_i in 0:lag) out[idx, lag_i + 1L] <- dplyr::lag(x_i, n = lag_i)
  }
  
  out
}

extract_rr_at_temp <- function(pred_obj, target_temp) {
  idx <- which.min(abs(pred_obj$predvar - target_temp))
  
  tibble(
    temperature = as.numeric(pred_obj$predvar[idx]),
    rr = as.numeric(pred_obj$allRRfit[idx]),
    low = as.numeric(pred_obj$allRRlow[idx]),
    high = as.numeric(pred_obj$allRRhigh[idx])
  )
}

get_cb_coef_vcov <- function(model, basis) {
  beta_all <- coef(model)
  vcov_all <- vcov(model)
  
  cb_terms <- grep("^cb", names(beta_all), value = TRUE)
  
  if (length(cb_terms) != ncol(basis)) {
    cb_terms <- intersect(names(beta_all), colnames(basis))
  }
  
  if (length(cb_terms) != ncol(basis)) {
    stop("Could not identify all cross-basis coefficients in the fitted model.")
  }
  
  list(
    coef = beta_all[cb_terms],
    vcov = vcov_all[cb_terms, cb_terms, drop = FALSE]
  )
}

### ---------- 3. Cross-basis specification ----------

model_dat <- make_analysis_data(dat, OUTCOME)

temp_knots <- quantile(model_dat[[TEMP_VAR]], probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
temp_bounds <- range(model_dat[[TEMP_VAR]], na.rm = TRUE)
lag_knots <- logknots(MAIN_LAG, nk = 3)

argvar_main <- list(fun = "ns", knots = temp_knots, Boundary.knots = temp_bounds)
arglag_main <- list(fun = "ns", knots = lag_knots, intercept = TRUE)

cb <- crossbasis(
  x = model_dat[[TEMP_VAR]],
  lag = MAIN_LAG,
  group = model_dat$pref,
  argvar = argvar_main,
  arglag = arglag_main
)

x_lag_all <- make_lag_matrix_by_group(
  x = model_dat[[TEMP_VAR]],
  group = model_dat$pref,
  lag = MAIN_LAG
)

fit_rows <- model_dat$ind > 0 &
  complete.cases(cb) &
  complete.cases(x_lag_all) &
  !is.na(model_dat$y)

if (!any(fit_rows)) stop("No usable rows after lag construction.")

### ---------- 4. Fit the one-stage model ----------

fit <- gnm(
  y ~ cb + holiday,
  eliminate = stratum,
  family = quasipoisson(),
  data = model_dat,
  subset = fit_rows
)

### ---------- 5. MMT and cumulative RR curve ----------

temp_fit <- model_dat[[TEMP_VAR]][fit_rows]

temp_p1 <- as.numeric(quantile(temp_fit, 0.01, na.rm = TRUE))
temp_p50 <- as.numeric(quantile(temp_fit, 0.50, na.rm = TRUE))
temp_p99 <- as.numeric(quantile(temp_fit, 0.99, na.rm = TRUE))

temp_grid_search <- sort(unique(c(
  seq(temp_p50, temp_p99, by = 0.1),
  temp_p50,
  temp_p99
)))

pred_search <- crosspred(
  cb,
  model = fit,
  at = temp_grid_search,
  cen = temp_p50
)

mmt <- temp_grid_search[which.min(pred_search$allRRfit)]
mmt_percentile <- mean(temp_fit <= mmt, na.rm = TRUE) * 100

temp_grid_plot <- seq(min(temp_fit), max(temp_fit), length.out = 250)

pred_final <- crosspred(
  cb,
  model = fit,
  at = temp_grid_plot,
  cen = mmt
)

rr_p1 <- extract_rr_at_temp(pred_final, temp_p1)
rr_p99 <- extract_rr_at_temp(pred_final, temp_p99)
cb_par <- get_cb_coef_vcov(fit, cb)

rr_curve <- tibble(
  temperature = as.numeric(pred_final$predvar),
  rr = as.numeric(pred_final$allRRfit),
  low = as.numeric(pred_final$allRRlow),
  high = as.numeric(pred_final$allRRhigh)
)

rr_summary <- tibble(
  outcome = OUTCOME,
  n_rows_fitted = sum(fit_rows),
  deaths_fitted = sum(model_dat$y[fit_rows]),
  mmt = mmt,
  mmt_percentile = mmt_percentile,
  p1_temperature = temp_p1,
  p1_rr = rr_p1$rr,
  p1_low = rr_p1$low,
  p1_high = rr_p1$high,
  p99_temperature = temp_p99,
  p99_rr = rr_p99$rr,
  p99_low = rr_p99$low,
  p99_high = rr_p99$high
)

write_csv(rr_summary, file.path(OUT_DIR, "example_rr_summary.csv"))

### ---------- 6. Plot example cumulative RR curve ----------

png(
  filename = file.path(OUT_DIR, "example_rr_curve.png"),
  width = 1600,
  height = 1200,
  res = 180
)

par(
  mar = c(5.2, 5.2, 1.5, 1.5),
  mgp = c(3.0, 0.8, 0),
  las = 1,
  bty = "l"
)

plot(
  rr_curve$temperature,
  rr_curve$rr,
  type = "n",
  xlab = "Mean temperature (°C)",
  ylab = "Cumulative relative risk",
  ylim = range(c(rr_curve$low, rr_curve$high), finite = TRUE),
  cex.lab = 1.25,
  cex.axis = 1.05
)

polygon(
  c(rr_curve$temperature, rev(rr_curve$temperature)),
  c(rr_curve$low, rev(rr_curve$high)),
  border = NA,
  col = adjustcolor("grey70", alpha.f = 0.45)
)

lines(
  rr_curve$temperature,
  rr_curve$rr,
  lwd = 2
)

abline(h = 1, lty = 2, col = "grey40")
abline(v = mmt, lty = 3, col = "grey30")

points(mmt, 1, pch = 19)

legend(
  "topright",
  legend = c(
    "Cumulative RR",
    "95% CI",
    paste0("MMT = ", round(mmt, 1), " °C")
  ),
  lty = c(1, NA, 3),
  lwd = c(2, NA, 1),
  pch = c(NA, 15, NA),
  pt.cex = c(NA, 1.5, NA),
  col = c("black", "grey70", "grey30"),
  bty = "n",
  cex = 0.95
)

dev.off()

### ---------- 7. Save model package for attributable-risk example ----------

model_package <- list(
  temp_var = TEMP_VAR,
  outcome = OUTCOME,
  main_lag = MAIN_LAG,
  temp_knots = temp_knots,
  temp_bounds = temp_bounds,
  lag_knots = lag_knots,
  argvar_main = argvar_main,
  arglag_main = arglag_main,
  mmt = mmt,
  mmt_percentile = mmt_percentile,
  temp_p1 = temp_p1,
  temp_p99 = temp_p99,
  cb_coef = cb_par$coef,
  cb_vcov = cb_par$vcov,
  n_fit_rows = sum(fit_rows),
  deaths_fitted = sum(model_dat$y[fit_rows])
)

saveRDS(
  model_package,
  file.path(OUT_DIR, "example_model_package.rds")
)

message(
  "Finished. MMT = ", round(mmt, 2),
  " C; RR figure, summary, and model package written to ",
  OUT_DIR, "/"
)