# 03_shapley_decomposition.R
# Exact five-factor Shapley decomposition of attributable-death magnitude and cause shares.
# Uses synthetic factor data only; outputs are not study estimates.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

### ---------- 1. Paths and factor definitions ----------

DATA_FILE <- file.path("data", "example_decomposition_data.csv")
OUT_DIR <- "output"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(DATA_FILE)) stop("Run this script from the repository root. Missing: ", DATA_FILE)

FACTOR_ORDER <- c("population", "aging", "mortality_rate", "temperature", "er")
FACTOR_LABELS <- c(
  population = "Population size",
  aging = "Age structure",
  mortality_rate = "Mortality rate",
  temperature = "Temperature-death distribution",
  er = "Exposure-response"
)
SHAPLEY_TOL <- 1e-8

make_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], make_permutations(x[-i]))))
}

STATE_GRID <- expand.grid(
  setNames(rep(list(0:1), length(FACTOR_ORDER)), FACTOR_ORDER),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(across(all_of(FACTOR_ORDER), as.integer))
STATE_GRID$key <- apply(STATE_GRID[FACTOR_ORDER], 1L, paste0, collapse = "")
STATE_KEYS <- STATE_GRID$key
SHAPLEY_PERMUTATIONS <- make_permutations(FACTOR_ORDER)

exact_shapley <- function(value_matrix) {
  value_matrix <- as.matrix(value_matrix)
  if (is.null(rownames(value_matrix))) stop("value_matrix must have 5-bit state keys as row names.")
  missing_states <- setdiff(STATE_KEYS, rownames(value_matrix))
  if (length(missing_states) > 0L) stop("Missing states: ", paste(missing_states, collapse = ", "))
  value_matrix <- value_matrix[STATE_KEYS, , drop = FALSE]

  contribution <- matrix(
    0, nrow = length(FACTOR_ORDER), ncol = ncol(value_matrix),
    dimnames = list(FACTOR_ORDER, colnames(value_matrix))
  )

  for (perm_i in seq_len(nrow(SHAPLEY_PERMUTATIONS))) {
    state <- setNames(rep(0L, length(FACTOR_ORDER)), FACTOR_ORDER)
    old_value <- value_matrix[paste0(state, collapse = ""), , drop = TRUE]
    for (factor_i in SHAPLEY_PERMUTATIONS[perm_i, ]) {
      state[factor_i] <- 1L
      new_value <- value_matrix[paste0(state, collapse = ""), , drop = TRUE]
      contribution[factor_i, ] <- contribution[factor_i, ] + new_value - old_value
      old_value <- new_value
    }
  }
  contribution / nrow(SHAPLEY_PERMUTATIONS)
}

### ---------- 2. Read and validate synthetic factor data ----------

dat <- read_csv(DATA_FILE, show_col_types = FALSE)
required <- c(
  "pref", "age_group", "cause", "burden_type",
  "P0", "P1", "S0", "S1", "M0", "M1",
  "AF_t0_er0", "AF_t1_er0", "AF_t0_er1", "AF_t1_er1"
)
missing_columns <- setdiff(required, names(dat))
if (length(missing_columns) > 0L) stop("Missing columns: ", paste(missing_columns, collapse = ", "))

select_af <- function(data, temperature_state, er_state) {
  if (temperature_state == 0L && er_state == 0L) return(data$AF_t0_er0)
  if (temperature_state == 1L && er_state == 0L) return(data$AF_t1_er0)
  if (temperature_state == 0L && er_state == 1L) return(data$AF_t0_er1)
  data$AF_t1_er1
}

### ---------- 3. Construct all 32 counterfactual states ----------

coalition_values <- lapply(seq_len(nrow(STATE_GRID)), function(i) {
  state <- STATE_GRID[i, ]
  p <- if (state$population == 1L) dat$P1 else dat$P0
  s <- if (state$aging == 1L) dat$S1 else dat$S0
  m <- if (state$mortality_rate == 1L) dat$M1 else dat$M0
  af <- select_af(dat, state$temperature, state$er)

  dat %>%
    transmute(
      cause,
      burden_type,
      an_per_year = p * s * m * af
    ) %>%
    group_by(burden_type, cause) %>%
    summarise(an_per_year = sum(an_per_year), .groups = "drop") %>%
    mutate(
      population = state$population,
      aging = state$aging,
      mortality_rate = state$mortality_rate,
      temperature = state$temperature,
      er = state$er,
      key = state$key
    )
}) %>% bind_rows()

### ---------- 4. Decompose change in attributable-death magnitude ----------

magnitude_results <- lapply(unique(coalition_values$burden_type), function(burden_i) {
  state_total <- coalition_values %>%
    filter(burden_type == burden_i) %>%
    group_by(key) %>%
    summarise(value = sum(an_per_year), .groups = "drop")

  value_matrix <- matrix(state_total$value, ncol = 1L, dimnames = list(state_total$key, "AN_per_year"))
  shap <- exact_shapley(value_matrix)
  delta <- value_matrix["11111", 1] - value_matrix["00000", 1]

  tibble(
    burden_type = burden_i,
    factor = rownames(shap),
    factor_label = unname(FACTOR_LABELS[rownames(shap)]),
    contribution_deaths_per_year = as.numeric(shap[, 1]),
    reference_deaths_per_year = as.numeric(value_matrix["00000", 1]),
    target_deaths_per_year = as.numeric(value_matrix["11111", 1]),
    net_change_deaths_per_year = as.numeric(delta),
    reconciliation_error = sum(shap[, 1]) - delta
  )
}) %>% bind_rows()

if (any(abs(magnitude_results$reconciliation_error) > SHAPLEY_TOL)) {
  stop("Magnitude Shapley decomposition failed the reconciliation check.")
}

write_csv(magnitude_results, file.path(OUT_DIR, "example_shapley_magnitude.csv"))

### ---------- 5. Decompose change in cause composition ----------

composition_results <- lapply(unique(coalition_values$burden_type), function(burden_i) {
  wide <- coalition_values %>%
    filter(burden_type == burden_i) %>%
    mutate(an_positive = pmax(an_per_year, 0)) %>%
    group_by(key) %>%
    mutate(share_percent = 100 * an_positive / sum(an_positive)) %>%
    ungroup() %>%
    select(key, cause, share_percent) %>%
    pivot_wider(names_from = cause, values_from = share_percent)

  share_matrix <- as.matrix(wide[, -1, drop = FALSE])
  rownames(share_matrix) <- wide$key
  shap <- exact_shapley(share_matrix)
  causes <- colnames(share_matrix)

  bind_rows(lapply(causes, function(cause_i) {
    delta <- share_matrix["11111", cause_i] - share_matrix["00000", cause_i]
    tibble(
      burden_type = burden_i,
      cause = cause_i,
      factor = rownames(shap),
      factor_label = unname(FACTOR_LABELS[rownames(shap)]),
      contribution_percentage_points = as.numeric(shap[, cause_i]),
      reference_share_percent = as.numeric(share_matrix["00000", cause_i]),
      target_share_percent = as.numeric(share_matrix["11111", cause_i]),
      net_change_percentage_points = as.numeric(delta),
      reconciliation_error = sum(shap[, cause_i]) - delta
    )
  }))
}) %>% bind_rows()

if (any(abs(composition_results$reconciliation_error) > SHAPLEY_TOL)) {
  stop("Cause-share Shapley decomposition failed the reconciliation check.")
}

write_csv(composition_results, file.path(OUT_DIR, "example_shapley_composition.csv"))

message("Finished. Exact five-factor Shapley outputs written to ", OUT_DIR, "/")
