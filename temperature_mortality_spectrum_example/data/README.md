# Synthetic example data

Both files in this folder were generated specifically for this repository and contain no records from the Japanese Vital Statistics data used in the study.

## `example_daily_data.csv`

Columns:

- `date`: calendar date
- `pref`: synthetic prefecture identifier
- `temperature`: synthetic daily mean temperature in degrees Celsius
- `deaths`: synthetic daily all-cause death count
- `holiday`: synthetic/example national-holiday indicator (0/1)

The file is used by `01_temperature_mortality_model.R` and `02_attributable_mortality.R`.

## `example_decomposition_data.csv`

Each row is a synthetic prefecture × age group × cause × burden-type cell. Columns contain reference-period (`0`) and target-period (`1`) values for:

- `P0`, `P1`: mean annual total population
- `S0`, `S1`: age-group population shares
- `M0`, `M1`: annual cause-specific mortality rates
- `AF_t0_er0`, `AF_t1_er0`, `AF_t0_er1`, `AF_t1_er1`: attributable fractions under the four combinations of reference/target temperature–death distribution (`t`) and reference/target exposure–response function (`er`)

The file is used only to demonstrate the decomposition mechanics. Values are synthetic and are not estimates from the study.
