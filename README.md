# Heat, urban green space and wearable-derived stress

R code for the analyses in *Urban green-space attributes modify associations between heat and wearable-derived stress*. The four modules use city-level aggregate data. Individual-level stress processing, green-space indicator construction and other individual-level analyses are not included. The supplied `example_data.rds` files contain synthetic data and do not reproduce the study estimates.

## Modules

| Script | Analysis | Main output |
|---|---|---|
| `01_DLNMs/DLNM.R` | City-specific quasi-binomial distributed lag non-linear models (DLNMs) | `first_stage.rds`, reduced coefficients and covariance matrices |
| `02_Meta_Regression/Meta_Regression.R` | Multivariate meta-regression and green-space effect modification | `meta_analysis.rds`, odds ratios, RORs and tests |
| `03_Subgroup_Analyses/Subgroup_Analyses.R` | Sex, age and bioclimatic-region analyses | Subgroup and regional estimates |
| `04_Future_Projections/Future_Projections.R` | Historical and projected heat-related burden under green-space adaptation scenarios | `projections.rds`, burden estimates and city classifications |

Shared settings and functions are defined in `DLNM.R` and loaded by subsequent scripts. Each example includes the fitted objects needed from earlier stages, so the modules can be run separately.

## Requirements

R 4.4 or later, with `dlnm`, `mixmeta` and `MASS`. The code was prepared with R 4.6.0, dlnm 2.4.10, mixmeta 1.2.2 and MASS 7.3-65.

```r
install.packages(c("dlnm", "mixmeta", "MASS"), repos = "https://cloud.r-project.org")
```

## Input data

Each module reads an RDS file containing a named list. The list includes `module` (integer 1–4), `synthetic` (logical), `provenance` and the objects below. Dates use R's `Date` class, identifiers and categories are character values, and measurements and counts are numeric. Missing values are recorded as `NA`.

| Module | Required objects |
|---|---|
| 1 | `city_day_outcomes`, `historical_city_temperature` |
| 2 | `first_stage`, `city_attributes` |
| 3 | `first_stage`, `city_attributes`, `city_contextual_covariates`, `subgroup_city_day_outcomes`, `historical_city_temperature` |
| 4 | `gvi_model`, `city_attributes`, `city_day_outcomes`, `historical_city_temperature`, `future_city_temperature`, `city_month_vegetation`, `city_month_cooling_efficiency`, `city_month_baseline_probabilities`, `temperature_support` |

### Data fields

- `city_day_outcomes`: `city_id,date,high_stress_count,non_high_stress_count,temperature,relative_humidity,wind_speed,precipitation,aqi,holiday`. One record per city and date. Units are °C, %, m/s and mm/day for temperature, relative humidity, wind speed and precipitation, respectively; `holiday` is 0 or 1.
- `subgroup_city_day_outcomes`: the same fields, with `stratum_type,stratum`. Sex groups are `Female` and `Male`; age groups are `18-34`, `35-49`, `50-64` and `65plus`.
- `city_attributes`: `city_id,region,cooling,gvi,accessibility,population_2020`. Cooling is signed air-temperature reduction capacity (°C); GVI is a fraction from 0 to 1; accessibility is the 15-minute accessibility index; population is the 2020 census resident count.
- `historical_city_temperature`: `city_id,date,temperature`. Complete April–September daily temperature records for the analytical and historical reference years.
- `future_city_temperature`: `city_id,date,temperature,ssp,climate_model`. Daily April–September temperatures for 2030–2100 under `ssp126`, `ssp245` and `ssp585`. Each `climate_model` identifies one GCM with one daily temperature series, and the same GCM set is used across the corresponding scenarios.
- `city_month_vegetation`: `city_id,month,ndvi`; monthly NDVI (−1 to 1) for April–September.
- `city_month_cooling_efficiency`: `city_id,month,cooling_efficiency`; signed cooling efficiency in °C per unit NDVI, with positive values indicating cooling.
- `city_month_baseline_probabilities`: `city_id,month,p0,reference_days,reference_user_days,high_stress_count,first_reference_date,last_reference_date`; city-month baseline probabilities and the counts used to calculate them.

`city_contextual_covariates` includes `city_id` and seven variables: `pct_age_0_14`, `pct_age_65_plus`, `pct_urban_builtup`, `education_years`, `pct_one_person_households`, `pct_nonlocal_registered` and `pct_agricultural_employment`. The extended adjustment adds `gdp_100m_cny` (CNY 100 million), `nighttime_light` (nW cm−2 sr−1) and `ac_per_100_households` (units per 100 households). Proportion variables are expressed as percentages.

Bioclimatic regions use the codes `humid_subtropical_tropical_forest`, `humid_temperate_forest`, `semi_arid_temperate_grassland`, `arid_temperate_desert` and `qinghai_tibetan_alpine`.

`first_stage` contains city-level `coefficients` (`city_id,b1,b2,b3`), corresponding full 3×3 `covariance` matrices, `climate` (`mean_temperature,temperature_range`), `basis`, `ledger` and `temperature_support`. `gvi_model` is the `continuous$gvi` result from Module 2. `temperature_support` records each city's fitted temperature range and observation dates (`city_id,historical_min,historical_max,model_days,first_model_date,last_model_date`).

Before aggregation, high-stress days are defined as days with stress strictly above each user's P90, calculated from valid observations in the study period. This requires at least five measurements per day and ten valid days per user. The same thresholds are used for subgroup data.

## Run the examples

From the repository root, run each script with its example input. Output directories must be empty and located outside the repository.

```sh
Rscript 01_DLNMs/DLNM.R --input 01_DLNMs/example_data.rds --output ../results_dlnm
Rscript 02_Meta_Regression/Meta_Regression.R --input 02_Meta_Regression/example_data.rds --output ../results_meta
Rscript 03_Subgroup_Analyses/Subgroup_Analyses.R --input 03_Subgroup_Analyses/example_data.rds --output ../results_subgroups
Rscript 04_Future_Projections/Future_Projections.R --input 04_Future_Projections/example_data.rds --output ../results_future
```

The examples cover 74 fictional cities. An input file can be inspected with:

```r
str(readRDS("01_DLNMs/example_data.rds"))
```

For study-data runs, set `synthetic = FALSE` and create input bundles with the same fields and fitted objects from the corresponding preceding analyses. For example, the `first_stage` object for Module 2 can be loaded from `../results_dlnm/first_stage.rds`. Save the input bundle as an RDS file and supply its path with `--input`.

## Analytical methods

The DLNM module uses city-specific models based on April–September observations with at least 10 users per city-day and 50 eligible city-days per city. Quasi-binomial models with a logit link relate high-stress counts to temperature using a natural spline with national P33 and P66 knots and a lag 0–5 natural spline with 3 degrees of freedom and an intercept. Models adjust for relative humidity and wind speed (natural splines, 3 df each), `log1p(precipitation)`, AQI, day of week, holidays, year-specific seasonality (`ns(doy,3):factor(year)`) and long-term trend (`ns(time,1)`). Lag histories are reset between city-years.

National P33, P50, P66, P90, P95 and P99 temperatures are calculated from complete April–September daily meteorological records for the included cities and analytical years, with equal weight per city-day. P50 is the reference temperature; P90, P95 and P99 are comparison temperatures. The temperature spline has a lower boundary at the minimum of those records and an upper boundary equal to the mean of the city-specific maximum temperatures among eligible first-stage records before lag-related exclusions. The same temperature basis is used in subsequent modules. City mean temperatures and temperature ranges are calculated from the complete meteorological records.

The meta-regression module pools reduced city coefficients and their full covariance matrices using `mixmeta`, with city mean temperature and temperature range as covariates. Random-effects covariance is unstructured. Linear and 3-df spline forms for continuous green-space modifiers are compared by ML likelihood-ratio tests; selected models are fitted by REML. Analyses include green-space tertiles, ratios of odds ratios (RORs) and attenuation estimates.

The subgroup module uses the common temperature basis and national green-space tertiles for sex and age analyses. Bioclimatic-region analyses use region-specific tertiles and a joint region-by-green-space interaction model. Contextual adjustment is examined using the seven- and ten-variable sets described above.

For future projections, city-month baseline probabilities are calculated from valid warm-season observations at temperatures at or below the common P50. The fitted odds ratios are converted to probabilities using `OR*p0/(1-p0+OR*p0)`, and only positive excess probability above baseline at temperatures above P50 contributes to the heat-related burden. Historical and future burdens are expressed as mean annual April–September high-stress user-days per 100,000 people. Adaptation scenarios increase GVI, enhance cooling or combine both changes.

For each SSP and adaptation scenario, burdens are calculated separately from each GCM's daily temperatures, averaged across reporting years and then equally across GCMs within each city. National estimates use fixed 2020 city populations. Parameter uncertainty is estimated from 1,000 joint draws of the fitted coefficients and covariance matrix, using the 2.5th and 97.5th percentiles. Outputs include burdens, adaptation offsets and gaps, and city classifications.

## Data availability

Huawei individual-level data are subject to privacy and data-use restrictions and are not included in this repository. Reproducing the study estimates requires access to the corresponding authorized data and other study inputs.
