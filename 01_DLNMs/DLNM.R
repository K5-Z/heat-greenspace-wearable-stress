#!/usr/bin/env Rscript

.dlnm_file <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  if (length(paths)) tail(paths, 1)[[1]] else {
    arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (!length(arg)) stop("Load this file with source() or Rscript")
    sub("^--file=", "", arg[1])
  }
})
.dlnm_root <- dirname(dirname(normalizePath(.dlnm_file, winslash = "/", mustWork = TRUE)))

# Common analysis definitions for city models, pooling and projection.
analysis_settings <- function() {
  list(
    study = list(future_years = c(2030L, 2100L), warm_months = 4:9,
      min_users = 10L, min_city_days = 50L,
      temperature_knot_probs = c(.33, .66), lag_max = 5L, lag_df = 3L,
      alpha = .05,
      ssps = c("ssp126", "ssp245", "ssp585")),
    choices = list(quantile_type = 7L, tertile_probs = c(1/3, 2/3),
      decade_starts = c(2030L, 2050L, 2070L, 2090L), decade_length = 10L),
    bootstrap = list(repetitions = 1000L, seed = 20261027L, batch_size = 20L),
    numerics = list(condition_warning = 1e12, glm_maxit = 100L,
      meta_maxiter = 1000L, meta_tolerance = 1e-7, small_cell_warning = 5L))
}

assert_columns <- function(data, fields, label = deparse(substitute(data))) {
  absent <- setdiff(fields, names(data))
  if (length(absent)) stop(label, " missing fields: ", paste(absent, collapse = ", "))
  invisible(TRUE)
}

assert_unique <- function(data, keys) {
  assert_columns(data, keys)
  if (anyNA(data[keys]) || anyDuplicated(data[keys]))
    stop("Missing or duplicate key: ", paste(keys, collapse = "+"))
  invisible(TRUE)
}

assert_finite <- function(x, label) {
  if (!is.numeric(x) || any(!is.finite(x))) stop(label, " must contain finite numeric values")
  invisible(TRUE)
}

validate_covariance <- function(v, n = nrow(v), positive_definite = FALSE, tolerance = 1e-10) {
  if (!is.matrix(v) || !identical(dim(v), as.integer(c(n, n))) || any(!is.finite(v)))
    stop("Invalid covariance dimension or non-finite entries")
  scale <- max(abs(v), .Machine$double.eps)
  if (max(abs(v - t(v))) > tolerance * scale) stop("Covariance is not symmetric")
  ev <- eigen(v, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -tolerance * scale) stop("Covariance is not positive semidefinite")
  if (positive_definite && min(ev) <= tolerance * scale) stop("Covariance is singular or numerically unstable")
  list(min_eigenvalue = min(ev), max_eigenvalue = max(ev),
       condition = if (min(ev) <= 0) Inf else max(ev) / min(ev))
}

validate_outcomes <- function(d, subgroup = FALSE) {
  fields <- c("city_id", "date", "high_stress_count", "non_high_stress_count",
              "temperature", "relative_humidity", "wind_speed", "precipitation", "aqi", "holiday")
  assert_columns(d, fields)
  if (subgroup) assert_columns(d, c("stratum_type", "stratum"))
  d$date <- as.Date(d$date)
  assert_unique(d, c("city_id", "date", if (subgroup) c("stratum_type", "stratum")))
  for (nm in c("high_stress_count", "non_high_stress_count")) {
    x <- d[[nm]]
    if (!is.numeric(x) || any(!is.na(x) & (!is.finite(x) | x < 0 | x != floor(x))))
      stop(nm, " must be nonnegative integer counts or NA")
  }
  for (nm in c("temperature", "relative_humidity", "wind_speed", "precipitation", "aqi", "holiday"))
    if (!is.numeric(d[[nm]]) || any(is.infinite(d[[nm]]))) stop("Invalid numeric field: ", nm)
  if (any(!is.na(d$precipitation) & d$precipitation < 0) ||
      any(!is.na(d$wind_speed) & d$wind_speed < 0) ||
      any(!is.na(d$relative_humidity) & (d$relative_humidity < 0 | d$relative_humidity > 100)) ||
      any(!is.na(d$holiday) & !d$holiday %in% c(0, 1))) stop("Out-of-domain environmental/calendar value")
  d
}

region_levels <- function(include_alpine = FALSE) {
  c("humid_subtropical_tropical_forest", "humid_temperate_forest",
    "semi_arid_temperate_grassland", "arid_temperate_desert",
    if (include_alpine) "qinghai_tibetan_alpine")
}

validate_attributes <- function(a) {
  assert_columns(a, c("city_id", "region", "cooling", "gvi", "accessibility", "population_2020"))
  assert_unique(a, "city_id")
  for (nm in c("cooling", "gvi", "accessibility", "population_2020")) assert_finite(a[[nm]], nm)
  if (any(a$gvi < 0 | a$gvi > 1)) stop("GVI must be a fraction in [0,1], not a percentage")
  if (any(a$population_2020 <= 0) || any(a$accessibility < 0)) stop("Invalid attribute domain")
  if (any(!a$region %in% region_levels(TRUE))) stop("Unrecognized bioclimatic region")
  a
}

capture_warnings <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
  })
  list(value = value, warnings = unique(messages))
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  invisible(path)
}

save_model <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path, version = 3)
  invisible(path)
}

check_temperature_basis <- function(spec, cfg) {
  k <- spec$argvar$knots; b <- spec$argvar$Boundary.knots
  q <- spec$quantiles
  eq <- function(x,y) isTRUE(all.equal(unname(x),unname(y),tolerance=1e-12))
  valid <- identical(spec$argvar$fun,"ns") && identical(spec$argvar$intercept,FALSE) &&
    length(k)==2L && length(b)==2L && all(is.finite(c(k,b))) &&
    b[1]<k[1] && k[1]<k[2] && k[2]<b[2] &&
    identical(names(q),c("P33","P50","P66","P90","P95","P99")) &&
    all(is.finite(q)) && all(diff(q)>=0) && eq(k,q[c("P33","P66")]) &&
    eq(spec$reference,q["P50"]) && eq(spec$contrasts,q[c("P90","P95","P99")]) &&
    identical(spec$arglag$fun,"ns") && isTRUE(spec$arglag$intercept) &&
    is.null(spec$arglag$knots) && isTRUE(spec$arglag$df==cfg$study$lag_df) &&
    isTRUE(spec$lag==cfg$study$lag_max)
  if (!valid) stop("Inconsistent shared temperature/lag basis; refit all stages together")
  lb <- do.call(dlnm::onebasis,c(list(x=0:spec$lag),spec$arglag))
  if (ncol(lb)!=3L || qr(lb)$rank!=3L) stop("Lag basis must have three degrees of freedom")
  src <- spec$source
  if ((length(src$years)==0L || anyNA(src$years) || anyDuplicated(src$years)) || !identical(src$months,cfg$study$warm_months) ||
      !eq(b[1],src$national_minimum) || !eq(b[2],mean(src$city_analysis_maxima$maximum)))
    stop("Missing or inconsistent temperature-source metadata")
  invisible(TRUE)
}

read_module_input <- function(path, module, fields) {
  if (!file.exists(path)) stop("Missing input RDS (no example fallback): ", path)
  x <- readRDS(path)
  if (!is.list(x) || !identical(x$module, as.integer(module)) ||
      !is.logical(x$synthetic) || length(x$synthetic) != 1L || is.na(x$synthetic))
    stop("Input requires module (integer) and synthetic (TRUE/FALSE) metadata")
  assert_columns(x, fields, "RDS bundle")
  forbidden <- c("user_id", "userid", "device_id", "deviceid", "latitude", "longitude", "gps",
                 "avgstress", "maxstress", "minstress", "laststress")
  for (d in Filter(is.data.frame, x)) {
    if (!nrow(d)) stop("Empty input table")
    if (any(tolower(names(d)) %in% forbidden)) stop("Individual records are outside this aggregate interface")
    if ("city_id" %in% names(d) && !is.character(d$city_id)) stop("city_id must be character")
  }
  message(if (x$synthetic) "Synthetic example input." else
    "User-supplied aggregate input.")
  x
}

module_arguments <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 4L || !identical(args[c(1, 3)], c("--input", "--output")))
    stop("Usage: Rscript <module.R> --input <bundle.rds> --output <new external directory>")
  list(input = args[2], output = args[4])
}

prepare_output <- function(path, root, input) {
  # Keep all generated data outside the nine-file public repository.
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  dest <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (tolower(dest) == tolower(root) || startsWith(tolower(dest), paste0(tolower(root), "/")))
    stop("Output directory must be outside the public repository")
  if (length(list.files(dest, all.files = TRUE, no.. = TRUE))) stop("Output directory must be empty")
  dest
}

finish_module <- function(result, input, path) {
  result$synthetic <- input$synthetic
  result$module <- input$module
  result$settings <- analysis_settings()
  result$software <- list(R = R.version.string,
    packages = vapply(c("dlnm", "mixmeta", "MASS"), function(p)
      if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1)))
  save_model(result, path)
  message("Completed module ", input$module, "; results: ", dirname(path))
  invisible(result)
}

validate_parameter_covariance <- function(beta, V) {
  validate_covariance(V, length(beta), positive_definite = TRUE)
  if (is.null(names(beta)) || anyDuplicated(names(beta)) ||
      !identical(rownames(V), names(beta)) || !identical(colnames(V), names(beta)))
    stop("Coefficient/covariance parameter names and order must match exactly")
  invisible(TRUE)
}

# Complete meteorological calendars are independent of outcome availability.
primary_temperature_calendar <- function(historical, city_ids, cfg) {
  assert_columns(historical,c("city_id","date","temperature"))
  historical$date <- as.Date(historical$date)
  assert_unique(historical,c("city_id","date"))
  year <- as.integer(format(historical$date,"%Y")); month <- as.integer(format(historical$date,"%m"))
  h <- historical[historical$city_id %in% city_ids & year %in% cfg$study$primary_years &
    month %in% cfg$study$warm_months,c("city_id","date","temperature")]
  assert_finite(h$temperature,"Complete primary warm-season meteorology")
  if (!setequal(h$city_id,city_ids)) stop("Missing meteorology for an included city")
  for (id in city_ids) for (yr in cfg$study$primary_years) {
    dates <- h$date[h$city_id==id & format(h$date,"%Y")==as.character(yr)]
    expected <- seq(as.Date(sprintf("%d-04-01",yr)),as.Date(sprintf("%d-09-30",yr)),by="day")
    if (!setequal(dates,expected)) stop("Incomplete primary meteorological calendar: ",id," / ",yr)
  }
  h[order(h$city_id,h$date),]
}

temperature_basis_spec <- function(calendar, analysis_records, cfg) {
  # Quantiles/lower boundary: one value per included city-date, full warm-season meteorological calendars for the included analytical years.
  # Upper boundary: equal-city mean of maxima in eligible first-stage records BEFORE
  # lag-history loss. Table 13 support remains the separate actually fitted-row maximum.
  assert_unique(calendar,c("city_id","date")); assert_finite(calendar$temperature,"National calendar")
  assert_unique(analysis_records,c("city_id","date"))
  assert_finite(analysis_records$temperature,"First-stage analysis temperatures")
  if (!setequal(calendar$city_id,analysis_records$city_id)) stop("Basis source city sets differ")
  q <- setNames(unname(quantile(calendar$temperature,c(.33,.5,.66,.9,.95,.99),
    type=cfg$choices$quantile_type)),c("P33","P50","P66","P90","P95","P99"))
  maxima <- aggregate(temperature~city_id,analysis_records,max)
  names(maxima)[2] <- "maximum"
  boundary <- c(min(calendar$temperature),mean(maxima$maximum))
  spec <- list(argvar=list(fun="ns",knots=unname(q[c("P33","P66")]),Boundary.knots=boundary,intercept=FALSE),
    arglag=list(fun="ns",df=cfg$study$lag_df,intercept=TRUE),lag=cfg$study$lag_max,
    quantiles=q,reference=unname(q["P50"]),contrasts=q[c("P90","P95","P99")],
    source=list(years=sort(unique(as.integer(format(calendar$date,"%Y")))),months=cfg$study$warm_months,
      city_ids=sort(unique(calendar$city_id)),n_city_days=nrow(calendar),quantile_type=cfg$choices$quantile_type,
      quantile_sample="Complete meteorological city-date calendar; equal city-day weight; independent of wearable validity",
      national_minimum=min(calendar$temperature),city_analysis_maxima=maxima,
      upper_boundary_sample="Arithmetic mean of city maxima in complete eligible first-stage outcome/covariate records before lag loss",
      climate_sample="Complete primary warm-season meteorological calendar"))
  check_temperature_basis(spec,cfg)
  spec
}

make_grouped_crossbasis <- function(calendar, spec) {
  assert_unique(calendar, c("city_id", "date"))
  year <- format(calendar$date, "%Y")
  group <- interaction(calendar$city_id, year, drop = TRUE)
  for (ix in split(seq_len(nrow(calendar)), group)) {
    if (length(ix) <= spec$lag || any(diff(as.integer(calendar$date[ix])) != 1L))
      stop("Crossbasis requires consecutive ordered daily calendars within city-year")
  }
  dlnm::crossbasis(calendar$temperature, lag = spec$lag,
                  argvar = spec$argvar, arglag = spec$arglag, group = group)
}

temperature_basis <- function(x, spec, center = spec$reference) {
  b <- do.call(dlnm::onebasis, c(list(x = x), spec$argvar))
  ref <- do.call(dlnm::onebasis, c(list(x = center), spec$argvar))
  sweep(unclass(b), 2, as.numeric(ref), "-")
}

reduce_city_crossbasis <- function(cb, model, city_reference) {
  # Explicit coefficients avoid dependence on a caller-local crossbasis name.
  ids <- grep("^cbv[0-9]+\\.l[0-9]+$", names(coef(model)))
  if (length(ids) != ncol(cb) || anyNA(coef(model)[ids])) stop("Non-estimable crossbasis coefficients")
  reduced <- dlnm::crossreduce(cb, type = "overall", coef = coef(model)[ids],
                              vcov = vcov(model)[ids, ids, drop = FALSE],
                              model.link = "logit", cen = city_reference)
  theta <- unname(coef(reduced)); s <- unname(vcov(reduced))
  if (length(theta) != 3L) stop("Expected three reduced temperature coefficients")
  validate_covariance(s, 3L, positive_definite = TRUE)
  list(coefficients = theta, covariance = s, reduction = reduced)
}

# Supplementary Note 5, exact primary covariates. No outcome imputation.
eligible_city_days <- function(d, cfg) {
  fields <- c("high_stress_count", "non_high_stress_count", "temperature",
              "relative_humidity", "wind_speed", "precipitation", "aqi", "holiday")
  yr <- as.integer(format(d$date, "%Y")); mo <- as.integer(format(d$date, "%m"))
  n <- d$high_stress_count + d$non_high_stress_count
  (is.null(cfg$study$primary_years) | yr %in% cfg$study$primary_years) & mo %in% cfg$study$warm_months &
    complete.cases(d[fields]) & !is.na(n) & n >= cfg$study$min_users
}

align_observed_temperature <- function(d, historical_temperature) {
  if (is.null(historical_temperature)) stop("Independent historical temperature calendar is required")
  h <- historical_temperature
  assert_columns(h, c("city_id", "date", "temperature"))
  h$date <- as.Date(h$date)
  assert_unique(h, c("city_id", "date"))
  if (!is.numeric(h$temperature) || any(is.infinite(h$temperature))) stop("Invalid historical temperature")
  key <- function(z) paste(z$city_id, z$date, sep = ":")
  observed <- h$temperature[match(key(d), key(h))]
  paired <- is.finite(d$temperature) & is.finite(observed)
  if (any(abs(d$temperature[paired] - observed[paired]) > 1e-8))
    stop("Outcome and historical temperature sources disagree on matching city-dates")
  d$temperature <- observed
  d
}

build_city_calendar <- function(d, cfg, historical_temperature) {
  if (length(unique(d$city_id)) != 1L) stop("One city required")
  d <- align_observed_temperature(d, historical_temperature)
  eligible <- eligible_city_days(d, cfg)
  d$eligible <- eligible
  years <- cfg$study$primary_years
  dates <- as.Date(unlist(lapply(years, function(y) {
    seq(as.Date(sprintf("%d-04-01", y)), as.Date(sprintf("%d-09-30", y)), by = "day")
  })), origin = "1970-01-01")
  cal <- merge(data.frame(date = dates), d, by = "date", all.x = TRUE, sort = TRUE)
  cal$city_id <- unique(d$city_id)
  # Exposure availability is independent of wearable-outcome availability.
  cal <- align_observed_temperature(cal, historical_temperature)
  cal$eligible[is.na(cal$eligible)] <- FALSE
  cal
}

fit_city_dlnm <- function(d, spec, cfg, historical_temperature) {
  cfg$study$primary_years <- spec$source$years
  d <- validate_outcomes(d)
  d <- align_observed_temperature(d, historical_temperature)
  city <- unique(d$city_id)
  if (length(city) != 1L) stop("fit_city_dlnm requires one city")
  eligible <- eligible_city_days(d, cfg)
  base <- list(city_id = city, eligible_days = sum(eligible),
               eligible_user_days = sum((d$high_stress_count + d$non_high_stress_count)[eligible]),
               model_days = 0L, model_user_days = 0, status = "ineligible", reason = "fewer_than_50_complete_city_days",
               warnings = character(), start = NA_character_, end = NA_character_)
  if (sum(eligible) < cfg$study$min_city_days) return(base)
  tryCatch({
    cal <- build_city_calendar(d, cfg, historical_temperature)
    cb_all <- make_grouped_crossbasis(cal, spec)
    ok <- cal$eligible & complete.cases(cb_all)
    dd <- cal[ok, , drop = FALSE]
    cb <- cb_all[ok, , drop = FALSE]
    # Keep original crossbasis attributes for crossreduce after row subsetting.
    for (nm in setdiff(names(attributes(cb_all)), c("dim", "dimnames"))) attr(cb, nm) <- attr(cb_all, nm)
    base$model_days <- nrow(dd)
    base$model_user_days <- sum(dd$high_stress_count + dd$non_high_stress_count)
    if (!nrow(dd)) stop("No complete lag histories")
    base$start <- as.character(min(dd$date)); base$end <- as.character(max(dd$date))
    dd$dow <- factor(as.POSIXlt(dd$date)$wday, levels = 0:6)
    dd$year <- factor(format(dd$date, "%Y"))
    dd$doy <- as.integer(format(dd$date, "%j"))
    dd$time <- as.numeric(dd$date - as.Date(sprintf("%d-01-01",min(cfg$study$primary_years))))
    f <- cbind(high_stress_count, non_high_stress_count) ~ cb +
      splines::ns(relative_humidity, df = 3) + splines::ns(wind_speed, df = 3) +
      log1p(precipitation) + aqi + dow + holiday +
      splines::ns(doy, df = 3):year + splines::ns(time, df = 1)
    mm <- model.matrix(f, dd)
    if (qr(mm)$rank < ncol(mm)) stop("First-stage design is rank deficient")
    fit <- capture_warnings(glm(f, data = dd, family = quasibinomial(link = "logit"),
                                na.action = na.fail, control = glm.control(maxit = cfg$numerics$glm_maxit)))
    base$warnings <- fit$warnings
    if (!isTRUE(fit$value$converged)) stop("First-stage GLM did not converge")
    if (fit$value$df.residual <= 0) stop("No residual degrees of freedom")
    city_ref <- unname(quantile(d$temperature[eligible], .5, type = cfg$choices$quantile_type))
    red <- reduce_city_crossbasis(cb, fit$value, city_ref)
    diag <- validate_covariance(red$covariance, 3L, positive_definite = TRUE)
    base$status <- "ok"; base$reason <- ""
    base$coefficients <- red$coefficients; base$covariance <- red$covariance
    base$model <- fit$value; base$reduction <- red$reduction
    base$temperature_support <- data.frame(city_id = city,
      historical_min = min(dd$temperature), historical_max = max(dd$temperature),
      model_days = nrow(dd), first_model_date = as.character(min(dd$date)),
      last_model_date = as.character(max(dd$date)))
    base$diagnostics <- c(diag, list(converged = TRUE, dispersion = summary(fit$value)$dispersion,
                                    design_condition = kappa(mm), iterations = fit$value$iter,
                                    lag_incomplete_days = sum(cal$eligible & !complete.cases(cb_all))))
    base
  }, error = function(e) {
    base$status <- "failed"; base$reason <- conditionMessage(e); base
  })
}

run_first_stage <- function(outcomes, cfg, output_dir, spec = NULL, historical_temperature = NULL) {
  if (is.null(historical_temperature)) stop("Independent historical temperature calendar is required")
  historical_temperature$date <- as.Date(historical_temperature$date)
  outcomes <- validate_outcomes(outcomes)
  outcomes <- align_observed_temperature(outcomes, historical_temperature)
  if (!is.null(spec)) cfg$study$primary_years <- spec$source$years
  eligible <- eligible_city_days(outcomes, cfg)
  nvalid <- table(outcomes$city_id[eligible])
  primary_ids <- names(nvalid)[nvalid >= cfg$study$min_city_days]
  if (!length(primary_ids)) {
    ids <- unique(outcomes$city_id)
    counts <- as.integer(table(factor(outcomes$city_id[eligible], levels = ids)))
    write_table(data.frame(city_id = ids, status = "ineligible", eligible_days = counts,
      reason = "fewer_than_50_complete_city_days"), file.path(output_dir, "city_fit_ledger.csv"))
    stop("No eligible cities; see city_fit_ledger.csv")
  }
  primary <- outcomes[eligible & outcomes$city_id %in% primary_ids, ]
  if (is.null(spec)) cfg$study$primary_years <- sort(unique(as.integer(format(primary$date,"%Y"))))
  calendar <- primary_temperature_calendar(historical_temperature, primary_ids, cfg)
  if (is.null(spec)) spec <- temperature_basis_spec(calendar, primary, cfg)
  # Supplied population-stratum specifications inherit the PRIMARY national sample;
  # do not recompute national quantiles, boundaries or lag specifications within a subgroup.
  check_temperature_basis(spec, cfg)
  fits <- lapply(split(outcomes, outcomes$city_id), function(d)
    fit_city_dlnm(d, spec, cfg, historical_temperature[historical_temperature$city_id == d$city_id[1], ]))
  ledger <- do.call(rbind, lapply(fits, function(x) data.frame(
    city_id = x$city_id, status = x$status, reason = x$reason,
    eligible_days = x$eligible_days, model_days = x$model_days,
    eligible_user_days = x$eligible_user_days, model_user_days = x$model_user_days,
    start = x$start, end = x$end, warnings = paste(x$warnings, collapse = " | "),
    dispersion = if (is.null(x$diagnostics)) NA_real_ else x$diagnostics$dispersion,
    covariance_condition = if (is.null(x$diagnostics)) NA_real_ else x$diagnostics$condition,
    design_condition = if (is.null(x$diagnostics)) NA_real_ else x$diagnostics$design_condition)))
  write_table(ledger, file.path(output_dir, "city_fit_ledger.csv"))
  save_model(spec, file.path(output_dir, "basis_specification.rds"))
  write_table(data.frame(percentile=names(spec$quantiles),temperature=unname(spec$quantiles)),
    file.path(output_dir,"national_temperature_quantiles.csv"))
  # Preserve every failure before stopping, rather than silently selecting successful cities.
  if (any(ledger$status == "failed")) stop("Eligible city models failed; see city_fit_ledger.csv")
  kept <- fits[vapply(fits, function(x) x$status == "ok", logical(1))]
  coefficients <- data.frame(city_id = names(kept), do.call(rbind, lapply(kept, `[[`, "coefficients")))
  names(coefficients) <- c("city_id", "b1", "b2", "b3")
  covariance <- do.call(rbind, lapply(kept, function(x) {
    z <- expand.grid(row = 1:3, col = 1:3)
    data.frame(city_id = x$city_id, z, covariance = as.vector(x$covariance))
  }))
  climate <- do.call(rbind, lapply(split(calendar, calendar$city_id), function(d)
    data.frame(city_id = d$city_id[1], mean_temperature = mean(d$temperature),
               temperature_range = diff(range(d$temperature)))))
  write_table(coefficients, file.path(output_dir, "first_stage_coefficients.csv"))
  write_table(covariance, file.path(output_dir, "first_stage_covariance_matrices.csv"))
  write_table(climate, file.path(output_dir, "city_climate.csv"))
  covariance_list <- lapply(kept, function(z) {
    v <- z$covariance; dimnames(v) <- list(c("b1", "b2", "b3"), c("b1", "b2", "b3")); v
  })
  result <- list(coefficients = coefficients, covariance = covariance_list, climate = climate,
                 basis = spec, ledger = ledger,
                 temperature_support = do.call(rbind, lapply(kept, `[[`, "temperature_support")))
  write_table(result$temperature_support, file.path(output_dir, "historical_model_temperature_support.csv"))
  save_model(result, file.path(output_dir, "first_stage.rds"))
  result
}

run_dlnm_module <- function(input, output, root = .dlnm_root) {
  x <- read_module_input(input, 1L, c("city_day_outcomes", "historical_city_temperature"))
  out <- prepare_output(output, root, x); cfg <- analysis_settings()
  z <- run_first_stage(x$city_day_outcomes, cfg, out,
                       historical_temperature = x$historical_city_temperature)
  stopifnot(nrow(z$coefficients) == length(z$covariance), all(z$ledger$status %in% c("ok", "ineligible")))
  finish_module(z, x, file.path(out, "first_stage.rds"))
}

if (sys.nframe() == 0L) {
  args <- module_arguments()
  run_dlnm_module(args$input, args$output)
}
