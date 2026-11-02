#!/usr/bin/env Rscript

.future_file <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  if (length(paths)) tail(paths, 1)[[1]] else {
    arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (!length(arg)) stop("Load this file with source() or Rscript")
    sub("^--file=", "", arg[1])
  }
})
.future_root <- dirname(dirname(normalizePath(.future_file, winslash = "/", mustWork = TRUE)))
source(file.path(.future_root, "02_Meta_Regression/Meta_Regression.R"), local = TRUE)

# Cumulative OR evaluated at each day's temperature and city GVI (Note 8).
surface_design <- function(fit, temperature, gvi, cfg) {
  if (fit$attribute != "gvi") stop("Projection surface must be the continuous GVI model")
  if (length(temperature) != length(gvi)) stop("Surface inputs have unequal lengths")
  assert_finite(temperature, "Projection temperature")
  boundary <- fit$basis$argvar$Boundary.knots
  out <- temperature < boundary[1] | temperature > boundary[2]
  # Natural-spline extrapolation retains linear tails.
  check_temperature_basis(fit$basis, cfg)
  design <- point_design(fit, continuous_newdata(fit, gvi), temperature)
  attr(design, "outside_temperature_support") <- out
  design
}

predict_surface <- function(fit, temperature, gvi, cfg, beta = fit$beta) {
  D <- surface_design(fit, temperature, gvi, cfg)
  if (length(beta) != ncol(D)) stop("Wrong surface parameter vector")
  data.frame(temperature = temperature, gvi = gvi, log_or = as.numeric(D %*% beta),
             or = exp(as.numeric(D %*% beta)), extrapolated = attr(D, "outside_temperature_support"))
}

odds_to_probability <- function(or, p0) {
  if (any(!is.finite(or) | or <= 0) || any(!is.finite(p0) | p0 < 0 | p0 > 1)) stop("Invalid OR or baseline probability")
  plogis(log(or) + qlogis(p0))
}

estimate_baseline_probabilities <- function(outcomes, cfg, basis) {
  check_temperature_basis(basis,cfg)
  # National P50 from the fitted primary meteorological calendar.
  reference <- basis$reference
  d <- validate_outcomes(outcomes)
  d$month <- as.integer(format(d$date, "%m"))
  n <- d$high_stress_count + d$non_high_stress_count
  cfg$study$primary_years <- NULL # The supplied baseline dataset already defines its observation period.
  keep <- eligible_city_days(d, cfg) & !is.na(d$temperature) & d$temperature <= reference
  d <- d[keep, ]; d$total <- d$high_stress_count + d$non_high_stress_count
  groups <- split(d, interaction(d$city_id, d$month, drop = TRUE))
  if (!length(groups)) return(data.frame(city_id = character(), month = integer(), p0 = numeric(),
    reference_days = integer(), reference_user_days = numeric(), high_stress_count = numeric(),
    first_reference_date = character(), last_reference_date = character()))
  result <- do.call(rbind, lapply(groups, function(z) {
    p <- sum(z$high_stress_count) / sum(z$total)
    data.frame(city_id = z$city_id[1], month = z$month[1], p0 = p,
                reference_days = nrow(z), reference_user_days = sum(z$total),
                high_stress_count = sum(z$high_stress_count),
                first_reference_date = as.character(min(z$date)), last_reference_date = as.character(max(z$date)))
  }))
  rownames(result) <- NULL
  result
}

daily_excess <- function(log_or, p0, temperature, reference) {
  if (length(reference)!=1L || !is.finite(reference)) stop("A fitted national P50 is required")
  if (any(!is.finite(log_or)) || any(!is.finite(p0) | p0 < 0 | p0 > 1)) stop("Invalid daily risk input")
  p <- plogis(qlogis(p0) + log_or)
  delta <- p - p0
  # Retain only positive heat-related excess probability.
  delta <- pmax(delta, 0)
  delta[temperature <= reference] <- 0
  1e5 * delta
}

adaptation_metrics <- function(existing, adapted, historical) {
  additional <- existing - historical
  avoided <- existing - adapted
  residual <- adapted - historical
  data.frame(future_burden = adapted, existing_burden = existing, historical_burden = historical,
              additional_burden = additional, avoided_burden = avoided,
              offset_percent = ifelse(additional > 0, 100 * avoided / additional, NA_real_),
              residual_additional_burden = residual, adaptation_gap = pmax(residual, 0))
}

classify_gap <- function(existing, combined, historical) {
  ifelse(existing <= historical, "no_excess_current",
          ifelse(combined <= historical, "offset_by_combined", "residual_after_combined"))
}

make_adaptation_targets <- function(attributes, vegetation, efficiency, cfg) {
  attributes <- validate_attributes(attributes)
  assert_columns(vegetation, c("city_id", "month", "ndvi"))
  assert_columns(efficiency, c("city_id", "month", "cooling_efficiency"))
  assert_unique(vegetation, c("city_id", "month")); assert_unique(efficiency, c("city_id", "month"))
  assert_finite(vegetation$ndvi, "NDVI"); assert_finite(efficiency$cooling_efficiency, "Cooling efficiency")
  if (any(vegetation$ndvi < -1 | vegetation$ndvi > 1)) stop("NDVI outside [-1,1]")
  grid <- expand.grid(city_id = attributes$city_id, month = cfg$study$warm_months, stringsAsFactors = FALSE)
  key <- function(x) paste(x$city_id, x$month, sep = ":")
  iv <- match(key(grid), key(vegetation)); ie <- match(key(grid), key(efficiency))
  if (anyNA(iv) || anyNA(ie)) stop("Missing city-month vegetation/cooling efficiency")
  a <- attributes[match(grid$city_id, attributes$city_id), ]
  grid$region <- a$region; grid$gvi_observed <- a$gvi
  grid$ndvi_observed <- vegetation$ndvi[iv]; grid$cooling_efficiency <- efficiency$cooling_efficiency[ie]
  gtarget <- tapply(attributes$gvi, attributes$region, quantile, probs = .95, type = cfg$choices$quantile_type)
  grid$gvi_target <- unname(gtarget[grid$region])
  grid$gvi_enhanced <- pmax(grid$gvi_observed, grid$gvi_target)
  g <- interaction(grid$region, grid$month, drop = TRUE)
  vtarget <- tapply(grid$ndvi_observed, g, quantile, probs = .95, type = cfg$choices$quantile_type)
  grid$ndvi_target <- unname(vtarget[as.character(g)])
  grid$vegetation_increment <- pmax(0, grid$ndvi_target - grid$ndvi_observed)
  grid$temperature_reduction <- grid$cooling_efficiency * grid$vegetation_increment
  grid
}

apply_adaptation <- function(temperature, gvi, temperature_reduction, gvi_enhanced,
                              scenario = c("existing", "gvi", "cooling", "combined")) {
  scenario <- match.arg(scenario)
  list(temperature = if (scenario %in% c("cooling", "combined")) temperature - temperature_reduction else temperature,
       gvi = if (scenario %in% c("gvi", "combined")) gvi_enhanced else gvi)
}

classify_pathway <- function(gvi_offset, cooling_offset, gap_class) {
  # Note 8: >=10 percentage-point advantage, unless both pathways offset >=50%.
  result <- rep("not_applicable", length(gap_class))
  eligible <- gap_class == "residual_after_combined"
  if (any(!is.finite(gvi_offset[eligible]) | !is.finite(cooling_offset[eligible])))
    stop("Residual-gap cities require finite single-pathway offsets")
  result[eligible] <- "no_single_pathway_priority"
  both <- gvi_offset >= 50 & cooling_offset >= 50
  result[eligible & !both & gvi_offset - cooling_offset >= 10] <- "gvi"
  result[eligible & !both & cooling_offset - gvi_offset >= 10] <- "cooling"
  result
}

# SI Table 13: equal-member percentages by SSP and decade, with explicit denominators.
temperature_exceedance <- function(future, historical_support, cfg, upper_boundary) {
  assert_columns(future, c("city_id", "date", "ssp", "climate_model", "temperature"))
  assert_unique(future, c("city_id", "date", "ssp", "climate_model"))
  assert_columns(historical_support, c("city_id", "historical_max"))
  assert_unique(historical_support, "city_id")
  assert_finite(historical_support$historical_max, "Historical model-row maximum")
  if (length(upper_boundary) != 1L || !is.finite(upper_boundary))
    stop("A finite upper boundary from the fitted basis is required")
  future$date <- as.Date(future$date)
  year <- as.integer(format(future$date, "%Y"))
  month <- as.integer(format(future$date, "%m"))
  future$period <- NA_character_
  for (start in cfg$choices$decade_starts)
    future$period[year >= start & year < start + cfg$choices$decade_length] <- paste0(start, "s")
  future <- future[month %in% cfg$study$warm_months & !is.na(future$period), ]
  if (!nrow(future)) stop("No projected warm-season days in requested decades")
  assert_finite(future$temperature, "Future temperature")
  future$historical_max <- historical_support$historical_max[match(future$city_id, historical_support$city_id)]
  groups <- split(future, interaction(future$ssp, future$period, drop = TRUE))
  do.call(rbind, lapply(groups, function(d) {
    make_row <- function(measure, threshold) {
      evaluable <- is.finite(threshold)
      members <- split(seq_len(nrow(d)), d$climate_model)
      percentages <- vapply(members, function(ix) {
        ok <- ix[evaluable[ix]]
        if (!length(ok)) NA_real_ else 100 * mean(d$temperature[ok] > threshold[ok])
      }, numeric(1))
      data.frame(ssp = d$ssp[1], period = d$period[1], measure = measure,
        n_members = length(members), n_cities = length(unique(d$city_id[evaluable])),
        n_city_day_members = sum(evaluable),
        n_exceeding_city_day_members = sum(d$temperature[evaluable] > threshold[evaluable]),
        n_days_without_reference = sum(!evaluable), percentage = mean(percentages),
        spline_upper_boundary = upper_boundary)
    }
    rbind(make_row("above_city_historical_maximum", d$historical_max),
      make_row("above_prediction_spline_boundary", rep(upper_boundary, nrow(d))))
  }))
}

# Supplementary Note 8: 1,000 joint parametric simulations; full model covariance.
draw_model_parameters <- function(fit, cfg) {
  validate_parameter_covariance(fit$beta, fit$V)
  b <- cfg$bootstrap$repetitions
  if (length(b) != 1L || !is.finite(b) || b != 1000L || b != floor(b))
    stop("Supplementary Note 8 specifies 1000 parametric simulations")
  set.seed(cfg$bootstrap$seed)
  draws <- MASS::mvrnorm(b, mu = fit$beta, Sigma = fit$V)
  colnames(draws) <- names(fit$beta)
  draws
}

validate_parameter_draws <- function(fit, draws) {
  if (!is.matrix(draws) || !is.numeric(draws) || nrow(draws) < 2L ||
      !identical(colnames(draws), names(fit$beta)) || any(!is.finite(draws)))
    stop("Joint draws must have exact fitted parameter order")
  draws
}

evaluate_exposure_batch <- function(prepared, draws, cfg) {
  eta <- prepared$D %*% t(draws)
  delta <- pmax(plogis(eta + qlogis(prepared$p0)) - prepared$p0, 0)
  delta[prepared$temperature <= prepared$reference, ] <- 0
  totals <- rowsum(1e5 * delta * prepared$annual_day_weight, prepared$group, reorder = FALSE)
  s <- prepared$summary
  key <- paste(s$city_id, s$ssp, s$period, sep = "|")
  levels <- sort(unique(key)); group <- match(key, levels)
  # Same equal-GCM averaging order as evaluate_exposure(), for every parameter draw.
  means <- rowsum(totals, group, reorder = TRUE) / tabulate(group, length(levels))
  list(keys = s[match(levels, key), c("city_id", "ssp", "period"), drop = FALSE], burden = means)
}

burden_metric_matrices <- function(existing, adapted, historical) {
  additional <- existing - historical
  avoided <- existing - adapted
  residual <- adapted - historical
  offset <- 100 * avoided / additional
  offset[additional <= 0] <- NA_real_
  list(future_burden = adapted, existing_burden = existing, historical_burden = historical,
    additional_burden = additional, avoided_burden = avoided, offset_percent = offset,
    residual_additional_burden = residual, adaptation_gap = pmax(residual, 0))
}

projection_batch <- function(prepared, draws, cfg, reference) {
  h <- evaluate_exposure_batch(prepared$historical, draws, cfg)
  f <- lapply(prepared$future, evaluate_exposure_batch, draws = draws, cfg = cfg)
  keys <- f$existing$keys
  for (z in f) if (!identical(z$keys, keys)) stop("Bootstrap scenario keys differ")
  hist <- h$burden[match(keys$city_id, h$keys$city_id), , drop = FALSE]
  w <- prepared$attributes$population_2020[match(keys$city_id, prepared$attributes$city_id)]
  groups <- split(seq_len(nrow(keys)), paste(keys$ssp, keys$period, sep = "|"))
  signature <- function(d) do.call(paste, c(d[c("city_id", "ssp", "period", "adaptation")], sep = "|"))
  reference_key <- signature(reference)
  result <- lapply(burden_metric_matrices(f$existing$burden, f$existing$burden, hist),
    function(x) matrix(NA_real_, nrow(reference), nrow(draws)))
  for (a in names(f)) {
    values <- burden_metric_matrices(f$existing$burden, f[[a]]$burden, hist)
    rows <- match(signature(cbind(keys, adaptation = a)), reference_key)
    for (nm in names(values)) result[[nm]][rows, ] <- values[[nm]]
    for (ix in groups) {
      average <- function(m) matrix(colSums(m[ix, , drop = FALSE] * w[ix]) / sum(w[ix]), nrow = 1)
      national <- burden_metric_matrices(average(f$existing$burden), average(f[[a]]$burden), average(hist))
      national$adaptation_gap <- average(values$adaptation_gap)
      k <- data.frame(city_id = "NATIONAL", ssp = keys$ssp[ix[1]], period = keys$period[ix[1]], adaptation = a)
      row <- match(signature(k), reference_key)
      for (nm in names(national)) result[[nm]][row, ] <- national[[nm]]
    }
  }
  result
}

propagate_uncertainty <- function(prepared, fit, cfg, draws = NULL) {
  if (is.null(draws)) draws <- draw_model_parameters(fit, cfg)
  draws <- validate_parameter_draws(fit, draws)
  point <- evaluate_prepared_projection(prepared, fit$beta, cfg)$metrics
  keys <- c("city_id", "ssp", "period", "adaptation")
  reference <- point[keys]
  metrics <- setdiff(names(point), keys)
  samples <- setNames(lapply(metrics, function(nm) matrix(NA_real_, nrow(point), nrow(draws))), metrics)
  batches <- split(seq_len(nrow(draws)), ceiling(seq_len(nrow(draws)) / cfg$bootstrap$batch_size))
  for (ix in batches) {
    if (ix[1] == 1L || (ix[1] - 1L) %% 200L == 0L)
      message("Joint parametric simulations: ", ix[1], "-", tail(ix, 1), "/", nrow(draws))
    values <- projection_batch(prepared, draws[ix, , drop = FALSE], cfg, reference)
    for (nm in metrics) samples[[nm]][, ix] <- values[[nm]]
  }
  limits <- do.call(rbind, lapply(metrics, function(nm) {
    m <- samples[[nm]]
    q <- t(apply(m, 1, function(z) if (all(!is.finite(z))) c(NA_real_, NA_real_) else
      quantile(z[is.finite(z)], c(.025, .975), names = FALSE, type = cfg$choices$quantile_type)))
    cbind(reference, metric = nm, lower = q[, 1], upper = q[, 2],
      valid_draws = rowSums(is.finite(m)), total_draws = nrow(draws))
  }))
  list(limits = limits, parameter_draws = draws,
    coverage = paste("Joint multivariate-normal draws from the selected GVI model's fixed coefficients and full vcov.",
      "Identical draws are used for historical climate and all future scenarios.",
      "Baseline probabilities, cooling efficiency, targets, climate members and 2020 populations are fixed.",
      "Only health-model fixed-parameter uncertainty is propagated.",
      "Offset intervals condition on positive additional burden; valid_draws identifies undefined ratios."))
}

validate_temperature_calendar <- function(d, cfg, historical = FALSE) {
  # Each climate_model ID represents one independent GCM with one daily series.
  fields <- c("city_id", "date", "temperature", if (!historical) c("ssp", "climate_model"))
  assert_columns(d, fields); d$date <- as.Date(d$date)
  assert_unique(d, c("city_id", "date", if (!historical) c("ssp", "climate_model")))
  assert_finite(d$temperature, "Daily temperature")
  d$year <- as.integer(format(d$date, "%Y")); d$month <- as.integer(format(d$date, "%m"))
  d <- d[d$month %in% cfg$study$warm_months, ]
  years <- if (historical) sort(unique(d$year)) else seq(cfg$study$future_years[1], cfg$study$future_years[2])
  d <- d[d$year %in% years, ]
  if (!nrow(d)) stop("Empty temperature analysis period")
  if (historical) { d$ssp <- "historical"; d$climate_model <- "observed" }
  else {
    if (!setequal(unique(d$ssp), cfg$study$ssps)) stop("Missing or unexpected SSP scenarios")
  }
  groups <- interaction(d$city_id, d$ssp, d$climate_model, drop = TRUE)
  for (ix in split(seq_len(nrow(d)), groups)) {
    counts <- table(factor(d$year[ix], levels = years))
    if (any(counts != 183L)) stop("Incomplete April-September calendar; no annualization of partial years")
  }
  # Each city must carry the same SSP/member set, preventing moving spatial denominators.
  membership <- split(paste(d$ssp, d$climate_model, sep = ":"), d$city_id)
  sets <- lapply(membership, function(x) sort(unique(x)))
  if (!all(vapply(sets, identical, logical(1), sets[[1]]))) stop("Climate member sets differ across cities")
  d
}

collapse_exposure_days <- function(d, cfg, historical = FALSE) {
  if (historical) { d$period <- "historical"; divisor <- length(unique(d$year)) }
  else {
    d$period <- NA_character_
    for (start in cfg$choices$decade_starts)
      d$period[d$year >= start & d$year < start + cfg$choices$decade_length] <- paste0(start, "s")
    d <- d[!is.na(d$period), ]; divisor <- cfg$choices$decade_length
  }
  # Exact duplicate exposure days can be collapsed without approximating T or OR.
  # Counts divided by number of years preserve mean ANNUAL burden.
  columns <- c("city_id", "ssp", "climate_model", "period", "month", "temperature")
  key <- do.call(paste, c(d[columns], sep = "|"))
  unique_key <- unique(key); index <- match(key, unique_key)
  out <- d[match(unique_key, key), columns]
  out$annual_day_weight <- tabulate(index, nbins = length(unique_key)) / divisor
  rownames(out) <- NULL
  out
}

prepare_projection <- function(fit, historical, future, attributes, baseline, targets, cfg) {
  attributes <- validate_attributes(attributes)
  historical <- validate_temperature_calendar(historical, cfg, TRUE)
  future <- validate_temperature_calendar(future, cfg, FALSE)
  if (!setequal(unique(historical$city_id), attributes$city_id) ||
      !setequal(unique(future$city_id), attributes$city_id)) stop("Projection city sets must match attributes exactly")
  assert_columns(baseline, c("city_id", "month", "p0", "reference_days", "reference_user_days"))
  assert_unique(baseline, c("city_id", "month"))
  if (any(!is.finite(baseline$p0) | baseline$p0 < 0 | baseline$p0 > 1) ||
      any(!is.finite(baseline$reference_days) | baseline$reference_days < 1L)) stop("Invalid baseline probabilities/provenance counts")
  if (any(!is.finite(baseline$reference_user_days) | baseline$reference_user_days < cfg$study$min_users * baseline$reference_days) ||
      any(baseline$reference_days != floor(baseline$reference_days)) ||
      any(baseline$reference_user_days != floor(baseline$reference_user_days))) stop("Invalid baseline reference counts")
  h <- collapse_exposure_days(historical, cfg, TRUE)
  f <- collapse_exposure_days(future, cfg, FALSE)
  key <- function(d) paste(d$city_id, d$month, sep = ":")
  prepare_one <- function(d, scenario) {
    ib <- match(key(d), key(baseline)); it <- match(key(d), key(targets))
    if (anyNA(ib)) stop("Missing low-temperature baseline cells: ", paste(unique(key(d)[is.na(ib)]), collapse = ", "))
    if (anyNA(it)) stop("Missing adaptation targets")
    t <- targets[it, ]
    change <- apply_adaptation(d$temperature, t$gvi_observed, t$temperature_reduction, t$gvi_enhanced, scenario)
    ci <- match(d$city_id, fit$data$city_id)
    if (anyNA(ci)) stop("Projection city missing from continuous model")
    D <- surface_design(fit, change$temperature, change$gvi, cfg)
    group_key <- paste(d$city_id, d$ssp, d$climate_model, d$period, sep = "|")
    groups <- unique(group_key); gid <- match(group_key, groups)
    summary <- d[match(groups, group_key), c("city_id", "ssp", "climate_model", "period")]
    list(D = D, p0 = baseline$p0[ib], temperature = change$temperature, reference = fit$basis$reference,
          annual_day_weight = d$annual_day_weight, group = gid, summary = summary,
          extrapolated_annual_days = sum(d$annual_day_weight[attr(D, "outside_temperature_support")]))
  }
  list(historical = prepare_one(h, "existing"),
        future = setNames(lapply(c("existing", "gvi", "cooling", "combined"), function(a) prepare_one(f, a)),
                           c("existing", "gvi", "cooling", "combined")), attributes = attributes)
}

evaluate_exposure <- function(prepared, beta, cfg) {
  burden <- daily_excess(as.numeric(prepared$D %*% beta), prepared$p0, prepared$temperature, prepared$reference) * prepared$annual_day_weight
  totals <- rowsum(burden, prepared$group, reorder = FALSE)[, 1]
  result <- prepared$summary; result$burden <- as.numeric(totals)
  groups <- split(seq_len(nrow(result)), paste(result$city_id, result$ssp, result$period, sep = "|"))
  do.call(rbind, lapply(groups, function(ix) {
    z <- result[ix, ]; out <- z[1, c("city_id", "ssp", "period"), drop = FALSE]
    # Equal GCM weight AFTER nonlinear daily risk conversion and annual/period aggregation.
    out$burden <- mean(z$burden); out
  }))
}

evaluate_prepared_projection <- function(prepared, beta, cfg, historical_beta = beta) {
  hist <- evaluate_exposure(prepared$historical, historical_beta, cfg)
  fut <- lapply(prepared$future, evaluate_exposure, beta = beta, cfg = cfg)
  # All scenarios use identical row keys; changes in temperature cannot change membership.
  keys <- c("city_id", "ssp", "period")
  for (z in fut) if (!identical(z[keys], fut$existing[keys])) stop("Scenario rows do not align")
  h <- hist$burden[match(fut$existing$city_id, hist$city_id)]
  if (anyNA(h)) stop("Historical baseline city missing")
  city <- do.call(rbind, lapply(names(fut), function(a)
    cbind(fut[[a]][keys], adaptation = a, adaptation_metrics(fut$existing$burden, fut[[a]]$burden, h))))
  attrs <- prepared$attributes
  ww <- attrs$population_2020[match(fut$existing$city_id, attrs$city_id)]
  groups <- split(seq_len(nrow(fut$existing)), paste(fut$existing$ssp, fut$existing$period, sep = "|"))
  national <- do.call(rbind, lapply(groups, function(ix) {
    existing <- weighted.mean(fut$existing$burden[ix], ww[ix]); history <- weighted.mean(h[ix], ww[ix])
    do.call(rbind, lapply(names(fut), function(a) {
      metrics <- adaptation_metrics(existing, weighted.mean(fut[[a]]$burden[ix], ww[ix]), history)
      # Population mean of city-specific positive gaps, distinct from signed residual.
      metrics$adaptation_gap <- weighted.mean(pmax(fut[[a]]$burden[ix] - h[ix], 0), ww[ix])
      cbind(data.frame(city_id = "NATIONAL", ssp = fut$existing$ssp[ix[1]], period = fut$existing$period[ix[1]], adaptation = a), metrics)
    }))
  }))
  classes <- fut$existing[keys]
  classes$gap_class <- classify_gap(fut$existing$burden, fut$combined$burden, h)
  additional <- fut$existing$burden - h
  ga <- ifelse(additional > 0, 100 * (fut$existing$burden - fut$gvi$burden) / additional, NA_real_)
  ca <- ifelse(additional > 0, 100 * (fut$existing$burden - fut$cooling$burden) / additional, NA_real_)
  classes$higher_contributing_pathway <- classify_pathway(ga, ca, classes$gap_class)
  classes$gvi_offset_percent <- ga; classes$cooling_offset_percent <- ca
  metrics <- rbind(city, national); rownames(metrics) <- NULL; rownames(classes) <- NULL
  list(metrics = metrics, classifications = classes)
}

run_future_projections <- function(fit, historical, future, attributes, outcomes, vegetation,
                                   efficiency, cfg, output_dir, historical_support) {
  check_temperature_basis(fit$basis, cfg)
  validate_parameter_covariance(fit$beta, fit$V)
  if (!setequal(attributes$city_id, fit$data$city_id))
    stop("Adaptation targets require the complete primary-model city universe")
  outcomes <- align_observed_temperature(validate_outcomes(outcomes), historical)
  baseline <- estimate_baseline_probabilities(outcomes, cfg, fit$basis)
  write_table(baseline, file.path(output_dir, "city_month_baseline_probabilities.csv"))
  grid <- expand.grid(city_id = attributes$city_id, month = cfg$study$warm_months,
                      stringsAsFactors = FALSE)
  key <- function(d) paste(d$city_id, d$month, sep = ":")
  missing <- grid[!key(grid) %in% key(baseline), , drop = FALSE]
  if (nrow(missing)) {
    write_table(missing, file.path(output_dir, "missing_baseline_cells.csv"))
    stop("Projection cannot estimate p0 for city-month cells: ", paste(key(missing), collapse = ", "))
  }
  targets <- make_adaptation_targets(attributes, vegetation, efficiency, cfg)
  write_table(targets, file.path(output_dir, "adaptation_targets.csv"))
  prepared <- prepare_projection(fit, historical, future, attributes, baseline, targets, cfg)
  point <- evaluate_prepared_projection(prepared, fit$beta, cfg)
  write_table(point$metrics, file.path(output_dir, "burden_metrics.csv"))
  write_table(point$classifications, file.path(output_dir, "city_adaptation_classification.csv"))
  write_table(temperature_exceedance(future, historical_support, cfg, fit$basis$argvar$Boundary.knots[2]),
              file.path(output_dir, "temperature_exceedance.csv"))
  uncertainty <- propagate_uncertainty(prepared, fit, cfg)
  write_table(uncertainty$limits, file.path(output_dir, "uncertainty_intervals.csv"))
  save_model(list(point = point, uncertainty = uncertainty, targets = targets, baseline = baseline),
              file.path(output_dir, "projections.rds"))
  invisible(list(point = point, uncertainty = uncertainty, targets = targets, baseline = baseline))
}

run_future_module <- function(input, output, root = .future_root) {
  fields <- c("gvi_model", "historical_city_temperature", "future_city_temperature", "city_attributes",
    "city_day_outcomes", "city_month_vegetation", "city_month_cooling_efficiency",
    "city_month_baseline_probabilities", "temperature_support")
  x <- read_module_input(input, 4L, fields)
  out <- prepare_output(output, root, x); cfg <- analysis_settings()
  aligned <- align_observed_temperature(validate_outcomes(x$city_day_outcomes), x$historical_city_temperature)
  p0 <- estimate_baseline_probabilities(aligned, cfg, x$gvi_model$basis)
  supplied <- x$city_month_baseline_probabilities
  assert_columns(supplied, c("city_id", "month", "p0", "reference_days", "reference_user_days"))
  assert_unique(supplied, c("city_id", "month"))
  key <- function(d) paste(d$city_id, d$month, sep=":")
  if (!setequal(key(p0), key(supplied))) stop("Supplied and reconstructed baseline city-month sets differ")
  j <- match(key(p0), key(supplied))
  for (nm in c("p0", "reference_days", "reference_user_days"))
    if (!isTRUE(all.equal(p0[[nm]], supplied[[nm]][j], tolerance=1e-10)))
      stop("Supplied baseline disagrees with reference-day counts: ", nm)
  z <- run_future_projections(x$gvi_model, x$historical_city_temperature, x$future_city_temperature,
    x$city_attributes, aligned, x$city_month_vegetation, x$city_month_cooling_efficiency,
    cfg, out, x$temperature_support)
  stopifnot(nrow(z$uncertainty$parameter_draws) == 1000L)
  finish_module(z, x, file.path(out, "projections.rds"))
}

if (sys.nframe() == 0L) {
  args <- module_arguments()
  run_future_module(args$input, args$output)
}
