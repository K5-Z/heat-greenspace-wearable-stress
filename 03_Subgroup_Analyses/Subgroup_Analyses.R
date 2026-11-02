#!/usr/bin/env Rscript

.subgroup_file <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  if (length(paths)) tail(paths, 1)[[1]] else {
    arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (!length(arg)) stop("Load this file with source() or Rscript")
    sub("^--file=", "", arg[1])
  }
})
.subgroup_root <- dirname(dirname(normalizePath(.subgroup_file, winslash = "/", mustWork = TRUE)))
source(file.path(.subgroup_root, "02_Meta_Regression/Meta_Regression.R"), local = TRUE)

run_population_analysis <- function(subgroup_outcomes, primary, attributes, cuts, cfg, output_dir, historical_temperature) {
  # P50/P95, analysis years and lag specification are inherited from the complete PRIMARY national calendar.
  check_temperature_basis(primary$basis,cfg)
  d <- validate_outcomes(subgroup_outcomes, subgroup = TRUE)
  expected <- list(sex = c("Female", "Male"), age = c("18-34", "35-49", "50-64", "65plus"))
  actual <- unique(paste(d$stratum_type, d$stratum, sep = ":"))
  required <- unlist(lapply(names(expected), function(t) paste(t, expected[[t]], sep = ":")))
  if (!setequal(actual, required)) stop("Population strata do not match the six specified groups")
  result <- list()
  for (type in names(expected)) for (stratum in expected[[type]]) {
    key <- paste(type, stratum, sep = "_")
    dd <- d[d$stratum_type == type & d$stratum == stratum, ]
    message("Population DLNM: ", key)
    first <- run_first_stage(dd, cfg, file.path(output_dir, key, "first_stage"), spec = primary$basis,
                             historical_temperature = historical_temperature)
    # Climate reference sample is frozen to primary cities; green cutoffs are never recomputed here.
    meta <- assemble_meta_data(first, attributes, climate = primary$climate)
    fits <- list()
    for (a in c("cooling", "gvi")) {
      fits[[a]] <- fit_categorical(meta, a, cuts, cfg)
      fits[[a]]$climate_reference <- c(mean_temperature = mean(primary$climate$mean_temperature),
                                      temperature_range = mean(primary$climate$temperature_range))
      summary <- summarize_categorical(fits[[a]], primary$basis$contrasts["P95"])
      for (nm in c("or", "contrasts", "tests")) write_table(summary[[nm]],
        file.path(output_dir, key, paste0(a, "_", nm, ".csv")))
      write_table(data.frame(group = names(fits[[a]]$group_counts), n_cities = as.integer(fits[[a]]$group_counts)),
        file.path(output_dir, key, paste0(a, "_counts.csv")))
    }
    result[[key]] <- list(first_stage = first, fits = fits)
    save_model(result[[key]], file.path(output_dir, key, "population_models.rds"))
  }
  result
}

# Seven census covariates (Note 6/Table 10C) versus ten (Note 7/Table 12B).
context_variables <- function(specification = c("census7", "extended10")) {
  specification <- match.arg(specification)
  census <- c("pct_age_0_14", "pct_age_65_plus", "pct_urban_builtup", "education_years",
               "pct_one_person_households", "pct_nonlocal_registered", "pct_agricultural_employment")
  c(census, if (specification == "extended10") c("gdp_100m_cny", "nighttime_light", "ac_per_100_households"))
}

add_context <- function(meta, context, specification, cfg) {
  vars <- context_variables(specification)
  assert_columns(context, c("city_id", vars)); assert_unique(context, "city_id")
  ix <- match(meta$data$city_id, context$city_id)
  if (anyNA(ix)) stop("Missing city contextual covariates; no automatic city exclusion")
  d <- context[ix, vars, drop = FALSE]
  centers <- scales <- numeric(length(vars)); names(centers) <- names(scales) <- vars
  for (v in vars) {
    assert_finite(d[[v]], v)
    centers[v] <- mean(d[[v]]); scales[v] <- sd(d[[v]])
    if (scales[v] <= 0) stop("Constant contextual covariate: ", v)
    meta$data[[paste0("context_z_", v)]] <- (d[[v]] - centers[v]) / scales[v]
  }
  list(meta = meta, terms = paste0("context_z_", vars),
       metadata = data.frame(variable = vars, center = centers, scale = scales, specification = specification))
}

run_contextual_adjustment <- function(meta, context, cuts, cfg, output_dir, primary_fits) {
  results <- list()
  for (a in green_attributes()) {
    primary <- primary_fits[[a]]
    base <- summarize_categorical(primary, meta$basis$contrasts)
    for (specification in c("census7", "extended10")) {
      extended <- add_context(meta, context, specification, cfg)
      fit <- fit_categorical(extended$meta, a, cuts, cfg, extended$terms)
      adjusted <- summarize_categorical(fit, meta$basis$contrasts)
      key <- paste(a, specification, sep = "_")
      results[[key]] <- list(primary = primary, adjusted = fit)
      for (nm in c("or", "contrasts", "tests")) {
        write_table(rbind(cbind(specification = "primary", base[[nm]]),
                          cbind(specification = specification, adjusted[[nm]])),
                    file.path(output_dir, paste0(key, "_", nm, ".csv")))
      }
      difference <- data.frame(temperature = base$contrasts$temperature,
        attenuation_primary = base$contrasts$attenuation,
        attenuation_adjusted = adjusted$contrasts$attenuation,
        attenuation_difference = adjusted$contrasts$attenuation - base$contrasts$attenuation,
        ror_difference = adjusted$contrasts$ror - base$contrasts$ror)
      write_table(difference, file.path(output_dir, paste0(key, "_difference.csv")))
      write_table(extended$metadata, file.path(output_dir, paste0(key, "_covariate_scaling.csv")))
    }
  }
  save_model(results, file.path(output_dir, "context_models.rds"))
  results
}

fit_bioclimatic <- function(meta, attribute, cutoffs, cfg, extra_terms = character()) {
  keep <- meta$data$region %in% region_levels()
  meta <- subset_meta(meta, keep)
  meta$data$region <- factor(meta$data$region, levels = region_levels())
  if (any(table(meta$data$region) == 0)) stop("All four bioclimatic regions are required")
  meta <- attach_green_groups(meta, attribute, cutoffs, "region")
  counts <- as.data.frame(table(region = meta$data$region, group = meta$data$green_group))
  if (any(counts$Freq == 0)) stop("Empty region-by-tertile cell")
  if (any(counts$Freq < cfg$numerics$small_cell_warning)) warning("Small region-by-tertile cell; inspect counts")
  rhs <- reformulate(c("region * green_group", "mean_temperature", "temperature_range", extra_terms))
  fit <- fit_mixmeta(meta, rhs, cfg)
  fit$attribute <- attribute; fit$group_counts <- counts
  # Joint region-by-green-space interaction model.
  fit$interaction_test <- coefficient_wald(fit, grep(":", fit$columns, value = TRUE))
  fit
}

summarize_bioclimatic <- function(fit, temperatures) {
  rows <- contrasts <- tests <- list()
  for (r in region_levels()) {
    d <- reference_data(fit, 3)
    d$region <- factor(r, levels = region_levels())
    d$green_group <- factor(c("Low", "Medium", "High"), levels = c("Low", "Medium", "High"))
    rows[[r]] <- do.call(rbind, lapply(1:3, function(i)
      cbind(region = r, group = as.character(d$green_group[i]), predict_or(fit, d[i, , drop = FALSE], temperatures))))
    contrasts[[r]] <- cbind(region = r, compare_high_low(fit, d[3, ], d[1, ], temperatures))
    x <- meta_model_matrix(fit, d)
    L <- rbind(reduced_design(x[2, ]) - reduced_design(x[1, ]),
               reduced_design(x[3, ]) - reduced_design(x[1, ]))
    tests[[r]] <- rbind(cbind(region = r, test = "overall_curve", wald_contrast(fit$beta, fit$V, L)),
                        cbind(region = r, test = "low_high_curve", compare_curves(fit, d[3, ], d[1, ])))
  }
  list(or = do.call(rbind, rows), contrasts = do.call(rbind, contrasts),
       tests = do.call(rbind, tests), interaction = fit$interaction_test)
}

run_bioclimatic_analysis <- function(meta, context, cfg, output_dir) {
  regional <- subset_meta(meta, meta$data$region %in% region_levels())
  cuts <- compute_cutoffs(regional$data, cfg, "region")
  write_table(cuts, file.path(output_dir, "regional_tertile_cutoffs.csv"))
  write_table(data.frame(city_id = meta$data$city_id, region = meta$data$region,
                         included = meta$data$region %in% region_levels()), file.path(output_dir, "region_eligibility.csv"))
  result <- list()
  adjusted <- add_context(regional, context, "census7", cfg)
  for (a in green_attributes()) {
    result[[a]] <- list(primary = fit_bioclimatic(regional, a, cuts, cfg))
    if (a %in% c("cooling", "gvi"))
      result[[a]]$census7 <- fit_bioclimatic(adjusted$meta, a, cuts, cfg, adjusted$terms)
    summaries <- lapply(result[[a]], summarize_bioclimatic, temperatures = meta$basis$contrasts[2])
    for (nm in c("or", "contrasts", "tests", "interaction"))
      write_table(do.call(rbind, lapply(names(summaries), function(s) cbind(specification = s, summaries[[s]][[nm]]))),
                  file.path(output_dir, paste0(a, "_", nm, ".csv")))
    write_table(result[[a]]$primary$group_counts, file.path(output_dir, paste0(a, "_counts.csv")))
    if (!is.null(summaries$census7)) {
      b <- summaries$primary$contrasts; z <- summaries$census7$contrasts
      write_table(data.frame(region = b$region, temperature = b$temperature,
                             attenuation_primary = b$attenuation, attenuation_adjusted = z$attenuation,
                             attenuation_difference = z$attenuation - b$attenuation),
                  file.path(output_dir, paste0(a, "_adjustment_difference.csv")))
    }
  }
  save_model(result, file.path(output_dir, "joint_bioclimatic_models.rds"))
  result
}

run_heterogeneity <- function(first, second, subgroup_outcomes, attributes, context, cfg, output_dir, historical_temperature) {
  population <- run_population_analysis(subgroup_outcomes, first, attributes,
                                        second$national_cutoffs, cfg, file.path(output_dir, "population"), historical_temperature)
  bioclimatic <- run_bioclimatic_analysis(second$meta, context, cfg, file.path(output_dir, "bioclimatic"))
  contextual <- run_contextual_adjustment(second$meta, context, second$national_cutoffs, cfg,
                                          file.path(output_dir, "contextual"), primary_fits = second$categorical)
  invisible(list(population = population, bioclimatic = bioclimatic, contextual = contextual))
}

run_subgroup_module <- function(input, output, root = .subgroup_root) {
  fields <- c("first_stage", "city_attributes", "city_contextual_covariates",
              "subgroup_city_day_outcomes", "historical_city_temperature")
  x <- read_module_input(input, 3L, fields)
  out <- prepare_output(output, root, x); cfg <- analysis_settings()
  meta <- assemble_meta_data(x$first_stage, x$city_attributes)
  cuts <- compute_cutoffs(meta$data, cfg, "national")
  primary_fits <- setNames(lapply(green_attributes(), function(a) fit_categorical(meta, a, cuts, cfg)), green_attributes())
  second <- list(meta = meta, national_cutoffs = cuts, categorical = primary_fits)
  z <- run_heterogeneity(x$first_stage, second, x$subgroup_city_day_outcomes,
    x$city_attributes, x$city_contextual_covariates, cfg, out, x$historical_city_temperature)
  stopifnot(length(z$population) == 6L)
  finish_module(z, x, file.path(out, "subgroup_analysis.rds"))
}

if (sys.nframe() == 0L) {
  args <- module_arguments()
  run_subgroup_module(args$input, args$output)
}
