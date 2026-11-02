#!/usr/bin/env Rscript

.meta_file <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  if (length(paths)) tail(paths, 1)[[1]] else {
    arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (!length(arg)) stop("Load this file with source() or Rscript")
    sub("^--file=", "", arg[1])
  }
})
.meta_root <- dirname(dirname(normalizePath(.meta_file, winslash = "/", mustWork = TRUE)))
source(file.path(.meta_root, "01_DLNMs/DLNM.R"), local = TRUE)

# Supplementary Note 5: multivariate pooling and effect modification.
assemble_meta_data <- function(first, attributes, climate = first$climate) {
  attributes <- validate_attributes(attributes)
  b <- first$coefficients; s <- first$covariance
  assert_columns(b, c("city_id", "b1", "b2", "b3")); assert_unique(b, "city_id")
  assert_unique(climate, "city_id")
  if (!is.list(s) || is.data.frame(s) || is.null(names(s)) || anyDuplicated(names(s)) ||
      !setequal(names(s), b$city_id)) stop("Coefficient/covariance city sets differ")
  if (any(!b$city_id %in% attributes$city_id) || any(!b$city_id %in% climate$city_id))
    stop("Missing city attributes or climate covariates")
  data <- cbind(b, attributes[match(b$city_id, attributes$city_id), setdiff(names(attributes), "city_id"), drop = FALSE],
                climate[match(b$city_id, climate$city_id), c("mean_temperature", "temperature_range"), drop = FALSE])
  rownames(data) <- NULL
  assert_finite(data$mean_temperature, "City mean temperature")
  assert_finite(data$temperature_range, "City temperature range")
  if (any(data$temperature_range <= 0)) stop("City temperature ranges must be positive")
  covs <- lapply(b$city_id, function(id) {
    v <- s[[id]]
    validate_covariance(v, 3L, positive_definite = TRUE)
    if (!identical(rownames(v), c("b1", "b2", "b3")) || !identical(colnames(v), c("b1", "b2", "b3")))
      stop("Each city covariance must use b1,b2,b3 row and column order")
    v
  })
  names(covs) <- b$city_id
  list(data = data, S = covs, basis = first$basis)
}

subset_meta <- function(meta, keep) {
  if (is.character(keep)) keep <- match(keep, meta$data$city_id)
  if (anyNA(keep)) stop("Unknown city in meta subset")
  list(data = meta$data[keep, , drop = FALSE], S = meta$S[keep], basis = meta$basis)
}

fit_mixmeta <- function(meta, rhs = ~ mean_temperature + temperature_range, cfg, method = "reml") {
  if (!method %in% c("ml", "reml")) stop("Only ML and REML are supported")
  check_temperature_basis(meta$basis, cfg)
  mf <- model.frame(rhs, data = meta$data, na.action = na.fail)
  x <- model.matrix(rhs, mf)
  if (qr(x)$rank != ncol(x)) stop("Second-stage design is rank deficient")
  if (nrow(x) <= ncol(x)) stop("Insufficient cities for second-stage model")
  for (i in seq_along(meta$S)) validate_covariance(meta$S[[i]], 3L, TRUE)
  f <- as.formula(paste("cbind(b1,b2,b3)", paste(deparse(rhs), collapse = " ")))
  captured <- capture_warnings(mixmeta::mixmeta(f, S = meta$S, data = meta$data,
    method = method, bscov = "unstr", na.action = na.fail,
    control = list(maxiter = cfg$numerics$meta_maxiter, reltol = cfg$numerics$meta_tolerance)))
  model <- captured$value
  if (!isTRUE(model$converged)) stop("mixmeta did not converge")
  if (length(coef(model)) != 3L * ncol(x)) stop("Unexpected mixmeta coefficient dimension")
  validate_parameter_covariance(coef(model), vcov(model))
  vd <- validate_covariance(vcov(model), length(coef(model)), TRUE)
  pd <- validate_covariance(model$Psi, 3L)
  flags <- character()
  if (vd$condition > cfg$numerics$condition_warning) flags <- c(flags, "Ill-conditioned fixed-effect covariance")
  if (pd$condition > cfg$numerics$condition_warning) flags <- c(flags, "Random-effects covariance near boundary; inspect diagnostics")
  for (flag in flags) warning(flag)
  list(model = model, rhs = terms(rhs, data = meta$data), xlevels = stats::.getXlevels(terms(rhs), mf),
       contrasts = attr(x, "contrasts"), columns = colnames(x),
       beta = coef(model), V = vcov(model), Psi = model$Psi, basis = meta$basis,
       climate_reference = c(mean_temperature = mean(meta$data$mean_temperature),
                             temperature_range = mean(meta$data$temperature_range)),
       data = meta$data, diagnostics = list(n_cities = nrow(x), n_predictors = ncol(x),
         converged = model$converged, method = method, logLik = as.numeric(logLik(model)),
         fixed_covariance = vd, random_covariance = pd, design_condition = kappa(x), warnings = unique(c(captured$warnings, flags))))
}

meta_model_matrix <- function(fit, newdata) {
  mf <- model.frame(fit$rhs, newdata, na.action = na.fail, xlev = fit$xlevels)
  x <- model.matrix(fit$rhs, mf, contrasts.arg = fit$contrasts)
  if (!identical(colnames(x), fit$columns)) stop("Prediction design columns differ from fitted model")
  x
}

reduced_design <- function(x) {
  # mixmeta stores vector coefficients by predictor, with 3 outcomes per predictor.
  kronecker(matrix(as.numeric(x), nrow = 1), diag(3))
}

point_design <- function(fit, newdata, temperatures) {
  x <- meta_model_matrix(fit, newdata)
  b <- temperature_basis(temperatures, fit$basis)
  if (nrow(x) == 1L) return(b %*% reduced_design(x[1, ]))
  if (nrow(x) != nrow(b)) stop("Temperature and modifier prediction rows must match")
  out <- matrix(NA_real_, nrow(x), length(fit$beta))
  for (j in seq_len(ncol(x))) out[, ((j - 1L) * 3L + 1L):(j * 3L)] <- b * x[, j]
  out
}

reference_data <- function(fit, n = 1L) {
  d <- fit$data[rep(1L, n), , drop = FALSE]
  for (nm in c("mean_temperature", "temperature_range")) d[[nm]] <- fit$climate_reference[[nm]]
  for (nm in grep("^context_z_", names(d), value = TRUE)) d[[nm]] <- 0
  d
}

# Contrasts use the full joint coefficient covariance.
linear_contrast <- function(beta, V, d, alpha = .05) {
  validate_covariance(V, length(beta))
  if (!is.null(names(beta))) validate_parameter_covariance(beta, V)
  if (length(d) != length(beta) || any(!is.finite(d)) || any(!is.finite(beta)))
    stop("Contrast and coefficient dimensions or values invalid")
  estimate <- sum(d * beta)
  variance <- as.numeric(crossprod(d, V %*% d))
  if (variance < 0) stop("Negative contrast variance")
  se <- sqrt(variance)
  if (se == 0 && estimate != 0) stop("Nonzero contrast has zero uncertainty")
  z <- if (se == 0) 0 else estimate / se
  critical <- qnorm(1 - alpha / 2)
  lo <- estimate - critical * se; hi <- estimate + critical * se
  data.frame(log_estimate = estimate, variance = variance, se = se, z = z,
             p_value = 2 * pnorm(-abs(z)), ratio = exp(estimate),
             ratio_low = exp(lo), ratio_high = exp(hi),
             attenuation = 100 * (1 - exp(estimate)),
             attenuation_low = 100 * (1 - exp(hi)),
             attenuation_high = 100 * (1 - exp(lo)))
}

wald_contrast <- function(beta, V, L) {
  L <- as.matrix(L)
  validate_covariance(V, length(beta))
  if (!is.null(names(beta))) validate_parameter_covariance(beta, V)
  if (ncol(L) != length(beta) || qr(L)$rank != nrow(L)) stop("Invalid or redundant Wald contrast")
  d <- as.numeric(L %*% beta); W <- L %*% V %*% t(L)
  validate_covariance(W, nrow(L), TRUE)
  statistic <- as.numeric(crossprod(d, solve(W, d)))
  data.frame(statistic = statistic, df = nrow(L), p_value = pchisq(statistic, nrow(L), lower.tail = FALSE))
}

coefficient_wald <- function(fit, columns) {
  j <- match(columns, fit$columns)
  if (anyNA(j) || !length(j)) stop("Unknown Wald coefficient columns")
  ids <- unlist(lapply(j, function(z) (3L * (z - 1L) + 1L):(3L * z)))
  wald_contrast(fit$beta, fit$V, diag(length(fit$beta))[ids, , drop = FALSE])
}

predict_or <- function(fit, newdata, temperatures, alpha = .05) {
  L <- point_design(fit, newdata, temperatures)
  out <- do.call(rbind, lapply(seq_len(nrow(L)), function(i) linear_contrast(fit$beta, fit$V, L[i, ], alpha)))
  names(out)[names(out) == "ratio"] <- "or"
  names(out)[names(out) == "ratio_low"] <- "or_low"
  names(out)[names(out) == "ratio_high"] <- "or_high"
  cbind(temperature = temperatures, out[c("or", "or_low", "or_high", "log_estimate", "se", "p_value")])
}

compare_high_low <- function(fit, high, low, temperatures, alpha = .05) {
  D <- point_design(fit, high, temperatures) - point_design(fit, low, temperatures)
  ans <- do.call(rbind, lapply(seq_len(nrow(D)), function(i) linear_contrast(fit$beta, fit$V, D[i, ], alpha)))
  names(ans)[names(ans) == "ratio"] <- "ror"
  names(ans)[names(ans) == "ratio_low"] <- "ror_low"
  names(ans)[names(ans) == "ratio_high"] <- "ror_high"
  cbind(temperature = temperatures, ans)
}

compare_curves <- function(fit, high, low) {
  L <- reduced_design(meta_model_matrix(fit, high)[1, ]) - reduced_design(meta_model_matrix(fit, low)[1, ])
  wald_contrast(fit$beta, fit$V, L)
}

green_attributes <- function() c("cooling", "gvi", "accessibility")

compute_cutoffs <- function(attributes, cfg, scope = c("national", "region")) {
  scope <- match.arg(scope)
  groups <- if (scope == "national") list(national = attributes) else split(attributes, attributes$region)
  do.call(rbind, lapply(names(groups), function(r) {
    d <- groups[[r]]
    do.call(rbind, lapply(green_attributes(), function(a) {
      assert_finite(d[[a]], a)
      q <- quantile(d[[a]], probs = cfg$choices$tertile_probs, type = cfg$choices$quantile_type, names = FALSE)
      if (q[1] >= q[2]) stop("Tied tertiles cannot identify three groups: ", r, " / ", a)
      data.frame(scope = scope, region = r, attribute = a, lower = q[1], upper = q[2], n_cities = nrow(d))
    }))
  }))
}

assign_tertiles <- function(values, lower, upper) {
  assert_finite(values, "Green-space values")
  if (any(!is.finite(c(lower, upper))) || any(lower >= upper)) stop("Invalid tertile boundaries")
  factor(ifelse(values <= lower, "Low", ifelse(values <= upper, "Medium", "High")),
         levels = c("Low", "Medium", "High"))
}

attach_green_groups <- function(meta, attribute, cutoffs, scope) {
  if (!all(cutoffs$scope == scope)) stop("National and regional tertile scopes cannot be mixed")
  cuts <- cutoffs[cutoffs$attribute == attribute, ]
  region <- if (scope == "national") rep("national", nrow(meta$data)) else as.character(meta$data$region)
  ix <- match(region, cuts$region)
  if (anyNA(ix)) stop("Missing green-space cutoff")
  meta$data$green_group <- assign_tertiles(meta$data[[attribute]], cuts$lower[ix], cuts$upper[ix])
  meta
}

fit_categorical <- function(meta, attribute, cutoffs, cfg, extra_terms = character()) {
  meta <- attach_green_groups(meta, attribute, cutoffs, "national")
  counts <- table(meta$data$green_group)
  if (any(counts == 0)) stop("Empty national green-space group")
  if (any(counts < cfg$numerics$small_cell_warning)) warning("Small national tertile cell")
  rhs <- reformulate(c("green_group", "mean_temperature", "temperature_range", extra_terms))
  fit <- fit_mixmeta(meta, rhs, cfg)
  fit$attribute <- attribute; fit$group_counts <- counts
  fit
}

summarize_categorical <- function(fit, temperatures, alpha = .05, grid = NULL) {
  data <- reference_data(fit, 3)
  data$green_group <- factor(c("Low", "Medium", "High"), levels = c("Low", "Medium", "High"))
  ors <- do.call(rbind, lapply(1:3, function(i) cbind(group = as.character(data$green_group[i]),
    predict_or(fit, data[i, , drop = FALSE], temperatures, alpha))))
  high <- data[3, , drop = FALSE]; low <- data[1, , drop = FALSE]
  tests <- rbind(cbind(test = "overall_curve", coefficient_wald(fit, c("green_groupMedium", "green_groupHigh"))),
                 cbind(test = "low_high_curve", compare_curves(fit, high, low)))
  curves <- if (is.null(grid)) NULL else do.call(rbind, lapply(1:3, function(i)
    cbind(group = as.character(data$green_group[i]), predict_or(fit, data[i, , drop = FALSE], grid, alpha))))
  list(or = ors, contrasts = compare_high_low(fit, high, low, temperatures, alpha), tests = tests, curves = curves)
}

standardize_modifier <- function(x, cfg) {
  assert_finite(x, "Continuous modifier")
  scale <- sd(x)
  if (scale <= 0) stop("Continuous modifier has zero variance")
  list(value = (x - mean(x)) / scale, center = mean(x), scale = scale, range = range(x))
}

fit_continuous <- function(meta, attribute, cfg) {
  std <- standardize_modifier(meta$data[[attribute]], cfg)
  meta$data$modifier_z <- std$value
  b <- splines::ns(std$value, df = 3, intercept = FALSE)
  for (j in 1:3) meta$data[[paste0("modifier_s", j)]] <- b[, j]
  linear_rhs <- ~ modifier_z + mean_temperature + temperature_range
  spline_rhs <- ~ modifier_s1 + modifier_s2 + modifier_s3 + mean_temperature + temperature_range
  ml_linear <- fit_mixmeta(meta, linear_rhs, cfg, "ml")
  ml_spline <- fit_mixmeta(meta, spline_rhs, cfg, "ml")
  lr <- 2 * (as.numeric(logLik(ml_spline$model)) - as.numeric(logLik(ml_linear$model)))
  if (lr < -1e-5) stop("Nested ML spline likelihood smaller than linear likelihood; inspect optimization")
  lr <- max(lr, 0)
  df <- length(ml_spline$beta) - length(ml_linear$beta)
  if (df != 6L) stop("Expected six additional fixed-effect parameters")
  p <- pchisq(lr, df = df, lower.tail = FALSE)
  chosen <- if (p < cfg$study$alpha) "spline" else "linear"
  final <- fit_mixmeta(meta, if (chosen == "linear") linear_rhs else spline_rhs, cfg, "reml")
  cols <- if (chosen == "linear") "modifier_z" else paste0("modifier_s", 1:3)
  final$attribute <- attribute
  final$modifier <- c(std[c("center", "scale", "range")], list(form = chosen,
    spline_knots = attr(b, "knots"), spline_boundary = attr(b, "Boundary.knots")))
  final$comparison <- data.frame(attribute = attribute, statistic = lr, df = df, p_nonlinearity = p,
                                 selected = chosen)
  final$modifier_test <- coefficient_wald(final, cols)
  final$ml_models <- list(linear = ml_linear, spline = ml_spline)
  final
}

continuous_newdata <- function(fit, values, climate = NULL) {
  if (any(!is.finite(values)) || any(values < fit$modifier$range[1] | values > fit$modifier$range[2]))
    stop("Continuous modifier predictions outside empirical support")
  d <- reference_data(fit, length(values))
  z <- (values - fit$modifier$center) / fit$modifier$scale
  d$modifier_z <- z
  b <- splines::ns(z, knots = fit$modifier$spline_knots,
                   Boundary.knots = fit$modifier$spline_boundary, intercept = FALSE)
  for (j in 1:3) d[[paste0("modifier_s", j)]] <- b[, j]
  if (!is.null(climate)) {
    assert_columns(climate, c("mean_temperature", "temperature_range"))
    if (nrow(climate) != length(values)) stop("Climate prediction rows do not match")
    d$mean_temperature <- climate$mean_temperature; d$temperature_range <- climate$temperature_range
  }
  d
}

run_meta_analysis <- function(first, attributes, cfg, output_dir) {
  # Inherit complete-calendar national P50/P90/P95/P99 and shared fitted basis from Module 1.
  check_temperature_basis(first$basis,cfg)
  meta <- assemble_meta_data(first, attributes)
  cuts <- compute_cutoffs(meta$data, cfg, "national")
  write_table(cuts, file.path(output_dir, "national_tertile_cutoffs.csv"))
  grid <- seq(first$basis$argvar$Boundary.knots[1], max(first$basis$argvar$Boundary.knots[2], first$basis$contrasts), length.out = 121)
  national <- fit_mixmeta(meta, cfg = cfg)
  write_table(predict_or(national, reference_data(national), first$basis$contrasts), file.path(output_dir, "national_contrasts.csv"))
  write_table(predict_or(national, reference_data(national), grid), file.path(output_dir, "national_curve.csv"))
  categorical <- continuous <- list()
  for (a in green_attributes()) {
    catfit <- fit_categorical(meta, a, cuts, cfg)
    categorical[[a]] <- catfit
    tables <- summarize_categorical(catfit, first$basis$contrasts, cfg$study$alpha, grid)
    for (nm in names(tables)) write_table(tables[[nm]], file.path(output_dir, paste0(a, "_categorical_", nm, ".csv")))
    continuous[[a]] <- fit_continuous(meta, a, cfg)
    write_table(continuous[[a]]$comparison, file.path(output_dir, paste0(a, "_functional_form.csv")))
    write_table(continuous[[a]]$modifier_test, file.path(output_dir, paste0(a, "_continuous_wald.csv")))
    save_model(catfit, file.path(output_dir, paste0(a, "_categorical_model.rds")))
    save_model(continuous[[a]], file.path(output_dir, paste0(a, "_continuous_model.rds")))
  }
  result <- list(meta = meta, national = national, categorical = categorical, continuous = continuous,
                 national_cutoffs = cuts)
  save_model(result, file.path(output_dir, "meta_analysis.rds"))
  result
}

run_meta_module <- function(input, output, root = .meta_root) {
  x <- read_module_input(input, 2L, c("first_stage", "city_attributes"))
  out <- prepare_output(output, root, x); cfg <- analysis_settings()
  z <- run_meta_analysis(x$first_stage, x$city_attributes, cfg, out)
  finish_module(z, x, file.path(out, "meta_analysis.rds"))
}

if (sys.nframe() == 0L) {
  args <- module_arguments()
  run_meta_module(args$input, args$output)
}
