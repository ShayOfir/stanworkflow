# wflow.R -- Bayesian modelling workflow S3 class
# Wraps the compile -> prior_predictive -> fit -> diagnose -> infer lifecycle.

# ---- Constructor -----------------------------------------------------------

#' Create a new wflow object
#'
#' A `wflow` encapsulates a complete Bayesian modelling workflow:
#' model code, priors, sampler settings, fitted models, data, and inference.
#'
#' @param model_name   character; name of a subfolder under `models_folder`
#'   containing `model_code.stan`, `priors.csv`, and optionally
#'   `diagnostics_params.csv`. `priors.csv` has two columns (prior name,
#'   value), one row per `prior_*` data variable; a vector prior is written in
#'   one cell with values separated by `;` (e.g. `"0; 0.5"`). Matrix priors are
#'   set in R (`wf$compiled$priors$prior_x <- matrix(...)`). Names the model
#'   does not declare give a warning and are ignored.
#' @param models_folder character; root folder for model definitions
#'   (default `"./models/"`).
#' @param prior_sampler_params_list named list of overrides for `sampler_params()`.
#' @param post_sampler_params_list named list of overrides for `sampler_params()`.
#'   NULL uses package defaults.
#' @param diagnostics_params optional character vector of the variables
#'   `run_diagnostics()` checks (e.g. `c("beta_0", "beta_1", "sigma")`).
#'   If NULL, read from `diagnostics_params.csv` (or `.txt`) in the model
#'   folder: one name per line or comma-separated; a header line such as
#'   `parameters`, quotes, blank lines and `#` comments are ignored. Names the
#'   model does not define give a warning and are ignored. If neither is given,
#'   the model's parameters and transformed parameters are checked.
#' @return an object of class `wflow`
#' @export
wflow <- function(model_name,
                  models_folder      = "./models/",
                  prior_sampler_params_list = NULL,
                  post_sampler_params_list = NULL,
                  diagnostics_params = NULL) {

  stopifnot(is.character(model_name), length(model_name) == 1L)

  model_dir <- file.path(models_folder, model_name)
  if (!dir.exists(model_dir))
    stop(sprintf("Model directory not found: '%s'", model_dir))

  # ---- load model code ----------------------------------------------------
  stan_file <- file.path(model_dir, "model_code.stan")
  if (!file.exists(stan_file))
    stop(sprintf("'model_code.stan' not found in '%s'", model_dir))

  # ---- load priors ---------------------------------------------------------
  priors_file <- file.path(model_dir, "priors.csv")
  priors_list <- if (file.exists(priors_file)) {
    # a cell may hold a vector as "0; 0.5; 1" (one value per group/coefficient)
    read.csv(priors_file, colClasses = "character") |>
      tibble::deframe() |>
      as.list() |>
      lapply(function(v) as.numeric(trimws(strsplit(v, ";", fixed = TRUE)[[1]])))
  } else {
    message(sprintf("No priors.csv found in '%s'; priors will need to be set manually.", model_dir))
    list()
  }

  # ---- load diagnostics params ---------------------------------------------
  # the argument wins; otherwise diagnostics_params.csv / .txt in the model folder
  if (is.null(diagnostics_params)) {
    diag_file <- file.path(model_dir, c("diagnostics_params.csv", "diagnostics_params.txt"))
    diag_file <- diag_file[file.exists(diag_file)][1]
    if (!is.na(diag_file)) diagnostics_params <- .read_name_list(diag_file)
  }

  # ---- compile -------------------------------------------------------------

  compiled <- stanmodel() |> compile_from_file(stan_file)
  # names the model does not declare are usually leftovers from a switched
  # prior family (e.g. beta(a, b) -> normal(mean, sd)): warn and ignore them
  unknown <- setdiff(names(priors_list), names(compiled$priors))
  if (length(unknown))
    warning(sprintf("priors.csv in '%s' has value(s) the model does not declare (ignored): %s",
                    model_dir, paste(unknown, collapse = ", ")), call. = FALSE)
  compiled$priors <- modifyList(compiled$priors, priors_list[setdiff(names(priors_list), unknown)])
  diagnostics_params <- .check_diagnostics_params(diagnostics_params, compiled, model_dir)

  # ---- sampler params ------------------------------------------------------
  # a prior predictive check needs fewer draws: 500 + 500 (as in regimen())
  prior_base <- modifyList(sampler_params(), list(iter_warmup = 500L, iter_sampling = 500L))
  prior_sp <- sampler_params(modifyList(prior_base, as.list(prior_sampler_params_list %||% list())))

  post_sp <- if (!is.null(post_sampler_params_list)) {
    modifyList(sampler_params(), post_sampler_params_list)
  } else {
    sampler_params()
  }

  structure(
    list(
      model_name         = model_name,         # name of the model (subfolder under models_folder) 
      models_folder      = models_folder,      # root folder for model codes and accesory data
      diagnostics_params = diagnostics_params, # optional vector of parameter names for diagnostics
      prior_sampler_params = prior_sp,      # sampler params used for prior predictive fit
      post_sampler_params = post_sp,        # sampler params used for posterior fit
      compiled           = compiled,        # stanmodel, not yet fitted
      prior_fit          = NULL,            # stanmodel after prior predictive fit
      post_fit           = NULL,            # stanmodel after posterior fit
      data               = NULL,            # data list: $data, $cond_df, $extra
      diagnostics        = NULL,            # output of diagnose()
      inference          = list(),          # named list of infer_contrasts() results
      prior_log          = list(),          # one entry per prior predictive run (see prior_log())
      tuned             = FALSE             # whether hyperpriors have been tuned
    ),
    class = "wflow"
  )
}

# names from a forgiving text/CSV file: one per line or comma/semicolon
# separated; quotes, blank lines, "#" comments and a header word are ignored
#' @keywords internal
.read_name_list <- function(file) {
  x <- readLines(file, warn = FALSE, encoding = "UTF-8")
  x <- sub("^\ufeff", "", x)                    # byte-order mark
  x <- sub("#.*$", "", x)                              # comments
  x <- unlist(strsplit(x, "[,;\t]"))
  x <- trimws(gsub("[\"']", "", x))
  x <- x[nzchar(x)]
  header <- c("parameter", "parameters", "param", "params", "variable", "variables", "name", "names")
  if (length(x) && tolower(x[1]) %in% header) x <- x[-1]
  unique(x)
}

# keep diagnostics names the model defines (any block; "beta[2]" -> "beta"); warn about the rest
#' @keywords internal
.check_diagnostics_params <- function(params, compiled, model_dir) {
  if (is.null(params)) return(NULL)
  v <- tryCatch(compiled$compiled_model$variables(), error = function(e) NULL)
  known <- unlist(lapply(v[c("parameters", "transformed_parameters", "generated_quantities")], names))
  if (is.null(known)) return(params)
  bad <- params[!sub("\\[.*$", "", params) %in% known]
  if (length(bad))
    warning(sprintf("diagnostics_params for '%s': not defined by the model (ignored): %s",
                    basename(model_dir), paste(bad, collapse = ", ")), call. = FALSE)
  keep <- setdiff(params, bad)
  if (length(keep)) keep else NULL
}

# ---- print ------------------------------------------------------------------

#' @exportS3Method
print.wflow <- function(x, ...) {
  cat("<wflow>\n")
  cat("  model         :", x$model_name, "\n")
  cat("  compiled      :", !is.null(x$compiled), "\n")
  cat("  prior_fit     :", !is.null(x$prior_fit), "\n")
  cat("  hyperpriors tuned:", x$tuned, "\n")
  cat("  post_fit      :", !is.null(x$post_fit), "\n")
  cat("  data attached :", !is.null(x$data), "\n")
  if (length(x$inference) > 0L)
    cat("  inference     :", paste(names(x$inference), collapse = ", "), "\n")
  invisible(x)
}

# ---- attach_data ------------------------------------------------------------

#' Attach data to a wflow object
#'
#' @param wf      a `wflow` object
#' @param data    data frame of observations (trial-level)
#' @param extra   named list of additional Stan data (matrices, scalars, etc.)
#' @param cond_df optional condition-level data frame for `*_cond_array` vars
#' @param outcome_var name of the outcome column in `data` (default `"y"`)
#' @param ... additional arguments passed to methods
#' @return the `wflow` with `$data` populated
#' @export
attach_data <- function(wf, ...) UseMethod("attach_data")

#' @rdname attach_data
#' @exportS3Method
attach_data.wflow <- function(wf,
                               data,
                               extra       = list(),
                               cond_df     = NULL,
                               outcome_var = "y",
                               ...) {
  wf$data <- list(
    data        = data,
    extra       = extra,
    cond_df     = cond_df,
    outcome_var = outcome_var
  )
  wf
}

# ---- prior_predictive -------------------------------------------------------

#' Run predictive checks
#'
#' Fits the model with `skip_likelihood = 1` and stores the result in
#' `$prior_fit`. Use `ppc(wf$prior_fit, ...)` for visual checks.
#'
#' @param wf  a `wflow` with data attached via `attach_data()`
#' @param prior_predictive logical; run prior predictive checks (default TRUE)
#' @param sampler_params a `sampler_params` object, containing the parameters for the prior predictive sampling
#' @param ... additional arguments passed to `ppc.stanmodel()`.
#' @return the `wflow` with `$prior_fit` populated
#' @export
run_predictive_checks <- function(wf,...) UseMethod("run_predictive_checks")

#' @rdname run_predictive_checks
#' @exportS3Method
run_predictive_checks.wflow <- function(wf, prior_predictive = TRUE, sampler_params = NULL,  ...) {
  if (prior_predictive) .prior_check(wf, sampler_params, ...) else .posterior_check(wf, ...)
}

# print each outcome's predictive-check results
.print_checks <- function(results, label) {
  for (i in seq_along(results)) {
    cat(sprintf("%s predictive check for outcome '%s':
", label, names(results)[i]))
    if (inherits(results[[i]], "ppc_results_grouped")) display(results[[i]])
    else print(results[[i]], plot = TRUE)
  }
}

# ppc() returns one result or a named list per outcome; always return the list
.as_outcome_list <- function(res, outcome_var) {
  if (inherits(res, c("ppc_results_grouped", "ppc_results")))
    res <- stats::setNames(list(res), outcome_var)
  res
}

#' @keywords internal
.prior_check <- function(wf, sampler_params = NULL, ...) {
  .check_data(wf)
  .check_compiled(wf)
  # the wflow's prior sampler settings; a `sampler_params` (or list) argument overrides them
  sp <- wf$prior_sampler_params %||% wf$compiled$sampler_params
  if (!is.null(sampler_params)) sp <- modifyList(sp, as.list(sampler_params))
  wf$prior_fit <- wf$compiled
  wf$prior_fit$sampler_params <- sampler_params(as.list(sp))

  message("Running prior predictive checks for '", wf$model_name, "'...")
  wf$prior_fit <- wf$prior_fit |>
    fit(
      data             = wf$data$data,
      extra            = wf$data$extra,
      cond_df          = wf$data$cond_df,
      outcome_var      = wf$data$outcome_var,
      prior_predictive = TRUE
    )
  wf$prior_ppc_results <- .as_outcome_list(ppc(wf$prior_fit, ...), wf$data$outcome_var)
  wf$prior_log <- c(wf$prior_log %||% list(), list(.prior_log_entry(wf)))
  .print_checks(wf$prior_ppc_results, "Prior")
  wf
}

#' @keywords internal
.posterior_check <- function(wf, ...) {
  message("Performing posterior predictive checks '", wf$model_name, "'.")
  if (is.null(wf$post_fit) || !inherits(wf$post_fit, "stanmodel"))
    stop("Posterior fit not found or not as `stanmodel`. Please run `fit_posterior()` first.")
  wf$posterior_ppc_results <- .as_outcome_list(ppc(wf$post_fit, ...), wf$data$outcome_var)
  .print_checks(wf$posterior_ppc_results, "Posterior")
  wf
}

# --- accept_priors -------------------------------------------------------------
#' Accept the priors after the prior predictive check
#'
#' Prior tuning is interactive: set the prior values, run `tune()`, look at the
#' prior predictive check, change the values, repeat. When the check looks
#' right, `accept_priors()` records that decision: the wflow (or the regimen's
#' wflow) is marked as tuned and the last `prior_log()` entry as accepted.
#' Fitting without accepted priors gives a warning. `set_as_tuned()` is the
#' earlier name of the same function.
#' @param wf a `wflow` or a `regimen`
#' @param ... unused
#' @return the object with accepted priors (a regimen also re-saves its prior file)
#' @export
accept_priors <- function(wf, ...) UseMethod("accept_priors")

#' @rdname accept_priors
#' @export
set_as_tuned <- function(wf, ...) accept_priors(wf, ...)

#' @rdname accept_priors
#' @exportS3Method
accept_priors.wflow <- function(wf, ...) {
  wf$tuned <- TRUE
  n <- length(wf$prior_log)
  if (n) wf$prior_log[[n]]$accepted <- TRUE   # the priors of the last prior check were accepted
  wf
}

# ---- prior tuning log -------------------------------------------------------

# One log entry per prior predictive run: prior values, compact PPC numbers and
# only the minimal plots (density overlay, mean-sd), per group when grouped.
#' @keywords internal
.prior_log_entry <- function(wf) {
  keep <- c("density_overlay", "mean_sd_2d")
  compact <- function(res) {
    if (inherits(res, "ppc_results_grouped"))
      return(list(summary = as.data.frame(res$summary),
                  plots   = lapply(res$group_results, function(g) g$plots[intersect(keep, names(g$plots))])))
    list(summary = as.data.frame(as.list(c(stats::setNames(res$bayes_pvalues, paste0("p_B_", names(res$bayes_pvalues))),
                                           stats::setNames(res$coverage, paste0("cov_", names(res$coverage))))),
                                 check.names = FALSE),
         plots   = res$plots[intersect(keep, names(res$plots))])
  }
  ppc <- lapply(wf$prior_ppc_results, compact)
  list(time     = Sys.time(),
       priors   = as.list(Filter(Negate(is.na), unclass(wf$prior_fit$priors))),
       accepted = FALSE,
       summary  = lapply(ppc, `[[`, "summary"),
       plots    = lapply(ppc, `[[`, "plots"))
}

#' Prior tuning log
#'
#' Each prior predictive run (`tune()`, `run_predictive_checks()`) records the
#' prior values used, the predictive-check numbers and the minimal plots
#' (density overlay and mean-SD). `set_as_tuned()` marks the last run as accepted.
#'
#' @param wf a `wflow` (or a `regimen`, whose `$wf` is used)
#' @param plots logical; also return the stored ggplot objects (default FALSE)
#' @return a data frame with one row per run, outcome and group: `run`, `time`,
#'   `accepted`, one column per prior (vector values pasted with ", "),
#'   `outcome` and the predictive-check columns. With `plots = TRUE`, a list
#'   with `table` and `plots` (`plots[[run]][[outcome]]`).
#' @export
prior_log <- function(wf, plots = FALSE) {
  if (inherits(wf, "regimen")) wf <- wf$wf
  log <- wf$prior_log %||% list()
  rows <- lapply(seq_along(log), function(i) {
    e  <- log[[i]]
    pv <- lapply(e$priors, function(v) if (length(v) > 1) paste(v, collapse = ", ") else v)
    head <- data.frame(run = i, time = e$time, accepted = e$accepted, pv, check.names = FALSE)
    dplyr::bind_rows(lapply(names(e$summary), function(o)
      cbind(head, outcome = o, e$summary[[o]], row.names = NULL)))
  })
  tab <- dplyr::bind_rows(rows)
  if (plots) list(table = tab, plots = lapply(log, `[[`, "plots")) else tab
}

# ---- fit_posterior ----------------------------------------------------------

#' Fit the posterior model
#'
#' @param wf  a `wflow` with data attached
#' @param sampler_params_list a `sampler_params` object, containing the parameters for the posterior sampling
#' @param ... additional arguments forwarded to `fit.stanmodel()`
#' @return the `wflow` with `$post_fit` populated
#' @export
fit_posterior <- function(wf, sampler_params_list, ...) UseMethod("fit_posterior")

#' @rdname fit_posterior
#' @exportS3Method
fit_posterior.wflow <- function(wf, sampler_params_list = NULL, ...) {
  .check_data(wf)

  message("Fitting posterior for '", wf$model_name, "'...")

  if (!wf$tuned)
    warning(.model_not_tunned_msg(wf$model_name))

  if (!is.null(sampler_params_list)) {
    wf$post_sampler_params <- modifyList(wf$post_sampler_params, sampler_params_list)
  }

  # Start from compiled (carries priors + sampler_params) rather than
  # prior_fit, so prior_predictive sampling doesn't inflate the chain count.
  wf$post_fit <- wf$compiled 
  wf$post_fit$sampler_params <- wf$post_sampler_params
  wf$post_fit <- wf$post_fit |>
    fit(
      data        = wf$data$data,
      extra       = wf$data$extra,
      cond_df     = wf$data$cond_df,
      outcome_var = wf$data$outcome_var,
      ...
    )
  wf
}

# ---- diagnose ---------------------------------------------------------------

#' Run convergence diagnostics on the posterior fit
#'
#' Wraps `diagnose.stanmodel()`. Results are stored in `$diagnostics`.
#'
#' @param wf   a `wflow` with `$post_fit` populated
#' @param params character vector of parameter names; NULL uses the params
#'   loaded from `diagnostics_params.csv` (or all non-yrep params if absent)
#' @param plot  logical; render traceplots (default FALSE)
#' @param ... additional arguments passed to methods
#' @return the `wflow` invisibly (diagnostics printed as side effect)
#' @export
run_diagnostics <- function(wf, ...) UseMethod("run_diagnostics")

#' @rdname run_diagnostics
#' @exportS3Method
run_diagnostics.wflow <- function(wf, params = NULL, plot = FALSE, ...) {
  .check_post_fit(wf)
  
  params_use <- params %||% wf$diagnostics_params
  wf$diagnostics <- wf$post_fit |> diagnose(params = params_use, plot = plot)
  print(wf$diagnostics)
  wf
}

# ---- read-only checks on a wflow --------------------------------------------

#' Predictive check and diagnostics of a wflow's fits (read-only)
#'
#' Compute a prior or posterior predictive check, or convergence diagnostics,
#' from the fits stored in a `wflow` (or a `regimen`) and return the result.
#' The object itself is not changed, so you can try other test statistics or
#' groupings without touching the stored results.
#' @param x a `wflow` or `regimen`
#' @param type `"posterior"` (default) or `"prior"`
#' @param ... passed to [ppc.stanmodel()] (e.g. `.by`, `stats`) or
#'   [diagnose.stanmodel()] (e.g. `params`, `plot`)
#' @return the `ppc()` / `diagnose()` result
#' @name wflow_checks
NULL

#' @rdname wflow_checks
#' @exportS3Method
ppc.wflow <- function(x, type = c("posterior", "prior"), ...) {
  ppc(.wflow_fit(x, match.arg(type)), ...)
}

#' @rdname wflow_checks
#' @param params parameters to check; default the wflow's `diagnostics_params`
#' @exportS3Method
diagnose.wflow <- function(x, params = NULL, ...) {
  diagnose(.wflow_fit(x, "posterior"), params = params %||% x$diagnostics_params, ...)
}

#' @keywords internal
.wflow_fit <- function(wf, type) {
  fit <- if (type == "prior") wf$prior_fit else wf$post_fit
  if (is.null(fit))
    stop(sprintf("No %s fit in this wflow: run %s first.", type,
                 if (type == "prior") "`tune()`" else "`fit_posterior()`"), call. = FALSE)
  fit
}

# inference from a fit that failed its convergence checks: repeat the warning
#' @keywords internal
.warn_failed_diagnostics <- function(wf) {
  d <- wf$diagnostics
  if (inherits(d, "stanworkflow_diagnostics") && !isTRUE(d$ok))
    warning(structure(class = c("stanworkflow_convergence_warning", "warning", "condition"),
      list(message = paste0("Inference from a fit that failed its convergence diagnostics (",
                            paste(d$problems, collapse = "; "), "). Interpret with caution."),
           call = NULL)))
  invisible(wf)
}

# ---- run_inference ------------------------------------------------------------------

#' Compute posterior contrasts and store results
#'
#' A thin wrapper around `infer_contrasts()` that stores output in
#' `$inference[[name]]`.
#'
#' @param wf       a `wflow` with `$post_fit` populated
#' @param contrasts a `contrast` object or named list of `contrast` objects
#' @param name     character; key under `$inference` to store the result
#'   (default `"main"`)
#' @param data     data frame to use for inference; NULL uses `$data$data`
#' @param extra      optional named list of additional Stan data passed to `get_stan_data()` (default NULL)
#' @param numeric_arrays optional character vector (or named list of column subsets) naming numeric array variables used by continuous contrasts (default NULL)
#' @param extra      optional named list of additional Stan data passed to `get_stan_data()` (default NULL)
#' @param numeric_arrays optional character vector (or named list of column subsets) naming numeric array variables used by continuous contrasts (default NULL)
#' @param ...      additional arguments forwarded to `infer_contrasts()`
#' @return the `wflow` with `$inference[[name]]` populated
#' @export
run_inference <- function(wf, ...) UseMethod("run_inference")

#' @rdname run_inference
#' @exportS3Method
run_inference.wflow <- function(wf,
                                 contrasts,
                                 name = "main",
                                 data = NULL,
                                 extra = NULL,
                                 numeric_arrays = NULL,
                                 ...) {

  .check_post_fit(wf)
  .warn_failed_diagnostics(wf)
  data_use <- data %||% wf$data$data
  extra_use <- extra %||% wf$data$extra
  wf$inference[[name]] <- infer_contrasts(
    model     = wf$post_fit,
    contrasts = contrasts,
    data      = data_use,
    extra     = extra_use,
    numeric_arrays = numeric_arrays,
    ...
  )
  wf
}

# ---- tune -------------------------------------------------------------------
#' Tune hyperpriors via prior predictive checks
#' 
#' This function runs prior predictive checks to help tune hyperpriors. It fits the model with `skip_likelihood = 1` and stores the result in `$prior_fit`. Use `ppc(wf$prior_fit, ...)` for visual checks.
#' @param wf  a `wflow` with data attached via `attach_data()` (or a `regimen`)
#' @param sampler_params a `sampler_params` object, containing the parameters for the prior predictive sampling
#' @param file optional file path to save the `wflow` object after tuning
#' @param ... additional arguments passed to `ppc.stanmodel()`. 
#' @return the `wflow` with `$prior_fit` populated
#' @export
tune <- function(wf, ...) UseMethod("tune")

#' @rdname tune
#' @exportS3Method
tune.wflow <- function(wf, sampler_params = NULL, file = NULL, ...) {
  .check_data(wf)
  wf <- wf |> run_predictive_checks(prior_predictive = TRUE, sampler_params = sampler_params, ...)
  if (!is.null(file)) {
    wf |> write(file)
  }
  wf
}

# ---- run ----------------------------------------------------------------

#' Run the full workflow in one call
#'
#' Convenience wrapper: prior predictive (optional) -> posterior fit ->
#' diagnostics -> inference.
#'
#' @param wf                   a `wflow` with data attached (or a `regimen`)
#' @param contrasts            contrast object(s) passed to `run_inference()`
#' @param sampler_params  optional list of sampler overrides for this run
#' @param posterior_predictive logical; run posterior predictive checks (default TRUE)
#' @param diagnostic_plot      logical; render traceplots (default FALSE)
#' @param diagnostic_params    list of parameters for diagnostics (default NULL)
#' @param file                 optional file path to save the `wflow` object after running
#' @param summary_fn           aggregation function passed to `infer_contrasts()` (default `mean`)
#' @param prob                 HPDI width passed to `infer_contrasts()` (default 0.95)
#' @param draws_method         `"descriptive"`, `"counterfactual"` (default) or `"parameter"`, passed to `infer_contrasts()`
#' @param gq_steps             grid resolution for counterfactual sweeps (default 5)
#' @param summary_fn           aggregation function passed to `infer_contrasts()` (default `mean`)
#' @param prob                 HPDI width passed to `infer_contrasts()` (default 0.95)
#' @param draws_method         `"descriptive"`, `"counterfactual"` (default) or `"parameter"`, passed to `infer_contrasts()`
#' @param gq_steps             grid resolution for counterfactual sweeps (default 5)
#' @param ...                  forwarded to `run_inference()`
#' @return the fully-populated `wflow`
#' @export
run <- function(wf, ...) UseMethod("run")

#' @rdname run
#' @exportS3Method
run.wflow <- function(wf,
                      contrasts,
                      sampler_params = NULL,
                      posterior_predictive = TRUE,
                      diagnostic_params = NULL,
                      diagnostic_plot   = FALSE,
                      summary_fn        = mean,
                      prob              = 0.95,
                      draws_method      = "counterfactual",
                      gq_steps          = 5,
                      file = NULL,
                           ...) {

  
  .check_data(wf)
  .check_compiled(wf)

  # tuning is not checked since it is checked in fit_posterior() and a warning is issued if not tuned
  # no point in double warnings
                             
  if (!is.null(sampler_params)) {
    wf$post_sampler_params <- modifyList(wf$post_sampler_params,
                                              sampler_params)
  }
  
  wf <- wf |> fit_posterior()

  if (posterior_predictive)
    wf <- wf |> run_predictive_checks(prior_predictive = FALSE, ...)
  
  wf <- wf |> run_diagnostics(plot = diagnostic_plot, params = diagnostic_params)

  if (!is.list(contrasts)) {
    contrasts <- list(main = contrasts)
  }

  if (draws_method == 'none'){
    message('Skipping inference step since draws_method is set to "none".')
  } else {
    for (co in names(contrasts)) {
       wf <- wf |> run_inference(
         contrasts = contrasts[[co]],
         name = co,
       draws_method = draws_method,
       prob = prob,
       summary_fn = summary_fn,
       gq_steps = gq_steps
     ) 
    }
  }
  if (!is.null(file)) {
    wf |> write(file)
  }

  wf
}

# ---- Write/ Read (fst-compressed RDS) ---------------------------------------

#' Write a wflow to disk using compressed serialization
#'
#' Saves a `wflow` object as an RDS file holding an `fst`-compressed
#' serialization of the object (the format of `JWileymisc::saveRDSfst()`, so
#' files written by earlier versions and by this version are interchangeable).
#'
#' @param x         a `wflow` object
#' @param path      file path (`.rds` or `.fst`); extension does not affect
#'   serialization format
#' @param dir_create logical; if `TRUE` (default) and directory does not exist,
#'   create it recursively
#' @param compress  integer; compression level (0-100). Higher = more compression
#'   but slower and more memory: ZSTD at 100 needs about 9 times the object's size
#'   in RAM, enough to kill R silently on a large posterior. Default `50`.
#' @param algorithm character; `"ZSTD"` (default) or `"LZ4"` (faster, less
#'   compression)
#' @param ... ignored (kept for backward compatibility)
#'
#' @details
#' The file is `saveRDS(fst::compress_fst(serialize(x, NULL)))`. If that
#' fails, the object is saved with plain `saveRDS()` (with a warning);
#' [load_wflow()] reads both.
#'
#' @return the `wflow` invisibly
#'
#' @seealso [load_wflow()]
#'
#' @examples
#' \dontrun{
#' wf <- wflow(...) |> fit_posterior()
#' wf |> write("./models/m1_posterior.rds")
#' wf |> write("./models/m1_posterior.rds", compress = 0, algorithm = "LZ4")
#' }
#'
#' @export
write <- function(x, path, dir_create = TRUE, ...) UseMethod("write")

#' @rdname write
#' @method write wflow
#' @export
write.wflow <- function(x,
                        path,
                        dir_create = TRUE,
                        compress   = 50,
                        algorithm  = "ZSTD",
                        ...) {
  stopifnot(
    inherits(x, "wflow"),
    is.character(path), length(path) == 1L,
    is.numeric(compress), compress >= 0L, compress <= 100L,
    is.character(algorithm), length(algorithm) == 1L
  )

  dir_path <- dirname(path)

  if (!dir.exists(dir_path)) {
    if (dir_create) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      message("Created directory: ", dir_path)
    } else {
      stop(sprintf(
        "Directory does not exist: '%s'. Set dir_create = TRUE to create it.",
        dir_path
      ))
    }
  }

  tryCatch(
    {
      .save_rds_fst(x, path, compression = compress, algorithm = algorithm)
      message(sprintf(
        "wflow saved to: %s (compress=%d, algorithm=%s)",
        path, compress, algorithm
      ))
    },
    error = function(e) {
      warning(sprintf(
        "Failed to save wflow: %s\nTrying fallback with base::saveRDS()",
        e$message
      ))
      saveRDS(x, path)
      message("Fallback: wflow saved using base::saveRDS() to: ", path)
    }
  )

  invisible(x)
}

#' Load a wflow from disk
#'
#' Reads a `wflow` saved with [write.wflow()]: either the `fst`-compressed
#' format (also written by `JWileymisc::saveRDSfst()`) or a plain RDS file.
#'
#' @param path file path to the wflow (`.rds` or `.fst`); format auto-detected
#'
#' @return a `wflow` object
#'
#' @seealso [write.wflow()]
#'
#' @examples
#' \dontrun{
#' wf <- load_wflow("./models/m1_posterior.rds")
#' print(wf)
#' }
#'
#' @export
load_wflow <- function(path) {
  stopifnot(is.character(path), length(path) == 1L)

  if (!file.exists(path))
    stop(sprintf("File not found: '%s'", path))

  wf <- .read_rds_fst(path)

  stopifnot(
    inherits(wf, "wflow"),
    !is.null(wf$model_name),
    !is.null(wf$compiled)
  )

  message("wflow loaded from: ", path)
  wf
}

# Same bytes as JWileymisc::saveRDSfst() / readRDSfst().
#' @keywords internal
.save_rds_fst <- function(x, path, compression = 50, algorithm = "ZSTD") {
  saveRDS(fst::compress_fst(serialize(x, NULL), compressor = algorithm,
                            compression = compression, hash = TRUE),
          file = path, version = 3)
}

#' @keywords internal
.read_rds_fst <- function(path) {
  x <- readRDS(path)
  if (is.raw(x)) unserialize(fst::decompress_fst(x)) else x
}

# ---- helpers ----------------------------------------------------------------

#' @keywords internal
.check_data <- function(wf) {
  if (is.null(wf$data))
    stop("No data attached. Call `attach_data()` first.")
}

#' @keywords internal
.check_compiled <- function(wf) {
  if (is.null(wf$compiled))
    stop("Model not compiled. Call `compile()` first.")
}


#' @keywords internal
.check_post_fit <- function(wf) {
  if (is.null(wf$post_fit) || !inherits(wf$post_fit, "stanmodel"))
    stop("Model not fitted or not a stanmodel. Call `fit_posterior()` (or `run_all()`) first.")
}



#' @keywords internal
.model_not_tunned_msg <- function(model_name) {
  sprintf(
    "Model '%s' has not been marked as tuned: its priors were not accepted. Run `tune()`, look at the prior predictive check, then call `accept_priors()`.\n",
    model_name
  )
}
