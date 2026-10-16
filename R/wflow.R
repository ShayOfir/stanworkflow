# wflow.R — Bayesian modelling workflow S3 class
# Wraps the compile -> prior_predictive -> fit -> diagnose -> infer lifecycle.

# ---- Constructor -----------------------------------------------------------

#' Create a new wflow object
#'
#' A `wflow` encapsulates a complete Bayesian modelling workflow:
#' model code, priors, sampler settings, fitted models, data, and inference.
#'
#' @param model_name   character; name of a subfolder under `models_folder`
#'   containing `model_code.stan`, `priors.csv`, and optionally
#'   `diagnostics_params.csv`.
#' @param models_folder character; root folder for model definitions
#'   (default `"./models/"`).
#' @param prior_sampler_params_list named list of overrides for `sampler_params()`.
#' @param post_sampler_params_list named list of overrides for `sampler_params()`.
#'   NULL uses package defaults.
#' @return an object of class `wflow`
#' @export
wflow <- function(model_name,
                  models_folder      = "./models/",
                  prior_sampler_params_list = NULL,
                  post_sampler_params_list = NULL) {

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
    read.csv(priors_file) |>
      tibble::deframe() |>
      as.list()
  } else {
    message(sprintf("No priors.csv found in '%s'; priors will need to be set manually.", model_dir))
    list()
  }

  # ---- load diagnostics params ---------------------------------------------
  diag_file <- file.path(model_dir, "diagnostics_params.csv")
  diagnostics_params <- if (file.exists(diag_file)) {
    read.csv(diag_file) |> dplyr::pull(parameters)
  } else {
    NULL
  }

  # ---- compile -------------------------------------------------------------

  compiled <- stanmodel() |> compile_from_file(stan_file)
  compiled$priors <- modifyList(compiled$priors, priors_list)

  # ---- sampler params ------------------------------------------------------
  prior_sp <- if (!is.null(prior_sampler_params_list)) {
    modifyList(sampler_params(), prior_sampler_params_list)
  } else {
    sampler_params()
  }

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
      tuned             = FALSE             # whether hyperpriors have been tuned
    ),
    class = "wflow"
  )
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
#' @return the `wflow` with `$data` populated
#' @export
attach_data <- function(wf, ...) UseMethod("attach_data")

#' @exportS3Method
attach_data.wflow <- function(wf,
                               data,
                               extra       = list(),
                               cond_df     = NULL,
                               outcome_var = "y") {
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
#' @param ppc logical; run `ppc()` on the prior fit (default FALSE)
#' @param ... additional arguments passed to `ppc.stanmodel()`. Ignored if `ppc = FALSE`.
#' @return the `wflow` with `$prior_fit` populated
#' @export
run_predictive_checks <- function(wf,...) UseMethod("run_predictive_checks")

#' @exportS3Method
run_predictive_checks.wflow <- function(wf, prior_predictive = TRUE, sampler_params = NULL,  ...) {

  if (prior_predictive){
      .check_data(wf)
      .check_compiled(wf)
    
      wf$prior_fit <- wf$compiled
      if (inherits(sampler_params,'sampler_params'))
        {
        wf$prior_fit$sampler_params <- wf$prior_fit$sampler_params |> modifyList(sampler_params)
        }
      else {
        warning('Prior sampler_params is not a `sampler_params` class.\n')
      }

      message("Running prior predictive checks for '", wf$model_name, "'...")
      wf$prior_fit <- wf$prior_fit |>
        fit(
          data             = wf$data$data,
          extra            = wf$data$extra,
          cond_df          = wf$data$cond_df,
          outcome_var      = wf$data$outcome_var,
          prior_predictive = TRUE,
        )


      wf$prior_ppc_results <- wf$prior_fit |> ppc(...)
      if (inherits(wf$prior_ppc_results, "ppc_results_grouped") || inherits(wf$prior_ppc_results, "ppc_results")) {
        wf$prior_ppc_results <- list(wf$prior_ppc_results)
        names(wf$prior_ppc_results) <- wf$data$outcome_var
      }
      for (i in seq_along(wf$prior_ppc_results)) {
        cat(sprintf("Prior predictive check for outcome '%s':\n", names(wf$prior_ppc_results)[i]))
        if (inherits(wf$prior_ppc_results[[i]], "ppc_results_grouped")) {
          wf$prior_ppc_results[[i]] |> display()
        } else {
          wf$prior_ppc_results[[i]] |> print(plot = TRUE)
        }
      }
  } else {
      message("Performing posterior predictive checks '", wf$model_name, "'.")
      if (is.null(wf$post_fit) || !inherits(wf$post_fit, "stanmodel")) {
        stop("Posterior fit not found or not as `stanmodel`. Please run `fit_posterior()` first.")
      }
      wf$posterior_ppc_results <- wf$post_fit |> ppc(...)
      if (inherits(wf$posterior_ppc_results, "ppc_results_grouped") || inherits(wf$posterior_ppc_results, "ppc_results")) {
        wf$posterior_ppc_results <- list(wf$posterior_ppc_results)
        names(wf$posterior_ppc_results) <- wf$data$outcome_var
      }
      for (i in seq_along(wf$posterior_ppc_results)) {
        cat(sprintf("Posterior predictive check for outcome '%s':\n", names(wf$posterior_ppc_results)[i]))
        if (inherits(wf$posterior_ppc_results[[i]], "ppc_results_grouped")) {
          wf$posterior_ppc_results[[i]] |> display()
        } else {
          wf$posterior_ppc_results[[i]] |> print(plot = TRUE)
        }
      }
  }

  
  wf
  }

# --- set_as_tuned -------------------------------------------------------------
#' Mark the wflow as having tuned hyperpriors
#' @param wf a `wflow` object
#' @return the `wflow` with `$tuned` set to TRUE
#' @export
set_as_tuned <- function(wf) UseMethod("set_as_tuned")

#' @exportS3Method
set_as_tuned.wflow <- function(wf) {
  wf$tuned <- TRUE
  wf
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
#' @return the `wflow` invisibly (diagnostics printed as side effect)
#' @export
run_diagnostics <- function(wf, ...) UseMethod("run_diagnostics")

#' @exportS3Method
run_diagnostics.wflow <- function(wf, params = NULL, plot = FALSE) {
  .check_post_fit(wf)
  
  params_use <- params %||% wf$diagnostics_params
  wf$diagnostics <- wf$post_fit |> diagnose(params = params_use, plot = plot)
  wf
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
#' @param ...      additional arguments forwarded to `infer_contrasts()`
#' @return the `wflow` with `$inference[[name]]` populated
#' @export
run_inference <- function(wf, ...) UseMethod("run_inference")

#' @exportS3Method
run_inference.wflow <- function(wf,
                                 contrasts,
                                 name = "main",
                                 data = NULL,
                                 extra = NULL,
                                 numeric_arrays = NULL,
                                 ...) {

  .check_post_fit(wf)
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
#' @param wf  a `wflow` with data attached via `attach_data()`
#' @param sampler_params a `sampler_params` object, containing the parameters for the prior predictive sampling
#' @param file optional file path to save the `wflow` object after tuning
#' @param ... additional arguments passed to `ppc.stanmodel()`. 
#' @return the `wflow` with `$prior_fit` populated
#' @export
tune <- function(wf, sampler_params = NULL, file = NULL, ...) UseMethod("tune")

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
#' @param wf                   a `wflow` with data attached
#' @param contrasts            contrast object(s) passed to `run_inference()`
#' @param sampler_params  optional list of sampler overrides for this run
#' @param posterior_predictive logical; run posterior predictive checks (default TRUE)
#' @param diagnostic_plot      logical; render traceplots (default FALSE)
#' @param diagnostic_params    list of parameters for diagnostics (default NULL)
#' @param file                 optional file path to save the `wflow` object after running
#' @param ...                  forwarded to `run_inference()`
#' @return the fully-populated `wflow`
#' @export
run <- function(wf, ...) UseMethod("run")

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

# ---- Write/ Read (MODIFIED with saveRDSfst) --------------------------------

#' Write a wflow to disk using compressed serialization
#'
#' Saves a `wflow` object to disk using `JWileymisc::saveRDSfst()` for
#' efficient compression. Significantly faster and more disk-efficient than
#' standard `saveRDS()`.
#'
#' @param x         a `wflow` object
#' @param path      file path (`.rds` or `.fst`); extension does not affect
#'   serialization format
#' @param dir_create logical; if `TRUE` (default) and directory does not exist,
#'   create it recursively
#' @param compress  integer; compression level (0-100). Higher = more compression
#'   but slower. Default `100` (maximum compression, recommended for wflow objects
#'   with large fitted models)
#' @param algorithm character; compression algorithm passed to `saveRDSfst()`.
#'   Options: `"ZSTD"` (default, fast), `"LZ4"` (very fast, less compression),
#'   `"ZSTD_SMALL"` (slow, smallest file)
#' @param ... additional arguments passed to `JWileymisc::saveRDSfst()`
#'   (e.g., `data_frame`, `version`, `old_format`)
#'
#' @details
#' Uses `JWileymisc::saveRDSfst()` instead of base `saveRDS()`:
#' - **Speed**: ~3-10× faster write/read for large objects
#' - **Disk**: ~30-50% smaller file sizes (especially for fitted Stan models)
#' - **Compatibility**: RDS format preserved; loads with `load_wflow()`
#'
#' Typical compression settings:
#' - `compress = 100, algorithm = "ZSTD"`: maximum compression, ~1-2s extra time
#' - `compress = 50, algorithm = "ZSTD"`: balanced (recommended default)
#' - `compress = 0, algorithm = "LZ4"`: fastest, still ~50% smaller than uncompressed
#'
#' @return the `wflow` invisibly
#'
#' @seealso [load_wflow()], [JWileymisc::saveRDSfst()]
#'
#' @examples
#' \dontrun{
#' wf <- wflow(...) |> fit_posterior()
#'
#' # Write with default maximum compression
#' wf |> write("./models/m1_posterior.rds")
#'
#' # Write with fast compression (smaller file)
#' wf |> write("./models/m1_posterior.rds", compress = 50)
#'
#' # Write uncompressed (fastest)
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
                        compress   = 100,
                        algorithm  = "ZSTD",
                        ...) {
  stopifnot(
    inherits(x, "wflow"),
    is.character(path), length(path) == 1L,
    is.numeric(compress), compress >= 0L, compress <= 100L,
    is.character(algorithm), length(algorithm) == 1L
  )

  dir_path <- dirname(path)

  # ---- MODIFIED: Create directory if needed ------
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

  # ---- MODIFIED: Use saveRDSfst with compression options ------
  tryCatch(
    {
      JWileymisc::saveRDSfst(
        x,
        file      = path,
        compress  = compress,
        algorithm = algorithm,
        ...
      )
      message(sprintf(
        "wflow saved to: %s (compress=%d, algorithm=%s)",
        path, compress, algorithm
      ))
    },
    error = function(e) {
      stop(sprintf(
        "Failed to save wflow: %s\nTrying fallback with base::saveRDS()",
        e$message
      ))
      # Fallback to base saveRDS if saveRDSfst fails
      saveRDS(x, path)
      message("Fallback: wflow saved using base::saveRDS() to: ", path)
    }
  )

  invisible(x)
}

#' Load a wflow from disk
#'
#' Deserializes a `wflow` object saved with [write.wflow()] using
#' `JWileymisc::readRDSfst()`. Handles both compressed (fst format) and
#' standard RDS files transparently.
#'
#' @param path file path to the wflow (`.rds` or `.fst`); format auto-detected
#'
#' @return a `wflow` object
#'
#' @details
#' Uses `JWileymisc::readRDSfst()` which:
#' - Auto-detects compression format (fst or RDS)
#' - ~3-10× faster than base `readRDS()` for compressed files
#' - Falls back gracefully to base `readRDS()` for standard `.rds` files
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

  # ---- MODIFIED: Try loadRDSfst first, fallback to readRDS ------
  wf <- tryCatch(
    {
      JWileymisc::readRDSfst(path)
    },
    error = function(e) {
      warning(sprintf(
        "readRDSfst() failed (%s); attempting fallback with base::readRDS()",
        e$message
      ))
      readRDS(path)
    }
  )

  stopifnot(
    inherits(wf, "wflow"),
    !is.null(wf$model_name),
    !is.null(wf$compiled)
  )

  message("wflow loaded from: ", path)
  wf
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
    "Model '%s' has not been marked as tuned. Consider running `tune()` first.\n When you are satisfied with the hyperpriors, call `set_as_tuned()` to mark the wflow as tuned.\n",
    model_name
  )
}