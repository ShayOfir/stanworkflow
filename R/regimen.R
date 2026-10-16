# regimen.R — OOP wrapper for the full Bayesian modelling regimen

# ---- Constructor ------------------------------------------------------------

#' Create a new regimen object
#'
#' Instantiates a complete Bayesian modelling regimen that orchestrates the
#' prepare → tune → fit → infer workflow. A regimen encapsulates data,
#' model configuration, contrasts, sampler parameters, an underlying wflow,
#' and inference results with explicit status tracking.
#'
#' @param model_name character scalar; name of a subfolder under `models_folder`
#'   containing `model_code.stan` and `priors.csv`.
#' @param data_list list with required element `data` (data frame of observations)
#'   and optional elements `extra` (named list of additional Stan data like
#'   matrices or scalars) and `cond_df` (optional condition-level data frame
#'   for `*_cond_array` variables).
#' @param contrasts named list of `contrast` objects specifying the comparisons
#'   to compute during inference. All elements must inherit from class `contrast`.
#' @param outcome_var character scalar; name of the outcome column in
#'   `data_list$data` (default `"y"`).
#' @param models_folder character scalar; root directory containing model
#'   subfolders (default `"./models/"`).
#' @param results_path character scalar; directory where prior, posterior, and
#'   inference results are saved (default `"./results/"`). Created if missing.
#' @param ppc_by character vector; column names for grouping variables in
#'   posterior predictive checks (default `NULL`, no grouping).
#' @param ppc_outcomes character vector or list of `outcome_spec` objects
#'   specifying which quantities to extract for PPC (default `NULL`).
#' @param prior_sampler_params named list of sampler parameter overrides for
#'   prior predictive sampling. Merged with package defaults via `modifyList()`.
#'   Standard keys: `iter_warmup`, `iter_sampling`, `refresh`, `adapt_delta`,
#'   `max_treedepth`, `thin`, `seed`, `parallel_chains`, `threads_per_chain`.
#' @param post_sampler_params named list of sampler parameter overrides for
#'   posterior sampling. Merged with package defaults via `modifyList()`.
#' @param draws_method character scalar; how to generate draws for contrasts.
#'   `"counterfactual"` (default) sweeps the contrast variable across a grid;
#'   `"descriptive"` computes contrasts directly from posterior draws.
#' @param gq_steps integer scalar; number of grid points for counterfactual
#'   sweeps (default `5L`). Ignored if `draws_method = "descriptive"`.
#' @param prob numeric scalar ∈ (0,1); width of highest posterior density
#'   intervals (HPDI) for credible intervals in summaries (default `0.95`).
#'
#' @return An object of class `regimen`, a list containing:
#'   - Identity: `model_name`, `models_folder`, `results_path`, `paths`
#'   - Data: `data_list`, `outcome_var`
#'   - Inference config: `contrasts`, `draws_method`, `gq_steps`, `prob`
#'   - PPC config: `ppc_by`, `ppc_outcomes`
#'   - Sampler params: `prior_sampler_params`, `post_sampler_params`
#'   - Runtime state: `wf` (wflow object), `inference` (results list)
#'   - Status: `status` (list with `prepared`, `tuned`, `fitted`, `inferred` flags)
#'
#' @seealso [run()], [prepare()], [tune()], [fit()], [infer()]
#'
#' @examples
#' \dontrun{
#' reg <- regimen(
#'   model_name = "m1",
#'   data_list = list(
#'     data = my_data,
#'     extra = list(N_roi = 10),
#'     cond_df = NULL
#'   ),
#'   contrasts = list(
#'     group_effect = contrast(contrast_var = "group", ...),
#'     interaction  = contrast(contrast_var = "group", two_way = "direction", ...)
#'   ),
#'   results_path = "./results_m1/"
#' )
#'
#' # Run full pipeline
#' reg <- run(reg)
#' print(reg$inference$group_effect)
#' }
#'
#' @export
regimen <- function(model_name,
                    data_list,
                    contrasts            = NULL,
                    outcome_var          = "y",
                    models_folder        = "./models/",
                    results_path         = "./results/",
                    ppc_by               = NULL,
                    ppc_outcomes         = NULL,
                    prior_sampler_params = list(iter_warmup   = 500L,
                                                iter_sampling = 500L,
                                                refresh       = 500L),
                    post_sampler_params  = list(iter_warmup   = 1000L,
                                                iter_sampling = 1000L,
                                                adapt_delta   = 0.99,
                                                max_treedepth = 15L,
                                                refresh       = 500L),
                    draws_method         = "counterfactual",
                    gq_steps             = 5L,
                    prob                 = 0.95) {

  stopifnot(is.character(model_name), length(model_name) == 1L)
  stopifnot(is.list(data_list), "data" %in% names(data_list))
  if (!is.null(contrasts)) stopifnot(is.list(contrasts), length(contrasts) >= 1L)
  if (!all(vapply(contrasts, inherits, logical(1L), "contrast")))
    stop("All elements of `contrasts` must be `contrast` objects.")
  if (is.null(names(contrasts)) || any(names(contrasts) == ""))
    stop("`contrasts` must be a fully named list.")

  draws_method <- match.arg(draws_method, c("counterfactual", "descriptive"))

  prior_sp <- modifyList(sampler_params(), prior_sampler_params)
  post_sp  <- modifyList(sampler_params(), post_sampler_params)

  .make_paths <- function(rp, mn) {
    list(
      prior     = file.path(rp, paste0(mn, "_prior.rds")),
      posterior = file.path(rp, paste0(mn, "_post.rds")),
      inference = file.path(rp, paste0(mn, "_inference.rds"))
    )
  }

  structure(
    list(
      model_name    = model_name,
      models_folder = models_folder,
      results_path  = results_path,
      paths         = .make_paths(results_path, model_name),
      data_list     = data_list,
      outcome_var   = outcome_var,
      contrasts     = contrasts,
      draws_method  = draws_method,
      gq_steps      = gq_steps,
      prob          = prob,
      ppc_by        = ppc_by,
      ppc_outcomes  = ppc_outcomes,
      prior_sampler_params = prior_sp,
      post_sampler_params  = post_sp,
      wf        = NULL,
      inference = list(),
      status = list(
        prepared = FALSE,
        tuned    = FALSE,
        fitted   = FALSE,
        inferred = FALSE
      )
    ),
    class = "regimen"
  )
}

# ---- print ------------------------------------------------------------------

#' Print a regimen summary
#'
#' Displays model identity, configuration, and current status of the regimen
#' pipeline.
#'
#' @param x a `regimen` object
#' @param ... additional arguments (unused)
#'
#' @return the `regimen` invisibly
#'
#' @export
print.regimen <- function(x, ...) {
  cat("<regimen>\n")
  cat("  model         :", x$model_name, "\n")
  cat("  models_folder :", x$models_folder, "\n")
  cat("  results_path  :", x$results_path, "\n")
  cat("  outcome_var   :", x$outcome_var, "\n")
  cat("  contrasts     :", paste(names(x$contrasts), collapse = ", "), "\n")
  cat("  draws_method  :", x$draws_method, "\n")
  cat("  gq_steps      :", x$gq_steps, "\n")
  cat("  prob          :", x$prob, "\n")
  cat("\n  Status:\n")
  cat("    prepared :", x$status$prepared, "\n")
  cat("    tuned    :", x$status$tuned, "\n")
  cat("    fitted   :", x$status$fitted, "\n")
  cat("    inferred :", x$status$inferred, "\n")
  if (length(x$inference) > 0L)
    cat("  inference keys:", paste(names(x$inference), collapse = ", "), "\n")
  invisible(x)
}

# ---- prepare ----------------------------------------------------------------

#' Prepare a regimen for analysis
#'
#' Compiles the Stan model and attaches data to the underlying `wflow` object.
#' This step is required before tuning, fitting, or inference.
#'
#' @param x a `regimen` object
#' @param ... additional arguments passed to [wflow()]
#'
#' @return the `regimen` with `$wf` populated and `$status$prepared = TRUE`
#'
#' @seealso [run()], [tune()], [fit()], [infer()]
#'
#' @export
prepare <- function(x, ...) UseMethod("prepare")

#' @rdname prepare
#' @export
prepare.regimen <- function(x, ...) {
  if (!dir.exists(x$results_path))
    dir.create(x$results_path, recursive = TRUE, showWarnings = FALSE)

  x$wf <- wflow(
    model_name                = x$model_name,
    models_folder             = x$models_folder,
    prior_sampler_params_list = as.list(x$prior_sampler_params),
    post_sampler_params_list  = as.list(x$post_sampler_params),
    ...
  ) |>
    attach_data(
      data        = x$data_list$data,
      extra       = x$data_list$extra   %||% list(),
      cond_df     = x$data_list$cond_df %||% NULL,
      outcome_var = x$outcome_var
    )

  x$status$prepared <- TRUE
  x
}

# ---- tune -------------------------------------------------------------------

#' Run prior predictive checks for a regimen
#'
#' Samples from the prior predictive distribution (likelihood skipped) and
#' stores the resulting `wflow` to disk at `$paths$prior`. Useful for
#' prior sensitivity checks and model debugging before fitting.
#'
#' @param x a `regimen` with `$status$prepared = TRUE`
#' @param save logical; if `TRUE` (default), saves the wflow to `$paths$prior`
#' @param ... additional arguments passed to [tune()]
#'
#' @return the `regimen` with `$wf$prior_fit` populated and `$status$tuned = TRUE`
#'
#' @seealso [run()], [fit()], [prepare()]
#'
#' @export
tune <- function(x, ...) UseMethod("tune")

#' @rdname tune
#' @export
tune.regimen <- function(x, save = TRUE, ...) {
  .check_prepared_regimen(x)

  x$wf <- x$wf |>
    tune(
      sampler_params = x$prior_sampler_params,
      outcomes       = x$ppc_outcomes,
      .by            = x$ppc_by,
      file           = if (save) x$paths$prior else NULL,
      ...
    )

  x$status$tuned <- TRUE
  x
}

# ---- fit --------------------------------------------------------------------

#' Fit the posterior for a regimen
#'
#' Runs posterior MCMC sampling. Optionally runs posterior predictive checks
#' and diagnostics. Saves the fitted `wflow` to `$paths$posterior`.
#'
#' @param x a `regimen` object (will be prepared automatically if needed)
#' @param posterior_predictive logical; if `TRUE` (default), run posterior
#'   predictive checks after fitting
#' @param diagnostic_plot logical; if `TRUE`, render traceplots and other
#'   diagnostic visualizations (default `FALSE`)
#' @param diagnostic_params character vector; names of parameters to include
#'   in diagnostics (default `NULL` = all scalar parameters)
#' @param save logical; if `TRUE` (default), save fitted wflow to
#'   `$paths$posterior`
#' @param load_prior logical; if `TRUE` (default) and prior wflow exists on
#'   disk, load it instead of re-preparing
#' @param ... additional arguments passed to [wflow()]
#'
#' @return the `regimen` with `$wf$post_fit` populated and `$status$fitted = TRUE`
#'
#' @seealso [run()], [prepare()], [tune()], [infer()]
#'
#' @export
fit <- function(x, ...) UseMethod("fit")

#' @rdname fit
#' @export
fit.regimen <- function(x,
                        posterior_predictive = TRUE,
                        diagnostic_plot      = FALSE,
                        diagnostic_params    = NULL,
                        save                 = TRUE,
                        load_prior           = TRUE,
                        ...) {
  if (load_prior && file.exists(x$paths$prior)) {
    message("Loading prior wflow from: ", x$paths$prior)
    x$wf <- load_wflow(x$paths$prior)
  } else {
    .check_prepared_regimen(x)
  }

  x$wf <- x$wf |>
    run(
      contrasts            = NULL,
      sampler_params       = as.list(x$post_sampler_params),
      posterior_predictive = posterior_predictive,
      diagnostic_params    = diagnostic_params,
      diagnostic_plot      = diagnostic_plot,
      outcomes             = x$ppc_outcomes,
      .by                  = x$ppc_by,
      draws_method         = "none",
      file                 = if (save) x$paths$posterior else NULL,
      ...
    )

  x$status$fitted   <- TRUE
  x$status$prepared <- TRUE
  x
}

# ---- infer ------------------------------------------------------------------

#' Run inference for all contrasts in a regimen
#'
#' Iterates over the named contrasts in `$contrasts`, calls [run_inference()]
#' for each, and collects results into `$inference`. Results are saved to
#' `$paths$inference` and the flag `$status$inferred` is set to `TRUE`.
#'
#' @param x a `regimen` with `$status$fitted = TRUE`
#' @param contrasts optional named list of `contrast` objects to override
#'   `$contrasts` (default `NULL` uses the regimen's contrasts)
#' @param save logical; if `TRUE` (default), save inference list to
#'   `$paths$inference`
#' @param load_post logical; if `TRUE` (default) and posterior wflow exists on
#'   disk, load it instead of using `$wf`
#' @param ... additional arguments passed to [run_inference()]
#'
#' @return the `regimen` with `$inference` populated and `$status$inferred = TRUE`
#'
#' @seealso [run()], [fit()], [prepare()], [get_inference()]
#'
#' @export
infer <- function(x, ...) UseMethod("infer")

#' @rdname infer
#' @export
infer.regimen <- function(x,
                          contrasts = NULL,
                          save      = TRUE,
                          load_post = TRUE,
                          ...) {
  contrasts_use <- contrasts %||% x$contrasts

  if (load_post && file.exists(x$paths$posterior)) {
    message("Loading posterior wflow from: ", x$paths$posterior)
    x$wf <- load_wflow(x$paths$posterior)
  } else {
    .check_fitted_regimen(x)
  }

  for (co_name in names(contrasts_use)) {
    message(sprintf("Running inference for contrast '%s'...", co_name))
    x$wf <- x$wf |>
      run_inference(
        contrasts    = contrasts_use[[co_name]],
        name         = co_name,
        draws_method = x$draws_method,
        gq_steps     = x$gq_steps,
        prob         = x$prob,
        ...
      )
    x$inference[[co_name]] <- x$wf$inference[[co_name]]
  }

  x$status$inferred <- TRUE

  if (save) {
    saveRDS(x$inference, x$paths$inference)
    message("Inference saved to: ", x$paths$inference)
  }

  x
}

# ---- quick_fit (NEW) ---------------------------------------------------------

#' Fit posterior quickly without checks (for mass-univariate workflows)
#'
#' Skips prior and posterior predictive checks and diagnostics to save time
#' and disk space. Only runs: prepare → fit. Useful for mass-univariate
#' applications where many models are fitted across ROIs/conditions.
#'
#' @param x a `regimen` object (will be prepared automatically if needed)
#' @param save logical; if `TRUE` (default), save fitted wflow to
#'   `$paths$posterior`
#' @param load_prior logical; if `TRUE` (default) and prior wflow exists on
#'   disk, load it instead of re-preparing. Note: prior predictive checks
#'   are skipped regardless (default `FALSE` since prior is typically not
#'   saved in quick_fit mode)
#' @param ... additional arguments passed to [wflow()]
#'
#' @details
#' Compared to `fit()`:
#' - Skips prior predictive checks (no `$prior_fit` saved)
#' - Skips posterior predictive checks
#' - Skips convergence diagnostics
#' - Does NOT run inference (`$inference` remains empty)
#'
#' This reduces computation time significantly and saves disk space. Use when:
#' - Fitting many univariate models (e.g., one per ROI)
#' - Prior checks already validated in a pilot study
#' - Focus is on point estimates rather than visual diagnostics
#'
#' After `quick_fit()`, call `infer()` separately if needed.
#'
#' @return the `regimen` with `$wf$post_fit` populated and `$status$fitted = TRUE`.
#'   Status flags `$status$tuned` and `posterior_predictive` results remain
#'   `FALSE` / `NULL`.
#'
#' @seealso [fit()], [run()], [infer()], [prepare()]
#'
#' @examples
#' \dontrun{
#' # Fit many models quickly
#' for (roi in rois) {
#'   reg <- regimen(..., contrasts = list(...)) |>
#'     quick_fit(save = TRUE)
#' }
#' }
#'
#' @export
quick_fit <- function(x, ...) UseMethod("quick_fit")

#' @rdname quick_fit
#' @export
quick_fit.regimen <- function(x,
                              save       = TRUE,
                              load_prior = FALSE,
                              ...) {
  # Prepare if not already done
  if (load_prior && file.exists(x$paths$prior)) {
    message("Loading prior wflow from: ", x$paths$prior)
    x$wf <- load_wflow(x$paths$prior)
  } else {
    .check_prepared_regimen(x)
    if (!x$status$prepared) {
      x <- x |> prepare(...)
    }
  }

  # ---- SKIP: prior predictive checks ------
  # (no tune() call)

  # ---- RUN: posterior fit only (no checks, no diagnostics) ------
  message("--- Quick fit (posterior only, skipping checks) ---")
  
  x$wf <- x$wf |>
    fit_posterior(
      sampler_params_list = as.list(x$post_sampler_params),
      ...
    )

  # ---- Mark as tuned to avoid warning in fit_posterior
  x$wf <- x$wf |> set_as_tuned()

  # ---- SKIP: posterior predictive checks ------
  # (no run_predictive_checks() for posterior)

  # ---- SKIP: diagnostics ------
  # (no run_diagnostics() call)

  # ---- Save posterior if requested ------
  if (save) {
    x$wf |> write(x$paths$posterior)
  }

  # ---- Update status ------
  x$status$fitted   <- TRUE
  x$status$prepared <- TRUE
  # NOTE: $status$tuned remains FALSE (prior checks were skipped)

  x
}

# ---- run.regimen (MODIFIED with quick_fit option) ----------------------------

#' Execute the full regimen pipeline (with quick_fit option)
#'
#' Orchestrates the complete workflow: prepare → tune → fit → infer.
#' Individual steps can be skipped via the `steps` argument. Each step is
#' idempotent if the corresponding wflow or inference file exists on disk.
#'
#' Alternatively, use `steps = "quick_fit"` to skip prior/posterior checks
#' and diagnostics (fast, disk-efficient for mass-univariate workflows).
#'
#' @param x a `regimen` object
#' @param steps character vector; subset of:
#'   - `c("prepare", "tune", "fit", "infer")` (default: all four)
#'   - `"quick_fit"` (special: prepare + fit, skip all checks)
#'   - Any custom subset of the above
#'   Controls which pipeline stages to execute.
#' @param save logical; if `TRUE` (default), intermediate results are saved to disk
#' @param ... additional arguments forwarded to individual step functions
#'   ([prepare()], [tune()], [fit()], [quick_fit()], [infer()])
#'
#' @return the populated `regimen` with all requested steps completed
#'
#' @seealso [prepare()], [tune()], [fit()], [quick_fit()], [infer()], [resume()]
#'
#' @examples
#' \dontrun{
#' # Run all steps normally
#' reg <- run(reg)
#'
#' # Run only prepare and fit, skip tune
#' reg <- run(reg, steps = c("prepare", "fit"))
#'
#' # Quick fit (no checks) - ideal for mass-univariate
#' reg <- run(reg, steps = "quick_fit")
#'
#' # Quick fit then infer
#' reg <- run(reg, steps = c("quick_fit", "infer"))
#'
#' # Run only inference with new contrasts
#' reg <- update_contrasts(reg, new_contrast_list) |>
#'   run(steps = "infer")
#' }
#'
#' @export
run <- function(x, ...) UseMethod("run")

#' @rdname run
#' @export
run.regimen <- function(x,
                        steps = c("prepare", "tune", "fit", "quick_fit", "infer"),
                        save  = TRUE,
                        ...) {
  

  valid_steps <- c("prepare", "tune", "fit", "quick_fit", "infer")
  unknown <- setdiff(steps, valid_steps)
  if (length(unknown) > 0L)
    stop(sprintf("Unknown step(s): %s. Valid: %s",
                 paste(unknown, collapse = ", "),
                 paste(valid_steps, collapse = ", ")))

  # ---- Check for conflicting steps ------
  # MODIFIED: Prevent simultaneous fit and quick_fit
  if ("fit" %in% steps && "quick_fit" %in% steps)
    stop("Cannot run both 'fit' and 'quick_fit' in same pipeline. Choose one.")

  if ("prepare" %in% steps) {
    message("--- Step: prepare (compile + attach data) ---")
    x <- x |> prepare(...)
  }

  if ("tune" %in% steps) {
    message("--- Step: tune (prior predictive) ---")
    x <- x |> tune(save = save, ...)
  }

  # ---- MODIFIED: Add quick_fit branch ------
  if ("quick_fit" %in% steps) {
    message("--- Step: quick_fit (posterior only, skipping checks) ---")
    x <- x |> quick_fit(save = save, ...)
  }

  if ("fit" %in% steps) {
    message("--- Step: fit (posterior + checks) ---")
    x <- x |> fit(save = save, ...)
  }

  if ("infer" %in% steps) {
    message("--- Step: infer ---")
    x <- x |> infer(save = save, ...)
  }

  x
}

#' Resume a regimen from the last saved checkpoint
#'
#' Inspects which checkpoint files exist on disk (`$paths`) and re-runs only
#' the necessary steps to reach completion. Useful for continuing after crashes,
#' network failures, or iterative development.
#'
#' @param x a `regimen` object (status flags and paths are checked)
#' @param ... additional arguments passed to [run()]
#'
#' @return the `regimen` with missing steps executed
#'
#' @details
#' Checkpoint priority:
#' 1. If `$paths$posterior` exists and `$status$fitted = FALSE`, load posterior
#'    and skip fit.
#' 2. If `$paths$prior` exists and `$status$tuned = FALSE`, load prior and skip
#'    tune.
#' 3. If `$paths$inference` exists and `$status$inferred = FALSE`, load inference
#'    and mark complete.
#'
#' @seealso [run()], [save_regimen()], [load_regimen()]
#'
#' @examples
#' \dontrun{
#' # Regimen crashed after fit; resume from checkpoint
#' reg <- load_regimen("./m1_regimen.rds") |>
#'   resume()
#' }
#'
#' @export
resume <- function(x, ...) UseMethod("resume")

#' @rdname resume
#' @export
resume.regimen <- function(x, ...) {
  steps_needed <- character()

  if (!x$status$prepared) {
    steps_needed <- c(steps_needed, "prepare")
  }

  if (!x$status$tuned && !file.exists(x$paths$prior)) {
    steps_needed <- c(steps_needed, "tune")
  } else if (file.exists(x$paths$prior) && !x$status$tuned) {
    message("Prior wflow found on disk — skipping tune, will load from: ",
            x$paths$prior)
    x$status$tuned <- TRUE
  }

  if (!x$status$fitted && !file.exists(x$paths$posterior)) {
    steps_needed <- c(steps_needed, "fit")
  } else if (file.exists(x$paths$posterior) && !x$status$fitted) {
    message("Posterior wflow found on disk — skipping fit, will load from: ",
            x$paths$posterior)
    x$status$fitted <- TRUE
  }

  if (!x$status$inferred && !file.exists(x$paths$inference)) {
    steps_needed <- c(steps_needed, "infer")
  } else if (file.exists(x$paths$inference) && !x$status$inferred) {
    message("Inference found on disk — loading from: ", x$paths$inference)
    x$inference       <- readRDS(x$paths$inference)
    x$status$inferred <- TRUE
  }

  if (length(steps_needed) == 0L) {
    message("Regimen already complete. Nothing to run.")
    return(x)
  }

  message("Resuming with steps: ", paste(steps_needed, collapse = " -> "))
  x |> run(steps = steps_needed, ...)
}

# ---- save / load ------------------------------------------------------------

#' Save a regimen to disk
#'
#' Serializes the entire regimen object (metadata, configuration, wflow, and
#' inference results) to an RDS file. Useful for checkpointing and resuming
#' analyses across sessions.
#'
#' @param x a `regimen` object
#' @param path character; file path for the `.rds` output
#'
#' @return the `regimen` invisibly
#'
#' @seealso [load_regimen()], [resume()]
#'
#' @export
save_regimen <- function(x, path) {
  stopifnot(inherits(x, "regimen"))
  dir_path <- dirname(path)
  if (!dir.exists(dir_path))
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
  message("regimen saved to: ", path)
  invisible(x)
}

#' Load a regimen from disk
#'
#' Deserializes a regimen object previously saved with [save_regimen()].
#'
#' @param path character; file path to the `.rds` regimen file
#'
#' @return a `regimen` object
#'
#' @seealso [save_regimen()], [resume()]
#'
#' @export
load_regimen <- function(path) {
  x <- readRDS(path)
  stopifnot(inherits(x, "regimen"))
  x
}

# ---- accessors --------------------------------------------------------------

#' Extract inference results for a named contrast
#'
#' Retrieves the `inference` object corresponding to a single named contrast
#' from a completed regimen.
#'
#' @param x a `regimen` object
#' @param contrast_name character scalar; key in `$inference` list
#'
#' @return an `inference` object, or `NULL` if inference not yet run for that
#'   contrast
#'
#' @seealso [infer()], [update_contrasts()]
#'
#' @export
get_inference <- function(x, contrast_name) {
  stopifnot(inherits(x, "regimen"))
  if (!contrast_name %in% names(x$inference))
    stop(sprintf("Contrast '%s' not found. Available: %s",
                 contrast_name,
                 paste(names(x$inference), collapse = ", ")))
  x$inference[[contrast_name]]
}

#' Update data in a regimen
#'
#' Replaces `$data_list` with new data and resets all status flags to `FALSE`.
#' Call [run()] or individual step functions to re-execute the pipeline.
#'
#' @param x a `regimen` object
#' @param data_list list with `data` (required) and optional `extra`, `cond_df`
#'   elements
#'
#' @return the `regimen` with updated data and reset status
#'
#' @details
#' Resetting status forces re-execution of all pipeline steps on the new data.
#' To update data *and* re-run only inference (keeping fitted model), use
#' [update_contrasts()] instead.
#'
#' @seealso [update_contrasts()], [run()]
#'
#' @export
update_data <- function(x, data_list) {
  stopifnot(inherits(x, "regimen"))
  x$data_list <- data_list
  x$wf        <- NULL
  x$inference <- list()
  x$status    <- list(prepared = FALSE, tuned = FALSE,
                      fitted = FALSE, inferred = FALSE)
  message("Data updated. Status reset — re-run `run()` or individual steps.")
  x
}

#' Update contrasts in a regimen
#'
#' Replaces `$contrasts` with new contrast specifications and resets only the
#' `$status$inferred` flag to `FALSE`. Keeps the fitted model intact and only
#' re-runs inference with the new contrasts.
#'
#' @param x a `regimen` object
#' @param contrasts named list of `contrast` objects; all elements must inherit
#'   from class `contrast`
#'
#' @return the `regimen` with updated contrasts and `$status$inferred = FALSE`
#'
#' @seealso [update_data()], [infer()]
#'
#' @examples
#' \dontrun{
#' # After fitting, explore additional contrasts without re-fitting
#' reg <- reg |>
#'   update_contrasts(list(new_contrast_1 = ..., new_contrast_2 = ...)) |>
#'   infer()
#' }
#'
#' @export
update_contrasts <- function(x, contrasts) {
  stopifnot(inherits(x, "regimen"))
  stopifnot(is.list(contrasts), !is.null(names(contrasts)))
  if (!all(vapply(contrasts, inherits, logical(1L), "contrast")))
    stop("All elements of `contrasts` must be `contrast` objects.")
  x$contrasts       <- contrasts
  x$inference       <- list()
  x$status$inferred <- FALSE
  message("Contrasts updated. Re-run `infer()` to regenerate inference.")
  x
}

# ---- helpers ----------------------------------------------------------------

#' @keywords internal
.check_prepared_regimen <- function(x) {
  if (!x$status$prepared || is.null(x$wf))
    stop("Regimen not prepared. Call `prepare()` or `run()` first.")
}

#' @keywords internal
.check_fitted_regimen <- function(x) {
  if (!x$status$fitted)
    stop("Regimen not fitted. Call `fit()` or `run()` first.")
}


# ---- Run workflow in steps  -----------------------------------
#
# Option 1: Run the entire regimen in one go
# reg <- reg |> run(save = TRUE)
#
# Option 2: Run steps individually
# reg <- reg |> prepare()
# reg <- reg |> tune(save = TRUE)
# reg <- reg |> fit(save = TRUE)
# reg <- reg |> infer(save = TRUE)

# Option 3: Resume from saved state (checks disk for prior/posterior/inference)
# reg <- reg |> resume()

# ---- Access results ---------------------------------------------------------

# Print regimen status
# reg

# Get inference for a specific contrast
# intact_minus_stroke_inf <- get_inference(reg, "intact_minus_stroke")
# dir_group_interaction_inf <- get_inference(reg, "dir_group_interaction")

# Access all inference results
# all_inference <- reg$inference

# ---- Update data or contrasts and re-run (optional) -------------------------

# Update contrasts and re-run inference only
# new_contrasts <- list(
#   my_new_contrast = contrast(...)
# )
# reg <- reg |> update_contrasts(new_contrasts) |> infer_regimen(save = TRUE)

# Update data and start from scratch
# new_data_list <- list(data = ..., extra = list(), cond_df = NULL)
# reg <- reg |> update_data(new_data_list) |> run_regimen(save = TRUE)