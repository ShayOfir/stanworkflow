# stanmodels.R

# --- Global variables ---
.default_sampler_params <- list(
  chains          = 4L,
  parallel_chains = 4L,
  iter_warmup     = 1000L,
  iter_sampling   = 1000L,
  adapt_delta     = 0.8,
  max_treedepth   = 10L,
  show_messages   = TRUE,
  show_exceptions = TRUE,
  refresh         = 100L
)

# ---- internal functions ----

#' @keywords internal
.get_variables <- function(draws) {
  posterior::variables(draws) |>
    gsub("\\[.*\\]", "", x = _) |>
    unique() |>
    setdiff("lp__")
}





#' @keywords internal
.extract_prior_names <- function(model_code) {
  matches <- regmatches(
    model_code,
    gregexpr("\\bprior_[A-Za-z0-9_]+", model_code)
  )[[1]]
  unique(matches)
}

#' @keywords internal
.parse_stan_data_block <- function(model_code) {
  block <- regmatches(
    model_code,
    regexpr("data\\s*\\{[^}]*\\}", model_code)
  )
  if (length(block) == 0) stop("No data block found in Stan model code.")

  lines <- strsplit(block, "\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines <- lines[!grepl("^data\\s*\\{|^\\}|^//", lines)]

  parsed <- lapply(lines, function(ln) {
    ln <- gsub(";.*$", "", ln)
    ln <- trimws(ln)
    if (!nzchar(ln)) return(NULL)

    # Extract the last word (variable name) from the declaration
    name <- trimws(gsub("^.*\\s+(\\w+)\\s*(\\[.*\\])?$", "\\1", ln))
    name <- gsub("[^A-Za-z0-9_]", "", name)

    base_type <- dplyr::case_when(
      grepl("^int",    ln) ~ "int",
      grepl("^real",   ln) ~ "real",
      grepl("^vector", ln) ~ "vector",
      grepl("^array",  ln) ~ "array",
      grepl("^matrix", ln) ~ "matrix",
      TRUE                 ~ "unknown"
    )
    list(name = name, type = base_type)
  })
  parsed <- Filter(Negate(is.null), parsed)
  names(parsed) <- vapply(parsed, `[[`, character(1), "name")
  parsed
}

# --- priors class ----

#' Create a new priors object
#' @param prior_names character vector of prior parameter names
#' @return an object of class `priors`
#' @importFrom stats setNames
#' @export
priors <- function(prior_names = character()) {
  structure(
    setNames(rep(list(NA), length(prior_names)), prior_names),
    class = "priors"
  )
}

#' @exportS3Method
print.priors <- function(x, ...) {
  cat("<priors>\n")
  for (nm in names(x)) {
    cat(sprintf("  %-30s= %s\n", nm, x[[nm]]))
  }
  invisible(x)
}

# --- sampler_params class ---

#' Create a new sampler_params object
#' @param params named list of sampler parameters and their values
#' @return an object of class `sampler_params`
#' @export
sampler_params <- function(params = .default_sampler_params) {
  structure(params, class = "sampler_params")
}

#' @exportS3Method
print.sampler_params <- function(x, ...) {
  cat("<sampler_params>\n")
  for (nm in names(x)) {
    cat(sprintf("  %-20s: %s\n", nm, x[[nm]]))
  }
  invisible(x)
}

# --- stanmodel class ---

#' Create a new stanmodel object
#' @return an object of class `stanmodel`
#' @export
stanmodel <- function() {
  structure(
    list(
      compiled_model = NULL,
      model_code     = NULL,
      priors         = NULL,
      sampler_params = sampler_params(),
      fitted_model   = NULL
    ),
    class = "stanmodel"
  )
}

#' Compile a stanmodel from a string of code
#' @param x a stanmodel object
#' @param model_code a string containing the Stan model code
#' @param ... additional arguments passed to `cmdstanr::cmdstan_model()`
#' @export
compile_from_string <- function(x, ...) UseMethod("compile_from_string")

#' @exportS3Method
compile_from_string.stanmodel <- function(x, model_code, ...) {
  stan_file <- cmdstanr::write_stan_file(model_code)
  compile_from_file(x, stan_file, ...)
}

#' Compile a stanmodel from a .stan file
#' @param x a stanmodel object
#' @param stan_file path to a .stan file
#' @param ... additional arguments passed to `cmdstanr::cmdstan_model()`
#' @export
compile_from_file <- function(x, ...) UseMethod("compile_from_file")

#' @exportS3Method
compile_from_file.stanmodel <- function(x, stan_file, ...) {
  model_code       <- paste(readLines(stan_file), collapse = "\n")
  x$model_code     <- model_code
  x$priors         <- priors(.extract_prior_names(model_code))
  x$compiled_model <- cmdstanr::cmdstan_model(stan_file = stan_file, ...)
  x
}

#' Build a Stan data list from a dataframe
#' @param x a stanmodel object
#' @param df a dataframe
#' @param extra named list of additional pre-computed values (e.g. matrices)
#' @param cond_df an optional dataframe of conditioning variables for *_cond_array vars
#' @return a named list; NA entries need to be filled by the user
#' @importFrom dplyr n_distinct
#' @export
get_stan_data <- function(x, ...) UseMethod("get_stan_data")

#' @exportS3Method
get_stan_data.stanmodel <- function(x,
                                    df,
                                    extra   = list(),
                                    cond_df = NULL) {
  vars   <- .parse_stan_data_block(x$model_code)
  vnames <- names(vars)
  result <- vector("list", length(vnames))
  names(result) <- vnames
  unfilled <- character()
  
  
  # ---- resolve cond_df columns into *_cond_array vars ----------------------
  cond_array_filled <- character()
  if (!is.null(cond_df)) {
    stopifnot(is.data.frame(cond_df))

    # N_cond_array
    if ("N_cond_array" %in% vnames) {
      result[["N_cond_array"]] <- nrow(cond_df)
      cond_array_filled <- c(cond_array_filled, "N_cond_array")
    }

    # For each column in cond_df, look for a matching *_cond_array variable
    # Matching rules mirror the main df logic:
    #   col               -> col_cond_array        (direct numeric/logical)
    #   col (factor/char) -> col_id_cond_array      (integer-encoded)
    for (col in names(cond_df)) {
      # direct match: col_cond_array
      direct_nm <- paste0(col, "_cond_array")
      if (direct_nm %in% vnames) {
        result[[direct_nm]] <- cond_df[[col]]
        cond_array_filled   <- c(cond_array_filled, direct_nm)
      }

      # id match: col_id_cond_array
      id_nm <- paste0(col, "_id_cond_array")
      if (id_nm %in% vnames) {
        result[[id_nm]]   <- as.integer(factor(cond_df[[col]]))
        cond_array_filled <- c(cond_array_filled, id_nm)
      }
    }
  }

# ---- main df resolution (unchanged logic) --------------------------------
  for (nm in vnames) {

    if (nm %in% cond_array_filled) next   # already filled above

    if (nm == "N") { result[["N"]] <- nrow(df); next }

    if (grepl("^n_", nm)) {
      col <- sub("^n_", "", nm)
      if (col %in% names(df)) result[[nm]] <- dplyr::n_distinct(df[[col]])
      else { result[[nm]] <- NA; unfilled <- c(unfilled, nm) }
      next
    }

    if (grepl("_id$", nm)) {
      col <- sub("_id$", "", nm)
      if (col %in% names(df)) result[[nm]] <- as.integer(factor(df[[col]]))
      else { result[[nm]] <- NA; unfilled <- c(unfilled, nm) }
      next
    }

    # if (nm %in% names(outcome)) { result[[nm]] <- df[[outcome[[nm]]]]; next }

    if (grepl("^prior_", nm)) { result[[nm]] <- NA_real_; next }

    if (nm %in% names(df))
       { result[[nm]] <- df[[nm]]; next} # for datafields

    result[[nm]] <- NA
    unfilled <- c(unfilled, nm)
  }

  result[names(extra)] <- extra

 
  unfilled <- setdiff(unfilled, names(extra))
 
  if (length(unfilled) > 0){
    if (!(length(unfilled) == 1 && (unfilled[1]) == 'skip_likelihood')){
      print("Unfilled variables in stan_data:")
      print(unfilled)    
      message(
        "The following variables could not be auto-filled and need to be supplied via `extra`:\n  ",
        paste(unfilled, collapse = ", ")
      )
    }
  }

  result
}

#' Fit a stanmodel to data
#' @param x a stanmodel object
#' @param data a dataframe containing the data to fit; columns should match those expected by the model code and specified in `stan_data`
#' @param outcome_var the name of the outcome variable in `data` (default: "y")
#' @param stan_data an optional named list of data to pass to Stan; if NULL, will be generated from `data` and `outcome_var` using `get_stan_data()`
#' @param priors an optional `priors` object; overrides `prior_*` fields in `stan_data`
#' @param sampler_params an optional `sampler_params` object
#' @param prior_predictive if TRUE, sets `skip_likelihood` to 1 in `stan_data` for prior predictive checks (if the model supports it)
#' @param extra a named list of additional pre-computed values to merge into `stan_data` (e.g. matrices)
#' @param cond_df an optional dataframe of conditioning variables for *_cond_array vars
#' @export
fit <- function(x, ...) UseMethod("fit")

#' @exportS3Method
fit.stanmodel <- function(x,
                          data,
                          stan_data = NULL,
                          priors         = NULL,
                          sampler_params = NULL,
                          prior_predictive = FALSE,
                          outcome_var = 'y',
                          extra = list(),
                          cond_df = NULL,
                          ...) {
  if (is.null(x$compiled_model)) {
    stop("Model has not been compiled yet. Run `compile_from_string()` or `compile_from_file()` first.")
  }


  # if the model contains skip_liklihood flag, set it to 1 for prior predictive checks
  # Prepare data

     
  stan_data <- get_stan_data(
    x,
    df      = data,
    extra = extra,
    cond_df = cond_df
  )
   

  
  if (is.null(sampler_params)) sampler_params <- x$sampler_params
  stopifnot(inherits(sampler_params, "sampler_params"))

  # 1. merge priors first
  if (is.null(priors)) priors <- x$priors
  if (!is.null(priors)) {
    stopifnot(inherits(priors, "priors"))
    # only merge non-NA prior values — don't clobber already-set fields
    prior_vals <- Filter(Negate(is.na), priors)
    stan_data <- modifyList(stan_data, prior_vals)
  }

  # 2. skip_likelihood after priors merge
  if ("skip_likelihood" %in% names(stan_data)) {
      stan_data$skip_likelihood <- as.integer(prior_predictive)
    } else {
    if (prior_predictive)
      warning("Model does not contain `skip_likelihood` variable; cannot perform prior predictive checks.")
  }  

  # 3. NA check after everything is merged
  
  na_fields <- names(Filter(anyNA, stan_data))
  if (length(na_fields) > 0) {
    stop(
      "stan_data contains NA values in: ",
      paste(na_fields, collapse = ", "),
      "\nSupply values via `priors` or `extra` before fitting."
    )
  }

  x$fitted_model <- do.call(
    x$compiled_model$sample,
    c(list(data = stan_data), sampler_params, list(...))
  )

  x$data <- data
  x$cond_df <- cond_df
  x$extra <- extra
  x$stan_data <- stan_data
  x$outcome_var <- outcome_var


  # Draw posterior samples:
  x$draws <- x$fitted_model$draws()

  x$variables <- .get_variables(x$draws)
  return(x)
}


#' gq - method for using generate_quantities() method in stanmodel
#'
#' @param x a `stanmodel` object
#' @param new_data dataset to use for generated quantities
#' @param n_draws optional integer; if set and less than the number of available
#'   posterior draws, a random subsample of `n_draws` draws is used for
#'   `generate_quantities()` instead of the full fit. Trades Monte Carlo
#'   precision for speed and output size. Default `NULL` uses all draws
#'   (unchanged behavior).
#' @param output_dir optional directory for CmdStan's GQ output CSVs. Default
#'   `NULL` uses a fresh temp directory that is deleted automatically after
#'   the draws are extracted (avoids leaking temp files across repeated calls).
#'   If supplied, the caller owns cleanup of `output_dir`.
#' @param cache_compiled logical; if `TRUE` (default), the standalone-GQ model
#'   compiled from `x$model_code` is memoized per session (keyed by
#'   `model_code`) so repeated `gq()` calls on the same model skip
#'   `write_stan_file()`/`cmdstan_model()` after the first call. Set `FALSE`
#'   to force the original per-call recompile behavior.
#' @return draws object containing the generated quantities
gq <- function(x, ...) UseMethod("gq")

# Session-level cache of compiled standalone-GQ models, keyed by model_code,
# so gq.stanmodel() only pays for write_stan_file()+cmdstan_model() once per
# unique model within a session instead of once per call. This is purely a
# memoization of the exact same compile call below — set cache_compiled =
# FALSE on gq() to bypass it entirely and get the original per-call recompile
# behavior back with no code edit.
.gq_model_cache <- new.env(parent = emptyenv())

#' @exportS3Method
gq.stanmodel <- function (x, new_data, n_draws = NULL, output_dir = NULL, cache_compiled = TRUE){
  mod_gqp <- if (cache_compiled) .gq_model_cache[[x$model_code]] else NULL
  if (is.null(mod_gqp)) {
    gqp <- cmdstanr::write_stan_file(x$model_code)
    mod_gqp <- cmdstanr::cmdstan_model(gqp)
    if (cache_compiled) .gq_model_cache[[x$model_code]] <- mod_gqp
  }

  mod_fit <- x$fitted_model
  if (!is.null(n_draws)) {
    draws_full  <- mod_fit$draws()
    n_available <- posterior::ndraws(draws_full)
    idx <- if (n_available <= n_draws) seq_len(n_available) else
      sample.int(n_available, n_draws)
    mod_fit <- posterior::subset_draws(draws_full, draw = idx)
  }

  pc <- x$sampler_params$parallel_chains %||% 4L

  own_dir <- is.null(output_dir)
  if (own_dir) {
    output_dir <- tempfile("gq_")
    dir.create(output_dir, recursive = TRUE)
  }

  fit_gq <- mod_gqp$generate_quantities(mod_fit, data = new_data, parallel_chains = pc, output_dir = output_dir)
  out <- fit_gq$draws()
  rm(fit_gq)
  if (own_dir) unlink(output_dir, recursive = TRUE)
  out
}


# ---- ppc helpers ------------------------------------------------------------

#' @keywords internal
.interpret_bayes_pvalue <- function(p, stat_name) {
  judgement <- dplyr::case_when(
    p < 0.05               ~ "severe misfit",
    p < 0.10               ~ "moderate misfit",
    p >= 0.10 & p <= 0.90  ~ "acceptable",
    p > 0.90 & p <= 0.95   ~ "moderate misfit",
    p > 0.95               ~ "severe misfit"
  )
  direction <- dplyr::case_when(
    stat_name == "mean"   & p < 0.10 ~ " — y_rep means below observed (prior/likelihood biased low)",
    stat_name == "mean"   & p > 0.90 ~ " — y_rep means above observed (prior/likelihood biased high)",
    stat_name == "sd"     & p < 0.10 ~ " — model under-dispersed relative to data",
    stat_name == "sd"     & p > 0.90 ~ " — model over-dispersed relative to data",
    stat_name == "median" & p < 0.10 ~ " — replicated medians below observed",
    stat_name == "median" & p > 0.90 ~ " — replicated medians above observed",
    stat_name == "mad"    & p < 0.10 ~ " — model robust spread too narrow",
    stat_name == "mad"    & p > 0.90 ~ " — model robust spread too wide",
    stat_name == "min"    & p < 0.10 ~ " — model does not capture lower tail",
    stat_name == "min"    & p > 0.90 ~ " — model lower tail too extreme",
    stat_name == "max"    & p < 0.10 ~ " — model does not capture upper tail",
    stat_name == "max"    & p > 0.90 ~ " — model upper tail too extreme",
    TRUE                             ~ ""
  )
  sprintf("p_B = %.3f [%s%s]", p, judgement, direction)
}

#' @keywords internal
.interpret_coverage <- function(empirical, nominal) {
  tol_lo <- nominal - 0.10
  tol_hi <- nominal + 0.05
  judgement <- dplyr::case_when(
    empirical < tol_lo ~ sprintf(
      "under-coverage (%.1f%% < nominal %.0f%%) — predictive intervals too narrow; model may be overconfident",
      empirical * 100, nominal * 100
    ),
    empirical > tol_hi ~ sprintf(
      "over-coverage (%.1f%% > nominal %.0f%%) — predictive intervals too wide; prior likely dominates likelihood",
      empirical * 100, nominal * 100
    ),
    TRUE ~ sprintf(
      "near-nominal (%.1f%% ≈ %.0f%%) — well-calibrated",
      empirical * 100, nominal * 100
    )
  )
  judgement
}

# ---- ppc --------------------------------------------------------------------



# ---- ppc helpers ------------------------------------------------------------

#' @keywords internal
.interpret_bayes_pvalue <- function(p, stat_name) {
  judgement <- dplyr::case_when(
    p < 0.05               ~ "severe misfit",
    p < 0.10               ~ "moderate misfit",
    p >= 0.10 & p <= 0.90  ~ "acceptable",
    p > 0.90 & p <= 0.95   ~ "moderate misfit",
    p > 0.95               ~ "severe misfit"
  )
  direction <- dplyr::case_when(
    stat_name == "mean"   & p < 0.10 ~ " — y_rep means below observed (prior/likelihood biased low)",
    stat_name == "mean"   & p > 0.90 ~ " — y_rep means above observed (prior/likelihood biased high)",
    stat_name == "sd"     & p < 0.10 ~ " — model under-dispersed relative to data",
    stat_name == "sd"     & p > 0.90 ~ " — model over-dispersed relative to data",
    stat_name == "median" & p < 0.10 ~ " — replicated medians below observed",
    stat_name == "median" & p > 0.90 ~ " — replicated medians above observed",
    stat_name == "mad"    & p < 0.10 ~ " — model robust spread too narrow",
    stat_name == "mad"    & p > 0.90 ~ " — model robust spread too wide",
    stat_name == "min"    & p < 0.10 ~ " — model does not capture lower tail",
    stat_name == "min"    & p > 0.90 ~ " — model lower tail too extreme",
    stat_name == "max"    & p < 0.10 ~ " — model does not capture upper tail",
    stat_name == "max"    & p > 0.90 ~ " — model upper tail too extreme",
    TRUE                             ~ ""
  )
  sprintf("p_B = %.3f [%s%s]", p, judgement, direction)
}

#' @keywords internal
.interpret_coverage <- function(empirical, nominal) {
  tol_lo <- nominal - 0.10
  tol_hi <- nominal + 0.05
  dplyr::case_when(
    empirical < tol_lo ~ sprintf(
      "under-coverage (%.1f%% < nominal %.0f%%) — predictive intervals too narrow; model may be overconfident",
      empirical * 100, nominal * 100
    ),
    empirical > tol_hi ~ sprintf(
      "over-coverage (%.1f%% > nominal %.0f%%) — predictive intervals too wide; prior likely dominates likelihood",
      empirical * 100, nominal * 100
    ),
    TRUE ~ sprintf(
      "near-nominal (%.1f%% ≈ %.0f%%) — well-calibrated",
      empirical * 100, nominal * 100
    )
  )
}

#' Run a single PPC for one (y_obs, y_rep_mat) pair
#'
#' @param y_obs      numeric vector of observed values
#' @param y_rep_mat  draws matrix (n_draws × n_obs)
#' @param obs_var    character label used in plot titles
#' @param n_draws    number of draws to use for density overlay
#' @return a `ppc_results` object
#' @keywords internal
.run_single_ppc <- function(y_obs, y_rep_mat, obs_var, n_draws = 50L) {

  n_available <- nrow(y_rep_mat)
  idx_draws   <- if (n_available <= n_draws) seq_len(n_available) else
    sample.int(n_available, n_draws)
  y_rep_sub <- y_rep_mat[idx_draws, , drop = FALSE]

  # ---- test statistics ----------------------------------------------------
  stats <- list(
    mean   = mean,
    sd     = stats::sd,
    median = stats::median,
    mad    = stats::mad,
    min    = min,
    max    = max
  )

  T_obs <- vapply(stats, function(f) f(y_obs), numeric(1))

  T_rep <- vapply(stats, function(f) {
    apply(y_rep_mat, 1L, f)
  }, numeric(nrow(y_rep_mat)))
  # T_rep is n_draws × n_stats after vapply; transpose to n_stats × n_draws
  if (is.matrix(T_rep)) T_rep <- t(T_rep)

  bayes_pvalues <- vapply(names(stats), function(s) {
    mean(T_rep[s, ] >= T_obs[[s]])
  }, numeric(1))

  # ---- coverage -----------------------------------------------------------
  q50 <- matrixStats::colQuantiles(y_rep_mat, probs = c(0.25, 0.75))
  q90 <- matrixStats::colQuantiles(y_rep_mat, probs = c(0.05, 0.95))
  coverage <- c(
    `50%` = mean(y_obs >= q50[, 1] & y_obs <= q50[, 2]),
    `90%` = mean(y_obs >= q90[, 1] & y_obs <= q90[, 2])
  )

  # ---- plots --------------------------------------------------------------
  title_base <- sprintf("— %s (Gabry et al.)", obs_var)

  p_density <- bayesplot::ppc_dens_overlay(y_obs, y_rep_sub) +
    ggplot2::labs(
      title    = sprintf("Density overlay — %s", obs_var),
      subtitle = "Dark line = observed y; light lines = y_rep draws."
    )

  stat_plots <- lapply(names(stats), function(s) {
    bayesplot::ppc_stat(y_obs, y_rep_mat, stat = s) +
      ggplot2::labs(title = sprintf("T = %s %s", s, title_base))
  })
  names(stat_plots) <- names(stats)

  p_mean_sd <- bayesplot::ppc_stat_2d(y_obs, y_rep_mat, stat = c("mean", "sd")) +
    ggplot2::labs(title = sprintf("Joint mean-sd %s", title_base))

  plots <- c(
    list(density_overlay = p_density),
    stat_plots,
    list(mean_sd_2d = p_mean_sd)
  )

  # ---- interpretation -----------------------------------------------------
  interpretation <- list(
    bayes_pvalues = mapply(.interpret_bayes_pvalue,
                           p         = bayes_pvalues,
                           stat_name = names(bayes_pvalues),
                           SIMPLIFY  = FALSE),
    coverage      = mapply(.interpret_coverage,
                           empirical = coverage,
                           nominal   = c(0.50, 0.90),
                           SIMPLIFY  = FALSE)
  )

  ppc_results(
    plots          = plots,
    coverage       = coverage,
    bayes_pvalues  = bayes_pvalues,
    interpretation = interpretation,
    outcomes       = obs_var
  )
}

#' Combine per-group ppc_results into faceted ggplot objects
#'
#' @param group_results named list of `ppc_results` objects (one per group key)
#' @param group_label_df data frame with columns `.grp_key` + `.by` variables (factored)
#' @param obs_var  character label for the outcome variable
#' @param by_vars  character vector of grouping variable names
#' @return named list of faceted ggplot objects, one per plot type
#' @keywords internal
.build_faceted_plots <- function(group_results, group_label_df, obs_var, by_vars) {

  stopifnot(
    is.list(group_results),
    length(group_results) > 0,
    is.data.frame(group_label_df),
    is.character(obs_var),
    is.character(by_vars)
  )

  plot_names <- names(group_results[[1]]$plots)
  if (is.null(plot_names) || length(plot_names) == 0)
    stop("group_results[[1]]$plots is empty or missing.")

  # For each plot type, extract, combine, and facet
  result_plots <- lapply(plot_names, function(pnm) {

    # Step 1: Build each group's plot and extract layer data
    plot_data_list <- lapply(names(group_results), function(gk) {
      if (!pnm %in% names(group_results[[gk]]$plots)) {
        stop(sprintf("Plot '%s' not found in group '%s'.", pnm, gk))
      }

      p <- group_results[[gk]]$plots[[pnm]]
      grp_row <- group_label_df[group_label_df$.grp_key == gk, by_vars, drop = FALSE][1L, ]

      built <- ggplot2::ggplot_build(p)

      # Attach group label columns to each layer's computed data
      lapply(built$data, function(layer_df) {
        for (bv in by_vars) {
          layer_df[[bv]] <- grp_row[[bv]]
        }
        layer_df
      })
    })

    # Step 2: Row-bind layer data across all groups
    n_layers <- length(plot_data_list[[1]])
    all_layer_data <- lapply(seq_len(n_layers), function(li) {
      do.call(rbind, lapply(plot_data_list, `[[`, li))
    })

    # Step 3: Use template's raw plot structure, inject combined layer data
    template <- group_results[[1]]$plots[[pnm]]
    p_out    <- template

    # Replace each layer's data with combined built data AND replace stat
    # with StatIdentity so no recomputation occurs
    for (li in seq_len(n_layers)) {
      p_out$layers[[li]]$data        <- all_layer_data[[li]]
      p_out$layers[[li]]$stat        <- ggplot2::StatIdentity$new()
      p_out$layers[[li]]$stat_params <- list()
    }

    # Step 4: Add faceting
    if (length(by_vars) == 2L) {
      p_out <- p_out +
        ggplot2::facet_grid(
          rows = ggplot2::vars(!!ggplot2::sym(by_vars[2L])),
          cols = ggplot2::vars(!!ggplot2::sym(by_vars[1L]))
        )
    } else if (length(by_vars) == 1L) {
      p_out <- p_out +
        ggplot2::facet_wrap(ggplot2::vars(!!ggplot2::sym(by_vars[1L])))
    }

    p_out +
      ggplot2::labs(title = sprintf("%s — %s", pnm, obs_var))
  })

  # Name the result list with plot names
  names(result_plots) <- plot_names
  result_plots
}





#' Prior and posterior predictive checks
#'
#' Compares replicated datasets `y_rep` from model draws to the observed
#' outcome vector using the full set of graphical and numerical checks
#' recommended by Gabry et al. (2019) and Gelman et al. (2020).
#'
#' **Checks performed:**
#' - *Density overlay* (`ppc_dens_overlay`): visual band of replicated
#'   densities vs. observed. Too narrow = overconfident model; too wide =
#'   prior dominates.
#' - *Marginal test statistics* (`ppc_stat`): mean, sd, median, mad, min, max.
#'   Bayesian p-value p_B = P(T(y_rep) ≥ T(y)). Values outside (0.05, 0.95)
#'   flag misfit in that dimension (Gelman BDA3 §6.3).
#' - *Joint mean–sd plot* (`ppc_stat_2d`): Gabry et al. primary visual
#'   recommendation. The observed point should sit inside the bulk of the
#'   cloud; peripheral placement signals systematic location–scale misfit.
#' - *Coverage calibration*: empirical coverage of 50% and 90% predictive
#'   intervals. Systematic under-coverage → too tight; over-coverage → too
#'   wide (Gelman et al. 2020 §4).
#'
#' @param x a fitted `stanmodel` object
#' @param n_draws number of y_rep draws for the density overlay (default 50)
#' @param data optional dataframe to use for filtering; if NULL, uses `x$data`
#' @param filter_expr an optional dplyr-style expression to filter rows of `data`
#' @param filter_fn an optional function of the form `function(y_obs, y_rep)` that returns a list with filtered `y_obs` and `y_rep`. If both `filter_expr` and `filter_fn` are provided, `filter_expr` takes precedence.
#' @param outcomes an optional named list of outcome specifications. Each element should be a list with `obs_var` (name of the observed variable in `x$stan_data`) and `rep_var` (name of the replicated variable in `x$draws`). If NULL, defaults to the single outcome used at fit time.
#' @return a named list:
#'   \describe{
#'     \item{`plots`}{named list of ggplot objects: `density_overlay`,
#'       `mean`, `sd`, `median`, `mad`, `min`, `max`, `mean_sd_2d`}
#'     \item{`coverage`}{named numeric of empirical coverage at 50% and 90%;
#'       see `$interpretation$coverage`}
#'     \item{`bayes_pvalues`}{named numeric of p_B per test statistic;
#'       see `$interpretation$bayes_pvalues`}
#'     \item{`interpretation`}{named list of human-readable strings for
#'       `coverage` and `bayes_pvalues`}
#'   }
#' @references
#'   Gabry et al. (2019). Visualization in Bayesian workflow. JRSS-A 182(2).
#'   Gelman et al. (2020). Bayesian Workflow. arXiv:2011.01808.
#' @importFrom bayesplot ppc_dens_overlay ppc_stat ppc_stat_2d
#' @importFrom matrixStats colQuantiles
#' @importFrom posterior subset_draws as_draws_matrix
#' @export
ppc <- function(x, ...) UseMethod("ppc")



# ---- ppc --------------------------------------------------------------------
#' @exportS3Method
ppc.stanmodel <- function(x,
                           n_draws     = 50L,
                           data        = NULL,
                           filter_expr = NULL,
                           filter_fn   = NULL,
                           outcomes    = NULL,
                           .by         = NULL,
                           .levels     = NULL,
                           ...) {

  data <- data %||% x$data

  # ---- validate .by -------------------------------------------------------
  if (!is.null(.by)) {
    stopifnot(is.character(.by), length(.by) <= 2L)
    if (is.null(data))
      stop("`.by` requires `data` to be present (x$data or passed via `data =`).")
    missing_by <- setdiff(.by, names(data))
    if (length(missing_by) > 0)
      stop(sprintf("`.by` variables not found in data: %s",
                   paste(missing_by, collapse = ", ")))
    if (!is.null(.levels)) {
      stopifnot(is.list(.levels))
      missing_lvl <- setdiff(.by, names(.levels))
      if (length(missing_lvl) > 0)
        stop(sprintf("`.levels` missing entries for: %s",
                     paste(missing_lvl, collapse = ", ")))
    }
  }

  # ---- basic guards -------------------------------------------------------
  if (is.null(outcomes)) {
    stan_outcome             <- x$outcome_var %||% "y"
    outcomes                 <- list()
    outcomes[[stan_outcome]] <- list(obs_var = stan_outcome, rep_var = "y_rep")
  }
  if (is.null(x$draws))     stop("No draws found. Run `fit()` first.")
  if (is.null(x$stan_data)) stop("No stan_data found. Run `fit()` first.")

  # ---- global filter (applies to all outcomes) ----------------------------
  global_filter_fn <- NULL
  if (!is.null(filter_fn)) {
    stopifnot(is.function(filter_fn))
    global_filter_fn <- filter_fn
  }
  fe <- substitute(filter_expr)
  if (!is.null(fe)) {
    if (is.null(data)) stop("`filter_expr` supplied but no `data` found.")
    keep_idx         <- which(eval(fe, envir = as.list(data), enclos = parent.frame()))
    global_filter_fn <- function(y_obs, y_rep) {
      list(y_obs = y_obs[keep_idx],
           y_rep = y_rep[, keep_idx, drop = FALSE])
    }
  }

  # ---- process each outcome -----------------------------------------------
  results <- lapply(outcomes, function(spec) {

    obs_var <- spec$obs_var
    rep_var <- spec$rep_var
    source  <- spec$source %||% "stan_data"

    # ------------------------------------------------------------------
    # rep_filter_expr: evaluated against x$stan_data to select COLUMNS
    # of y_rep_mat (draw-level indices).
    #
    # Use case: mu draws have one column per condition in cond_df
    # (including response==0 conditions), while mu_y_obs only covers
    # response==1 conditions. rep_filter_expr selects the matching
    # subset of mu columns so ncol(y_rep_mat) == length(y_obs).
    #
    # Example:
    #   rep_filter_expr = rlang::expr(response_area_cond_array == 1L)
    # evaluated against x$stan_data, which contains
    # response_area_cond_array as a cond_array-length integer vector.
    # ------------------------------------------------------------------
    rep_col_indices <- NULL
    if (!is.null(spec$rep_filter_expr)) {
      rep_col_indices <- which(
        rlang::eval_tidy(spec$rep_filter_expr,
                         data = as.list(x$stan_data))
      )
      if (length(rep_col_indices) == 0L)
        stop(sprintf("outcome '%s': `rep_filter_expr` removed all draw columns.", obs_var))
    }

    # ---- resolve y_obs --------------------------------------------------
    if (source == "cond_df") {
      if (is.null(x$cond_df))
        stop(sprintf("outcome '%s': source='cond_df' but no cond_df stored.", obs_var))
      if (!obs_var %in% names(x$stan_data))
        stop(sprintf("obs_var '%s' not found in stan_data.", obs_var))
      cond_id <- x$stan_data[["cond_id"]]
      if (is.null(cond_id)) stop("stan_data$cond_id not found.")
      y_obs <- x$stan_data[[obs_var]][cond_id]
    } else {
      if (!obs_var %in% names(x$stan_data))
        stop(sprintf("obs_var '%s' not found in stan_data.", obs_var))
      y_obs <- x$stan_data[[obs_var]]
    }

    if (!is.numeric(y_obs))
      stop(sprintf("Resolved obs_var '%s' must be numeric.", obs_var))
    if (!rep_var %in% x$variables)
      stop(sprintf("rep_var '%s' not found in draws.", rep_var))

    # ---- y_rep matrix (n_draws x n_obs) ---------------------------------
    y_rep_mat <- posterior::subset_draws(x$draws, variable = rep_var) |>
      posterior::as_draws_matrix()

    # ---- apply rep_filter_expr: subset columns of y_rep_mat -------------
    # This aligns y_rep columns with y_obs when the rep variable covers
    # a superset of conditions (e.g. mu covers all cond_df rows including
    # response==0, but mu_y_obs only covers response==1 rows).
    if (!is.null(rep_col_indices)) {
      y_rep_mat <- y_rep_mat[, rep_col_indices, drop = FALSE]
    }

    if (ncol(y_rep_mat) != length(y_obs))
      stop(sprintf(paste0(
        "%s: after rep_filter_expr, y_rep has %d columns but y_obs has %d elements.\n",
        "  Ensure rep_filter_expr selects exactly the conditions present in stan_data[['%s']]."),
        obs_var, ncol(y_rep_mat), length(y_obs), obs_var
      ))

    # ---- apply global filter --------------------------------------------
    if (!is.null(global_filter_fn)) {
      filtered  <- global_filter_fn(y_obs, y_rep_mat)
      y_obs     <- filtered$y_obs
      y_rep_mat <- filtered$y_rep
      if (length(y_obs) == 0L || ncol(y_rep_mat) == 0L)
        stop(sprintf("%s: global filter removed all observations.", obs_var))
      if (length(y_obs) != ncol(y_rep_mat))
        stop(sprintf("%s: global filter_fn returned mismatched y_obs/y_rep lengths.", obs_var))
    }

    # =========================================================================
    # No grouping -> single PPC
    # =========================================================================
    if (is.null(.by)) {
      return(.run_single_ppc(y_obs, y_rep_mat, obs_var, n_draws))
    }

    # =========================================================================
    # Grouped PPC
    # data aligns 1:1 with y_obs: for extensive, data has N rows matching
    # trial-expanded cond obs; for intensive, data has N rows and y_obs
    # is also N (both are trial-level, response==1 only, no filter needed).
    # =========================================================================
    if (length(y_obs) != nrow(data))
      stop(sprintf(paste0(
        "%s: y_obs has %d elements but data has %d rows.\n",
        "  For intensive outcomes where y_obs is trial-level (stan_data$%s),\n",
        "  data must also be trial-level with the same N.\n",
        "  If mu draws cover more conditions than y_obs, use `rep_filter_expr`."),
        obs_var, length(y_obs), nrow(data), obs_var
      ))

    # ---- resolve group labels -------------------------------------------
    group_df <- data[, .by, drop = FALSE]

    for (bv in .by) {
      col <- group_df[[bv]]
      group_df[[bv]] <- if (!is.null(.levels) && bv %in% names(.levels)) {
        lvls <- .levels[[bv]]
        if (is.factor(col)) as.character(lvls[as.integer(col)])
        else                as.character(lvls[col])
      } else if (is.factor(col)) {
        as.character(col)
      } else {
        as.character(col)
      }
    }

    group_df$.grp_key <- do.call(paste, c(group_df[.by], sep = "__"))
    group_keys        <- unique(group_df$.grp_key)
    group_label_df    <- unique(group_df)

    level_order <- lapply(.by, function(bv) {
      col <- data[[bv]]
      if (!is.null(.levels) && bv %in% names(.levels))
        return(as.character(.levels[[bv]]))
      if (is.factor(col)) return(levels(col))
      unique(as.character(col))
    })
    names(level_order) <- .by

    for (bv in .by)
      group_label_df[[bv]] <- factor(group_label_df[[bv]], levels = level_order[[bv]])

    # ---- per-group PPC --------------------------------------------------
    group_results <- lapply(
      stats::setNames(group_keys, group_keys),
      function(gk) {
        idx_g   <- which(group_df$.grp_key == gk)
        y_obs_g <- y_obs[idx_g]
        y_rep_g <- y_rep_mat[, idx_g, drop = FALSE]
        grp_row <- group_label_df[group_label_df$.grp_key == gk, .by, drop = FALSE][1L, ]
        lbl     <- paste(paste(.by, unlist(grp_row), sep = "="), collapse = " | ")
        .run_single_ppc(y_obs_g, y_rep_g,
                        sprintf("%s [%s]", obs_var, lbl),
                        n_draws)
      }
    )

    # ---- faceted plots --------------------------------------------------
    #faceted <- .build_faceted_plots(group_results, group_label_df, obs_var, .by)

    # ---- summary table --------------------------------------------------
    summary_tbl <- purrr::imap_dfr(group_results, function(res, gk) {
      grp_row <- group_label_df[group_label_df$.grp_key == gk, .by, drop = FALSE][1L, ]
      bp  <- as.list(res$bayes_pvalues);  names(bp)  <- paste0("p_B_", names(bp))
      cov <- as.list(res$coverage);       names(cov) <- paste0("cov_", names(cov))
      dplyr::bind_cols(grp_row, as.data.frame(bp), as.data.frame(cov))
    })


    ppc_results_grouped(
      by            = .by,
      outcomes      = outcomes,
      group_results = group_results,
      summary       = summary_tbl
    )
  })

  if (length(results) == 1L) results[[1L]] else results
}




#' @exportS3Method
print.stanmodel <- function(x, ...) {
  cat("<stanmodel>\n")
  cat("  compiled  :", !is.null(x$compiled_model), "\n")
  cat("  fitted    :", !is.null(x$fitted_model), "\n")
  cat("  code lines:", if (!is.null(x$model_code)) length(strsplit(x$model_code, "\n")[[1]]) else 0L, "\n")
  cat("  priors    :\n")
  if (!is.null(x$priors)) print(x$priors) else cat("    <none>\n")
  invisible(x)
}


#' Create a ppc_results_grouped object
#' @keywords internal
ppc_results_grouped <- function(by, outcomes, group_results, summary) {
  structure(
    list(
      by = by,
      outcomes = outcomes,
      group_results = group_results,
      summary = summary
    ),
    class = "ppc_results_grouped"
  )
}


# ---- ppc_results class ------------------------------------------------------

#' Create a ppc_results object
#' @keywords internal
ppc_results <- function(plots, coverage, bayes_pvalues, interpretation, outcomes) {
  structure(
    list(
      plots          = plots,
      coverage       = coverage,
      bayes_pvalues  = bayes_pvalues,
      interpretation = interpretation,
      outcomes       = outcomes
    ),
    class = "ppc_results"
  )
}

#' Flatten y_rep draws into a long dataframe with factor labels
#'
#' @param x a fitted `stanmodel` obbject
#' @param data      the original data frame used to fit the model (nrow = N)
#' @param feature_names   character vector of column names in `data` to attach
#' @param enforce_features logical (= TRUE) whether to require that all `feature_names` are present in `data`; if FALSE, only the features that exist in `data` will be included
#' @param preserve_.draws logical (= FALSE) whether to keep the .chain, .iteration, and .draw columns in the output; if FALSE, these are dropped for a cleaner format
#' @return a tibble with columns: `outcome_type` (observed vs. posterior), `outcome` (numeric), and one column per feature in `feature_names`
#' @importFrom posterior as_draws_df subset_draws
#' @importFrom dplyr select mutate
#' @importFrom tidyr pivot_longer
#' @export
grouped_ppc <- function(x, ...) UseMethod("grouped_ppc")

#' @exportS3Method
grouped_ppc.stanmodel <- function(x, 
                                  data,
                                  feature_names,
                                  enforce_features = TRUE,
                                  preserve_.draws = FALSE) {

  fnames <- feature_names
  df <- data
  outcome_name <- x$outcome_var %||% "y"

  # ---- validate -----------------------------------------------------------
  stopifnot(inherits(x, "stanmodel"))
  stopifnot(is.data.frame(df))
  stopifnot(is.character(fnames))
  stopifnot(is.character(outcome_name))
  if (enforce_features)
    { 
      stopifnot(all(c(fnames, outcome_name) %in% names(df)))
    }
    else {
      # fnames will include only feature names that are present in the data
      fnames <- intersect(fnames, names(df))
    }

  if (is.null(x$draws)) {
    stop("No draws found in the model. Run `fit()` first.")
  }

  draws <- x$draws
  if (!"y_rep" %in% x$variables) {
    stop("No `y_rep` variable found in draws. Ensure your model code includes a generated quantity named `y_rep`.")
  }

  # ---- extract y_rep as draws_df ------------------------------------------
  # result: S rows (draws) × N columns (y_rep[1]...y_rep[N]) + draw metadata
  y_rep_df <- draws |>
    posterior::subset_draws(variable = "y_rep") |>
    posterior::as_draws_df()                       # S × (N + 3 metadata cols)

  n_obs   <- ncol(y_rep_df) - 3L  # subtract .chain, .iteration, .draw
  n_draws <- nrow(y_rep_df)

  if (n_obs != nrow(df)) {
    stop(sprintf(
      "y_rep has %d observations but df has %d rows — they must match.",
      n_obs, nrow(df)
    ))
  }


  # ---- pivot y_rep to long ------------------------------------------------
  y_rep_long <- y_rep_df |>
    tidyr::pivot_longer(
      cols      = starts_with("y_rep"),
      names_to  = ".obs_idx",
      values_to = "y_rep"
    ) |>
    dplyr::mutate(
      .obs_idx = as.integer(stringr::str_extract(.obs_idx, "\\d+"))
    )
  
   # define cols_to_tile here — outcome + features to carry through
    cols_to_tile <- c(outcome_name, fnames)   # <-- add this line

  # ---- repeat factor + outcome columns to match long format ---------------
  # pivot_longer produces rows ordered: all obs for draw 1, all obs for draw 2...
  # so simply repeating df n_draws times in the same order aligns perfectly
  factors_repeated <- df |>
    dplyr::select(dplyr::all_of(cols_to_tile)) |>
    dplyr::slice(rep(seq_len(nrow(df)), times = n_draws))

  # ---- bind columns directly — no join needed -----------------------------
  out <- dplyr::bind_cols(y_rep_long, factors_repeated) |>
    dplyr::select(
      .chain, .iteration, .draw,
      dplyr::all_of(cols_to_tile),
      y_rep
    ) |>
    # dplyr::select(-".obs_idx") |>  # drop the temporary obs index
    tidyr::pivot_longer(
      cols      = dplyr::all_of(c(outcome_name, "y_rep")),
      names_to  = "outcome_type",
      values_to = "outcome"
    ) |>
    dplyr::mutate(
      outcome_type = dplyr::recode(
        outcome_type,
        !!outcome_name := "observed",
        "y_rep"        := "posterior"
      )) |>
      dplyr::mutate(
        outcome_type = factor(outcome_type, levels = c("observed", "posterior")))
    if (!preserve_.draws)    
          out <- out |> select(-c(".draw", ".iteration", ".chain"))
  
    
        return(out)
  
}

#' Print a ppc_results object
#'
#' Displays a formatted textual summary of predictive check results.
#' Optionally renders all plots.
#'
#' @param x a `ppc_results` object
#' @param plot logical; if TRUE, all plots are rendered (default FALSE)
#' @param ... unused
#' @exportS3Method
print.ppc_results <- function(x, plot = FALSE, plot_minimal = TRUE, ...) {
  cat("<ppc_results>\n")

  # ---- Coverage -----------------------------------------------------------
  cat("\nInterval coverage calibration:\n")
  cat(sprintf("  %-6s  empirical  interpretation\n", "level"))
  cat(sprintf("  %-6s  ---------  --------------\n", "-----"))
  for (nm in names(x$coverage)) {
    cat(sprintf(
      "  %-6s  %-9.3f  %s\n",
      nm,
      x$coverage[[nm]],
      x$interpretation$coverage[[nm]]
    ))
  }

  # ---- Bayesian p-values --------------------------------------------------
  cat("\nBayesian p-values:\n")
  cat(sprintf("  %-8s  %s\n", "stat", "result"))
  cat(sprintf("  %-8s  %s\n", "-------", "------"))
  for (nm in names(x$bayes_pvalues)) {
    cat(sprintf(
      "  %-8s  %s\n",
      nm,
      x$interpretation$bayes_pvalues[[nm]]
    ))
  }

  # ---- Available plots ----------------------------------------------------
  cat("\nAvailable plots: $plots$",
      paste(names(x$plots), collapse = ", $plots$"), "\n")
  cat("  Use print(x, plots = TRUE) to render all, or x$plots$<name> individually.\n")

  # ---- Optional plot rendering --------------------------------------------
  if (plot) {
    for (nm in names(x$plots)) {
      if (!plot_minimal || (plot_minimal && (nm %in% c("density_overlay","mean_sd_2d")))) {
        print(x$plots[[nm]] + ggplot2::theme_minimal())
      } 
    }
  }

  invisible(x)
}


#' Print a ppc_results_grouped object
#'
#' Displays a formatted textual summary of predictive check results for each group,
#' then creates a grid of faceted plots (density_overlay and mean_sd_2d).
#'
#' @param x a `ppc_results_grouped` object
#' @param ... unused
#' @importFrom patchwork wrap_plots
#' @export
display <- function(x, ...) UseMethod("display")


#' @exportS3Method
display.ppc_results_grouped <- function(x, ...) {

  cat("<ppc_results_grouped>\n")
  cat(sprintf("Grouped by: %s\n", paste(x$by, collapse = ", ")))
  cat(sprintf("Number of groups: %d\n\n", length(x$group_results)))

  # ---- Print summary table ------------------------------------------------
  cat("Summary table:\n")
  print(x$summary)

  # ---- Print each group's results -----------------------------------------
  cat("\n--- Per-group diagnostics ---\n")
  for (gk in names(x$group_results)) {
    cat(sprintf("\n%s:\n", gk))
    print(x$group_results[[gk]], plot = FALSE)
  }

  # ---- Build faceted plots ------------------------------------------------
  
  # Extract density_overlay plots
  density_plots <- lapply(names(x$group_results), function(gk) {
    x$group_results[[gk]]$plots$density_overlay
  })
  names(density_plots) <- names(x$group_results)

  # Extract mean_sd_2d plots
  mean_sd_plots <- lapply(names(x$group_results), function(gk) {
    x$group_results[[gk]]$plots$mean_sd_2d
  })
  names(mean_sd_plots) <- names(x$group_results)

  # Determine grid layout based on number of .by variables
  n_groups <- length(x$group_results)

  if (length(x$by) == 1L) {
    # Single grouping variable: arrange in column
    nrow_grid <- n_groups
    ncol_grid <- 1L
  } else if (length(x$by) == 2L) {
    # Two grouping variables: first .by = columns, second .by = rows
    # Extract unique levels for each grouping variable
    by_data <- x$summary[, x$by, drop = FALSE]
    n_levels_by1 <- n_distinct(by_data[[x$by[1L]]])
    n_levels_by2 <- n_distinct(by_data[[x$by[2L]]])
    ncol_grid <- n_levels_by1
    nrow_grid <- n_levels_by2
  } else {
    # Fallback: single row
    nrow_grid <- 1L
    ncol_grid <- n_groups
  }

  # Create density_overlay grid
  cat("\n--- Density Overlay Grid ---\n")
  g_density <- patchwork::wrap_plots(
    density_plots,
    nrow = nrow_grid,
    ncol = ncol_grid
  ) +
    patchwork::plot_annotation(
      title = "Posterior Predictive Check — Density Overlay"
    )

  print(g_density)

  # Create mean_sd_2d grid
  cat("\n--- Mean-SD 2D Grid ---\n")
  g_mean_sd <- patchwork::wrap_plots(
    mean_sd_plots,
    nrow = nrow_grid,
    ncol = ncol_grid
  ) +
    patchwork::plot_annotation(
      title = "Posterior Predictive Check — Mean vs. SD"
    )
  print(g_mean_sd)
  

 invisible(x)  # ← Changed from return(list(...))
}

#' Convergence diagnostics for fitted stanmodel
#' Checks R-hat, effective sample size (ESS), divergences, treedepth, and energy diagnostics. Optionally renders traceplots.
#' @param x a fitted `stanmodel` object
#' @param params optional character vector of parameter names to check; if NULL, checks all
#' @param plot logical; if TRUE, renders traceplots for specified parameters (default FALSE)
#' @return a data frame summarizing convergence diagnostics for each parameter
#' @importFrom posterior summarise_draws subset_draws
#' @importFrom bayesplot mcmc_trace
#' @importFrom dplyr filter
#' @importFrom bayesplot mcmc_trace mcmc_rank_overlay mcmc_rank_hist
#' @export
diagnose <- function(x, ...) UseMethod("diagnose")


#' @exportS3Method
diagnose.stanmodel <- function(x, params = NULL, plot = FALSE) {
  if (!inherits(x, "stanmodel")) {
    stop("`x` must be a stanmodel object.")
  }

  if (is.null(x$fitted_model)) {
    stop("Model has not been fitted yet. Run `fit()` first.")
  } else {
    fit <- x$fitted_model
  }


    diag  <- fit$diagnostic_summary()
    div   <- diag$num_divergent
    td    <- diag$num_max_treedepth
    ebfmi <- diag$ebfmi

  
  cat("\n====================\n Convergence Summary\n====================\n")
  
  # default to scalar parameters only for plots — exclude y_rep and log_lik
  plot_params <- if (!is.null(params)) params else {
    x$variables[!grepl("^y_rep|^log_lik", x$variables)]
  }

  # Check if specified params exist in draws
  missing_params <- setdiff(plot_params, x$variables)
  if (length(missing_params) > 0) {
    warning(sprintf(
      "The following parameters were not found in the draws and will be ignored: %s",
      paste(missing_params, collapse = ", ")
    ))
    plot_params <- setdiff(plot_params, missing_params)
  }

  # Extract draws
  draws <- x$draws |> posterior::subset_draws(variable = plot_params) # filter to specified params if provided
  
  # Posterior package summary
  summ <- posterior::summarise_draws(draws)
  summ |> dplyr::select(variable, mean, sd, rhat, ess_bulk, ess_tail) |> print()
  
  # Rhat check
  cat("\n--- R-hat diagnostics ---\n")
  bad_rhat <- summ$variable[summ$rhat > 1.01]
  if (length(bad_rhat) == 0) {
    cat("All R-hat < 1.01 (good)\n")
  } else {
    cat("Parameters with R-hat > 1.01:\n")
    print(bad_rhat)
  }
  
  # ESS check
  cat("\n--- Effective Sample Size (ESS) ---\n")
  low_ess <- summ$variable[summ$ess_bulk < 400]
  if (length(low_ess) == 0) {
    cat("All ESS_bulk > 400 (good)\n")
  } else {
    cat("Parameters with ESS_bulk < 400:\n")
    print(low_ess)
  }
  
  # Divergences
  cat("\n--- Divergences ---\n")
  cat(sprintf("Total divergent transitions: %d\n", sum(div)))
  if (sum(div) > 0) {
    cat("Per chain:\n")
    for (i in seq_along(div)) cat(sprintf("  chain %d: %d\n", i, div[i]))
    cat("Action: increase `adapt_delta` toward 0.99 or reparameterise (non-centred)\n")
  }
  
  # Treedepth
  cat("\n--- Max treedepth ---\n")
  cat("Transitions hitting max treedepth:", td, "\n")
  
  # Energy Bayesian fraction of missing information (E-BFMI)
  cat("\n--- Energy diagnostics (E-BFMI) ---\n")
  print(ebfmi)
  if (any(ebfmi < 0.3)) {
    cat("Warning: E-BFMI < 0.3 indicates poor exploration\n")
  }
  
  # Optional traceplots + trankplots
  if (plot) {
    cat("\n--- Traceplots ---\n")
    g_trace <- bayesplot::mcmc_trace(draws) +
      ggplot2::labs(
        title    = "Traceplots",
        subtitle = "Chains should mix well (hairy caterpillar). Drifts or sticking indicate poor convergence."
      )
    print(g_trace)

    # ---- McElreath trankplots (Statistical Rethinking §9.5) ---------------
    # Draws are ranked across all chains combined; rank histograms per chain
    # are plotted. Interpretation:
    #   Uniform histogram → chain explored the full posterior (good mixing)
    #   U-shape           → chain stuck at tails (poor mixing)
    #   Humped shape      → chain stuck in centre (poor mixing)
    #   Chains diverge    → non-stationarity or poor initialisation
    cat("\n--- Trankplots (McElreath SR2 §9.5) ---\n")
    g_trank_overlay <- bayesplot::mcmc_rank_overlay(draws) +
      ggplot2::labs(
        title    = "Trankplots — rank overlay",
        subtitle = "Uniform, overlapping lines across chains = good mixing (McElreath SR2 §9.5)."
      )
    print(g_trank_overlay)

    g_trank_hist <- bayesplot::mcmc_rank_hist(draws) +
      ggplot2::labs(
        title    = "Trankplots — rank histograms",
        subtitle = "Uniform histograms per chain = good. U-shape or hump = poor mixing."
      )
    print(g_trank_hist)
  }
  
  invisible(summ)
}

# --- posterior parameters
#' Extract posterior draws for a specific parameter, attaching factor labels
#' @param x a fitted `stanmodel` object
#' @param parameter the name of the parameter to extract (e.g. "beta")
#' @param factors a named list mapping factor names to their levels, e.g. list
#'  condition = c("control", "treatment"), time = c("pre", "post"))
#' @param .widen_factor optional character name of a factor to pivot wider, creating one column per level of that factor
#' @return a tibble with columns for the parameter value and each factor, where factor levels are attached according to the indices in the parameter name.
#' @importFrom posterior as_draws_df subset_draws
#' @importFrom dplyr mutate select across
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_match
#' @importFrom purrr map map_chr
#' @importFrom rlang expr
#' @importFrom rlang !!
#' @export
ppc_by_parameter <- function(x, ...) UseMethod("ppc_by_parameter")

#' @exportS3Method
ppc_by_parameter.stanmodel <- function(x, parameter, factors, .widen_factor = NULL) {
  
  # check if draws is present
  if (is.null(x$draws)) {
    stop("No draws found in the model. Run `fit()` first.")
  }

  # check if `parameter` is present in draws
  if (!parameter %in% x$variables) {
    stop(sprintf("Parameter '%s' not found in draws.", parameter))
  }

  draws <- x$draws |>
          posterior::subset_draws(variable = parameter) |>
          posterior::as_draws_df() 
  
  # check if the number of factors in the draws df matches the length of the provided factors list
  n_factors_in_draws <- length(grep("\\[", names(draws)))

  # compute the product of the lengths of the factors to get the expected number of indices
  if (length(factors) == 0) {
    return(draws |> tibble())
  } else {
    n_expected_indices <- prod(sapply(factors, length))
  }

  if (n_factors_in_draws != n_expected_indices) {
    stop(sprintf(
      "Number of factors in draws (%d) does not match length of provided factors list (%d).",
      n_factors_in_draws, n_expected_indices
    ))
  }


  
  # build regex to capture all indices: [i,j,...] -> groups i, j, ...
  n_factors  <- length(factors)
  idx_pattern <- paste0(
    "\\[",
    paste(rep("(\\d+)", n_factors), collapse = ","),
    "\\]"
  )
  factor_names <- names(factors)

  result <- draws |>
    tidyr::pivot_longer(
      cols      = matches("\\["),
      names_to  = "parameter",
      values_to = "value"
    ) |>
    mutate(
      # extract indices as a list-column of integer vectors
      indices = map(
        str_match(parameter, idx_pattern)[, -1, drop = FALSE] |>
          (\(m) split(as.integer(m), row(m)))(),
        identity
      )
    ) 
   
    
   for (i in seq_along(factors)) {
      nm   <- names(factors)[[i]]
      lvls <- factors[[i]]

     # check that all indices for this factor are within the levels provided
     max_idx <- max(purrr::map_int(result$indices, \(idx) idx[[i]]))
      if (max_idx > length(lvls)) {
        stop(sprintf(
          "Factor '%s' (position %d) has %d levels but Stan index %d was found. Ensure `factors` are ordered to match Stan's index dimensions.",
          nm, i, length(lvls), max_idx
        ))
      }
      result <- dplyr::mutate(result, !!nm := purrr::map_chr(indices, \(idx) lvls[[idx[[i]]]]))
   }

  
  
    result <- 
      result |>
      dplyr::mutate(dplyr::across(
        dplyr::all_of(factor_names),
        \(col) factor(col, levels = factors[[dplyr::cur_column()]])
      )) |>
      dplyr::select(-parameter, -indices)
  
    if (!is.null(.widen_factor)) {
      # reshape to wide format with one column per factor combination
      id_factors <- setdiff(factor_names, .widen_factor)
      result <- result |>
        tidyr::pivot_wider(
          id_cols = c(".chain", ".iteration", ".draw", dplyr::all_of(id_factors)  ),
          names_from = dplyr::all_of(.widen_factor),
          values_from = "value"
        )
    }
    
  return(result)

}

#' Scenario analysis from posterior parameter draws
#'
#' Extracts wide-format posterior draws via `ppc_by_parameter()`, classifies
#' each draw into one of four 2-way scenarios, computes scenario probabilities,
#' and optionally plots the result.
#'
#' @param x              a fitted `stanmodel` object
#' @param parameter      character; name of the Stan parameter to extract (e.g. "beta")
#' @param factors        named list mapping factor names to their levels, passed to
#'                       `ppc_by_parameter()` (e.g. `list(direction = c("FLEX","EXT"), finger = c("T","I","M","R","L"))`)
#' @param scenario_factor character; name of the factor in `factors` to use for constructing
#'                       2-way scenarios — its first two levels become the outcome columns
#'                       (e.g. "direction" with levels c("FLEX", "EXT"))
#' @param scenario_levels optional character vector of length 2; custom outcome levels to use
#'                       instead of the first two levels of `scenario_factor`. If NULL, 
#'                       defaults to the first two levels.
#' @param scenario_labels optional character vector of length 2; human-readable labels for
#'                       the scenario levels (e.g. c("Flexion", "Extension")). If NULL,
#'                       defaults to `scenario_levels`.
#' @param effect_positive logical; if TRUE, positive values count as an effect; if FALSE,
#'                       negative values count as an effect (default TRUE)
#' @param plot            logical; if TRUE, render the 2-way scenario scatter plot (default FALSE)
#' @param facet_formula   formula or NULL; passed to `.plot_2way_scenarios()` for faceting
#'                       (e.g. `~ finger`). If NULL and one grouping variable exists,
#'                       a formula is constructed automatically.
#' @param title           optional character; plot title
#'
#' @return A list with:
#'   \describe{
#'     \item{`draws`}{wide-format draws tibble from `ppc_by_parameter()` with one column
#'       per level of `scenario_factor`}
#'     \item{`scenario_draws`}{draws with an additional `scenario` column classifying each
#'       draw into one of four 2-way scenarios}
#'     \item{`scenario_grid`}{tibble of scenario probabilities per grouping variable cell;
#'       format depends on `wide` parameter}
#'     \item{`plot`}{ggplot object (only present when `plot = TRUE`)}
#'   }
#'
#' @export
ppc_scenarios <- function(x, ...) UseMethod("ppc_scenarios")

#' @exportS3Method
ppc_scenarios.stanmodel <- function(x,
                                    parameter,
                                    factors,
                                    scenario_factor,
                                    scenario_levels         = NULL,
                                    scenario_labels         = NULL,
                                    effect_positive = TRUE,
                                    plot           = FALSE,
                                    facet_formula  = NULL,
                                    multiple_plots = NULL,
                                    title          = NULL) {

  stopifnot(
    is.character(parameter), length(parameter) == 1L,
    is.list(factors), !is.null(names(factors))
  )

  if (is.null(scenario_levels)){
    stopifnot(scenario_factor %in% names(factors) && length(factors[[scenario_factor]]) >= 2L)
    col1 <- factors[[scenario_factor]][1]
    col2 <- factors[[scenario_factor]][2]
    if (length(factors[[scenario_factor]]) > 2L) {
      warning(sprintf(
        "More than 2 levels found for scenario_factor '%s'. Only the first two levels will be used: '%s' and '%s'.",
        scenario_factor, col1, col2
      ))
    }
  }

  if (is.null(scenario_labels)){
    lv1 <- col1
    lv2 <- col2
  } else {
    stopifnot(length(scenario_labels) == 2L)
    lv1 <- scenario_labels[1]
    lv2 <- scenario_labels[2]
  }

  # 1. Extract wide draws — one column per level of widen_factor
  draws_wide <- ppc_by_parameter(x,
                                 parameter    = parameter,
                                 factors      = factors,
                                 .widen_factor = scenario_factor)
  
  
  # 2. Classify draws into 2-way scenarios
  scenario_draws <- draws_wide |>
    .define_scenarios(col1          = col1,
                      col2          = col2,
                      lv1           = lv1,
                      lv2           = lv2,
                      effect_positive = effect_positive)

  # 3. Compute scenario grid
  # by_vars defaults to all factor columns except widen_factor and draw metadata

  by_vars <- setdiff(names(factors), scenario_factor)
  

  scenario_grid <- scenario_draws |>
    .compute_scenario_grid(by_vars = by_vars, wide = TRUE)

  # 4. Auto-construct facet_formula from by_vars if not supplied
  if (is.null(facet_formula)) {
    facet_formula <- if (length(by_vars) == 1L) {
      stats::as.formula(paste("~", by_vars[1]))
    } else if (length(by_vars) == 2L) {
      stats::as.formula(paste(by_vars[2], "~", by_vars[1]))
    } else {
      NULL
    }
  }

  # 5. Optionally plot
  out <- list(
    draws          = draws_wide,
    scenario_draws = scenario_draws,
    scenario_grid  = scenario_grid
  )

  if (is.null(multiple_plots) || !multiple_plots %in% names(factors)) {
    out$plot <- .plot_2way_scenarios(
        scenario_draws = scenario_draws,
        col1           = col1,
        col2           = col2,
        lv1            = lv1,
        lv2            = lv2,
        facet_formula  = facet_formula,
        title          = title
      )
    } else 
    {
    plot_levels <- scenario_draws |> pull(!!sym(multiple_plots)) |> unique()
    out$plot <- list()
    for (plot_lvl in plot_levels){    
      scenario_draws_filtered <-
        scenario_draws |> filter(!!sym(multiple_plots) == plot_lvl)
      out$plot[[plot_lvl]] <- .plot_2way_scenarios(
          scenario_draws = scenario_draws_filtered,
          col1           = col1,
          col2           = col2,
          lv1            = lv1,
          lv2            = lv2,
          facet_formula  = facet_formula,
          title          = paste0(title, " - ", plot_lvl)
        )
    }
  }
  if (plot) print(out$plot)
  

  out
  }