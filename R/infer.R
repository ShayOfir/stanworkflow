#infer.R

# ---- outcome_spec S3 class --------------------------------------------------

#' Create an outcome specification for posterior contrast extraction
#'
#' Decouples *what to extract from draws* from *how to compare groups*.
#' Pass to `contrast(..., outcome = outcome_spec(...))`.
#'
#' @param gq_var   character; name of the GQ variable to extract from draws.
#'   Standard values: "y_rep" (default), "mu_marginal", "p_respond", "mu_y".
#' @param filter_fn optional function(draws_long_tibble) -> tibble applied after
#'   extraction but before contrasting. Use to drop sentinel rows, transform
#'   the outcome, etc. The tibble has columns: .chain, .iteration, .draw,
#'   .obs_idx, y, <feature cols>.
#'   Default NULL (no filter).
#' @param transform_fn optional function(numeric vector) -> numeric vector
#'   applied element-wise to the extracted values (e.g. exp for log-scale models).
#'   Default NULL (identity).
#' @param label optional human-readable label for this quantity (used in plots
#'   and summaries). Default NULL (derived from gq_var).
#' @return an object of class `outcome_spec`
#' @export
outcome_spec <- function(gq_var      = "y_rep",
                         filter_fn   = NULL,
                         transform_fn = NULL,
                         label       = NULL) {
  stopifnot(is.character(gq_var), length(gq_var) == 1L)
  if (!is.null(filter_fn))    stopifnot(is.function(filter_fn))
  if (!is.null(transform_fn)) stopifnot(is.function(transform_fn))

  structure(
    list(
      gq_var       = gq_var,
      filter_fn    = filter_fn,
      transform_fn = transform_fn,
      label        = label %||% gq_var
    ),
    class = "outcome_spec"
  )
}

#' @exportS3Method
print.outcome_spec <- function(x, ...) {
  cat("<outcome_spec>\n")
  cat("  gq_var      :", x$gq_var, "\n")
  cat("  transform   :", if (is.null(x$transform_fn)) "none" else "custom fn", "\n")
  cat("  filter      :", if (is.null(x$filter_fn))    "none" else "custom fn", "\n")
  cat("  label       :", x$label, "\n")
  invisible(x)
}

# ---- built-in outcome_spec constructors ------------------------------------

#' Standard y_rep outcome (default; backwards-compatible)
#' @param filter_sentinel if TRUE, drops rows where y <= sentinel (default FALSE)
#' @param sentinel_value  numeric sentinel to drop (default 0.001)
#' @export
os_yrep <- function(filter_sentinel = FALSE, sentinel_value = 0.001) {
  filter_fn <- if (filter_sentinel) {
    sv <- sentinel_value
    function(df) dplyr::filter(df, y > sv)
  } else {
    NULL
  }
  outcome_spec(gq_var = "y_rep", filter_fn = filter_fn, label = "y_rep")
}

#' Marginal expected outcome: E[y] = P(respond) * E y respond
#' Primary inferential target for hurdle models (Bürkner 2017, Woody et al. 2021)
#' @export
os_marginal <- function() {
  outcome_spec(gq_var = "mu_marginal", label = "E[y] (marginal)")
}

#' Extensive margin: P(respond)
#' @export
os_extensive <- function() {
  outcome_spec(gq_var = "p_respond", label = "P(respond)")
}

#' Intensive margin: E y respond
#' @export
os_intensive <- function() {
  outcome_spec(gq_var = "mu_y", label = "E y respond")
}

#' Compound hurdle outcome_spec: returns all three margins
#' infer_contrasts() detects this class and runs all three in one pass.
#' @export
os_hurdle <- function() {
  structure(
    list(
      marginal  = os_marginal(),
      extensive = os_extensive(),
      intensive = os_intensive(),
      label     = "hurdle (marginal + extensive + intensive)"
    ),
    class = c("outcome_spec_hurdle", "outcome_spec")
  )
}

#' @exportS3Method
print.outcome_spec_hurdle <- function(x, ...) {
  cat("<outcome_spec_hurdle>\n")
  cat("  margins: marginal (mu_marginal), extensive (p_respond),",
      "intensive (mu_y)\n")
  invisible(x)
}



# ---- contrast S3 class ------------------------------------------------------




#' Create a contrast specification object
#'
#' @param contrast_var  character scalar — the variable whose levels are contrasted
#' @param grouping_vars character vector of variables to group by
#' @param contrast_levels for categorical: length-2 character vector c("hi","lo").
#'   For continuous with quantiles: length-2 numeric vector of probabilities c(hi_q, lo_q).
#'   For continuous with absolute thresholds: length-2 numeric vector c(hi_threshold, lo_threshold).
#'   Default NULL uses first two factor levels (categorical) or c(0.75, 0.25) (continuous).
#' @param effect_positive logical; if TRUE, higher values of contrast_var are considered an effect in the "correct direction". Default TRUE.
#' @param continuous logical; if TRUE treat contrast_var as continuous. Default FALSE.
#' @param absolute_thresholds logical; if TRUE and continuous = TRUE, contrast_levels are
#'   treated as absolute cut-off values rather than quantile probabilities. Default FALSE.
#' @param filter_draws a quoted or unquoted dplyr filter expression to subset the
#'   draws before computing the contrast. E.g. quote(group == "stroke") or
#'   quote(CSTLL > 0.5). Default NULL (no filtering).
#' @param two_way character scalar or NULL; if not NULL, specifies the variable whose
#'   two levels define the two outcomes (e.g. "direction"). The contrast is then
#'   computed as a 4-scenario analysis: affects both, first only, second only, neither.
#'   Requires `continuous = FALSE`. Default NULL.
#' @param two_way_levels optional character vector of exactly 2 levels from `two_way`
#'   variable to use as the two outcomes. If NULL, uses first 2 levels automatically.
#' @return an object of class `contrast`
#' @export
# ---- contrast() validation: remove the continuous + two_way guard ----------
contrast <- function(
                     contrast_var,
                     grouping_vars       = NULL,
                     contrast_levels     = NULL,
                     effect_positive      = TRUE,
                     continuous          = FALSE,
                     absolute_thresholds = FALSE,
                     filter_draws        = NULL,
                     two_way             = NULL,
                     two_way_levels      = NULL,
                     outcome             = NULL) {
  if (!is.null(grouping_vars)) {
    stopifnot(is.character(grouping_vars), length(grouping_vars) >= 1)
  }
  stopifnot(is.character(contrast_var), length(contrast_var) == 1)

  if (continuous) {
    if (!is.null(contrast_levels)) {
      stopifnot(is.numeric(contrast_levels), length(contrast_levels) == 2)
      if (!absolute_thresholds) {
        stopifnot(all(contrast_levels >= 0 & contrast_levels <= 1))
      }
    }
  } else {
    if (absolute_thresholds)
      stop("absolute_thresholds = TRUE requires continuous = TRUE.")
    stopifnot(is.null(contrast_levels) ||
              (is.character(contrast_levels) && length(contrast_levels) == 2))
  }

  # <<< CHANGED: removed the guard that blocked continuous + two_way
  # two_way + continuous is now supported in counterfactual mode via
  # .apply_contrast_split(), which handles both cases uniformly.

  if (!is.null(outcome)) {
    stopifnot(inherits(outcome, "outcome_spec"))
  }

  structure(
    list(
      grouping_vars       = grouping_vars,
      contrast_var        = contrast_var,
      contrast_levels     = contrast_levels,
      effect_positive     = effect_positive,
      continuous          = continuous,
      absolute_thresholds = absolute_thresholds,
      filter_draws        = filter_draws,
      two_way             = two_way,
      two_way_levels      = two_way_levels,
      outcome             = outcome,
      # internal numeric_array-related meta-data:
      .array_source       = NULL,
      .array_cols         = NULL

    ),
    class = "contrast"
  )
}

#' Generate a human-readable contrast label
#'
#' @param x a `contrast` object
#'
#' @return character scalar with contrast label
#' @export
get_label <- function(x) {
  UseMethod("get_label")
}

#' @exportS3Method
get_label.contrast <- function(x) {
  # <<< ADDED: two_way label branch
  if (!is.null(x$two_way)) {
    lvls <- x$two_way_levels %||% "<auto>"
    return(sprintf("2-way scenarios on '%s' [%s]", x$two_way, paste(lvls, collapse = " / ")))
  }
  # <<< END ADDITION
  if (is.null(x$contrast_levels)) {
    if (x$continuous) {
      sprintf("75%% - 25%% (default quantiles)")
    } else {
      sprintf("<auto: first two levels>")
    }
  } else {
    if (x$continuous) {
      if (x$absolute_thresholds) {
        hi_thresh <- x$contrast_levels[1]
        lo_thresh <- x$contrast_levels[2]
        sprintf("(≤%g) - (>%g)", hi_thresh, lo_thresh)
      } else {
        hi_q <- x$contrast_levels[1]
        lo_q <- x$contrast_levels[2]
        sprintf("%d%% - %d%%", round(hi_q * 100), round(lo_q * 100))
      }
    } else {
      sprintf("%s - %s", x$contrast_levels[1], x$contrast_levels[2])
    }
  }
}

#' @exportS3Method
print.contrast <- function(x, ...) {
  cat("<contrast>\n")
  # <<< CHANGED: two_way is now a character variable name, not a logical
  if (!is.null(x$two_way)) {
    cat("  type           : 2-way scenarios\n")
    cat("  two_way        :", x$two_way, "\n")
    cat("  two_way_levels :", paste(x$two_way_levels %||% "<auto>", collapse = ", "), "\n")
    cat("  contrast_var   :", x$contrast_var, "\n")
    cat("  contrast_levels:", if (is.null(x$contrast_levels)) "<default>" else paste(x$contrast_levels, collapse = ", "), "\n")
    cat("  continuous     :", x$continuous, "\n")
    cat("  absolute_thresholds:", x$absolute_thresholds, "\n")
  } else {
  # <<< END CHANGE
    cat("  contrast_var   :", x$contrast_var, "\n")
    if (x$continuous) {
      if (x$absolute_thresholds) {
        thresholds <- if (!is.null(x$contrast_levels)) x$contrast_levels else c(NA, NA)
        cat("  type           : continuous (absolute thresholds)\n")
        cat("  hi threshold   : <=", thresholds[1], "\n")
        cat("  lo threshold   :  >", thresholds[2], "\n")
      } else {
        probs <- if (!is.null(x$contrast_levels)) x$contrast_levels else c(0.75, 0.25)
        cat("  type           : continuous (quantile split)\n")
        cat("  hi quantile    :", probs[1], "\n")
        cat("  lo quantile    :", probs[2], "\n")
      }
    } else {
      lev <- if (!is.null(x$contrast_levels))
               paste(x$contrast_levels, collapse = " - ")
             else
               "<auto: first two levels>"
      cat("  type           : categorical\n")
      cat("  levels         :", lev, "\n")
    }
  }
  cat("  grouping_vars  :", 
      if (is.null(x$grouping_vars)) "<none>" else paste(x$grouping_vars, collapse = ", "), "\n")
  if (!is.null(x$filter_draws)) {
    cat("  filter_draws   :", deparse(x$filter_draws), "\n")

  }
   cat("  outcome        :",
      if (is.null(x$outcome)) "y_rep (default)" else x$outcome$label, "\n")
  invisible(x)
}
#' @keywords internal
.restore_training_dimensions <- function(stan_data, original_stan_data, 
                                         new_data, data,
                                         cond_df_use, cond_df) {
  # Restore n_* scalars from training data
  n_vars <- grep("^n_", names(stan_data), value = TRUE)
  for (nv in n_vars) {
    if (!is.null(original_stan_data[[nv]]) &&
        !identical(stan_data[[nv]], original_stan_data[[nv]])) {
      stan_data[[nv]] <- original_stan_data[[nv]]
    }
  }

  # Restore *_id index arrays: re-encode using original factor levels
  # so that e.g. coil_id = 1 still means "fig8" as in the fitted model.
  id_vars <- grep("_id$", names(stan_data), value = TRUE)
  for (iv in id_vars) {
    col <- sub("_id$", "", iv)
    if (col %in% names(new_data) && col %in% names(data)) {
      original_levels <- levels(factor(data[[col]]))
      stan_data[[iv]] <- as.integer(factor(new_data[[col]],
                                           levels = original_levels))
    }
  }

  # Restore *_id_cond_array index arrays from condition-level data:
  # re-encode using original factor levels from training cond_df
  id_cond_vars <- grep("_id_cond_array$", names(stan_data), value = TRUE)
  for (iv in id_cond_vars) {
    col <- sub("_id_cond_array$", "", iv)
    if (col %in% names(cond_df_use) && col %in% names(cond_df)) {
      original_levels <- levels(factor(cond_df[[col]]))
      stan_data[[iv]] <- as.integer(factor(cond_df_use[[col]],
                                           levels = original_levels))
    }
  }

  stan_data
}

# ---- inference S3 class =====================================================

#' Create an inference result object
#'
#' Wraps the output of infer_contrasts() with metadata to enable consistent
#' printing and plotting across 1-way and 2-way contrast structures.
#'
#' @param result list output from infer_contrasts()
#' @param contrast the contrast object used to generate result
#' @return object of class "inference"
#' @keywords internal
as_inference <- function(result, contrast) {
  structure(
    result,
    class = c("inference", class(result))
  )
}

#' @exportS3Method
print.inference <- function(x, ...) {
  
  # Detect contrast type
  has_two_way    <- !is.null(x$contrast$two_way)
  has_scenarios  <- !is.null(x$scenario_grid)
  has_grid       <- !is.null(x$grid_summary)
  
  cat("\n<inference>\n")
  cat("  Type: ", if (has_two_way) "2-way (scenarios)" else "1-way", "\n")
  cat("  Contrast variable:", x$contrast$contrast_var, "\n")
  cat("  Grouping vars:", 
      if (is.null(x$contrast$grouping_vars)) "<none>" 
      else paste(x$contrast$grouping_vars, collapse = ", "), "\n")
  
  if (has_two_way) {
    cat("  Two-way outcome:", x$contrast$two_way, "\n")
    cat("  Outcome levels:", paste(x$lv1, x$lv2, sep = " / "), "\n\n")
    
    # Print scenario probabilities
    cat("Scenario Probabilities:\n")
    print(x$scenario_grid)
    
    # Generate and display 2-way plot
    cat("\nGenerating 2-way scenario plot...\n")
    p <- .plot_inference_2way(x)
    print(p)
    
  } else if (has_grid) {
    # Counterfactual continuous
    cat("  Counterfactual grid sweep (continuous)\n")
    cat("  Grid points:", length(unique(x$grid_summary$contrast_var_value)), "\n\n")
    
    cat("Grid Summary:\n")
    print(x$grid_summary)
    
    # Generate and display continuous grid plot
    cat("\nGenerating counterfactual grid plot...\n")
    p <- .plot_inference_grid(x)
    print(p)
    
  } else {
    # 1-way descriptive or categorical counterfactual
    cat("  Method: ", 
        if (!is.null(x$grid_posterior)) "counterfactual (categorical)" 
        else "descriptive", "\n\n")
    
    if (!is.null(x$summary)) {
      cat("Summary:\n")
      print(x$summary)
    }
    
    # Generate and display 1-way plot
    cat("\nGenerating 1-way contrast plot...\n")
    p <- .plot_inference_1way(x)
    print(p)
  }
  
  invisible(x)
}

# ---- inference plotting helpers =============================================

#' Plot 1-way contrast (descriptive or categorical counterfactual)
#'
#' Generates a density plot of posterior contrasts, faceted by grouping variables.
#'
#' @param x inference object
#' @keywords internal
.plot_inference_1way <- function(x) {
  
  posterior_data <- if (!is.null(x$posterior)) x$posterior 
                   else if (!is.null(x$grid_posterior)) x$grid_posterior
                   else stop("No posterior data found in inference object")
  
  grouping_vars <- x$contrast$grouping_vars %||% character()
  
  # Determine faceting strategy
  facet_formula <- if (length(grouping_vars) == 0) {
    NULL
  } else if (length(grouping_vars) == 1) {
    stats::as.formula(paste(grouping_vars[1], "~ ."))
  } else if (length(grouping_vars) == 2) {
    stats::as.formula(paste(grouping_vars[2], "~", grouping_vars[1]))
  } else {
    stats::as.formula(paste(grouping_vars[1], "~", grouping_vars[2]))
  }
  
  p <- ggplot2::ggplot(
    posterior_data,
    ggplot2::aes(x = contrast, fill = .data[[grouping_vars[1] %||% "contrast"]])
  ) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = sprintf("Posterior distribution: %s - %s",
                      x$contrast$contrast_levels[2] %||% "hi",
                      x$contrast$contrast_levels[1] %||% "lo"),
      x = "Contrast",
      y = "Density",
      fill = grouping_vars[1] %||% ""
    ) +
    ggplot2::theme(legend.position = "bottom")
  
  if (!is.null(facet_formula)) {
    p <- p + ggplot2::facet_grid(facet_formula)
  }
  
  p
}

#' Plot 2-way scenario contrast (generic scatter plot)
#'
#' Low-level plotting function for 2-way scenario draws. Plots one outcome
#' against another, colored by scenario, with optional faceting.
#'
#' @param scenario_draws tibble with columns: col1, col2, scenario, and optional faceting columns
#' @param col1 character; name of column for x-axis (e.g., "contrast_FLEX")
#' @param col2 character; name of column for y-axis (e.g., "contrast_EXT")
#' @param lv1 character; label for outcome 1 (e.g., "Flexion")
#' @param lv2 character; label for outcome 2 (e.g., "Extension")
#' @param facet_formula formula object or NULL; faceting specification.
#'   Pass as `~ var1` for single facet, `var2 ~ var1` for grid, or NULL for no faceting.
#' @param title optional character; plot title
#'
#' @return ggplot object
#' @keywords internal
.plot_2way_scenarios <- function(scenario_draws,
                                 col1,
                                 col2,
                                 lv1,
                                 lv2,
                                 facet_formula = NULL,
                                 title = NULL) {
  
  p <- ggplot2::ggplot(
    scenario_draws |> dplyr::mutate(scenario = factor(scenario)),
    ggplot2::aes(x = .data[[col1]], y = .data[[col2]], color = scenario)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title %||% sprintf("2-way scenario analysis: %s vs %s", lv1, lv2),
      x = paste(lv1),
      y = paste(lv2),
      color = "Scenario"
    ) +
    ggplot2::theme(legend.position = "bottom")
  
  if (!is.null(facet_formula)) {
    p <- p + ggplot2::facet_grid(facet_formula)
  }
  
  p
}

#' Plot 2-way scenario contrast (inference-aware wrapper)
#'
#' Extracts scenario draws and metadata from inference object, determines
#' appropriate faceting, and delegates to .plot_2way_scenarios() for rendering.
#'
#' @param x inference object with contrast, scenario_draws, lv1, lv2 fields
#' @keywords internal
.plot_inference_2way <- function(x) {
  
  outer_vars <- x$contrast$grouping_vars %||% character()
  
  lv1 <- x$lv1
  lv2 <- x$lv2
  col1 <- paste0("contrast_", lv1)
  col2 <- paste0("contrast_", lv2)
  
  # Determine faceting based on outer_vars and grid coverage
  facet_formula <- if (length(outer_vars) == 0 && length(unique(x$scenario_draws$contrast_var_value)) == 1) {
    NULL
  } else if (length(outer_vars) == 1 && length(unique(x$scenario_draws$contrast_var_value)) == 1) {
    stats::as.formula(paste(" ~", outer_vars[1]))
  } else if (length(outer_vars) == 2 && length(unique(x$scenario_draws$contrast_var_value)) == 1) {
    stats::as.formula(paste(outer_vars[2], "~", outer_vars[1]))
  } else if (length(outer_vars) == 1 && length(unique(x$scenario_draws$contrast_var_value)) > 1) {
    stats::as.formula(paste("contrast_var_value ~", outer_vars[1]))
  } else if (length(outer_vars) == 2 && length(unique(x$scenario_draws$contrast_var_value)) > 1) {
    stats::as.formula(paste("contrast_var_value ~", outer_vars[1], "+", outer_vars[2]))
  } else {
    NULL
  }
  
  # Delegate to generic plotting function
  .plot_2way_scenarios(
    scenario_draws = x$scenario_draws,
    col1 = col1,
    col2 = col2,
    lv1 = paste0('Contrast ', lv1),
    lv2 = paste0('Contrast ', lv2),
    facet_formula = facet_formula
  )
}

# #' Plot 2-way scenario contrast
# #'
# #' Generates a scatter plot of scenario draws (or bar plot of scenario probs)
# #' across two outcome dimensions, faceted by grouping variables.
# #'
# #' @param x inference object
# #' @keywords internal
# .plot_inference_2way <- function(x) {
  
#   outer_vars <- x$contrast$grouping_vars %||% character()

  
#   lv1 <- x$lv1
#   lv2 <- x$lv2
#   col1 <- paste0("contrast_", lv1)
#   col2 <- paste0("contrast_", lv2)
  
#   # Determine faceting
#   facet_formula <- if (length(outer_vars) == 0 && length(unique(x$scenario_draws$contrast_var_value)) == 1) {
#     NULL
#   } else if (length(outer_vars) == 1 && length(unique(x$scenario_draws$contrast_var_value)) == 1) {
#     stats::as.formula(paste(" ~", outer_vars[1]))
#   } else if (length(outer_vars) == 2 && length(unique(x$scenario_draws$contrast_var_value)) == 1) {
#     stats::as.formula(paste(outer_vars[2], "~", outer_vars[1]))
#   } else if (length(outer_vars) == 1 && length(unique(x$scenario_draws$contrast_var_value)) > 1) {
#     stats::as.formula(paste("contrast_var_value ~", outer_vars[1]))
#   } else if (length(outer_vars) == 2 && length(unique(x$scenario_draws$contrast_var_value)) > 1) {
#     stats::as.formula(paste("contrast_var_value ~", outer_vars[1], " * ", outer_vars[2]))
#   } else {
#     NULL
#   }
  
#   p <- ggplot2::ggplot(
#     x$scenario_draws |> dplyr::mutate(scenario = factor(scenario)),
#     ggplot2::aes(x = .data[[col1]], y = .data[[col2]], color = scenario)
#   ) +
#     ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#     ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#     ggplot2::geom_point(alpha = 0.5) +
#     ggplot2::theme_bw() +
#     ggplot2::labs(
#       title = sprintf("2-way scenario analysis: %s vs %s", lv1, lv2),
#       x = paste("Contrast", lv1),
#       y = paste("Contrast", lv2),
#       color = "Scenario"
#     ) +
#     ggplot2::theme(legend.position = "bottom")
  
#   if (!is.null(facet_formula)) {
#     p <- p + ggplot2::facet_grid(facet_formula)
#   }
  
#   p
# }

#' Plot counterfactual grid sweep (continuous contrast_var)
#'
#' Generates a line plot of scenario probabilities or contrast means
#' across the contrast variable grid.
#'
#' @param x inference object
#' @keywords internal
.plot_inference_grid <- function(x) {
  
  outer_vars <- setdiff(x$contrast$grouping_vars %||% character(),
                       x$contrast$two_way)
  
  # For scenario_grid (2-way continuous)
  if (!is.null(x$scenario_grid)) {
    facet_formula <- if (length(outer_vars) == 0) {
      NULL
    } else if (length(outer_vars) == 1) {
      stats::as.formula(paste("~", outer_vars[1]))
    } else {
      stats::as.formula(paste(outer_vars[2], "~", outer_vars[1]))
    }
    
    p <- ggplot2::ggplot(
      x$scenario_grid,
      ggplot2::aes(x = contrast_var_value, y = probability, color = scenario)
    ) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Scenario probabilities across contrast variable",
        x = x$contrast$contrast_var,
        y = "Probability",
        color = "Scenario"
      ) +
      ggplot2::theme(legend.position = "bottom")
    
    if (!is.null(facet_formula)) {
      p <- p + ggplot2::facet_wrap(facet_formula)
    }
    
  } else if (!is.null(x$grid_summary)) {
    # For grid_summary (1-way continuous)
    facet_formula <- if (length(outer_vars) == 0) {
      NULL
    } else if (length(outer_vars) == 1) {
      stats::as.formula(paste("~", outer_vars[1]))
    } else {
      stats::as.formula(paste(outer_vars[2], "~", outer_vars[1]))
    }
    
    p <- ggplot2::ggplot(
      x$grid_summary,
      ggplot2::aes(x = contrast_var_value, y = mean)
    ) +
      ggplot2::geom_line(color = "blue", linewidth = 1) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lo, ymax = hi),
        alpha = 0.2, fill = "blue"
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Posterior contrast across contrast variable",
        x = x$contrast$contrast_var,
        y = "Mean contrast",
        subtitle = "Band shows 95% credible interval"
      ) +
      ggplot2::theme(legend.position = "bottom")
    
    if (!is.null(facet_formula)) {
      p <- p + ggplot2::facet_wrap(facet_formula)
    }
  }
  
  p
}

# >>> NEW: helper to parse a virtual numeric_array column name
# Returns list(array = "roi", col = "CST_L") or NULL if not a virtual name.
#' @keywords internal
.parse_virtual_contrast_var <- function(contrast_var) {
  m <- regmatches(contrast_var, regexec("^([A-Za-z_][A-Za-z0-9_.]*?)\\[(.+)\\]$", contrast_var))[[1]]
  if (length(m) == 3L) list(array = m[2], col = m[3]) else NULL
}
# >>> END NEW


# ---- internal: simulate_at

#' Simulate posterior predictions at a fixed covariate profile
#'
#' @param model a fitted `stanmodel` object
#' @param value the value to assign to `contrast_var`
#' @param contrast_var name of the variable to set
#' @param outcome_var name of the outcome variable to extract
#' @param data synthetic dataframe (trial-level)
#' @param extra named list passed through to `get_stan_data()` 
#' @param cond_df optional condition-level dataframe
#' @keywords internal


.simulate_at <- function(model,
                         value,
                         contrast_var,
                         outcome_var = 'y_rep',
                         data = NULL,
                         extra = NULL,
                         cond_df = NULL,
                         n_draws = NULL) {
  if (is.null(data))   data   <- model$data
  if (is.null(cond_df)) cond_df <- model$cond_df
  if (is.null(extra))  extra  <- model$extra

  new_data <- data

  # >>> NEW: intercept virtual numeric_array contrast vars.
  # "roi[CST_L]" means: set extra$roi[, "CST_L"] = value for all rows,
  # zero out all other columns in extra$roi (counterfactual = isolate one column).
  # Also set the virtual column in new_data so downstream feature-column
  # extraction in .to_long_draws() finds it.
  parsed <- .parse_virtual_contrast_var(contrast_var)

  if (!is.null(parsed)) {
    arr_name <- parsed$array
    col_name <- parsed$col

    if (is.null(extra[[arr_name]]))
      stop(sprintf(".simulate_at: numeric_array '%s' not found in extra.", arr_name))
    if (!col_name %in% colnames(extra[[arr_name]]))
      stop(sprintf(".simulate_at: column '%s' not found in numeric_array '%s'.", col_name, arr_name))

    # Modify a copy of the matrix — zero all columns, set target column to value
    arr_cf <- extra[[arr_name]]
    arr_cf[, ]         <- 0          # zero all ROIs
    arr_cf[, col_name] <- value      # set the single counterfactual ROI
    extra[[arr_name]]  <- arr_cf     # replace in local extra copy

    # Keep virtual column in new_data for feature attachment in .to_long_draws()
    new_data[[contrast_var]] <- value

  } else {
    # >>> END NEW — original path for non-virtual contrast vars
    if (!is.null(contrast_var) && !is.null(value))
      new_data[[contrast_var]] <- value
  }

  # cond_df counterfactual (only applies to non-virtual vars)
  cond_df_use <- cond_df
  if (!is.null(cond_df_use) && !is.null(contrast_var) && !is.null(value) && is.null(parsed)) {
    cond_df_use[[contrast_var]] <- value
  }

  stan_data <- get_stan_data(
    model,
    df     = new_data,
    extra  = extra,       # <<< now contains modified arr_cf
    cond_df = cond_df_use
  )

  stan_data <- .restore_training_dimensions(
    stan_data, model$stan_data, new_data, data, cond_df_use, cond_df
  )

  prior_vals <- Filter(Negate(is.na), model$priors)
  stan_data  <- modifyList(stan_data, prior_vals)
  stan_data$skip_likelihood <- 0L

  na_fields <- names(Filter(anyNA, stan_data))
  if (length(na_fields) > 0)
    stop("stan_data has NA in: ", paste(na_fields, collapse = ", "))
  
  gq_draws <- model |> gq(new_data = stan_data, n_draws = n_draws)
  gq_draws |> subset_draws(outcome_var) |> as_draws_matrix()
}

#' Extract and reshape draws for a single outcome_spec
#'
#' Core extraction worker. Replaces the hard-coded y_rep path in
#' .get_yrep_draws. Called once per outcome_spec in the compound hurdle case.
#'
#' @param model   fitted stanmodel
#' @param data    data frame (nrow = N)
#' @param spec    an `outcome_spec` object
#' @param vars    character vector of feature columns to attach
#' @keywords internal
.extract_spec_draws <- function(model, data, spec, vars) {

  gq_var <- spec$gq_var

  if (!gq_var %in% model$variables) {
    stop(sprintf(
      "outcome_spec requests GQ variable '%s' but it was not found in draws.\n",
      gq_var,
      "Available variables: ", paste(model$variables, collapse = ", ")
    ))
  }

  fnames <- intersect(vars, names(data))
  n_obs  <- nrow(data)

  # Extract named GQ variable as draws_df (S draws × N cols + 3 metadata)
  draws_df <- posterior::subset_draws(model$draws, variable = gq_var) |>
    posterior::as_draws_df()

  n_draws <- nrow(draws_df)

  # Pivot to long: one row per (draw, observation)
  long_df <- draws_df |>
    tidyr::pivot_longer(
      cols      = dplyr::starts_with(gq_var),
      names_to  = ".obs_idx",
      values_to = "y"
    ) |>
    dplyr::mutate(
      .obs_idx = as.integer(stringr::str_extract(.obs_idx, "\\d+"))
    )

  # Apply transform_fn element-wise if present
  if (!is.null(spec$transform_fn)) {
    long_df <- dplyr::mutate(long_df, y = spec$transform_fn(y))
  }

  # Attach feature columns (tile data frame to match long format)
  if (length(fnames) > 0L) {
    features_rep <- data |>
      dplyr::select(dplyr::all_of(fnames)) |>
      dplyr::slice(rep(seq_len(n_obs), times = n_draws))

    long_df <- dplyr::bind_cols(long_df, features_rep)
  }

  out <- long_df |>
    dplyr::select(
      .chain, .iteration, .draw, .obs_idx,
      dplyr::all_of(fnames),
      y
    )

  # Apply filter_fn if present (operates on the full long tibble)
  if (!is.null(spec$filter_fn)) {
    out <- spec$filter_fn(out)
    if (nrow(out) == 0L)
      stop(sprintf(
        "outcome_spec filter_fn removed all rows for variable '%s'.", gq_var
      ))
  }

  out
}


# ---- internal: draw extraction ----------------------------------------------

#' Build the long yrep tibble, using either the descriptive or counterfactual
#' method.
#'
#' "descriptive"    — uses grouped_ppc() on observed draws (original path).
#' "counterfactual" — runs generate_quantities() at each point of a
#'                    contrast_var grid; returns a tibble with an extra
#'                    column `contrast_var_value` identifying the grid point.
#'                    The baseline (min of contrast_var) is included as
#'                    grid point 1 and used internally for differencing.
#'
#' @param model        fitted `stanmodel` object
#' @param data         data frame used for fitting (nrow = N)
#' @param vars         character vector of feature columns to attach
#' @param method       "descriptive" (default) or "counterfactual"  # <<< CHANGED from grouped_ppc/generated_quantities
#' @param gq_steps     integer; number of grid steps for counterfactual
#'                     (default 30, ignored for descriptive)
#' @param contrast_var character; name of the variable to sweep in the grid
#'                     (required when method = "counterfactual")
#' @param gq_n_draws   optional integer; if set, subsample this many posterior
#'                     draws before each `generate_quantities()` call in the
#'                     counterfactual path (see `gq.stanmodel()`). Ignored for
#'                     method = "descriptive".
#' @return tibble; descriptive path: .draw, .iteration, .chain, <vars>, y
#'         counterfactual path: same + contrast_var_value (numeric grid point)
#' @keywords internal

.get_yrep_draws <- function(model,
                           data,
                           vars,
                           method           = c("descriptive", "counterfactual"),
                           gq_steps         = 5,
                           contrast_var     = NULL,
                           contrast_levels  = NULL,
                           continous        = FALSE,
                           spec             = NULL,
                           gq_n_draws       = NULL)
{
  
  method <- match.arg(method)
  vars <- unique(vars)

  # Resolve spec: NULL -> legacy y_rep behaviour (backwards compatible)
  if (is.null(spec)) spec <- os_yrep()

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }

  # ---- descriptive path ---------------------------------------------------
  if (method == "descriptive") {
    yrep <- .extract_spec_draws(model, data, spec, vars)

    # Rename .obs_idx -> drop, harmonise column names for downstream use
    yrep <- yrep |> dplyr::select(-.obs_idx)

    return(yrep)
  }

  # ---- counterfactual path ------------------------------------------------
  # Counterfactual always extracts y_rep (the full mixture draw) from GQ.
  # After differencing, the result is in the scale of the GQ variable.
  # If spec requests a different variable, we override gq_var for the CF path
  # to use the specified variable instead of y_rep.
  # NOTE: p_respond and mu_y are deterministic given parameters, so CF
  # on them is well-defined and numerically cleaner than on y_rep.

  if (is.null(contrast_var)) {
    stop("`contrast_var` must be supplied when `method = 'counterfactual'`.")
  }
 # >>> CHANGED: for virtual numeric_array columns (e.g. "roi[CST_L]"),
  # the grid source is extra[[arr_name]][, col_name], NOT data[[contrast_var]].
  # The injected data column is a placeholder and may be constant.
  # For normal columns, behaviour is unchanged.
  parsed_cv <- .parse_virtual_contrast_var(contrast_var)

  if (!is.null(parsed_cv)) {
    # Virtual column: source grid from extra matrix column
    arr  <- model$extra[[parsed_cv$array]]
    if (is.null(arr))
      stop(sprintf(
        ".get_yrep_draws: numeric_array '%s' not found in model$extra.",
        parsed_cv$array
      ))
    if (!parsed_cv$col %in% colnames(arr))
      stop(sprintf(
        ".get_yrep_draws: column '%s' not found in numeric_array '%s'.",
        parsed_cv$col, parsed_cv$array
      ))
    x <- arr[, parsed_cv$col]
  } else {
    # Normal column: must exist in data
    if (!contrast_var %in% names(data))
      stop(sprintf("`contrast_var` '%s' was not found in `data`.", contrast_var))
    x <- data[[contrast_var]]
  }
  # >>> END CHANGE

  gq_var_cf <- spec$gq_var   # <<< CHANGED from hard-coded "y_rep"

  
  cvar_grid <- NULL

  # <<< CHANGED: categorical / dichotomous variables use observed levels
  if (is.factor(x) || is.character(x)) {
    if (is.null(contrast_levels)) {
      stop(sprintf(
        "contrast_levels must be specified for categorical contrast_var '%s'.",
        contrast_var
      ))
      } else {
       for (j in seq_along(contrast_levels)) {
        for (k in seq_along(contrast_levels[[j]])) {
         if (!contrast_levels[[j]][k] %in% unique(x)) {
          stop(sprintf(
            "contrast_levels[%d][%d] = '%s' is not a level of contrast_var '%s'.",
            j, k, contrast_levels[[j]][k], contrast_var
          ))
        } else {
          cvar_grid <- c(cvar_grid, contrast_levels[[j]][k])
        }
      }
    }
    cvar_grid <- unique(cvar_grid)
  }

    
} else {

    # <<< CHANGED: continuous variables still use a numeric grid
    xr <- range(x, na.rm = TRUE)
    if (isTRUE(all.equal(xr[1], xr[2]))) {
      cvar_grid <- xr[1]
    } else {
      cvar_grid <- seq(xr[1], xr[2], length.out = gq_steps)
    }
  }

  if (length(cvar_grid) < 2L) {
    stop(sprintf(
      "Contrast variable '%s' must have at least two distinct values for counterfactual contrasts.",
      contrast_var
    ))
  }



  feat_cols <- intersect(vars, names(data))
  feat_cols <- setdiff(feat_cols, contrast_var)

  .to_long_draws <- function(draws_obj, data_df, feature_cols) {
   
    draws_df <- posterior::subset_draws(draws_obj, variable = gq_var_cf) |>
      posterior::as_draws_df()

    n_obs <- ncol(draws_df) - 3L
    if (n_obs != nrow(data_df)) {
      stop(sprintf(
        "%s has %d observations but `data` has %d rows.",
        gq_var_cf, n_obs, nrow(data_df)
      ))
    }

    long_df <- draws_df |>
      tidyr::pivot_longer(
        cols      = starts_with(gq_var_cf),
        names_to  = ".obs_idx",
        values_to = "y"
      ) |>
      dplyr::mutate(
        .obs_idx = as.integer(stringr::str_extract(.obs_idx, "\\d+"))
      )

    if (length(feature_cols) > 0L) {
      feature_df <- data_df |>
        dplyr::select(dplyr::all_of(feature_cols)) |>
        dplyr::slice(rep(seq_len(nrow(data_df)), times = nrow(draws_df)))

      long_df <- dplyr::bind_cols(long_df, feature_df)
    }

    long_df |>
      dplyr::select(
        .chain, .iteration, .draw, .obs_idx, dplyr::all_of(feature_cols), y
      )
  }
  
 
  

  gq_contrasts <- vector("list", length(cvar_grid))

  # <<< NEW: distinguish continuous (return all grid points) from
  # categorical (return specified contrast_levels pair)
  
    if (continous || is.null(contrast_levels)) {
    # ---- continuous: simulate at each grid point, convert to long tibble --
    for (i in seq_along(cvar_grid)) {
      cat(sprintf("Simulating at %s = %s...\n", contrast_var, cvar_grid[i]))
      gq_i <- model |> .simulate_at(
        contrast_var = contrast_var,
        value        = cvar_grid[i],
        outcome_var  = gq_var_cf,
        n_draws      = gq_n_draws
      )

      # <<< FIX: gq_i is a draws_matrix; convert and pivot to long before mutating
      gq_contrasts[[i]] <- gq_i |>
        posterior::as_draws_df() |>                      # draws_matrix -> draws_df
        .to_long_draws(data, feat_cols) |>               # pivot to obs-level long tibble
        dplyr::mutate(
          contrast_var_value = cvar_grid[i],
          baseline_var_value = cvar_grid[1]
        )
    }

  } else {
    # ---- categorical: simulate at each level, convert, then difference ----

    gq_levels <- vector("list", length(cvar_grid))
    for (i in seq_along(cvar_grid)) {
      cat(sprintf("Simulating at %s = %s...\n", contrast_var, cvar_grid[i]))
      gq_levels[[i]] <- model |> .simulate_at(
        contrast_var = contrast_var,
        value        = cvar_grid[i],
        outcome_var  = gq_var_cf,
        n_draws      = gq_n_draws
      ) |>
        posterior::as_draws_df() |>                      # <<< FIX: same conversion
        .to_long_draws(data, feat_cols)
    }

    for (i in seq_along(contrast_levels)) {
      baseline_val <- contrast_levels[[i]][1]
      index_val    <- contrast_levels[[i]][2]
      baseline_idx <- which(cvar_grid == baseline_val)
      index_idx    <- which(cvar_grid == index_val)

      gq_contrasts[[i]] <- dplyr::left_join(
        gq_levels[[index_idx]],
        gq_levels[[baseline_idx]],
        by     = c(".chain", ".iteration", ".draw", ".obs_idx", feat_cols),
        suffix = c("_index", "_baseline")
      ) |>
        dplyr::mutate(
          y                  = y_index - y_baseline,
          contrast_var_value = index_val,
          baseline_var_value = baseline_val,
          contrast_name      = sprintf("%s - %s", index_val, baseline_val)
        ) |>
        dplyr::select(-dplyr::starts_with("y_"))
    }
  }

  out <- dplyr::bind_rows(gq_contrasts)
  out |>
    dplyr::select(
      .chain, .iteration, .draw, .obs_idx,
      dplyr::all_of(feat_cols),
      y, contrast_var_value, baseline_var_value
    )
} 

 
# # ---- internal: one contrast -------------------------------------------------

#' @keywords internal
.resolve_contrast_levels <- function(co, yrep) {
 
  if (co$continuous) return(NULL)

  col  <- yrep[[co$contrast_var]]
  lvls <- if (is.factor(col)) levels(col) else sort(unique(col))

  if (!is.null(co$contrast_levels)) return(co$contrast_levels)

  if (length(lvls) < 2)
    stop(sprintf("contrast_var '%s' has fewer than 2 levels.", co$contrast_var))

  warning(sprintf(
    "contrast_levels not specified for '%s'; using '%s' - '%s' by default.",
    co$contrast_var, lvls[1], lvls[2]
  ))
  lvls[1:2]
}


# shared helper that applies the hi/lo split for both continuous
# and categorical contrast_var, returning a list with:
#   - yrep_filtered : rows restricted to hi and lo groups
#   - contrast_var_sym : column name to pivot on ("hi"/"lo" or original levels)
#   - hi_lbl / lo_lbl : character labels used in the contrast
#' @keywords internal
.apply_contrast_split <- function(co, yrep) {

  if (co$continuous) {
    # ---- continuous: quantile or absolute threshold split -----------------
    if (co$absolute_thresholds) {
      if (is.null(co$contrast_levels))
        stop("contrast_levels must be specified when absolute_thresholds = TRUE.")
      hi_threshold <- co$contrast_levels[2]   # values <= this are "hi"
      lo_threshold <- co$contrast_levels[1]   # values >  this are "lo"
    } else {
      probs        <- if (!is.null(co$contrast_levels)) co$contrast_levels else c(0.75, 0.25)
      hi_threshold <- quantile(yrep[[co$contrast_var]], probs[2], na.rm = TRUE)
      lo_threshold <- quantile(yrep[[co$contrast_var]], probs[1], na.rm = TRUE)
    }

    yrep_filtered <- yrep |>
      dplyr::mutate(
        .contrast_group = dplyr::case_when(
          .data[[co$contrast_var]] >= hi_threshold ~ "hi",
          .data[[co$contrast_var]] <  lo_threshold ~ "lo",
          TRUE                                     ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(.contrast_group))

    list(
      yrep_filtered    = yrep_filtered,
      contrast_var_sym = ".contrast_group",
      hi_lbl           = "hi",
      lo_lbl           = "lo"
    )

  } else {
    # ---- categorical: filter to the two resolved levels -------------------
    lvls <- .resolve_contrast_levels(co, yrep)   # <<< reuses existing logic
    hi_lbl <- lvls[2]
    lo_lbl <- lvls[1]

    yrep_filtered <- yrep |>
      dplyr::filter(.data[[co$contrast_var]] %in% lvls)

    list(
      yrep_filtered    = yrep_filtered,
      contrast_var_sym = co$contrast_var,
      hi_lbl           = hi_lbl,
      lo_lbl           = lo_lbl
    )
  }
}

# helper: expand a single contrast referencing a numeric_array into
# one contrast per column. The virtual column name is "<array>[<colname>]"
# (or simply the colname if unambiguous) and is injected into data on-the-fly
# inside .simulate_at() / .get_yrep_draws().

#' @keywords internal
# >>> CHANGED: .expand_numeric_array_contrast accepts optional col subset
.expand_numeric_array_contrast <- function(co, model) {
  arr_name <- co$contrast_var
  if (is.null(model$numeric_arrays) || !arr_name %in% names(model$numeric_arrays))
    return(list(co))

  col_names <- model$numeric_arrays[[arr_name]]
  if (length(col_names) == 0L)
    stop(sprintf("numeric_array '%s' has no registered column names.", arr_name))

  # >>> NEW: honour column subset stored on model$numeric_array_cols[[arr_name]]
  # This is set by infer_contrasts() from the numeric_arrays argument when
  # the user passes a named list like list(roi = c("CST_L", "M1_L")).
  subset_cols <- model$numeric_array_cols[[arr_name]]
  if (!is.null(subset_cols)) {
    unknown <- setdiff(subset_cols, col_names)
    if (length(unknown) > 0)
      stop(sprintf(
        "numeric_array '%s': requested columns not found: %s\n  Available: %s",
        arr_name,
        paste(unknown, collapse = ", "),
        paste(col_names, collapse = ", ")
      ))
    col_names <- subset_cols
  }
  # >>> END NEW

  lapply(col_names, function(cn) {
    new_co <- co
    new_co$contrast_var  <- sprintf("%s[%s]", arr_name, cn)
    new_co$.array_source <- arr_name
    new_co$.array_col    <- cn
    new_co
  })
}
# >>> END CHANGE



# ---- infer_contrasts --------------------------------------------------------
#' Compute posterior predictive contrasts from draws
#'
#' Extracts and contrasts posterior predictive distributions across specified
#' grouping and contrast variables. Supports multiple inferential targets via
#' `outcome_spec` objects, including compound hurdle models. Returns different
#' result structures depending on contrast type and draws method.
#'
#' @param model        a fitted `stanmodel` object
#' @param contrasts    a single `contrast` object or a list of `contrast` objects
#' @param data         data frame used for fitting; if `NULL` uses `model$data` (default NULL)
#' @param summary_fn   function to aggregate y within each draw × condition cell (default `mean`)
#' @param prob         HPDI width for credible intervals (default 0.95)
#' @param plot         logical; if TRUE, render diagnostic plots (default FALSE)
#' @param draws_method "descriptive" (default) uses observed draws grouped by features.
#'   "counterfactual" sweeps `contrast_var` across its observed range via
#'   `generate_quantities()` and computes contrasts against the baseline (minimum value).
#' @param gq_steps     integer; grid resolution for counterfactual sweeps (default 5,
#'   ignored for descriptive method)
#' @param gq_n_draws   optional integer; if set, subsample this many posterior draws
#'   before each `generate_quantities()` call in the counterfactual path — trades
#'   Monte Carlo precision for speed and output size. Default `NULL` uses all draws
#'   (unchanged behavior). Ignored for `draws_method = "descriptive"`, which never
#'   calls `generate_quantities()`.
#'
#' @return
#' **Single contrast, non-hurdle outcome:**
#'   - `descriptive` + one-way (no two_way): list with `summary` (tibble, one row per grouping cell)
#'     and `posterior` (tibble of full draws)
#'   - `descriptive` + two-way: list with `scenario_probs` (tibble of scenario probabilities),
#'     `scenario_draws` (per-draw scenario classifications), `lv1`/`lv2` (outcome levels),
#'     `hi_lbl`/`lo_lbl` (contrast split labels), and `contrast` (the contrast object)
#'   - `counterfactual` + one-way (continuous): list with `grid_summary` (one row per grid point),
#'     `grid_posterior` (full draws per grid point), and `contrast`
#'   - `counterfactual` + one-way (categorical): list with `posterior`, `summary`, and `contrast`
#'   - `counterfactual` + two-way: list with `scenario_grid` (scenario probs per grid point),
#'     `scenario_draws` (per-draw classifications), `lv1`/`lv2`, `baseline_val`, and `contrast`
#'
#' **Single contrast with hurdle outcome** (via `os_hurdle()`):
#'   - list with named elements `marginal`, `extensive`, `intensive`, each containing
#'     the single-contrast result structure appropriate to the draws method
#'   - plus `contrast` (the original contrast object)
#'
#' **Multiple contrasts:**
#'   - list of result objects (one per input contrast), with the structure described above
#'   - if any input contrast uses `os_hurdle()`, the corresponding output element
#'     is itself a list of three margin-specific results
#'
#' @details
#' **Outcome specifications** allow flexible choice of inferential target:
#' - `os_yrep()`: full posterior predictive draws (default, backwards-compatible)
#' - `os_marginal()`: E(y) for hurdle models (primary inferential target)
#' - `os_extensive()`: P(respond), the extensive margin
#' - `os_intensive()`: E(y | respond), the intensive margin
#' - `os_hurdle()`: all three margins in one pass
#'
#' **Draws methods:**
#' - Descriptive: groups observed posterior predictive draws by feature columns
#'   and applies the contrast operation (hi - lo grouping, two-way scenario classification, etc.)
#' - Counterfactual: simulates posterior predictive at each grid point of `contrast_var`,
#'   then differences vs. baseline; produces grid-indexed results showing how contrasts
#'   vary across the range of the contrast variable
#'
#' **Two-way contrasts** (via `co$two_way`) classify each draw into one of four scenarios:
#' affects both outcomes, first only, second only, or neither. Supported in both
#' descriptive and counterfactual modes, for both categorical and continuous contrast variables.
#'
#' @export
infer_contrasts <- function(model,
                            contrasts,
                            data           = NULL,
                            extra          = NULL,
                            numeric_arrays = NULL,
                            summary_fn     = mean,
                            prob           = 0.95,
                            plot           = FALSE,
                            draws_method   = c("descriptive", "counterfactual"),
                            gq_steps       = 2,
                            gq_n_draws     = NULL) {
  
  draws_method <- match.arg(draws_method)
  if (inherits(contrasts, "contrast")) contrasts <- list(contrasts)
  stopifnot(all(vapply(contrasts, inherits, logical(1), "contrast")))

  if (is.null(data))  data  <- model$data
  if (is.null(extra)) extra <- model$extra

  # ---- register numeric_arrays + optional column subsets ------------------
  # numeric_arrays can be:
  #   character vector: c("roi", "smatt")       -> use all columns
  #   named list:       list(roi = c("CST_L", "M1_L"), smatt = NULL)
  #                     -> subset roi to two cols, use all smatt cols (NULL = all)
  # >>> NEW: normalise to named list internally
  if (is.character(numeric_arrays)) {
    numeric_arrays <- stats::setNames(
      vector("list", length(numeric_arrays)),
      numeric_arrays
    )  # list(roi = NULL, smatt = NULL) — NULL means all columns
  }

  if (is.null(model$numeric_arrays))     model$numeric_arrays     <- list()
  if (is.null(model$numeric_array_cols)) model$numeric_array_cols <- list()

  for (obj_name in names(numeric_arrays)) {
    obj <- extra[[obj_name]]
    if (is.null(obj) || (!is.matrix(obj) && !is.data.frame(obj))) {
      warning(sprintf(
        "numeric_array '%s' not found in extra or not a matrix/data.frame. Skipping.",
        obj_name
      ))
      next
    }
    col_nms <- colnames(obj)
    if (is.null(col_nms)) {
      warning(sprintf("numeric_array '%s' has no column names. Skipping.", obj_name))
      next
    }
    model$numeric_arrays[[obj_name]] <- col_nms   # full registry

    # >>> NEW: store requested subset (NULL = all, checked later in .expand_*)
    requested <- numeric_arrays[[obj_name]]
    model$numeric_array_cols[[obj_name]] <- requested  # NULL or char vector
    # >>> END NEW

    cat(sprintf(
      "  numeric_array '%s': %d total cols, using %s\n",
      obj_name,
      length(col_nms),
      if (is.null(requested))
        sprintf("all (%d)", length(col_nms))
      else
        sprintf("%d: %s", length(requested), paste(requested, collapse = ", "))
    ))
  }                             
  
  # ---- numeric_array expansion + virtual column injection ------------------
  expanded_array_contrasts <- list()
  for (co in contrasts) {
    sub <- .expand_numeric_array_contrast(co, model)
    expanded_array_contrasts <- c(expanded_array_contrasts, sub)
    
    
    for (sub_co_indx in 1:length(sub)) {
      sub_co <- sub[[sub_co_indx]]
      if (!is.null(sub_co$.array_source)) {
        vcol <- sub_co$contrast_var
        arr  <- extra[[sub_co$.array_source]]
        if (is.null(arr)) { warning(sprintf("numeric_array '%s' not found in extra.", sub_co$.array_source)); next }
        if (nrow(arr) != nrow(data)) stop(sprintf("numeric_array '%s' has %d rows but data has %d rows.", sub_co$.array_source, nrow(arr), nrow(data)))
        cn <- sub_co$.array_col
        if (!cn %in% colnames(arr)) { warning(sprintf("Column '%s' not found in numeric_array '%s'.", cn, sub_co$.array_source)); next }
        if (!vcol %in% names(data)) data[[vcol]] <- arr[, cn]
      }
    }
  }
  
  contrasts <- expanded_array_contrasts

  # ---- resolve outcome specs; expand hurdle contrasts ----------------------
  .resolve_spec <- function(co) if (is.null(co$outcome)) os_yrep() else co$outcome

  expanded_contrasts <- list()
  hurdle_map         <- list()

  for (i in seq_along(contrasts)) {
    co   <- contrasts[[i]]
    spec <- .resolve_spec(co)
    if (inherits(spec, "outcome_spec_hurdle")) {
      for (margin_name in c("marginal", "extensive", "intensive")) {
        sub_co         <- co
        sub_co$outcome <- spec[[margin_name]]
        expanded_contrasts <- c(expanded_contrasts, list(sub_co))
        hurdle_map[[length(expanded_contrasts)]] <- list(original_idx = i, margin = margin_name)
      }
    } else {
      expanded_contrasts <- c(expanded_contrasts, list(co))
      hurdle_map[[length(expanded_contrasts)]] <- list(original_idx = i, margin = NULL)
    }
  }
 

  # >>> CHANGED: collect all_needed_vars AFTER expansion so virtual column
  # names like "roi[CST_L]" are included from the expanded contrast objects.
  all_vars_in_expr <- function(expr) {
    if (rlang::is_symbol(expr)) return(as.character(expr))
    if (rlang::is_call(expr))   return(unlist(lapply(as.list(expr)[-1], all_vars_in_expr)))
    character()
  }

  .extract_vars <- function(co) {
    vars <- c(co$grouping_vars, co$contrast_var)
    if (!is.null(co$two_way))      vars <- c(vars, co$two_way)
    if (!is.null(co$filter_draws)) vars <- unique(c(vars, all_vars_in_expr(co$filter_draws)))
    vars
  }
  all_needed_vars <- unique(unlist(lapply(expanded_contrasts, .extract_vars), use.names = FALSE))
  # >>> END CHANGE

  contrast_levels_list <- lapply(expanded_contrasts, function(co) co$contrast_levels)
  contrast_continuous  <- vapply(expanded_contrasts, function(co) co$continuous, logical(1))

  if (draws_method == "counterfactual" && length(expanded_contrasts) > 1 &&
      length(unique(contrast_continuous)) > 1) {
    stop("All contrasts must be either continuous or categorical when draws_method = 'counterfactual'.")
  }

  # >>> CHANGED: group by (gq_var, contrast_var) instead of gq_var alone.
  # This ensures:
  #   - hurdle margins (same contrast_var, different gq_var) still group
  #     separately by gq_var as before — no regression.
  #   - numeric_array-expanded contrasts (same gq_var = "y_rep", different
  #     contrast_var per ROI) each get their OWN .get_yrep_draws() call.
  #   - genuinely unrelated contrasts with different contrast_var also work.
  specs_by_group <- split(
    seq_along(expanded_contrasts),
    vapply(expanded_contrasts, function(co)
      paste0(.resolve_spec(co)$gq_var, "|||", co$contrast_var),
      character(1)
    )
  )
  # >>> END CHANGE

  # ---- run per (gq_var, contrast_var) group --------------------------------
  flat_results <- vector("list", length(expanded_contrasts))

  #[1:3] for debugging



  for (group_key in names(specs_by_group)) {          # >>> CHANGED: renamed loop var
    cat(sprintf("Processing contrast group: %s [%d/%d]\n", group_key, which(names(specs_by_group) == group_key), length(specs_by_group)))
    idx_group <- specs_by_group[[group_key]]           # >>> CHANGED
    spec_rep  <- .resolve_spec(expanded_contrasts[[idx_group[1]]])

    yrep <- .get_yrep_draws(
      model           = model,
      data            = data,
      vars            = all_needed_vars,
      method          = draws_method,
      gq_steps        = gq_steps,
      contrast_var    = if (draws_method == "counterfactual")
                          expanded_contrasts[[idx_group[1]]]$contrast_var
                        else NULL,
      contrast_levels = if (draws_method == "counterfactual" &&
                            !expanded_contrasts[[idx_group[1]]]$continuous)
                          contrast_levels_list[idx_group]
                        else NULL,
      continous       = any(contrast_continuous[idx_group]),  # >>> CHANGED: scope to group
      spec            = spec_rep,
      gq_n_draws      = gq_n_draws
    )

    for (i in idx_group) {
      co <- expanded_contrasts[[i]]
      yrep_i <- if (draws_method == "counterfactual" &&
                    "contrast_name" %in% names(yrep) &&
                    !is.null(co$contrast_levels)) {
        expected_name <- sprintf("%s - %s", co$contrast_levels[2], co$contrast_levels[1])
        dplyr::filter(yrep, contrast_name == expected_name)
      } else {
        yrep
      }

      flat_results[[i]] <- if (draws_method == "counterfactual") {
        res <- .compute_counterfactual_contrast(yrep_i, co, summary_fn, prob)
        res$contrast <- co
        res
      } else if (!is.null(co$two_way)) {
        .compute_two_way_contrast(yrep_i, co, summary_fn, prob)
      } else {
        .compute_one_contrast(yrep_i, co, summary_fn, prob)
      }
    }
  }

  # ---- reassemble: collapse hurdle sub-contrasts; keep others flat ---------
  # >>> CHANGED: hurdle_map$original_idx now maps into the post-expansion
  # `contrasts` list (which may be larger than the user's original input if
  # numeric_array expansion ran). Non-hurdle results are assigned directly.
  n_orig  <- length(contrasts)
  results <- vector("list", n_orig)

  for (flat_i in seq_along(flat_results)) {
    hm     <- hurdle_map[[flat_i]]
    orig_i <- hm$original_idx
    margin <- hm$margin
    if (is.null(margin)) {
      results[[orig_i]] <- flat_results[[flat_i]]
    } else {
      if (is.null(results[[orig_i]])) results[[orig_i]] <- list()
      results[[orig_i]][[margin]] <- flat_results[[flat_i]]
    }
  }

 # >>> CHANGED: attach $contrast consistently for hurdle results here,
  # not in a separate loop. Remove the old separate loop entirely.
  for (i in seq_along(results)) {
    is_hurdle <- any(vapply(hurdle_map,
                            function(h) h$original_idx == i && !is.null(h$margin),
                            logical(1)))
    if (is_hurdle) results[[i]]$contrast <- contrasts[[i]]
  }
  # >>> END CHANGE

   # >>> CHANGED: wrap results into inference objects correctly.
  # - Non-hurdle: wrap flat result directly.
  # - Hurdle: wrap each margin separately, keyed by margin name.
  #   Exclude $contrast from the lapply to avoid wrapping it as a margin.
  .wrap_result <- function(res, co) {
    is_hurdle <- is.list(res) && all(c("marginal", "extensive", "intensive") %in% names(res))
    if (is_hurdle) {
      wrapped <- lapply(
        c("marginal", "extensive", "intensive"),
        function(margin_name) as_inference(res[[margin_name]], co)
      )
      names(wrapped) <- c("marginal", "extensive", "intensive")
      wrapped$contrast <- co
      wrapped
    } else {
      as_inference(res, co)
    }
  }

  if (length(results) == 1L) {
    .wrap_result(results[[1L]], contrasts[[1L]])
  } else {
    lapply(seq_along(results), function(i) .wrap_result(results[[i]], contrasts[[i]]))
  }
  # >>> END CHANGE

}

#' @keywords internal
.define_scenarios <- function(draw_wide, col1, col2, lv1, lv2, effect_positive = TRUE) {

  if (effect_positive) {
    draw_scenarios <- draw_wide |> dplyr::mutate(
      .affect_1 = .data[[col1]] > 0,
      .affect_2 = .data[[col2]] > 0)
  } else {
    draw_scenarios <- draw_wide |> dplyr::mutate(
      .affect_1 = .data[[col1]] < 0,
      .affect_2 = .data[[col2]] < 0) 
  } 

  draw_scenarios <- draw_scenarios |> dplyr::mutate(
    scenario  = dplyr::case_when(
      .affect_1 &  .affect_2  ~ "Affects both",
      .affect_1 & !.affect_2  ~ paste0("Affects ", lv1, " only"),
      !.affect_1 &  .affect_2 ~ paste0("Affects ", lv2, " only"),
      TRUE                    ~ "Affects neither"
    )
  ) |>
  dplyr::select(-dplyr::starts_with(".affect_"))
  
  return(draw_scenarios)
}


#' @keywords internal
.compute_scenario_grid <- function(draw_scenarios, by_vars, wide = FALSE) {
  out <- draw_scenarios |>
    dplyr::summarise(
      n = dplyr::n(),
      .by = dplyr::all_of(c(by_vars,'scenario'))
    ) |>
    dplyr::mutate(
      probability = n / sum(n),
      .by = dplyr::all_of(by_vars)
    ) |>
    dplyr::select(-n) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(by_vars)))

  if (wide) {
    out <- out |>
      tidyr::pivot_wider(
        names_from  = scenario,
        values_from = probability,
        values_fill = 0
      )
  }
  out
}

# ---- counterfactual contrast ------------------------------------------------

#' Compute scenario probabilities across a contrast_var grid (counterfactual)
#'
#' Called automatically by infer_contrasts() when draws_method = "counterfactual".
#' yrep already contains pre-differenced y values (grid point minus baseline)
#' and a `contrast_var_value` column from .get_yrep_draws().
#'
#' For each grid point and each posterior draw, summarises y per outer grouping
#' cell and (when co$two_way is set) per outcome level, then classifies into
#' 4 scenarios exactly as in .compute_two_way_contrast().
#'
#' @param yrep       long tibble from .get_yrep_draws(method = "counterfactual")
#' @param co         a `contrast` object
#' @param summary_fn aggregation function (default mean)
#' @param prob       HPDI width (unused here, kept for interface parity)
#' @return list with:
#'   - scenario_grid  : tibble of scenario probs per grid point (+ outer grouping)
#'   - scenario_draws : per-draw classifications per grid point
#'   - contrast_var   : name of the swept variable
#'   - lv1, lv2        : two-way outcome levels used in the contrast
#'   - baseline_val   : baseline value of contrast_var used for differencing
#'   - params         : list of additional parameters passed to the contrast computation
#' @keywords internal  
# ---- .compute_counterfactual_contrast ---------------------------------------
# Changes marked with # <<< CHANGED and # <<< NEW
.compute_counterfactual_contrast <- function(yrep, co, summary_fn, prob) {

  if (!is.null(co$filter_draws)) {
    yrep <- dplyr::filter(yrep, !!co$filter_draws)
    if (nrow(yrep) == 0L)
      stop("filter_draws removed all rows. Check your filter expression.")
  }

  outer_vars <- if (is.null(co$grouping_vars)) character() else co$grouping_vars

  # ---- 2. two-way branch --------------------------------------------------
  # <<< CHANGED: removed the continuous guard; now supported via
  # .apply_contrast_split() for both categorical and continuous contrast_var.
  if (!is.null(co$two_way)) {

    two_way_col <- co$two_way
    outcome_col <- yrep[[two_way_col]]
    all_lvls    <- if (is.factor(outcome_col)) levels(outcome_col)
                   else sort(unique(outcome_col))

    lvls <- if (!is.null(co$two_way_levels)) {
      co$two_way_levels
    } else {
      if (length(all_lvls) < 2L)
        stop(sprintf("two_way variable '%s' has fewer than 2 levels.", two_way_col))
      if (length(all_lvls) > 2L)
        warning(sprintf(
          "two_way variable '%s' has %d levels; using first two: '%s', '%s'.",
          two_way_col, length(all_lvls), all_lvls[1], all_lvls[2]
        ))
      all_lvls[1:2]
    }

    lv1 <- lvls[1]
    lv2 <- lvls[2]

    yrep_two <- yrep |>
      dplyr::filter(.data[[two_way_col]] %in% lvls)

     # ... existing code above: lv1, lv2, yrep_two assignment ...

    if (co$continuous) {
      # <<< REPLACED: old hi/lo quantile split discarded middle grid points
      # and produced a single contrast per draw regardless of gq_steps.
      # New approach: difference each grid point vs baseline, classify
      # scenario per grid point per draw -> scenario_draws has contrast_var_value.

      baseline_val <- min(yrep_two$contrast_var_value, na.rm = TRUE)

      yrep_base <- yrep_two |>
        dplyr::filter(contrast_var_value == baseline_val) |>
        dplyr::select(.draw, .iteration, .chain, .obs_idx,
                      dplyr::all_of(outer_vars),
                      dplyr::all_of(two_way_col),
                      y_baseline = y)

      yrep_diff <- yrep_two |>
        dplyr::filter(contrast_var_value != baseline_val) |>
        dplyr::left_join(
          yrep_base,
          by = c(".draw", ".iteration", ".chain", ".obs_idx",
                 outer_vars, two_way_col)
        ) |>
        dplyr::mutate(y = y - y_baseline) |>
        dplyr::select(-y_baseline)

      draw_by <- c(".draw", ".iteration", ".chain",
                   "contrast_var_value", outer_vars, two_way_col)

      draw_summary <- yrep_diff |>
        dplyr::summarise(y = summary_fn(y), .by = dplyr::all_of(draw_by))

      col1 <- paste0("contrast_", lv1)
      col2 <- paste0("contrast_", lv2)

      draw_wide <- draw_summary |>
        tidyr::pivot_wider(
          names_from   = dplyr::all_of(two_way_col),
          values_from  = y,
          names_prefix = "contrast_"
        )
      
      draw_scenarios <- draw_wide |> .define_scenarios(col1, col2, lv1, lv2, effect_positive = co$effect_positive)

      by_vars <- c("contrast_var_value", outer_vars, "scenario")

      scenario_grid <- draw_scenarios |>
        dplyr::summarise(n = dplyr::n(), .by = dplyr::all_of(by_vars)) |>
        dplyr::mutate(
          probability = n / sum(n),
          .by = dplyr::all_of(c("contrast_var_value", outer_vars))
        ) |>
        dplyr::select(-n) |>
        dplyr::arrange(dplyr::across(dplyr::all_of(by_vars)))

      return(list(
        scenario_grid  = scenario_grid,
        scenario_draws = draw_scenarios,
        contrast       = co,
        lv1            = lv1,
        lv2            = lv2,
        baseline_val   = baseline_val,
        params     = list(summary_fn = summary_fn,
                               prob       = prob,
                               draws_method = 'counterfactual')
      ))
      # <<< END REPLACEMENT

    } else {
      # ---- categorical two-way (unchanged) --------------------------------
      # ... existing code below unchanged ...--
      draw_by <- c(".draw", ".iteration", ".chain", "contrast_var_value",
                   outer_vars, two_way_col)

      draw_summary <- yrep_two |>
        dplyr::summarise(y = summary_fn(y), .by = dplyr::all_of(draw_by))

      draw_wide <- draw_summary |>
        tidyr::pivot_wider(
          names_from   = dplyr::all_of(two_way_col),
          values_from  = y,
          names_prefix = "contrast_"
        )

      col1 <- paste0("contrast_", lv1)
      col2 <- paste0("contrast_", lv2)

      draw_scenarios <- draw_wide |> .define_scenarios(col1, col2, lv1, lv2, effect_positive = co$effect_positive)
      
  

      by_vars <- c("contrast_var_value", outer_vars, "scenario")

      scenario_grid <- draw_scenarios |> .compute_scenario_grid()

      scenario_grid <- draw_scenarios |>
        dplyr::summarise(n = dplyr::n(), .by = dplyr::all_of(by_vars)) |>
        dplyr::mutate(
          probability = n / sum(n),
          .by = dplyr::all_of(c("contrast_var_value", outer_vars))
        ) |>
        dplyr::select(-n) |>
        dplyr::arrange(dplyr::across(dplyr::all_of(by_vars)))

      return(list(
        scenario_grid  = scenario_grid,
        scenario_draws = draw_scenarios,
        contrast       = co,
        lv1            = lv1,
        lv2            = lv2,
        params     = list(summary_fn = summary_fn,
                               prob       = prob,
                               draws_method = 'counterfactual')

      ))
    }
  }

  # ---- 3. one-way branch --------------------------------------------------
  # <<< CHANGED: for continuous contrast_var, exploit the full gq_steps grid
  # (contrast_var_value column) rather than collapsing to a single hi-lo pair.
  # This yields a posterior distribution of contrasts *per grid point*, returned
  # as `grid_posterior` + `grid_summary` analogous to scenario_grid in two-way.
  # For categorical (pre-differenced, contrast_name present), behaviour unchanged.

  if (co$continuous) {
    # ---- continuous one-way: full grid summary ---------------------------
    # yrep has contrast_var_value (grid point) and y = GQ at that point.
    # baseline_var_value holds the minimum (grid point 1 value).
    # y is already differenced: GQ(grid_pt) - GQ(baseline) from .get_yrep_draws.
    # <<< NEW: summarise per draw × contrast_var_value × outer_vars
    draw_by <- c(".draw", ".iteration", ".chain",
                 "contrast_var_value", "baseline_var_value", outer_vars)

    draw_summary <- yrep |>
      dplyr::summarise(y = summary_fn(y), .by = dplyr::all_of(draw_by))

    alpha <- (1 - prob) / 2

    # Per-grid-point posterior summary (mean + interval)
    grid_summary <- draw_summary |>
      dplyr::summarise(
        mean = mean(y),
        lo   = stats::quantile(y, probs = alpha,     names = FALSE, type = 1),
        hi   = stats::quantile(y, probs = 1 - alpha, names = FALSE, type = 1),
        P_positive = mean(y > 0),
        P_negative = mean(y < 0),
        .by  = dplyr::all_of(c("contrast_var_value", "baseline_var_value",
                                outer_vars))
      )

    # Full posterior draws (for density plots or further processing)
    grid_posterior <- draw_summary |>
      dplyr::rename(contrast = y) |>
      dplyr::mutate(
        contrast_label = paste0(
          co$contrast_var, "=", contrast_var_value,
          " - ", baseline_var_value
        )
      )

    return(list(
      grid_summary   = grid_summary,    # <<< NEW: one row per grid point
      grid_posterior = grid_posterior,  # <<< NEW: full draws per grid point
      contrast       = co
    ))
    # <<< END NEW continuous one-way block

  } else {
    # ---- categorical one-way (unchanged) ----------------------------------
    draw_by <- c(".draw", ".iteration", ".chain", "contrast_var_value",
                 "baseline_var_value", outer_vars)

    draw_summary <- yrep |>
      dplyr::summarise(y = summary_fn(y), .by = dplyr::all_of(draw_by))

    posterior <- draw_summary |>
      dplyr::mutate(
        contrast = y,
        contrast_label = paste0(
          co$contrast_var, "=", contrast_var_value,
          " - ", baseline_var_value
        )
      ) |>
      dplyr::select(
        .draw, .iteration, .chain,
        contrast_var_value, baseline_var_value,
        dplyr::all_of(outer_vars),
        contrast, contrast_label
      )

    alpha <- (1 - prob) / 2

    summary_tbl <- posterior |>
      dplyr::summarise(
        mean = mean(contrast),
        lo   = stats::quantile(contrast, probs = alpha,     names = FALSE, type = 1),
        hi   = stats::quantile(contrast, probs = 1 - alpha, names = FALSE, type = 1),
        .by  = dplyr::all_of(c("contrast_label", outer_vars))
      )

    list(
      posterior = posterior,
      summary   = summary_tbl,
      contrast  = co
    )
  }
}
# # ---- counterfactual contrast ------------------------------------------------


# ---- plot helpers -----------------------------------------------------------

#' @keywords internal
.plot_two_way_scenarios <- function(res) {
  # CHANGED: guard setdiff when grouping_vars is NULL
  outer_vars <- setdiff(res$contrast$grouping_vars %||% character(),
                        res$contrast$two_way)

  p <- ggplot2::ggplot(
    res$scenario_probs,
    ggplot2::aes(x = scenario, y = probability, fill = scenario)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::labs(
      title    = get_label(res$contrast),
      subtitle = sprintf("'%s' contrast: %s — %s",
                         res$contrast$contrast_var, res$lv1, res$lv2),
      x = "Scenario",
      y = "Posterior probability"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  if (length(outer_vars) == 1L) {
    p <- p + ggplot2::facet_wrap(~ .data[[outer_vars[1]]])
  } else if (length(outer_vars) >= 2L) {
    p <- p + ggplot2::facet_grid(
      .data[[outer_vars[1]]] ~ .data[[outer_vars[2]]]
    )
  }
  p
}

#' @keywords internal
.plot_counterfactual_scenarios <- function(res) {
  # CHANGED: guard setdiff when grouping_vars is NULL
  outer_vars <- if (!is.null(res$contrast$two_way))
    setdiff(res$contrast$grouping_vars %||% character(), res$contrast$two_way)
  else
    res$contrast$grouping_vars %||% character()

  p <- ggplot2::ggplot(
    res$scenario_grid,
    ggplot2::aes(x = contrast_var_value, y = probability, color = scenario)
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::labs(
      title    = get_label(res$contrast),
      subtitle = sprintf("Counterfactual sweep of '%s' (baseline = min observed)",
                         res$contrast_var),
      x     = res$contrast_var,
      y     = "Posterior probability",
      color = "Scenario"
    ) +
    ggplot2::theme_minimal()

  if (length(outer_vars) == 1L) {
    p <- p + ggplot2::facet_wrap(~ .data[[outer_vars[1]]])
  } else if (length(outer_vars) >= 2L) {
    p <- p + ggplot2::facet_grid(
      .data[[outer_vars[1]]] ~ .data[[outer_vars[2]]]
    )
  }
  p
}


#' @keywords internal
.compute_one_contrast <- function(yrep, co, summary_fn, prob) {
  
  cell_vars <- co$grouping_vars

  if (!is.null(co$filter_draws)) {
    yrep <- dplyr::filter(yrep, !!co$filter_draws)
    if (nrow(yrep) == 0)
      stop("filter_draws removed all rows from the draws data frame. Check your filter expression.")
  }

  if (co$continuous) {
    hi_lbl <- "hi"
    lo_lbl <- "lo"

    if (co$absolute_thresholds) {
      if (is.null(co$contrast_levels))
        stop("contrast_levels must be specified when absolute_thresholds = TRUE.")
      hi_threshold <- co$contrast_levels[2]
      lo_threshold <- co$contrast_levels[1]
    } else {
      probs        <- if (!is.null(co$contrast_levels)) co$contrast_levels else c(0.75, 0.25)
      hi_threshold <- quantile(yrep[[co$contrast_var]], probs[2], na.rm = TRUE)
      lo_threshold <- quantile(yrep[[co$contrast_var]], probs[1], na.rm = TRUE)
    }

    yrep_filtered <- yrep |>
      dplyr::mutate(
        .contrast_group = dplyr::case_when(
          .data[[co$contrast_var]] >= hi_threshold ~ "hi",
          .data[[co$contrast_var]] <  lo_threshold ~ "lo",
          TRUE ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(.contrast_group))

    contrast_var_sym <- ".contrast_group"

  } else {
    lvls   <- .resolve_contrast_levels(co, yrep)
    hi_lbl <- lvls[2]
    lo_lbl <- lvls[1]

    yrep_filtered    <- yrep |>
      dplyr::filter(.data[[co$contrast_var]] %in% lvls)

    contrast_var_sym <- co$contrast_var
  }

  # CHANGED: guard .by when cell_vars is NULL
  by_draw_vars <- c(".draw", ".iteration", ".chain",
                    cell_vars,              # NULL is dropped silently by c()
                    contrast_var_sym)

  draw_summary <- yrep_filtered |>
    dplyr::summarise(
      y = summary_fn(y),
      .by = dplyr::all_of(by_draw_vars)
    )

  draw_wide <- draw_summary |>
    tidyr::pivot_wider(
      names_from  = dplyr::all_of(contrast_var_sym),
      values_from = y
    )

  draw_contrasts <- draw_wide |>
    dplyr::mutate(
      contrast       = .data[[hi_lbl]] - .data[[lo_lbl]],
      contrast_label = sprintf("%s - %s", hi_lbl, lo_lbl)
    )
  cat(sprintf('Computing contrast: %s - %s\n', hi_lbl, lo_lbl))

  contrast_summary <- draw_contrasts |>
    dplyr::summarise(
      contrast_mean = mean(contrast),
      contrast_sd   = sd(contrast),
      hpdi_lo       = rethinking::HPDI(contrast, prob = prob)[[1]],
      hpdi_hi       = rethinking::HPDI(contrast, prob = prob)[[2]],
      P_positive    = mean(contrast > 0),
      P_negative    = mean(contrast < 0),
      .by           = dplyr::all_of(cell_vars)
    ) |>
    dplyr::mutate(
      contrast_label = sprintf("%s - %s", hi_lbl, lo_lbl),
      prob           = prob
    )

  list(
    summary   = contrast_summary,
    posterior = draw_contrasts |>
      dplyr::select(dplyr::all_of(c(".draw", ".iteration", ".chain",
                                    cell_vars, "contrast", "contrast_label")))
  )
}


# ---- two-way scenario contrast ----------------------------------------------

#' @keywords internal
.compute_two_way_contrast <- function(yrep, co, summary_fn, prob) {

  # ---- 1. apply filter_draws if present -----------------------------------
  if (!is.null(co$filter_draws)) {
    yrep <- dplyr::filter(yrep, !!co$filter_draws)
    if (nrow(yrep) == 0)
      stop("filter_draws removed all rows. Check your filter expression.")
  }

  # ---- 2. resolve the two outcome levels from the two_way variable --------
  two_way_col <- co$two_way
  outcome_col <- yrep[[two_way_col]]
  all_lvls    <- if (is.factor(outcome_col)) levels(outcome_col) else sort(unique(outcome_col))

  lvls <- if (!is.null(co$two_way_levels)) {
    co$two_way_levels
  } else {
    if (length(all_lvls) < 2)
      stop(sprintf("two_way variable '%s' has fewer than 2 levels.", two_way_col))
    if (length(all_lvls) > 2)
      warning(sprintf(
        "two_way variable '%s' has %d levels; using first two: '%s', '%s'. Set `two_way_levels` to choose explicitly.",
        two_way_col, length(all_lvls), all_lvls[1], all_lvls[2]
      ))
    all_lvls[1:2]
  }

  lv1 <- lvls[1]
  lv2 <- lvls[2]

  # ---- 3. restrict to the two outcome levels before splitting -------------
  yrep_two <- yrep |>
    dplyr::filter(.data[[two_way_col]] %in% lvls)

  # ---- 4. apply hi/lo contrast split (continuous OR categorical) ----------
  # <<< CHANGED: was calling .resolve_contrast_levels() directly and filtering
  # manually; now delegates to .apply_contrast_split() which handles both cases
  split      <- .apply_contrast_split(co, yrep_two)
  yrep_split <- split$yrep_filtered
  cvar_sym   <- split$contrast_var_sym
  hi_lbl     <- split$hi_lbl
  lo_lbl     <- split$lo_lbl
  # <<< END CHANGE

# CHANGED: outer_vars must handle NULL grouping_vars
  outer_vars <- setdiff(co$grouping_vars %||% character(), two_way_col)

  # ---- 5. per-draw per-outcome-level mean summary -------------------------
  draw_summary <- yrep_split |>
    dplyr::summarise(
      y = summary_fn(y),
      .by = c(".draw", ".iteration", ".chain",
              dplyr::all_of(c(outer_vars, cvar_sym, two_way_col)))
    )

  # ---- 6. pivot contrast split wide (hi vs lo) ----------------------------
  draw_wide <- draw_summary |>
    tidyr::pivot_wider(
      names_from  = dplyr::all_of(cvar_sym),
      values_from = y
    ) |>
    dplyr::mutate(
      contrast = .data[[hi_lbl]] - .data[[lo_lbl]]
    )

  # ---- 7. pivot two_way levels wide (one contrast col per outcome level) --
  draw_wide2 <- draw_wide |>
    dplyr::select(dplyr::all_of(c(
      ".draw", ".iteration", ".chain", outer_vars, two_way_col, "contrast"
    ))) |>
    tidyr::pivot_wider(
      names_from   = dplyr::all_of(two_way_col),
      values_from  = contrast,
      names_prefix = "contrast_"
    )

  col1 <- paste0("contrast_", lv1)
  col2 <- paste0("contrast_", lv2)

  # ---- 8. classify each draw into one of 4 scenarios ----------------------
  draw_scenarios <- draw_wide2 |>
    .define_scenarios(col1, col2, lv1, lv2, effect_positive = co$effect_positive)

  # ---- 9. posterior probabilities per scenario (and outer grouping) -------
  by_vars <- if (length(outer_vars) > 0) c(outer_vars, "scenario") else "scenario"

  scenario_probs <- draw_scenarios |>
    dplyr::summarise(
      n   = dplyr::n(),
      .by = dplyr::all_of(by_vars)
    ) |>
    dplyr::mutate(
      probability = n / sum(n),
      .by          = if (length(outer_vars) > 0) dplyr::all_of(outer_vars) else NULL
    ) |>
    dplyr::select(-n) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(by_vars)))

  list(
    scenario_probs = scenario_probs,
    scenario_draws = draw_scenarios,
    lv1            = lv1,
    lv2            = lv2,
    hi_lbl         = hi_lbl,   # <<< ADDED: expose split labels for print/plot
    lo_lbl         = lo_lbl
  )
}

#' @keywords internal  # <<< ADDED: helper plot for two_way results
.plot_two_way_scenarios <- function(res) {
  outer_vars <- setdiff(res$contrast$grouping_vars, res$contrast$two_way)

  p <- ggplot2::ggplot(
    res$scenario_probs,
    ggplot2::aes(x = scenario, y = probability, fill = scenario)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::labs(
      title    = get_label(res$contrast),
      subtitle = sprintf("'%s' contrast: %s — %s",
                         res$contrast$contrast_var,
                         res$lv1, res$lv2),
      x = "Scenario",
      y = "Posterior probability"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  if (length(outer_vars) == 1) {
    p <- p + ggplot2::facet_wrap(~ .data[[outer_vars[1]]])
  } else if (length(outer_vars) >= 2) {
    p <- p + ggplot2::facet_grid(
      .data[[outer_vars[1]]] ~ .data[[outer_vars[2]]]
    )
  }
  p
}
# <<< END ADDITION: .plot_two_way_scenarios


  # ...existing code...