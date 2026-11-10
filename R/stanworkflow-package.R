#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom lifecycle deprecated
#' @importFrom stats quantile sd
#' @importFrom utils modifyList read.csv
#' @importFrom rlang .data := sym
#' @importFrom dplyr pull starts_with matches
## usethis namespace: end
NULL

# Column names used via non-standard evaluation (dplyr/ggplot2).
utils::globalVariables(c(
  ".chain", ".contrast_group", ".draw", ".iteration", ".obs_idx",
  "baseline_var_value", "category", "contrast_label", "contrast_name",
  "contrast_var_value", "ess_bulk", "ess_tail", "hi", "indices", "lo", "n",
  "obs", "outcome_type", "parameters", "pareto_k", "probability", "rhat",
  "scenario", "value", ".w", "variable", "y", "y_baseline", "y_index", "y_rep"
))
