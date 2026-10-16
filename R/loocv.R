#loocv.R

# ---- loo_results class ------------------------------------------------------
#' Create a loo_results object
#' @keywords internal
loo_results <- function(loo_obj, summary, pareto_k, interpretation) {
  structure(
    list(
      loo            = loo_obj,
      summary        = summary,
      pareto_k       = pareto_k,
      interpretation = interpretation
    ),
    class = "loo_results"
  )
}

#' @exportS3Method
print.loo_results <- function(x, plot = FALSE, ...) {
  cat("<loo_results>\n")

  cat("\n--- Scalar metrics ---\n")
  cat(sprintf("  %-10s  %10s  %10s\n", "metric", "estimate", "se"))
  cat(sprintf("  %-10s  %10s  %10s\n", "------", "--------", "--"))
  for (i in seq_len(nrow(x$summary))) {
    cat(sprintf("  %-10s  %10.2f  %10.2f\n",
                x$summary$metric[i],
                x$summary$estimate[i],
                x$summary$se[i]))
  }

  cat("\n--- Interpretation ---\n")
  for (nm in names(x$interpretation)) {
    cat(sprintf("  [%s]\n    %s\n", nm, x$interpretation[[nm]]))
  }

  k_counts <- table(x$pareto_k$category)
  cat("\n--- Pareto-k summary (McElreath SR2 §7.5) ---\n")
  for (nm in names(k_counts)) {
    cat(sprintf("  %-26s: %d\n", nm, k_counts[[nm]]))
  }

  if (plot) print(x$pareto_k |> .plot_pareto_k())

  invisible(x)
}

# ---- loo_cv helpers ---------------------------------------------------------


#' @keywords internal
.interpret_loo_metrics <- function(loo_obj, n_params) {
  est       <- loo_obj$estimates
  elpd      <- est["elpd_loo", "Estimate"]
  elpd_se   <- est["elpd_loo", "SE"]
  p_loo_val <- est["p_loo",    "Estimate"]
  looic_val <- est["looic",    "Estimate"]

  # ---- elpd_loo -----------------------------------------------------------
  # Gelman et al. (2020) §7: elpd_loo approximates expected log pointwise
  # predictive density for new data; interpretation is relative across models.
  # A difference > 4 SE between models is considered practically meaningful.
  elpd_interp <- sprintf(
    "elpd_loo = %.2f (SE %.2f)\n    Higher is better. Meaningful cross-model differences are >> 4 * SE.\n    looic = %.2f (deviance scale: lower is better).",
    elpd, elpd_se, looic_val
  )

  # ---- p_loo --------------------------------------------------------------
  # McElreath SR2 §7.4: p_loo is the effective number of parameters —
  # analogous to the penalty term in WAIC. p_loo >> nominal parameter count
  # signals overfitting, influential observations, or model misspecification
  # (common with heavy-tailed or poorly-identified parameters).
  p_ratio <- p_loo_val / max(n_params, 1)
  p_flag  <- if (p_ratio > 2)
    "WARNING: p_loo >> nominal parameter count — possible overfitting or\n    misspecification (McElreath SR2 §7.4, Gelman 2020 §7)."
  else if (p_ratio > 1)
    "Slightly elevated relative to nominal; inspect high-k observations."
  else
    "Plausible relative to model size — penalty consistent with complexity."

  p_loo_interp <- sprintf(
    "p_loo = %.2f  (nominal scalar params ≈ %d, ratio %.1f)\n    %s",
    p_loo_val, n_params, p_ratio, p_flag
  )

  list(elpd_loo = elpd_interp, p_loo = p_loo_interp)
}


#' @keywords internal
.categorise_pareto_k <- function(k_vals) {
  # McElreath SR2 §7.5 & Vehtari et al. (2017):
  #   k < 0.5              → PSIS estimate reliable
  #   0.5 ≤ k < 0.7        → OK; some caution warranted
  #   0.7 ≤ k < 1.0        → bad; high-influence observation; IS unreliable
  #   k ≥ 1.0              → very bad; PSIS fails; moments may not exist
  dplyr::tibble(obs = seq_along(k_vals), pareto_k = k_vals) |>
    dplyr::mutate(
      category = dplyr::case_when(
        pareto_k <  0.5 ~ "good (k < 0.5)",
        pareto_k <  0.7 ~ "ok (0.5 \u2264 k < 0.7)",
        pareto_k <  1.0 ~ "bad (0.7 \u2264 k < 1.0)",
        TRUE            ~ "very bad (k \u2265 1.0)"
      ),
      category = factor(category, levels = c(
        "good (k < 0.5)", "ok (0.5 \u2264 k < 0.7)",
        "bad (0.7 \u2264 k < 1.0)", "very bad (k \u2265 1.0)"
      ))
    )
}



#' @keywords internal
.interpret_pareto_k <- function(pareto_k_tbl) {
  k       <- pareto_k_tbl$pareto_k
  n_bad   <- sum(k >= 0.7)
  n_vbad  <- sum(k >= 1.0)

  # McElreath SR2 §7.5: high-k observations are "influential" — the posterior
  # would shift substantially if they were removed. They are not necessarily
  # outliers in the residual sense but deserve substantive scrutiny.
  dplyr::case_when(
    n_vbad > 0 ~ sprintf(
      "%d obs with k \u2265 1.0: PSIS approximation unreliable for these points.\n    PSIS-LOO elpd estimate may be biased. Options (Gelman 2020 §7):\n      (a) refit with moment_matching = TRUE\n      (b) exact LOO for the %d flagged observations\n      (c) investigate whether these are genuine outliers warranting\n          a more robust likelihood (e.g. Student-t).",
      n_vbad, n_vbad
    ),
    n_bad > 0 ~ sprintf(
      "%d obs with k \u2265 0.7: high-influence observations detected.\n    PSIS estimate less reliable; consider moment_matching = TRUE.\n    Substantively: inspect these observations — are they outliers,\n    leverage points, or evidence of model misspecification?\n    (McElreath SR2 §7.5; Gelman 2020 §7).",
      n_bad
    ),
    TRUE ~ sprintf(
      "All Pareto-k < 0.7 (%d good, %d ok).\n    PSIS-LOO reliable; no influential observations flagged.",
      sum(k < 0.5), sum(k >= 0.5 & k < 0.7)
    )
  )
}


#' @keywords internal
.plot_pareto_k <- function(pareto_k_tbl) {
  k_colours <- c(
    "good (k < 0.5)"         = "#2166ac",
    "ok (0.5 \u2264 k < 0.7)" = "#92c5de",
    "bad (0.7 \u2264 k < 1.0)" = "#f4a582",
    "very bad (k \u2265 1.0)"  = "#d6604d"
  )
  ggplot2::ggplot(
    pareto_k_tbl,
    ggplot2::aes(x = obs, y = pareto_k, colour = category)
  ) +
    ggplot2::geom_point(size = 1.4, alpha = 0.75) +
    ggplot2::geom_hline(
      yintercept = c(0.5, 0.7, 1.0),
      linetype   = "dashed",
      colour     = "grey50",
      linewidth  = 0.4
    ) +
    ggplot2::scale_colour_manual(values = k_colours, name = "Pareto-k") +
    ggplot2::labs(
      title    = "PSIS-LOO Pareto-k diagnostics",
      subtitle = "k < 0.5 reliable; 0.5\u20130.7 OK; > 0.7 high influence; > 1.0 PSIS fails\n(McElreath SR2 \u00a77.5; Vehtari et al. 2017)",
      x        = "Observation index",
      y        = "Pareto k"
    ) +
    ggplot2::theme_bw()
}


#' LOO-CV diagnostics for a fitted stanmodel
#'
#' Computes Pareto-smoothed importance sampling LOO-CV (PSIS-LOO) via the
#' `loo` package. Requires `log_lik` in `generated quantities` (length N).
#'
#' **Conceptual grounding:**
#'
#' *Gelman et al. (2020) §7 — Model comparison and stacking:*
#' LOO-CV is the preferred comparison tool in the Bayesian workflow because it
#' estimates out-of-sample predictive accuracy without refitting the model. The
#' workflow recommendation is to first check each model's absolute fit via
#' `ppc()`, then use LOO to compare models only after both pass their
#' predictive checks. A model that wins LOO but fails its PPC is not to be
#' trusted. `loo_compare()` differences should be interpreted relative to their
#' SE: a difference < 2 SE is not decisive.
#'
#' *McElreath SR2 §7 — Ulysses' Compass:*
#' Information criteria penalise complexity to balance underfitting and
#' overfitting (the bias–variance tradeoff McElreath frames as "Scylla and
#' Charybdis"). PSIS-LOO approximates LOOCV and is preferable to WAIC when
#' Pareto-k diagnostics are examined alongside it. `p_loo` is the "effective
#' number of parameters" — McElreath's penalty term. When `p_loo` >> the
#' nominal parameter count the model is either overfitting or misspecified
#' (SR2 §7.5). The Pareto-k per-observation diagnostic directly reveals
#' *which* observations drive complexity; McElreath recommends inspecting
#' high-k points substantively rather than just statistically.
#'
#' **Metrics:**
#' - *elpd_loo*: expected log pointwise predictive density; higher = better
#'   out-of-sample accuracy (SR2 §7.3, Gelman §7).
#' - *p_loo*: effective parameters. Compare to nominal count (Gelman §7,
#'   SR2 §7.5).
#' - *looic*: -2 × elpd_loo, deviance scale. Lower = better. Use for
#'   `loo_compare()` across models fit to the same data (SR2 §7.4).
#' - *Pareto-k*: per-observation stability. k < 0.5 good; 0.5–0.7 OK;
#'   0.7–1.0 problematic (high-influence observation); > 1.0 PSIS fails.
#'   McElreath SR2 §7.5: high-k points are not outliers to discard — they are
#'   the most *informative* observations and deserve substantive scrutiny.
#'
#' @param x a fitted `stanmodel` object
#' @param moment_matching logical; apply moment matching for k > 0.7 obs
#'   (requires `loo` >= 2.4.0). Gelman §7 recommends this before discarding
#'   the LOO estimate. Default FALSE.
#' @param plot logical; render Pareto-k plot. Default FALSE.
#' @param r_eff logical; compute relative ESS correction for MCMC draws
#'   (recommended; improves Pareto-k accuracy). Default TRUE.
#' @return invisibly, a named list:
#'   \describe{
#'     \item{`loo`}{raw `loo` object from `loo::loo()`}
#'     \item{`summary`}{tibble of elpd_loo, p_loo, looic with SE}
#'     \item{`pareto_k`}{tibble of per-observation k values with category}
#'     \item{`interpretation`}{named list of human-readable strings}
#'   }
#' @references
#'   Gelman et al. (2020). Bayesian Workflow. arXiv:2011.01808 §7.
#'   McElreath (2020). Statistical Rethinking 2e. CRC Press. §7.
#'   Vehtari, Gelman & Gabry (2017). Practical Bayesian model evaluation
#'     using LOO-CV and WAIC. *Stat Comput* 27(5).
#' @importFrom posterior subset_draws as_draws_matrix
#' @importFrom dplyr tibble mutate case_when
#' @importFrom ggplot2 ggplot aes geom_point geom_hline scale_colour_manual labs theme_bw annotate
#' @export
loo_cv <- function(x, ...) UseMethod("loo_cv")

#' @exportS3Method
loo_cv.stanmodel <- function(x,
                              moment_matching = FALSE,
                              plot            = FALSE,
                              r_eff           = TRUE) {

  if (is.null(x$draws))        stop("No draws found. Run `fit()` first.")
  if (is.null(x$fitted_model)) stop("No fitted_model found. Run `fit()` first.")
  if (!"log_lik" %in% x$variables) {
    stop(
      "No `log_lik` variable found in draws.\n",
      "Add to your Stan model's generated quantities block:\n\n",
      "  vector[N] log_lik;\n",
      "  for (n in 1:N)\n",
      "    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);\n"
    )
  }

  # ---- extract log_lik [S x N] --------------------------------------------
  log_lik_mat <- posterior::subset_draws(x$draws, variable = "log_lik") |>
    posterior::as_draws_matrix()                                  # [S × N]

  # ---- relative ESS correction (Vehtari et al. 2017) ---------------------
  # Corrects for the fact that MCMC draws are autocorrelated, not i.i.d.
  r_eff_vals <- if (r_eff) {
    n_chains   <- x$fitted_model$num_chains()
    n_iter     <- x$fitted_model$metadata()$iter_sampling
    loo::relative_eff(
      exp(log_lik_mat),
      chain_id = rep(seq_len(n_chains), each = n_iter)
    )
  } else {
    rep(1, ncol(log_lik_mat))
  }

  # ---- PSIS-LOO (McElreath SR2 §7.5, Gelman §7) --------------------------
  loo_obj <- loo::loo(log_lik_mat, r_eff = r_eff_vals)

  # Moment matching: Gelman §7 — before concluding LOO is unreliable for
  # high-k observations, attempt to stabilise via moment matching.
  if (moment_matching) {
    loo_obj <- loo::loo(
      x$fitted_model,
      moment_match = TRUE,
      r_eff        = r_eff_vals
    )
  }

  # ---- scalar summary tibble ----------------------------------------------
  est  <- loo_obj$estimates
  summ <- dplyr::tibble(
    metric   = rownames(est),
    estimate = est[, "Estimate"],
    se       = est[, "SE"]
  )

  # ---- per-observation Pareto-k tibble ------------------------------------
  # McElreath SR2 §7.5: treat high-k points as substantively interesting,
  # not as errors to remove.
  k_vals   <- loo_obj$diagnostics$pareto_k
  pareto_k <- dplyr::tibble(obs = seq_along(k_vals), pareto_k = k_vals) |>
    dplyr::mutate(
      category = dplyr::case_when(
        pareto_k < 0.5 ~ "good (< 0.5)",
        pareto_k < 0.7 ~ "ok (0.5–0.7)",
        pareto_k < 1.0 ~ "bad (0.7–1.0)",
        TRUE           ~ "very bad (> 1.0)"
      ),
      category = factor(category, levels = c(
        "good (< 0.5)", "ok (0.5–0.7)", "bad (0.7–1.0)", "very bad (> 1.0)"
      ))
    )

  # ---- interpretation -----------------------------------------------------
  n_bad      <- sum(k_vals >= 0.7)
  n_very_bad <- sum(k_vals >= 1.0)
  p_loo_val  <- est["p_loo",  "Estimate"]
  elpd_val   <- est["elpd_loo", "Estimate"]
  elpd_se    <- est["elpd_loo", "SE"]
  looic_val  <- est["looic", "Estimate"]
  looic_se   <- est["looic", "SE"]

    # Nominal scalar parameter count — exclude observation-level auxiliaries
  # Must count *scalar* parameters from draws, not top-level variable names
  n_params <- posterior::subset_draws(
    x$draws,
    variable = x$variables[!grepl("^y_rep$|^log_lik$", x$variables)]
  ) |>
    posterior::as_draws_matrix() |>
    ncol()

  interp <- list(

    elpd_loo = sprintf(
      paste0(
        "elpd_loo = %.2f (SE %.2f).\n",
        "  Higher = better out-of-sample predictive accuracy (McElreath SR2 §7.3).\n",
        "  Interpret differences between models relative to SE: |Δelpd| < 2×SE is not decisive\n",
        "  (Gelman 2020 §7). Use loo::loo_compare() for formal comparison."
      ),
      elpd_val, elpd_se
    ),

    looic = sprintf(
      paste0(
        "looic = %.2f (SE %.2f).\n",
        "  Deviance scale: -2 × elpd_loo. Lower = better.\n",
        "  Comparable to AIC/DIC but estimates true out-of-sample deviance\n",
        "  (McElreath SR2 §7.4). Only meaningful in comparison across models\n",
        "  fit to identical data."
      ),
      looic_val, looic_se
    ),

    p_loo = {
      ratio <- p_loo_val / max(n_params, 1)
      severity <- dplyr::case_when(
        ratio > 3 ~
          "SEVERE: p_loo >> nominal params — strong evidence of overfitting or\n  model misspecification. The posterior is too flexible relative to the data\n  (McElreath SR2 §7.5, Gelman 2020 §7). Revisit the likelihood or priors.",
        ratio > 2 ~
          "WARNING: p_loo moderately exceeds nominal params — possible misspecification.\n  Run ppc() to check where the model fails (Gelman 2020 §4 + §7).",
        TRUE ~
          "p_loo plausible relative to nominal parameter count — no obvious overfitting."
      )
      sprintf(
        "p_loo = %.2f (nominal params ≈ %d, ratio = %.1f).\n  %s",
        p_loo_val, n_params, ratio, severity
      )
    },

    pareto_k = dplyr::case_when(
      n_very_bad > 0 ~ sprintf(
        paste0(
          "%d observation(s) with k > 1.0: PSIS approximation has failed for these points.\n",
          "  %d obs have k > 0.7 total.\n",
          "  Action (Gelman 2020 §7): (1) try moment_matching = TRUE; (2) refit without\n",
          "  those observations for exact LOO; (3) inspect them substantively —\n",
          "  McElreath SR2 §7.5 notes high-k points reveal where the model is most\n",
          "  sensitive to individual data and may signal model misspecification."
        ),
        n_very_bad, n_bad
      ),
      n_bad > 0 ~ sprintf(
        paste0(
          "%d observation(s) with k > 0.7 (high influence on posterior).\n",
          "  PSIS estimate is unreliable for these points.\n",
          "  Action: try moment_matching = TRUE, or inspect the high-k observations\n",
          "  substantively (McElreath SR2 §7.5). These are the observations the model\n",
          "  learns from most — they deserve scientific scrutiny, not automatic removal."
        ),
        n_bad
      ),
      TRUE ~ sprintf(
        paste0(
          "All Pareto-k < 0.7 (%d < 0.5; %d in 0.5–0.7).\n",
          "  PSIS-LOO approximation is reliable (Vehtari et al. 2017).\n",
          "  The LOO estimate can be trusted for model comparison (Gelman 2020 §7)."
        ),
        sum(k_vals < 0.5), sum(k_vals >= 0.5 & k_vals < 0.7)
      )
    )
  )

  # ---- console report ------------------------------------------------------
  cat("<loo_cv>\n")

  cat("\n--- Scalar metrics ---\n")
  cat(sprintf("  %-12s %10s %10s\n", "metric", "estimate", "se"))
  cat(sprintf("  %-12s %10s %10s\n", "------",  "--------", "--"))
  for (i in seq_len(nrow(summ))) {
    cat(sprintf("  %-12s %10.2f %10.2f\n",
                summ$metric[i], summ$estimate[i], summ$se[i]))
  }

  cat("\n--- Pareto-k summary ---\n")
  k_counts <- table(pareto_k$category)
  for (nm in names(k_counts)) cat(sprintf("  %-22s : %d\n", nm, k_counts[[nm]]))

  cat("\n--- Interpretation (Gelman 2020 §7; McElreath SR2 §7) ---\n")
  for (nm in names(interp)) {
    cat(sprintf("\n[%s]\n  %s\n", nm, interp[[nm]]))
  }

  # ---- optional Pareto-k plot (Gabry et al. 2019 recommendation) ----------
  if (plot) {
    k_colours <- c(
      "good (< 0.5)"     = "#2166ac",
      "ok (0.5–0.7)"     = "#92c5de",
      "bad (0.7–1.0)"    = "#f4a582",
      "very bad (> 1.0)" = "#d6604d"
    )
    g <- ggplot2::ggplot(
      pareto_k,
      ggplot2::aes(x = obs, y = pareto_k, colour = category)
    ) +
      ggplot2::geom_point(size = 1.4, alpha = 0.75) +
      ggplot2::geom_hline(
        yintercept = c(0.5, 0.7, 1.0),
        linetype   = "dashed",
        colour     = "grey40"
      ) +
      ggplot2::annotate(
        "text", x = 1, y = c(0.52, 0.72, 1.02),
        label  = c("k = 0.5 (good/ok)", "k = 0.7 (ok/bad)", "k = 1.0 (PSIS fails)"),
        hjust  = 0, size = 3, colour = "grey30"
      ) +
      ggplot2::scale_colour_manual(values = k_colours, name = "Pareto-k") +
      ggplot2::labs(
        title    = "PSIS-LOO Pareto-k diagnostics",
        subtitle = paste0(
          "k < 0.5: reliable  |  0.5–0.7: ok  |  > 0.7: high influence  |  > 1.0: PSIS fails\n",
          "High-k points = most informative observations (McElreath SR2 §7.5) — inspect substantively."
        ),
        x = "Observation index",
        y = "Pareto k"
      ) +
      ggplot2::theme_bw()
    print(g)
  }

  invisible(list(
    loo            = loo_obj,
    summary        = summ,
    pareto_k       = pareto_k,
    interpretation = interp
  ))
}