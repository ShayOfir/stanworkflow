# Counterfactual contrasts on the shared linear model (helper-model.R).
# mu_rep is linear in x, so a contrast between two x grid points is
# b[direction] * dx in every draw; mu0_rep does not depend on x, so its
# contrast is exactly 0 whenever baseline and grid point use the same draw;
# a contrast between two groups is g[k] - g[j] in every draw.

cf_contrast <- function(gq_var, ...) {
  contrast(contrast_var = "x", continuous = TRUE, two_way = "direction",
           two_way_levels = c("FLEX", "EXT"), outcome = outcome_spec(gq_var = gq_var), ...)
}

cf_run <- function(m, co, ...) {
  quiet(infer_contrasts(m, co, draws_method = "counterfactual", gq_steps = 2, ...))
}

test_that("counterfactual contrasts pair the same draws when gq_n_draws is set", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  sd0 <- cf_run(m, cf_contrast("mu0_rep"), gq_n_draws = 50)$scenario_draws
  expect_equal(nrow(sd0), 50L)
  expect_true(all(sd0$contrast_FLEX == 0) && all(sd0$contrast_EXT == 0))
})

test_that("two-way counterfactual contrast on mu_rep equals b * dx draw by draw", {
  skip_if_no_cmdstan()
  m  <- tm_fit()
  b  <- tm_draws(m, "b")
  dx <- diff(range(m$data$x))
  res <- cf_run(m, cf_contrast("mu_rep"))
  sd1 <- dplyr::arrange(res$scenario_draws, .draw)
  expect_equal(nrow(sd1), posterior::ndraws(m$draws))
  expect_lt(max(abs(sd1$contrast_FLEX - b[["b[1]"]][sd1$.draw] * dx),
                abs(sd1$contrast_EXT  - b[["b[2]"]][sd1$.draw] * dx)), 1e-4)

  # scenarios follow the signs of the two contrasts; probabilities sum to 1
  both <- sd1$contrast_FLEX > 0 & sd1$contrast_EXT > 0
  expect_identical(sd1$scenario == "Affects both", both)
  expect_equal(sum(res$scenario_grid$probability), 1)
  expect_equal(res$baseline_val, min(m$data$x))
  expect_output(print(res), "2-way")
})

test_that("finger-paper pattern: filter_draws quosure, grouping_vars, effect_positive = FALSE, custom grid data", {
  skip_if_no_cmdstan()
  m  <- tm_fit()
  b  <- tm_draws(m, "b")
  # as in ii_tract_counterfactual.R: data whose x range is the two loads of interest
  lohi   <- c(0.2, 0.8)
  d_grid <- m$data
  d_grid$x <- c(lohi, rep(lohi[1], nrow(d_grid) - 2))
  co <- cf_contrast("mu_rep", filter_draws = rlang::quo(group != "c"),
                    grouping_vars = "group", effect_positive = FALSE)
  res <- quiet(infer_contrasts(m, contrasts = co, data = d_grid,
                               draws_method = "counterfactual", gq_steps = 2))
  sd <- res$scenario_draws
  expect_setequal(unique(sd$group), c("s1", "s2"))
  expect_equal(nrow(sd), 2L * posterior::ndraws(m$draws))
  expect_lt(max(abs(sd$contrast_FLEX - b[["b[1]"]][sd$.draw] * diff(lohi)),
                abs(sd$contrast_EXT  - b[["b[2]"]][sd$.draw] * diff(lohi))), 1e-4)
  # effect_positive = FALSE: "affects" means the contrast is negative
  expect_identical(sd$scenario == "Affects both", sd$contrast_FLEX < 0 & sd$contrast_EXT < 0)
  expect_identical(sd$scenario == "Affects EXT only", sd$contrast_FLEX >= 0 & sd$contrast_EXT < 0)
  p <- dplyr::summarise(res$scenario_grid, p = sum(probability), .by = group)
  expect_equal(p$p, c(1, 1))
})

test_that("counterfactual GQ runs on the supplied data, with the fitted factor coding", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  b <- tm_draws(m, "b")
  # every row set to EXT (2nd level), as a factor with one level: the contrast must be
  # b[EXT] * dx for every row, so the per-draw mean is b[2] * dx
  d_ext <- m$data
  d_ext$direction <- factor("EXT")
  co <- contrast(contrast_var = "x", continuous = TRUE, outcome = outcome_spec(gq_var = "mu_rep"))
  r <- quiet(infer_contrasts(m, co, data = d_ext, draws_method = "counterfactual", gq_steps = 2))
  gp <- r$grid_posterior
  expect_lt(max(abs(gp$contrast - b[["b[2]"]][gp$.draw] * diff(range(m$data$x)))), 1e-4)
})

test_that("counterfactual on y_rep warns and points to mu_rep", {
  skip_if_no_cmdstan()
  expect_warning(cf_run(tm_fit(), cf_contrast("y_rep"), gq_n_draws = 10), "mu_rep")
})

test_that("one-way continuous counterfactual differences each grid point against the baseline", {
  skip_if_no_cmdstan()
  m  <- tm_fit()
  b  <- tm_draws(m, "b")
  xr <- range(m$data$x)
  co <- contrast(contrast_var = "x", continuous = TRUE, grouping_vars = "direction",
                 outcome = outcome_spec(gq_var = "mu_rep"))
  r <- quiet(infer_contrasts(m, co, draws_method = "counterfactual", gq_steps = 3))
  gp <- r$grid_posterior |> dplyr::filter(direction == "FLEX")
  expect_equal(nrow(gp), 2L * posterior::ndraws(m$draws))   # grid points 2 and 3, no baseline row
  expect_lt(max(abs(gp$contrast - b[["b[1]"]][gp$.draw] * (gp$contrast_var_value - xr[1]))), 1e-4)
  expect_equal(nrow(r$grid_summary), 4L)                       # 2 grid points x 2 directions
  expect_true(all(r$grid_summary$lo <= r$grid_summary$mean & r$grid_summary$mean <= r$grid_summary$hi))
  expect_output(print(r), "Grid Summary")
})

test_that("categorical counterfactual contrasts equal the group effect difference, one result per contrast", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  g <- tm_draws(m, "g")
  cs <- list(contrast("group", contrast_levels = c("c", "s1"), outcome = outcome_spec(gq_var = "mu_rep")),
             contrast("group", contrast_levels = c("c", "s2"), outcome = outcome_spec(gq_var = "mu_rep")))
  r <- cf_run(m, cs)
  expect_length(r, 2L)
  for (k in 1:2) {
    p <- r[[k]]$posterior
    expect_equal(nrow(p), posterior::ndraws(m$draws))
    expect_equal(unique(p$contrast_label), sprintf("group=s%d - c", k))
    expect_lt(max(abs(p$contrast - (g[[sprintf("g[%d]", k + 1)]] - g[["g[1]"]])[p$.draw])), 1e-4)
    expect_equal(nrow(r[[k]]$summary), 1L)
    expect_output(print(r[[k]]), "counterfactual (categorical)", fixed = TRUE)
  }
})

test_that("categorical counterfactual requires contrast_levels", {
  skip_if_no_cmdstan()
  expect_error(cf_run(tm_fit(), contrast("group", outcome = outcome_spec(gq_var = "mu_rep"))))
})

test_that("observation index is read from the brackets, not the variable name", {
  d <- posterior::draws_df(`mu0_rep[1]` = 1:2, `mu0_rep[2]` = 3:4, `mu0_rep[3]` = 5:6)
  m <- list(draws = d, variables = "mu0_rep")
  out <- stanworkflow:::.extract_spec_draws(m, data.frame(g = 1:3),
                                            outcome_spec(gq_var = "mu0_rep"), "g")
  expect_equal(out$.obs_idx, rep(1:3, 2))
  expect_equal(out$g, out$.obs_idx)
})
