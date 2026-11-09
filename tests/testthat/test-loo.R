test_that("loo_cv() matches loo::loo() computed directly", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  res <- suppressWarnings(quiet(loo_cv(m)))
  expect_s3_class(res$loo, "psis_loo")
  expect_true(all(c("elpd_loo", "p_loo", "looic") %in% res$summary$metric))
  expect_equal(nrow(res$pareto_k), 48)
  expect_s3_class(res$pareto_k$category, "factor")
  expect_named(res$interpretation, c("elpd_loo", "looic", "p_loo", "pareto_k"))

  ll <- posterior::as_draws_matrix(posterior::subset_draws(m$draws, "log_lik"))
  ref <- suppressWarnings(loo::loo(ll, r_eff = loo::relative_eff(exp(ll), chain_id = rep(1:2, each = 200))))
  expect_equal(unname(res$summary$estimate[res$summary$metric == "elpd_loo"]),
               unname(ref$estimates["elpd_loo", "Estimate"]), tolerance = 1e-8)
})

test_that("loo_cv() errors without log_lik or without draws", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  m$variables <- setdiff(m$variables, "log_lik")
  expect_error(loo_cv(m), "log_lik")
  m <- tm_fit(); m$draws <- NULL
  expect_error(loo_cv(m), "Run `fit()` first", fixed = TRUE)
})
