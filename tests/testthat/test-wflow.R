new_wf <- function() {
  wflow("tm", models_folder = tm_models, post_sampler_params_list = tm_sampler,
        prior_sampler_params_list = tm_sampler)
}

# fitted (untuned -> warns once, then tuned) wflow, shared across tests
.wf <- new.env()
fitted_wf <- function() {
  if (is.null(.wf$wf))
    .wf$wf <- new_wf() |> attach_data(tm_data()) |> set_as_tuned() |> fit_posterior() |> quiet()
  .wf$wf
}

test_that("wflow() errors on missing model dir or missing model_code.stan", {
  root <- tempfile("models"); dir.create(file.path(root, "empty"), recursive = TRUE)
  expect_error(wflow("nope", models_folder = root), "Model directory not found")
  expect_error(wflow("empty", models_folder = root), "model_code.stan")
})

test_that("wflow() loads priors, sampler overrides and prints the model name", {
  skip_if_no_cmdstan()
  wf <- new_wf()
  expect_s3_class(wf, "wflow")
  expect_equal(wf$compiled$priors$prior_b_sd, 2)
  expect_equal(wf$post_sampler_params$iter_sampling, 200)
  expect_equal(wf$prior_sampler_params$chains, 2)
  expect_false(wf$tuned)
  expect_match(paste(capture.output(print(wf)), collapse = "\n"), "tm")
})

test_that("guards error without data or post_fit", {
  skip_if_no_cmdstan()
  wf <- new_wf()
  expect_error(fit_posterior(wf), "No data attached")
  expect_error(run_inference(wf, contrasts = list()), "not fitted")
  expect_error(run_diagnostics(wf), "not fitted")
})

test_that("fit_posterior() warns when untuned and not once tuned", {
  skip_if_no_cmdstan()
  wf <- new_wf() |> attach_data(tm_data())
  short <- list(iter_sampling = 20L, iter_warmup = 50L)
  expect_warning(quiet(fit_posterior(wf, sampler_params_list = short)), "not been marked as tuned")
  expect_no_warning(quiet(fit_posterior(set_as_tuned(wf), sampler_params_list = short)),
                    message = "tuned")
  wf <- fitted_wf()
  expect_true(wf$tuned)
  expect_s3_class(wf$post_fit, "stanmodel")
  expect_s3_class(wf$post_fit$draws, "draws")
  expect_equal(posterior::ndraws(wf$post_fit$draws), 400)
})

test_that("run_diagnostics() stores diagnostics on a fitted wflow", {
  skip_if_no_cmdstan()
  wf <- quiet(run_diagnostics(fitted_wf()))
  expect_s3_class(wf$diagnostics, "tbl_df")
  expect_true("rhat" %in% names(wf$diagnostics))
})

test_that("tune() fits the prior, runs ppc, and saves a loadable wflow", {
  skip_if_no_cmdstan()
  f <- tempfile(fileext = ".rds")
  sp <- sampler_params(modifyList(sampler_params(), tm_sampler))
  wf <- new_wf() |> attach_data(tm_data()) |>
    tune(sampler_params = sp, file = f, .by = c("group", "direction")) |> quiet()
  expect_s3_class(wf$prior_fit, "stanmodel")
  expect_equal(wf$prior_fit$stan_data$skip_likelihood, 1)
  expect_type(wf$prior_ppc_results, "list")
  expect_named(wf$prior_ppc_results, "y")
  expect_s3_class(wf$prior_ppc_results[[1]], "ppc_results_grouped")
  expect_true(file.exists(f))
  wf2 <- quiet(load_wflow(f))
  expect_identical(wf2$prior_fit$draws, wf$prior_fit$draws)
})

test_that("write()/load_wflow() round trip a fitted wflow", {
  skip_if_no_cmdstan()
  wf <- fitted_wf()
  f <- file.path(tempfile("dir"), "wf.rds")  # dir is created
  quiet(write(wf, f))
  wf2 <- quiet(load_wflow(f))
  expect_s3_class(wf2, "wflow")
  expect_identical(wf2$post_fit$draws, wf$post_fit$draws)
  expect_identical(wf2$model_name, "tm")

  g <- tempfile(fileext = ".rds"); saveRDS(wf, g)
  expect_identical(quiet(load_wflow(g))$post_fit$draws, wf$post_fit$draws)
})

test_that("load_wflow() rejects missing files and non-wflow objects", {
  expect_error(load_wflow(tempfile(fileext = ".rds")), "File not found")
  f <- tempfile(fileext = ".rds"); saveRDS(list(a = 1), f)
  expect_error(load_wflow(f))
})

test_that("files written by JWileymisc::saveRDSfst() load and re-save byte-identically", {
  fx <- test_path("fixtures", "wflow_jwileymisc_format.rds")
  expect_true(is.raw(readRDS(fx)))
  expected <- structure(list(model_name = "old_format", compiled = list(model_code = "// fixture"),
                             post_fit = list(draws = data.frame(a = c(0.5, 1.5))), inference = list()),
                        class = "wflow")
  wf <- quiet(load_wflow(fx))
  expect_identical(wf, expected)
  out <- tempfile(fileext = ".rds")
  stanworkflow:::.save_rds_fst(wf, out)
  expect_identical(readBin(fx, "raw", 1e4), readBin(out, "raw", 1e4))
})

test_that("ppc() on a fitted model returns ppc_results / grouped results and plots", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  res <- ppc(m)
  expect_s3_class(res, "ppc_results")
  expect_true(all(vapply(res$plots, inherits, logical(1), "ggplot")))
  expect_no_error(quiet(print(res, plot = TRUE)))

  g <- ppc(m, .by = c("group", "direction"))
  expect_s3_class(g, "ppc_results_grouped")
  expect_s3_class(g$summary, "data.frame")
  expect_equal(nrow(g$summary), 6)
  expect_true(all(c("group", "direction") %in% names(g$summary)))
  expect_no_error(quiet(display(g)))
})

test_that("diagnose() returns the draws summary", {
  skip_if_no_cmdstan()
  s <- quiet(diagnose(tm_fit()))
  expect_s3_class(s, "tbl_df")
  expect_true(all(c("variable", "rhat", "ess_bulk") %in% names(s)))
  expect_false(any(grepl("^(y_rep|log_lik)", s$variable)))
  expect_true(all(s$rhat < 1.05))
})

test_that("ppc_by_parameter() labels a vector parameter by factor level", {
  skip_if_no_cmdstan()
  r <- suppressWarnings(ppc_by_parameter(tm_fit(), "b", factors = list(direction = c("FLEX", "EXT"))))
  expect_s3_class(r$direction, "factor")
  expect_equal(levels(r$direction), c("FLEX", "EXT"))
  expect_equal(nrow(r), 400 * 2)
  expect_true(all(c("value", "direction") %in% names(r)))
  w <- suppressWarnings(ppc_by_parameter(tm_fit(), "b", factors = list(direction = c("FLEX", "EXT")),
                        .widen_factor = "direction"))
  expect_equal(nrow(w), 400)
  expect_true(all(c("FLEX", "EXT") %in% names(w)))
})
