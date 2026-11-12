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
  expect_s3_class(wf$diagnostics, "stanworkflow_diagnostics")
  expect_true("rhat" %in% names(wf$diagnostics$summary))
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

test_that("files written by JWileymisc::saveRDSfst() load and re-save round-trip", {
  fx <- test_path("fixtures", "wflow_jwileymisc_format.rds")
  expect_true(is.raw(readRDS(fx)))
  expected <- structure(list(model_name = "old_format", compiled = list(model_code = "// fixture"),
                             post_fit = list(draws = data.frame(a = c(0.5, 1.5))), inference = list()),
                        class = "wflow")
  wf <- quiet(load_wflow(fx))
  expect_identical(wf, expected)
  out <- tempfile(fileext = ".rds")
  stanworkflow:::.save_rds_fst(wf, out)
  expect_identical(quiet(load_wflow(out)), expected)  # compressed bytes vary by zlib/platform
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

  g1 <- ppc(m, .by = "group")                       # single .by keeps its column name
  expect_equal(names(g1$summary)[1], "group")
})

test_that("diagnose() returns the draws summary", {
  skip_if_no_cmdstan()
  s <- quiet(diagnose(tm_fit()))$summary
  expect_s3_class(s, "tbl_df")
  expect_true(all(c("variable", "rhat", "ess_bulk") %in% names(s)))
  expect_false(any(grepl("^(y_rep|log_lik)", s$variable)))
  expect_true(all(s$rhat < 1.05))
})

test_that("parameter_draws() labels a vector parameter by factor level", {
  skip_if_no_cmdstan()
  r <- suppressWarnings(parameter_draws(tm_fit(), "b", factors = list(direction = c("FLEX", "EXT"))))
  expect_s3_class(r$direction, "factor")
  expect_equal(levels(r$direction), c("FLEX", "EXT"))
  expect_equal(nrow(r), 400 * 2)
  expect_true(all(c("value", "direction") %in% names(r)))
  w <- suppressWarnings(parameter_draws(tm_fit(), "b", factors = list(direction = c("FLEX", "EXT")),
                        .widen_factor = "direction"))
  expect_equal(nrow(w), 400)
  expect_true(all(c("FLEX", "EXT") %in% names(w)))
})

test_that("wflow() warns about priors.csv values the model does not declare, and ignores them", {
  skip_if_no_cmdstan()
  root <- tempfile("models"); dir.create(file.path(root, "tm"), recursive = TRUE)
  file.copy(file.path(tm_models, "tm", "model_code.stan"), file.path(root, "tm"))
  write.csv(data.frame(name = c("prior_b_sd", "prior_b_alpha"), value = c(2, 1)),
            file.path(root, "tm", "priors.csv"), row.names = FALSE)
  expect_warning(wf <- wflow("tm", models_folder = root), "does not declare.*prior_b_alpha")
  expect_named(wf$compiled$priors, "prior_b_sd")
  expect_equal(wf$compiled$priors$prior_b_sd, 2)
})

test_that("prior_log() records each prior check; set_as_tuned() marks the last one accepted", {
  skip_if_no_cmdstan()
  sp <- sampler_params(modifyList(sampler_params(), tm_sampler))
  wf <- new_wf() |> attach_data(tm_data())
  wf <- quiet(tune(wf, sampler_params = sp))
  wf$compiled$priors$prior_b_sd <- 0.5                      # user changes a prior and re-checks
  wf <- quiet(tune(wf, sampler_params = sp, .by = "group"))
  wf <- set_as_tuned(wf)

  tab <- prior_log(wf)
  expect_equal(unique(tab$run), 1:2)
  expect_equal(tab$prior_b_sd[tab$run == 1], 2)
  expect_equal(unique(tab$prior_b_sd[tab$run == 2]), 0.5)
  expect_equal(as.vector(tapply(tab$accepted, tab$run, unique)), c(FALSE, TRUE))
  expect_true(all(c("outcome", "p_B_mean", "cov_50%") %in% names(tab)))
  expect_equal(sum(tab$run == 2), nlevels(tm_data()$group))  # one row per group when grouped

  lg <- prior_log(wf, plots = TRUE)
  expect_named(lg$plots[[1]]$y, c("density_overlay", "mean_sd_2d"))
  expect_s3_class(lg$plots[[2]]$y[[1]]$density_overlay, "ggplot")
})

test_that("priors.csv: a vector prior is written in one cell, separated by ';'", {
  skip_if_no_cmdstan()
  root <- tempfile("models"); dir.create(file.path(root, "vp"), recursive = TRUE)
  writeLines("data { int<lower=1> K; vector[K] prior_loc; real<lower=0> prior_sd; }
parameters { vector[K] b; } model { b ~ normal(prior_loc, prior_sd); }",
             file.path(root, "vp", "model_code.stan"))
  write.csv(data.frame(name = c("prior_loc", "prior_sd"), value = c("0; 0.5;1", "2")),
            file.path(root, "vp", "priors.csv"), row.names = FALSE)
  wf <- wflow("vp", models_folder = root)
  expect_equal(wf$compiled$priors$prior_loc, c(0, 0.5, 1))
  expect_equal(wf$compiled$priors$prior_sd, 2)
})

test_that("tune() uses the wflow's prior sampler settings; an argument overrides them", {
  skip_if_no_cmdstan()
  wf <- wflow("tm", models_folder = tm_models,
              prior_sampler_params_list = list(chains = 2L, parallel_chains = 2L,
                                               iter_warmup = 100L, iter_sampling = 50L))
  wf <- attach_data(wf, tm_data())
  expect_no_warning(wf1 <- quiet(tune(wf)), message = "sampler_params")
  expect_equal(posterior::ndraws(wf1$prior_fit$draws), 100L)          # 2 chains x 50
  wf2 <- quiet(tune(wf, sampler_params = list(iter_sampling = 30L)))
  expect_equal(posterior::ndraws(wf2$prior_fit$draws), 60L)
})

test_that("prior check reports coverage neutrally; accept_priors() and set_as_tuned() are the same", {
  skip_if_no_cmdstan()
  sp <- sampler_params(modifyList(sampler_params(), tm_sampler))
  wf <- quiet(tune(new_wf() |> attach_data(tm_data()), sampler_params = sp))
  txt <- unlist(wf$prior_ppc_results$y$interpretation$coverage)
  expect_match(txt, "^prior predictive (50|90)% intervals cover")
  expect_false(any(grepl("likelihood|calibrated", txt)))
  expect_true(accept_priors(wf)$tuned)
  expect_identical(set_as_tuned(wf), accept_priors(wf))
})

test_that("ppc(stats =) adds custom statistics; prior checks show p-values without labels", {
  skip_if_no_cmdstan()
  wf  <- fitted_wf()
  top <- function(y) mean(y > 2)
  p <- ppc(wf, stats = list(prop_gt2 = top))                 # read-only method on a wflow
  expect_true("prop_gt2" %in% names(p$bayes_pvalues))
  expect_s3_class(p$plots$prop_gt2, "ggplot")
  expect_match(p$interpretation$bayes_pvalues$prop_gt2, "[", fixed = TRUE)   # posterior: labelled
  g <- ppc(wf, .by = "group", stats = list(prop_gt2 = top))
  expect_true("p_B_prop_gt2" %in% names(g$summary))

  sp <- sampler_params(modifyList(sampler_params(), tm_sampler))
  wp <- quiet(tune(new_wf() |> attach_data(tm_data()), sampler_params = sp))
  txt <- unlist(wp$prior_ppc_results$y$interpretation$bayes_pvalues)
  expect_true(all(grepl("^p_B = [0-9.]+$", txt)))                 # prior: numbers only
})

test_that("ppc()/diagnose() on a wflow are read-only", {
  skip_if_no_cmdstan()
  wf <- fitted_wf()
  before <- wf
  d <- quiet(diagnose(wf, params = "sigma"))
  expect_output(print(d), "Divergences")
  expect_s3_class(ppc(wf, .by = "group"), "ppc_results_grouped")
  expect_identical(wf, before)
  expect_error(ppc(new_wf(), type = "prior"), "run `tune")
})

test_that("inference on a fit that failed diagnostics repeats the warning (does not stop)", {
  skip_if_no_cmdstan()
  wf <- quiet(run_diagnostics(fitted_wf()))
  wf$diagnostics$ok <- FALSE; wf$diagnostics$problems <- "12 divergent transitions"
  co <- contrast("x", continuous = TRUE, parameter = "b", factors = "direction", delta = "range")
  expect_warning(capture.output(suppressMessages(wf <- run_inference(wf, co, draws_method = "parameter"))),
                 "12 divergent", class = "stanworkflow_convergence_warning")
  expect_s3_class(wf$inference$main, "inference")
})

test_that("diagnostics_params: forgiving file formats, argument override, unknown names warned", {
  f <- tempfile(fileext = ".csv")
  rd <- function(lines) { writeLines(lines, f); stanworkflow:::.read_name_list(f) }
  expect_equal(rd(c('"parameters"', '"a"', '"sigma"')), c("a", "sigma"))   # old format
  expect_equal(rd(c("Parameter", "a", "sigma")), c("a", "sigma"))           # any header word
  expect_equal(rd(c("a", "sigma")), c("a", "sigma"))                        # no header
  expect_equal(rd("a, sigma; b"), c("a", "sigma", "b"))                     # one line
  expect_equal(rd(c("# check these", "a  # intercepts", "", "b[2]")), c("a", "b[2]"))

  skip_if_no_cmdstan()
  root <- tempfile("models"); dir.create(file.path(root, "tm"), recursive = TRUE)
  file.copy(file.path(tm_models, "tm", c("model_code.stan", "priors.csv")), file.path(root, "tm"))
  writeLines(c("variables", "a", "b[2]", "sigam"), file.path(root, "tm", "diagnostics_params.txt"))
  expect_warning(wf <- wflow("tm", models_folder = root), "not defined by the model.*sigam")
  expect_equal(wf$diagnostics_params, c("a", "b[2]"))
  wf <- wflow("tm", models_folder = root, diagnostics_params = c("sigma", "g"))   # argument wins
  expect_equal(wf$diagnostics_params, c("sigma", "g"))
})

test_that("wflow() prior check defaults to 500 + 500 draws; overrides still apply", {
  skip_if_no_cmdstan()
  wf <- wflow("tm", models_folder = tm_models)
  expect_equal(wf$prior_sampler_params$iter_sampling, 500L)
  expect_equal(wf$prior_sampler_params$iter_warmup, 500L)
  expect_equal(wf$post_sampler_params$iter_sampling, 1000L)
  expect_equal(wflow("tm", models_folder = tm_models,
                     prior_sampler_params_list = list(iter_sampling = 50L))$prior_sampler_params$iter_sampling, 50L)
})
