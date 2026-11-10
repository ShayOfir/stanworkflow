# regimen lifecycle: construct -> run -> results, checkpoints, resume, quick_fit,
# serialisation. The usage patterns of the finger paper (regimen() without
# contrasts + run(steps = c("prepare", "quick_fit")), run(c("prepare", "tune")),
# load_wflow() of the checkpoints) are tested explicitly.

# the paper's two-way counterfactual contrast on mu_rep
paper_contrast <- function() {
  contrast(contrast_var = "x", continuous = TRUE, two_way = "direction",
           two_way_levels = c("FLEX", "EXT"), effect_positive = FALSE,
           outcome = outcome_spec(gq_var = "mu_rep", label = "expected y"))
}

tm_regimen <- function(contrasts = NULL, results_path = tempfile("reg_"), ...) {
  regimen(model_name = "tm", data_list = list(data = tm_data()), contrasts = contrasts,
          models_folder = tm_models, results_path = results_path,
          prior_sampler_params = tm_sampler, post_sampler_params = tm_sampler,
          gq_steps = 2L, ...)
}

not_tuned <- "not been marked as tuned"

# ---- no Stan -------------------------------------------------------------------

test_that("regimen() validates input and starts unprepared", {
  expect_error(regimen("tm", data_list = list(x = 1)))
  expect_error(regimen("tm", data_list = list(data = tm_data()), contrasts = list(a = 1)), "contrast")
  expect_error(regimen("tm", data_list = list(data = tm_data()), contrasts = list(paper_contrast())),
               "named")
  expect_error(regimen("tm", data_list = list(data = tm_data()), draws_method = "nope"))
  reg <- tm_regimen(list(cf = paper_contrast()), results_path = "res_dir")
  expect_s3_class(reg, "regimen")
  expect_false(any(unlist(reg$status)))
  expect_equal(reg$paths$posterior, file.path("res_dir", "tm_post.rds"))
  expect_equal(reg$post_sampler_params$iter_sampling, 200L)
  expect_equal(reg$post_sampler_params$adapt_delta, 0.8)   # package default, not regimen()'s 0.99
  expect_s3_class(reg$post_sampler_params, "sampler_params")
  expect_output(print(reg), "contrasts     : cf")
  expect_output(print(tm_regimen()), "regimen")             # no contrasts
})

test_that("update_contrasts(), update_data(), get_inference() and save/load_regimen()", {
  reg <- tm_regimen(list(cf = paper_contrast()))
  reg$status <- list(prepared = TRUE, tuned = TRUE, fitted = TRUE, inferred = TRUE)
  reg$inference <- list(cf = "result")
  expect_equal(get_inference(reg, "cf"), "result")
  expect_error(get_inference(reg, "zz"), "Available: cf")

  r2 <- suppressMessages(update_contrasts(reg, list(other = paper_contrast())))
  expect_named(r2$contrasts, "other")
  expect_length(r2$inference, 0L)
  expect_true(r2$status$fitted && !r2$status$inferred)
  expect_error(update_contrasts(reg, list(a = 1)), "contrast")

  r3 <- suppressMessages(update_data(reg, list(data = tm_data()[1:10, ])))
  expect_false(any(unlist(r3$status)))
  expect_null(r3$wf)
  expect_equal(nrow(r3$data_list$data), 10L)

  f <- file.path(tempfile(), "sub", "reg.rds")              # directory is created
  expect_message(save_regimen(reg, f), "saved")
  expect_identical(load_regimen(f), reg)
  saveRDS(list(a = 1), f2 <- tempfile(fileext = ".rds"))
  expect_error(load_regimen(f2))
})

test_that("run() rejects unknown and conflicting steps; steps check their prerequisites", {
  reg <- tm_regimen()
  expect_error(run(reg, steps = "bogus"), "Unknown step")
  expect_error(run(reg, steps = c("fit", "quick_fit")), "Cannot run both")
  expect_error(tune(reg), "not prepared")
  expect_error(suppressMessages(infer(reg)), "not fitted")
})

# ---- CmdStan -------------------------------------------------------------------

test_that("full lifecycle: prepare -> tune -> fit -> infer, checkpoints, results", {
  skip_if_no_cmdstan()
  tm_compiled()   # compile once (shared executable)
  reg <- tm_regimen(list(cf = paper_contrast()), ppc_by = c("group", "direction"))
  expect_no_warning(quiet(reg <- run(reg, steps = c("prepare", "tune", "fit", "infer"))),
                    message = not_tuned)                       # tune() counts as acceptance

  expect_true(all(unlist(reg$status)))
  expect_true(all(file.exists(unlist(reg$paths))))
  expect_output(print(reg), "inference keys: cf")

  # the inference matches the closed form (contrast = b * dx; affects = negative)
  inf <- get_inference(reg, "cf")
  expect_s3_class(inf, "inference")
  b  <- tm_draws(reg$wf$post_fit, "b")
  sd <- inf$scenario_draws
  expect_lt(max(abs(sd$contrast_FLEX - b[["b[1]"]][sd$.draw] * diff(range(tm_data()$x)))), 1e-4)
  expect_identical(readRDS(reg$paths$inference), reg$inference)

  # checkpoints load back as wflows with the pieces each step adds
  pr <- suppressMessages(load_wflow(reg$paths$prior))
  expect_s3_class(pr$prior_fit, "stanmodel")
  expect_equal(pr$prior_fit$stan_data$skip_likelihood, 1L)
  expect_equal(pr$compiled$priors$prior_b_sd, 2)              # from priors.csv
  expect_s3_class(pr$prior_ppc_results[[1]], "ppc_results_grouped")
  po <- suppressMessages(load_wflow(reg$paths$posterior))
  expect_equal(po$post_fit$stan_data$skip_likelihood, 0L)
  expect_equal(posterior::ndraws(po$post_fit$draws), 400L)
  expect_identical(po$post_fit$draws, reg$wf$post_fit$draws)
  expect_s3_class(po$posterior_ppc_results[[1]], "ppc_results_grouped")
  expect_true(po$tuned)                                       # fit after tune: no "not tuned" warning
  expect_true(all(c("variable", "rhat") %in% names(po$diagnostics)))

  # resume() on a fresh regimen object finds every checkpoint and refits nothing
  t_post <- file.mtime(reg$paths$posterior)
  fresh <- tm_regimen(list(cf = paper_contrast()), results_path = reg$results_path)
  quiet(fresh <- resume(fresh))
  expect_true(all(unlist(fresh$status)))
  expect_identical(fresh$inference, reg$inference)
  expect_equal(file.mtime(reg$paths$posterior), t_post)

  # with the inference checkpoint gone, resume() re-runs infer only, from the saved posterior
  unlink(reg$paths$inference)
  quiet(again <- resume(tm_regimen(list(cf = paper_contrast()), results_path = reg$results_path)))
  expect_true(file.exists(reg$paths$inference))
  expect_equal(file.mtime(reg$paths$posterior), t_post)
  expect_equal(again$inference$cf$scenario_draws, reg$inference$cf$scenario_draws)

  # complete regimen: nothing to do
  expect_message(resume(again), "already complete")
})

test_that("finger-paper pattern: regimen without contrasts, run(prepare, quick_fit), save = FALSE", {
  skip_if_no_cmdstan()
  tm_compiled()
  rp  <- tempfile("reg_")
  reg <- tm_regimen(results_path = rp, ppc_by = c("group", "direction"))
  expect_warning(quiet(reg <- run(reg, steps = c("prepare", "quick_fit"), save = FALSE)), not_tuned)
  expect_true(reg$status$prepared && reg$status$fitted)
  expect_false(reg$status$tuned || reg$status$inferred)
  expect_length(reg$inference, 0L)
  expect_length(list.files(rp), 0L)                 # nothing written
  sm <- reg$wf$post_fit                             # what the paper uses
  expect_s3_class(sm, "stanmodel")
  expect_equal(posterior::ndraws(sm$draws), 400L)
  expect_null(reg$wf$posterior_ppc_results)          # quick_fit skips checks
  expect_null(reg$wf$diagnostics)

  # infer() without contrasts is a no-op that still marks the regimen inferred
  quiet(r0 <- infer(reg, save = FALSE, load_post = FALSE))
  expect_length(r0$inference, 0L)
  expect_true(r0$status$inferred)
})

test_that("run(steps = 'quick_fit') prepares, saves the posterior, and infer() picks it up", {
  skip_if_no_cmdstan()
  tm_compiled()
  reg <- tm_regimen()
  expect_warning(quiet(reg <- run(reg, steps = "quick_fit")), not_tuned)
  expect_true(file.exists(reg$paths$posterior))
  expect_false(file.exists(reg$paths$prior))
  wf <- suppressMessages(load_wflow(reg$paths$posterior))
  expect_true(wf$tuned)                               # quick_fit marks the wflow tuned
  expect_identical(wf$post_fit$draws, reg$wf$post_fit$draws)

  # explore contrasts later without refitting (update_contrasts -> infer)
  co <- contrast("group", contrast_levels = c("c", "s2"), grouping_vars = "direction",
                 outcome = outcome_spec(gq_var = "mu_rep"))
  reg <- suppressMessages(update_contrasts(reg, list(grp = co)))
  quiet(reg <- infer(reg))
  g <- tm_draws(reg$wf$post_fit, "g")
  p <- get_inference(reg, "grp")$posterior
  expect_lt(max(abs(p$contrast - (g[["g[3]"]] - g[["g[1]"]])[p$.draw])), 1e-4)
  expect_identical(readRDS(reg$paths$inference)$grp$posterior, p)
})
