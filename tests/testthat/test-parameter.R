# Parameter route (draws_method = "parameter") and the parameter-draw helpers
# parameter_draws() / parameter_scenarios() and their aliases, on the shared linear
# model (helper-model.R): y = a[direction] + g[group] + b[direction] * x + e.

dirs   <- c("FLEX", "EXT")
groups <- c("c", "s1", "s2")

hpdi <- function(v, p = 0.95) as.numeric(coda::HPDinterval(coda::as.mcmc(v), prob = p))

par_run <- function(m, co, ...) quiet(infer_contrasts(m, co, draws_method = "parameter", ...))

# ---- no Stan ---------------------------------------------------------------------

test_that("contrast() validates the parameter-route fields; delta needs a numeric predictor", {
  expect_error(contrast("x", parameter = c("a", "b")))
  expect_error(contrast("x", delta = "sd"))
  expect_error(contrast("x", delta = c(1, 2)))
  co <- contrast("x", continuous = TRUE, parameter = "b", factors = "direction", delta = "iqr")
  expect_equal(co$link, "identity")
  expect_output(print(co), "parameter      : b")
  d <- data.frame(x = c(0, 1, 2, 3, 4), g = letters[1:5])
  expect_equal(stanworkflow:::.delta_value(co, d), 2)                          # q75 - q25
  expect_equal(stanworkflow:::.delta_value(modifyList(co, list(delta = "q10q90")), d), 3.2)
  expect_equal(stanworkflow:::.delta_value(modifyList(co, list(delta = "2sd")), d), 2 * sd(d$x))
  expect_equal(stanworkflow:::.delta_value(modifyList(co, list(delta = "range")), d), 4)
  expect_equal(stanworkflow:::.delta_value(modifyList(co, list(delta = 0.5)), d), 0.5)
  expect_null(stanworkflow:::.delta_value(modifyList(co, list(delta = NULL)), d))
  expect_error(stanworkflow:::.delta_value(contrast("g", delta = "iqr"), d), "continuous = TRUE")
  expect_match(stanworkflow:::.scale_label(co), "IQR of x")
  expect_equal(stanworkflow:::.scale_label(contrast("x")), "parameter")
})

test_that("old names ppc_by_parameter()/ppc_scenarios() are deprecated aliases; grouped_ppc() = predictive_draws()", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  rlang::local_options(lifecycle_verbosity = "warning")
  expect_warning(old <- ppc_by_parameter(m, "b", list(direction = dirs)), "parameter_draws")
  expect_identical(old, parameter_draws(m, "b", list(direction = dirs)))
  expect_warning(olds <- ppc_scenarios(m, "b", list(direction = dirs), scenario_factor = "direction"),
                 "parameter_scenarios")
  expect_identical(olds[1:3], parameter_scenarios(m, "b", list(direction = dirs),
                                                  scenario_factor = "direction")[1:3])
  expect_identical(grouped_ppc(m, data = tm_data(), feature_names = "group"),
                   predictive_draws(m, data = tm_data(), feature_names = "group"))
})

# ---- CmdStan ---------------------------------------------------------------------

test_that("parameter_draws() labels indices", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  b <- tm_draws(m, "b")
  long <- parameter_draws(m, "b", list(direction = dirs))
  expect_equal(nrow(long), 2L * posterior::ndraws(m$draws))
  expect_equal(levels(long$direction), dirs)
  ext <- dplyr::arrange(dplyr::filter(long, direction == "EXT"), .draw)
  expect_equal(ext$value, b[["b[2]"]])
  wide <- parameter_draws(m, "g", list(group = groups), .widen_factor = "group")
  expect_equal(names(wide), c(".chain", ".iteration", ".draw", groups))
  expect_error(parameter_draws(m, "b", list(direction = c("FLEX", "EXT", "X"))), "does not match")
  expect_error(parameter_draws(m, "nope", list()), "not found")
  expect_equal(names(parameter_draws(m, "sigma", list())), c("sigma", ".chain", ".iteration", ".draw"))
})

test_that("parameter_scenarios() classifies draws by sign", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  b <- tm_draws(m, "b")
  s <- parameter_scenarios(m, "b", list(direction = dirs), scenario_factor = "direction")
  sd <- dplyr::arrange(s$scenario_draws, .draw)
  expect_identical(sd$scenario == "Affects both", b[["b[1]"]] > 0 & b[["b[2]"]] > 0)
  expect_equal(sum(unlist(s$scenario_grid)), 1)
  expect_s3_class(s$plot, "ggplot")
  # scenario_levels picks (and orders) the two outcome columns
  s2 <- parameter_scenarios(m, "b", list(direction = dirs), scenario_factor = "direction",
                      scenario_levels = c("EXT", "FLEX"), effect_positive = FALSE)
  expect_identical(s2$scenario_draws$scenario == "Affects EXT only",
                   s2$scenario_draws$EXT < 0 & s2$scenario_draws$FLEX >= 0)
})

test_that("parameter route: per-level summaries match a hand computation", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  b <- tm_draws(m, "b")
  r <- par_run(m, contrast("x", continuous = TRUE, parameter = "b", factors = "direction"))
  expect_s3_class(r, "inference")
  expect_equal(r$scale, "parameter")
  expect_null(r$delta)
  ls <- r$level_summary
  expect_equal(as.character(ls$direction), dirs)
  v <- b[["b[2]"]]
  expect_equal(ls$median[2], median(v))
  expect_equal(c(ls$hpdi_lo[2], ls$hpdi_hi[2]), hpdi(v))
  expect_equal(ls$pd[2], max(mean(v > 0), mean(v < 0)))
  expect_equal(ls$P_positive[2], mean(v > 0))
  # continuous: the contrast is the per-level parameter itself
  expect_equal(r$summary$median, ls$median)
  expect_equal(r$contrast$grouping_vars, "direction")
  # numeric delta multiplies; prob is used for the HPDI
  r2 <- par_run(m, contrast("x", continuous = TRUE, parameter = "b", factors = list(direction = dirs),
                            delta = 0.5), prob = 0.8)
  expect_equal(r2$level_summary$median[2], median(0.5 * v))
  expect_equal(c(r2$level_summary$hpdi_lo[2], r2$level_summary$hpdi_hi[2]), hpdi(0.5 * v, 0.8))
  expect_equal(r2$scale, sprintf("parameter %s 0.5", intToUtf8(0xd7)))
})

test_that("parameter route: differences between levels, draw by draw (one pair, all pairs, filter)", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  b <- tm_draws(m, "b")
  g <- tm_draws(m, "g")
  r <- par_run(m, contrast("direction", contrast_levels = c("EXT", "FLEX"),
                           parameter = "b", factors = "direction"))
  p <- dplyr::arrange(r$posterior, .draw)
  expect_equal(p$contrast, b[["b[1]"]] - b[["b[2]"]])                # FLEX - EXT
  expect_equal(unique(p$contrast_label), "FLEX - EXT")
  expect_equal(r$summary$median, median(b[["b[1]"]] - b[["b[2]"]]))
  expect_equal(nrow(r$level_summary), 2L)

  # all pairs when contrast_levels is NULL
  ra <- par_run(m, contrast("group", parameter = "g", factors = list(group = groups)))
  expect_equal(ra$summary$contrast_label, c("s1 - c", "s2 - c", "s2 - s1"))
  s2s1 <- dplyr::arrange(dplyr::filter(ra$posterior, contrast_label == "s2 - s1"), .draw)
  expect_equal(s2s1$contrast, g[["g[3]"]] - g[["g[2]"]])
  expect_equal(ra$contrast$grouping_vars, "contrast_label")

  # filter_draws on a factor; a contrast_var that is not a factor errors clearly
  rf <- par_run(m, contrast("x", continuous = TRUE, parameter = "b", factors = "direction",
                            filter_draws = rlang::quo(direction == "EXT")))
  expect_equal(as.character(unique(rf$posterior$direction)), "EXT")
  expect_error(par_run(m, contrast("x", parameter = "b", factors = "direction")), "not one of")
  expect_error(par_run(m, contrast("x", continuous = TRUE)), "needs `parameter`")
})

test_that("identity: parameter x range(x) equals the counterfactual mu_rep contrast draw by draw", {
  skip_if_no_cmdstan()
  m  <- tm_fit()
  ep <- FALSE   # sign convention must carry over as well
  cf <- quiet(infer_contrasts(m, contrast("x", continuous = TRUE, two_way = "direction",
                                          two_way_levels = dirs, effect_positive = ep,
                                          outcome = outcome_spec(gq_var = "mu_rep")),
                              draws_method = "counterfactual", gq_steps = 2))
  pr <- par_run(m, contrast("x", continuous = TRUE, two_way = "direction", two_way_levels = dirs,
                            effect_positive = ep, parameter = "b", factors = "direction",
                            delta = "range"))
  expect_equal(pr$delta, diff(range(m$data$x)))
  expect_match(pr$scale, "range of x")
  a <- dplyr::arrange(cf$scenario_draws, .draw)
  z <- dplyr::arrange(pr$scenario_draws, .draw)
  expect_equal(z$.draw, a$.draw)
  expect_lt(max(abs(z$contrast_FLEX - a$contrast_FLEX), abs(z$contrast_EXT - a$contrast_EXT)), 1e-4)
  expect_identical(z$scenario, a$scenario)
  expect_equal(pr$scenario_grid$probability,
               dplyr::arrange(cf$scenario_grid, scenario)$probability[
                 match(pr$scenario_grid$scenario, sort(cf$scenario_grid$scenario))])
  expect_output(print(pr), "Scale: parameter")

  # one-way: per-direction counterfactual grid posterior vs parameter posterior
  cf1 <- quiet(infer_contrasts(m, contrast("x", continuous = TRUE, grouping_vars = "direction",
                                           outcome = outcome_spec(gq_var = "mu_rep")),
                               draws_method = "counterfactual", gq_steps = 2))
  pr1 <- par_run(m, contrast("x", continuous = TRUE, parameter = "b", factors = "direction",
                             delta = "range"))
  j <- dplyr::inner_join(cf1$grid_posterior, pr1$posterior, by = c(".draw", "direction"))
  expect_equal(nrow(j), 2L * posterior::ndraws(m$draws))
  expect_lt(max(abs(j$contrast.x - j$contrast.y)), 1e-4)
})

test_that("averaging over factors: equal and data weights match hand computations", {
  skip_if_no_cmdstan()
  m  <- tm_fit()
  bg <- tm_draws(m, "bg")
  cell <- function(d, k) bg[[sprintf("bg[%d,%d]", d, k)]]
  # unbalanced data, so data weights differ from equal weights
  d  <- tm_data()
  d  <- d[-which(d$direction == "FLEX" & d$group == "s2")[1:6], ]
  d  <- d[-which(d$direction == "EXT" & d$group == "c")[1:3], ]
  n  <- table(d$direction, d$group)      # FLEX/EXT x c/s1/s2
  bf <- c("direction", "group")
  ap <- function(co) dplyr::arrange(co$posterior, .draw, direction)

  # NULL grouping_vars: every factor a cell, the same as naming all of them
  r0 <- par_run(m, contrast("x", continuous = TRUE, parameter = "bg", factors = bf), data = d)
  expect_identical(r0, par_run(m, contrast("x", continuous = TRUE, grouping_vars = bf,
                                           parameter = "bg", factors = bf), data = d))
  expect_equal(nrow(r0$level_summary), 6L)

  re <- par_run(m, contrast("x", continuous = TRUE, grouping_vars = "direction",
                            parameter = "bg", factors = bf), data = d)
  expect_equal(nrow(re$level_summary), 2L)
  expect_match(re$scale, "averaged over group (equal weights)", fixed = TRUE)
  expect_equal(re$contrast$grouping_vars, "direction")
  pe <- ap(re)
  expect_equal(pe$contrast[pe$direction == "EXT"], (cell(2, 1) + cell(2, 2) + cell(2, 3)) / 3)

  rd <- par_run(m, contrast("x", continuous = TRUE, grouping_vars = "direction",
                            parameter = "bg", factors = bf, weights = "data"), data = d)
  expect_match(rd$scale, "(data weights)", fixed = TRUE)
  pd <- ap(rd)
  w  <- n["FLEX", ] / sum(n["FLEX", ])
  expect_equal(pd$contrast[pd$direction == "FLEX"],
               w[["c"]] * cell(1, 1) + w[["s1"]] * cell(1, 2) + w[["s2"]] * cell(1, 3))
  expect_equal(rd$level_summary$median[1], median(pd$contrast[pd$direction == "FLEX"]))

  # categorical contrast_var stays a cell: FLEX - EXT of the group averages
  rc <- par_run(m, contrast("direction", contrast_levels = c("EXT", "FLEX"),
                            grouping_vars = "direction", parameter = "bg", factors = bf), data = d)
  pc <- dplyr::arrange(rc$posterior, .draw)
  expect_equal(pc$contrast, (cell(1, 1) + cell(1, 2) + cell(1, 3) -
                               cell(2, 1) - cell(2, 2) - cell(2, 3)) / 3)

  # two_way kept as a cell, everything else averaged
  rt <- par_run(m, contrast("x", continuous = TRUE, two_way = "direction", grouping_vars = "direction",
                            parameter = "bg", factors = bf), data = d)
  expect_equal(sum(rt$scenario_grid$probability), 1)

  # errors: unknown grouping_vars; data weights need data columns
  expect_error(par_run(m, contrast("x", continuous = TRUE, grouping_vars = "finger",
                                   parameter = "bg", factors = bf)), "not in names")
  expect_error(par_run(m, contrast("x", continuous = TRUE, grouping_vars = "direction", weights = "data",
                                   parameter = "bg", factors = list(direction = dirs, grp = groups))),
               "not a column")
})

test_that("identity: data-weighted average of bg x range(x) equals the counterfactual mu_bg_rep contrast", {
  skip_if_no_cmdstan()
  # counterfactual GQ runs on the fitted data, so refit (no recompile) on unbalanced
  # data to make data weights differ from equal weights
  d <- tm_data()
  d <- d[-which(d$direction == "FLEX" & d$group == "s2")[1:6], ]
  m <- tm_compiled() |> fit(data = d)
  cf <- quiet(infer_contrasts(m, contrast("x", continuous = TRUE, grouping_vars = "direction",
                                          outcome = outcome_spec(gq_var = "mu_bg_rep")),
                              data = d, draws_method = "counterfactual", gq_steps = 2))
  pr <- par_run(m, contrast("x", continuous = TRUE, grouping_vars = "direction", weights = "data",
                            parameter = "bg", factors = c("direction", "group"), delta = "range"),
                data = d)
  j <- dplyr::inner_join(cf$grid_posterior, pr$posterior, by = c(".draw", "direction"))
  expect_equal(nrow(j), 2L * posterior::ndraws(m$draws))
  expect_lt(max(abs(j$contrast.x - j$contrast.y)), 1e-4)
  # equal weights do not reproduce it on unbalanced data
  pe <- par_run(m, contrast("x", continuous = TRUE, grouping_vars = "direction",
                            parameter = "bg", factors = c("direction", "group"), delta = "range"),
                data = d)
  je <- dplyr::inner_join(cf$grid_posterior, pe$posterior, by = c(".draw", "direction"))
  expect_gt(max(abs(je$contrast.x - je$contrast.y)), 1e-3)
})

test_that("non-identity link gives a one-time message pointing to the counterfactual route", {
  skip_if_no_cmdstan()
  m  <- tm_fit()
  co <- contrast("x", continuous = TRUE, parameter = "b", factors = "direction", link = "logit")
  rlang::local_options(rlib_message_verbosity = "default")
  rlang::reset_message_verbosity("stanworkflow_parameter_link")
  expect_message(infer_contrasts(m, co, draws_method = "parameter"), "counterfactual")
  expect_no_message(infer_contrasts(m, co, draws_method = "parameter"))
  expect_no_message(infer_contrasts(m, contrast("x", continuous = TRUE, parameter = "b",
                                                factors = "direction"), draws_method = "parameter"))
})

test_that("regimen(draws_method = 'parameter') stores, saves and prints parameter results", {
  skip_if_no_cmdstan()
  tm_compiled()
  reg <- regimen(model_name = "tm", data_list = list(data = tm_data()),
                 contrasts = list(
                   slope = contrast("x", continuous = TRUE, two_way = "direction",
                                    parameter = "b", factors = "direction", delta = "2sd"),
                   dir   = contrast("direction", contrast_levels = c("EXT", "FLEX"),
                                    parameter = "b", factors = "direction")),
                 models_folder = tm_models, results_path = tempfile("reg_"),
                 prior_sampler_params = tm_sampler, post_sampler_params = tm_sampler,
                 draws_method = "parameter")
  expect_output(print(reg), "draws_method  : parameter")
  quiet(suppressWarnings(reg <- run(reg, steps = c("quick_fit", "infer"))))
  expect_true(reg$status$inferred)
  sl <- get_inference(reg, "slope")
  expect_s3_class(sl, "inference")
  expect_equal(sl$delta, 2 * sd(reg$data_list$data$x))
  expect_match(sl$scale, "2 SD of x")
  expect_false(is.null(sl$scenario_grid))
  expect_equal(unique(get_inference(reg, "dir")$posterior$contrast_label), "FLEX - EXT")
  expect_identical(readRDS(reg$paths$inference)$dir$summary, get_inference(reg, "dir")$summary)
  expect_output(print(sl), "Scale: parameter .*2 SD of x")
  expect_output(print(get_inference(reg, "dir")), "Method:  parameter")
})
