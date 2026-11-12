# ---- helpers ----------------------------------------------------------------

make_minimal_model_code <- function() {
  "
data {
  int<lower=0> N;
  int<lower=0> n_group;
  array[N] int<lower=1, upper=n_group> group_id;
  vector[N] y;
}
parameters { real mu; real<lower=0> sigma; }
model { y ~ normal(mu, sigma); }
generated quantities { vector[N] y_rep; for (i in 1:N) y_rep[i] = normal_rng(mu, sigma); }
"
}

make_minimal_df <- function(n = 10) {
  data.frame(
    group = factor(rep(c("A", "B"), each = n / 2)),
    y     = rnorm(n)
  )
}

make_mock_draws <- function(n = 10) {
  arr <- array(
    rnorm(n * 4 * 3),
    dim      = c(n, 4, 3),
    dimnames = list(NULL, NULL, c("mu", "sigma", "y_rep[1]"))
  )
  class(arr) <- c("draws_array", "draws", "array")
  arr
}

make_compiled_mock <- function(n = 10) {
  fitted_mock <- list(draws = function(...) make_mock_draws(n),
                      return_codes = function() 0L,
                      num_chains = function() 4L,
                      output_files = function(...) character(4))
  list(sample    = function(...) fitted_mock,
       variables = function() list(data = list(N = NULL, n_group = NULL,
                                               group_id = NULL, y = NULL)))
}

# a real CmdStanModel without compiling: enough for $variables()
uncompiled_model <- function(code) {
  m                <- stanmodel()
  m$model_code     <- code
  m$compiled_model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)
  m
}

make_compiled_stanmodel <- function(n = 10) {
  m                <- stanmodel()
  m$model_code     <- make_minimal_model_code()
  m$priors         <- priors(character())
  m$compiled_model <- make_compiled_mock(n)
  m
}

fast_sampler <- sampler_params(list(
  chains          = 1L,
  parallel_chains = 1L,
  iter_warmup     = 100L,
  iter_sampling   = 100L,
  adapt_delta     = 0.8,
  max_treedepth   = 10L,
  show_messages   = FALSE,
  show_exceptions = FALSE,
  refresh         = 0L
))
fast2 <- modifyList(fast_sampler, list(chains = 2L, parallel_chains = 2L))

# ---- stanmodel constructor --------------------------------------------------

test_that("stanmodel creates correct structure", {
  m <- stanmodel()
  expect_s3_class(m, "stanmodel")
  expect_null(m$compiled_model)
  expect_null(m$model_code)
  expect_null(m$priors)
  expect_null(m$fitted_model)
  expect_s3_class(m$sampler_params, "sampler_params")
})

# ---- compile_from_string ----------------------------------------------------

test_that("compile_from_string populates model fields", {
  skip_if_not(cmdstanr::cmdstan_version() > "0", "CmdStan not available")

  stan_code <- "
  data {
    int<lower=1> N;
    vector[N] y;
    real prior_mu;
    real<lower=0> prior_sigma;
  }
  parameters { real mu; real<lower=0> sigma; }
  model { mu ~ normal(prior_mu, prior_sigma); y ~ normal(mu, sigma); }
  "

  m <- stanmodel() |> compile_from_string(stan_code)

  expect_false(is.null(m$compiled_model))
  expect_equal(m$model_code, stan_code)
  expect_s3_class(m$priors, "priors")
  expect_named(m$priors, c("prior_mu", "prior_sigma"))
})

# ---- sampler_params ---------------------------------------------------------

test_that("sampler_params creates correct structure", {
  sp <- sampler_params()
  expect_s3_class(sp, "sampler_params")
  expect_named(sp, names(.default_sampler_params))
})

test_that("sampler_params accepts custom params", {
  sp <- sampler_params(modifyList(
    .default_sampler_params,
    list(chains = 2L, adapt_delta = 0.95)
  ))
  expect_equal(sp$chains, 2L)
  expect_equal(sp$adapt_delta, 0.95)
})

# ---- priors -----------------------------------------------------------------

test_that("priors creates correct structure", {
  p <- priors(c("prior_alpha", "prior_sigma"))
  expect_s3_class(p, "priors")
  expect_named(p, c("prior_alpha", "prior_sigma"))
  expect_true(all(sapply(p, is.na)))
})

test_that("priors works with empty input", {
  p <- priors()
  expect_s3_class(p, "priors")
  expect_length(p, 0)
})

# ---- data-block names (stanc via $variables()) -------------------------------

test_that("data names come from stanc: comments with braces, multi-declarations, multi-line", {
  skip_if_no_cmdstan()
  code <- "
functions { real f(real x) { return x; } }
data {
  int<lower=1> N;            // rows {braces in a comment}
  int<lower=1> n_movie, n_rater;
  array[N] int<lower=1, upper=n_movie>
     movie_id;
  vector[N] y;
  real<lower=0> prior_scale;  // prior_fake in a comment
  int<lower=0, upper=1> skip_likelihood;
}
transformed data { real prior_not_data = 1; }
parameters { real mu; }
model { mu ~ normal(0, prior_scale); }
"
  m <- uncompiled_model(code)
  expect_equal(stanworkflow:::.data_names(m),
               c("N", "n_movie", "n_rater", "movie_id", "y", "prior_scale", "skip_likelihood"))
  m <- compile_from_file(stanmodel(), cmdstanr::write_stan_file(code), compile = FALSE)
  expect_named(m$priors, "prior_scale")
})

test_that("data names: empty data block, and a clear error without a model object", {
  skip_if_no_cmdstan()
  expect_length(stanworkflow:::.data_names(uncompiled_model("parameters { real mu; } model { mu ~ normal(0, 1); }")), 0)
  m <- stanmodel(); m$model_code <- make_minimal_model_code()
  expect_error(get_stan_data(m, make_minimal_df()), "compile_from_string")
})

# ---- get_stan_data ----------------------------------------------------------

test_that("get_stan_data() auto-fills N, n_group, group_id, and outcome", {
  skip_if_no_cmdstan()
  m            <- uncompiled_model(make_minimal_model_code())
  df           <- make_minimal_df(10)

  sd <- suppressMessages(
    get_stan_data(m, df)
  )

  expect_equal(sd$N, nrow(df))
  expect_equal(sd$n_group, dplyr::n_distinct(df$group))
  expect_equal(sd$group_id, as.integer(factor(df$group)))
  expect_equal(sd$y, df$y)
})

test_that("get_stan_data() messages on unfilled variables", {
  skip_if_no_cmdstan()
  m            <- uncompiled_model(make_minimal_model_code())
  df           <- make_minimal_df(10)["group"]

  expect_message(get_stan_data(m, df), regexp = "y")
})

test_that("get_stan_data() includes skip_likelihood as NA", {
  skip_if_no_cmdstan()
  sdata <- suppressMessages(get_stan_data(tm_compiled(), df = tm_data()))

  expect_true("skip_likelihood" %in% names(sdata))
  expect_true(is.na(sdata$skip_likelihood))
})

test_that("get_stan_data() codes a non-factor column locale-independently and says so", {
  skip_if_no_cmdstan()
  m  <- uncompiled_model(make_minimal_model_code())
  df <- data.frame(group = c("b", "B", "a", "b"), y = 1:4)
  expect_warning(sd <- suppressMessages(get_stan_data(m, df)), "level order: B, a, b", class = "stanworkflow_id_warning")
  expect_equal(sd$group_id, c(3L, 1L, 2L, 3L))
  expect_equal(sd$n_group, 3L)
})

test_that("get_stan_data() warns about empty factor levels and keeps factor order", {
  skip_if_no_cmdstan()
  m  <- uncompiled_model(make_minimal_model_code())
  df <- data.frame(group = factor(c("z", "x", "z"), levels = c("z", "y", "x")), y = 1:3)
  expect_warning(sd <- suppressMessages(get_stan_data(m, df)), "empty factor level.*: y")
  expect_equal(sd$group_id, c(1L, 2L, 1L))
  expect_equal(sd$n_group, 2L)
  df$group <- droplevels(df$group)
  expect_no_warning(suppressMessages(get_stan_data(m, df)))
})

test_that("get_stan_data() reports unfilled variables once, as a message only", {
  skip_if_no_cmdstan()
  m <- uncompiled_model(make_minimal_model_code())
  df <- data.frame(group = factor(c("A", "B")))
  expect_output(expect_message(get_stan_data(m, df), "y"), NA)
})

# ---- fit: guards ------------------------------------------------------------

test_that("fit() errors when model is not compiled", {
  m <- stanmodel()
  expect_error(fit(m, data = make_minimal_df()), regexp = "compiled")
})

test_that("fit() errors when stan_data contains NA fields", {
  n  <- 10
  m  <- make_compiled_stanmodel(n)
  df <- make_minimal_df(n)

  sd <- list(
    N = n, n_group = 2L,
    group_id = as.integer(factor(df$group)),
    y = NA
  )

  expect_error(
    fit(m, data = df, stan_data = sd, sampler_params = fast_sampler),
    regexp = "NA"
  )
})

test_that("fit() errors when priors argument is not a priors object", {
  n  <- 10
  m  <- make_compiled_stanmodel(n)
  df <- make_minimal_df(n)

  sd <- list(N = n, n_group = 2L,
             group_id = as.integer(factor(df$group)), y = df$y)

  expect_error(
    fit(m, data = df, stan_data = sd,
        priors = list(prior_mu = 0), sampler_params = fast_sampler)
  )
})

# ---- fit: data paths --------------------------------------------------------

test_that("fit() accepts a pre-built stan_data list", {
  n  <- 10
  m  <- make_compiled_stanmodel(n)
  df <- make_minimal_df(n)

  sd <- list(
    N = n, n_group = 2L,
    group_id = as.integer(factor(df$group)),
    y = df$y
  )

  result <- fit(m, data = df, stan_data = sd, sampler_params = fast_sampler)

  expect_s3_class(result, "stanmodel")
  expect_false(is.null(result$fitted_model))
  expect_equal(result$stan_data, sd)
})

test_that("fit() generates stan_data from a dataframe when stan_data = NULL", {
  n  <- 10
  m  <- make_compiled_stanmodel(n)
  df <- make_minimal_df(n)

  result <- suppressMessages(
    fit(m, data = df, outcome_var = "value",
        stan_data = NULL, sampler_params = fast_sampler)
  )

  expect_s3_class(result, "stanmodel")
  expect_equal(result$stan_data$N, n)
  expect_equal(result$outcome_var, "value")
  expect_equal(result$data, df)
})

test_that("fit() stores outcome_var and data on the model", {
  n  <- 10
  m  <- make_compiled_stanmodel(n)
  df <- make_minimal_df(n)

  sd <- list(N = n, n_group = 2L,
             group_id = as.integer(factor(df$group)), y = df$y)

  result <- fit(m, data = df, outcome_var = "value",
                stan_data = sd, sampler_params = fast_sampler)

  expect_equal(result$outcome_var, "value")
  expect_equal(result$data, df)
})

# ---- fit: prior_predictive --------------------------------------------------

test_that("fit() sets skip_likelihood = 1 when prior_predictive = TRUE", {
  n  <- 10
  m  <- make_compiled_stanmodel(n)
  df <- make_minimal_df(n)

  sd <- list(N = n, n_group = 2L,
             group_id = as.integer(factor(df$group)),
             y = df$y, skip_likelihood = 0L)

  result <- fit(m, data = df, stan_data = sd,
                prior_predictive = TRUE, sampler_params = fast_sampler)

  expect_equal(result$stan_data$skip_likelihood, 1L)
})

test_that("fit() warns when prior_predictive = TRUE but skip_likelihood absent", {
  n  <- 10
  m  <- make_compiled_stanmodel(n)
  df <- make_minimal_df(n)

  sd <- list(N = n, n_group = 2L,
             group_id = as.integer(factor(df$group)),
             y = df$y)

  expect_warning(
    fit(m, data = df, stan_data = sd,
        prior_predictive = TRUE, sampler_params = fast_sampler),
    regexp = "skip_likelihood"
  )
})

test_that("fit() runs with skip_likelihood = 0L (CmdStan integration)", {
  skip_if_no_cmdstan()
  m     <- tm_compiled()
  df    <- tm_data()
  sdata <- suppressMessages(get_stan_data(m, df = df))
  sdata$skip_likelihood <- 0L

  result <- fit(m, data = df, stan_data = sdata, sampler_params = fast_sampler)
  expect_false(is.null(result$fitted_model))
})

# ---- fit: priors merge ------------------------------------------------------

test_that("fit() merges priors object into stan_data", {
  n  <- 10
  m  <- make_compiled_stanmodel(n)
  df <- make_minimal_df(n)

  sd      <- list(N = n, n_group = 2L,
                  group_id = as.integer(factor(df$group)),
                  y = df$y, prior_mu = NA_real_)
  p       <- priors("prior_mu")
  p$prior_mu <- 0

  result <- fit(m, data = df, stan_data = sd,
                priors = p, sampler_params = fast_sampler)

  expect_equal(result$stan_data$prior_mu, 0)
})

# ---- fit: return structure --------------------------------------------------

test_that("fit() populates draws and variables, excludes lp__", {
  n  <- 10
  m  <- make_compiled_stanmodel(n)
  df <- make_minimal_df(n)

  sd <- list(N = n, n_group = 2L,
             group_id = as.integer(factor(df$group)), y = df$y)

  result <- fit(m, data = df, stan_data = sd, sampler_params = fast_sampler)

  expect_false(is.null(result$draws))
  expect_true(is.character(result$variables))
  expect_false("lp__" %in% result$variables)
})

# ---- gq (CmdStan integration) ------------------------------------------------

.fit_gq_test_model <- function() tm_fit()

test_that("gq() with n_draws = NULL returns all available draws", {
  skip_if_no_cmdstan()

  m <- .fit_gq_test_model()
  n_available <- posterior::ndraws(m$fitted_model$draws())

  out <- gq(m, new_data = m$stan_data)

  expect_equal(posterior::ndraws(out), n_available)
})

test_that("gq() with n_draws = k < available subsamples to exactly k draws", {
  skip_if_no_cmdstan()

  m <- .fit_gq_test_model()

  out <- gq(m, new_data = m$stan_data, n_draws = 20L)

  expect_equal(posterior::ndraws(out), 20L)
})

test_that("gq() cleans up its own temp output_dir after the call", {
  skip_if_no_cmdstan()

  m <- .fit_gq_test_model()

  before <- list.files(tempdir(), pattern = "^gq_")
  gq(m, new_data = m$stan_data, n_draws = 20L)
  after  <- list.files(tempdir(), pattern = "^gq_")

  expect_equal(after, before)
})

test_that("gq() with cache_compiled = TRUE returns correct draws across repeated calls", {
  skip_if_no_cmdstan()

  m <- .fit_gq_test_model()

  out1 <- gq(m, new_data = m$stan_data, n_draws = 10L, cache_compiled = TRUE)
  expect_equal(posterior::ndraws(out1), 10L)

  # second call should reuse the memoized compiled model and still work
  out2 <- gq(m, new_data = m$stan_data, n_draws = 10L, cache_compiled = TRUE)
  expect_equal(posterior::ndraws(out2), 10L)
})

test_that("gq() with cache_compiled = FALSE reproduces original per-call recompile behavior", {
  skip_if_no_cmdstan()

  m <- .fit_gq_test_model()

  out <- gq(m, new_data = m$stan_data, cache_compiled = FALSE)

  expect_equal(posterior::ndraws(out), posterior::ndraws(m$fitted_model$draws()))
})

test_that("fit() keeps all chains when set.seed() repeats within a minute", {
  skip_if_no_cmdstan()
  m <- tm_compiled(); df <- tm_data()
  sd <- suppressMessages(get_stan_data(m, df = df)); sd$skip_likelihood <- 0L
  for (i in 1:3) { set.seed(1); f <- fit(m, data = df, stan_data = sd, sampler_params = fast2) }
  expect_equal(posterior::nchains(f$draws), 2L)
})

test_that(".check_chains() fails loudly on a failed chain", {
  expect_error(stanworkflow:::.check_chains(list(return_codes = function() c(0L, 1L))), "chain(s) 2", fixed = TRUE)
  expect_silent(stanworkflow:::.check_chains(list(return_codes = function() c(0L, 0L))))
})

test_that(".complete_fit() rebuilds a fit whose chain cmdstanr wrongly marked failed", {
  skip_if_no_cmdstan()
  f <- tm_fit()$fitted_model
  all_files <- f$output_files()
  # simulate the race: chain 2 exited 0 with a complete CSV, but output_files() skips it
  lossy <- list(return_codes = function() c(0L, 0L), num_chains = function() 2L,
                output_files = function(include_failed = FALSE)
                  if (include_failed) all_files else all_files[1])
  g <- stanworkflow:::.complete_fit(lossy)
  expect_equal(posterior::nchains(g$draws()), 2L)
  expect_s3_class(g, "CmdStanMCMC")                          # still usable by gq()/loo
  expect_identical(stanworkflow:::.complete_fit(f), f)       # complete fit passes through
  lossy_gone <- lossy
  lossy_gone$output_files <- function(include_failed = FALSE)
    if (include_failed) c(all_files[1], "nope.csv") else all_files[1]
  expect_error(stanworkflow:::.complete_fit(lossy_gone), "missing")
})
