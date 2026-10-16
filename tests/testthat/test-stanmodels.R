# # ---- helpers ----------------------------------------------------------------

# make_minimal_model_code <- function() {
#   "
# data {
#   int<lower=0> N;
#   int<lower=0> n_group;
#   array[N] int<lower=1, upper=n_group> group_id;
#   vector[N] y;
# }
# parameters { real mu; real<lower=0> sigma; }
# model { y ~ normal(mu, sigma); }
# generated quantities { vector[N] y_rep; for (i in 1:N) y_rep[i] = normal_rng(mu, sigma); }
# "
# }

# make_minimal_df <- function(n = 10) {
#   data.frame(
#     group = rep(c("A", "B"), each = n / 2),
#     value = rnorm(n)
#   )
# }


# test_that("stanmodel creates correct structure", {
#   m <- stanmodel()                          # was: new_stan_model()
#   expect_s3_class(m, "stanmodel")
#   expect_null(m$compiled_model)
#   expect_null(m$model_code)
#   expect_null(m$priors)
#   expect_null(m$fitted_model)
#   expect_s3_class(m$sampler_params, "sampler_params")
# })

# test_that("compile_from_string populates model fields", {
#   skip_if_not(cmdstanr::cmdstan_version() > "0", "CmdStan not available")

#   stan_code <- "
#   data {
#     int<lower=1> N;
#     vector[N] y;
#     real prior_mu;
#     real<lower=0> prior_sigma;
#   }
#   parameters { real mu; real<lower=0> sigma; }
#   model { mu ~ normal(prior_mu, prior_sigma); y ~ normal(mu, sigma); }
#   "

#   m <- stanmodel()                          # was: nstan_model()
#   m <- compile_from_string(m, stan_code)

#   expect_false(is.null(m$compiled_model))
#   expect_equal(m$model_code, stan_code)
#   expect_s3_class(m$priors, "priors")
#   expect_named(m$priors, c("prior_mu", "prior_sigma"))
# })

# test_that("sampler_params creates correct structure", {
#   sp <- sampler_params()                                         # was: new_sampler_params()
#   expect_s3_class(sp, "sampler_params")
#   expect_named(sp, names(.default_sampler_params))
# })

# test_that("sampler_params accepts custom params", {
#   sp <- sampler_params(modifyList(                              # was: new_sampler_params()
#     .default_sampler_params,
#     list(chains = 2L, adapt_delta = 0.95)
#   ))
#   expect_equal(sp$chains, 2L)
#   expect_equal(sp$adapt_delta, 0.95)
# })

# test_that("priors creates correct structure", {
#   p <- priors(c("prior_alpha", "prior_sigma"))           # was: new_priors
#   expect_s3_class(p, "priors")
#   expect_named(p, c("prior_alpha", "prior_sigma"))
#   expect_true(all(sapply(p, is.na)))
# })

# test_that("priors works with empty input", {
#   p <- priors()                                          # was: new_priors
#   expect_s3_class(p, "priors")
#   expect_length(p, 0)
# })


# test_that(".parse_stan_data_block extracts correct variable names", {
#   code <- "
#   data {
#     int<lower=1> N;
#     int<lower=1> n_group;
#     array[N] int<lower=1, upper=n_group> group_id;
#     vector[N] y;
#     real          prior_mu_mean;
#     real<lower=0> prior_mu_sd;
#     real<lower=0> prior_sigma_rate;
#     int<lower=0, upper=1> skip_likelihood;
#   }
#   "
#   parsed <- stanworkflow:::.parse_stan_data_block(code)
#   expect_named(
#     parsed,
#     c("N", "n_group", "group_id", "y",
#       "prior_mu_mean", "prior_mu_sd", "prior_sigma_rate", "skip_likelihood")
#   )
# })

# test_that(".parse_stan_data_block handles matrix declarations", {
#   code <- "
#   data {
#     int<lower=1> N;
#     int<lower=1> P;
#     matrix[N, P] X;
#   }
#   "
#   parsed <- stanworkflow:::.parse_stan_data_block(code)
#   expect_named(parsed, c("N", "P", "X"))
#   expect_equal(parsed$X$type, "matrix")
# })

# test_that(".parse_stan_data_block handles 2D array declarations", {
#   code <- "
#   data {
#     int<lower=1> N;
#     int<lower=1> P;
#     int<lower=1> Q;
#     array[N, P] real X;
#     array[N, P, Q] int Y;
#   }
#   "
#   parsed <- stanworkflow:::.parse_stan_data_block(code)
#   expect_named(parsed, c("N", "P", "Q", "X", "Y"))
#   expect_equal(parsed$X$type, "array")
#   expect_equal(parsed$Y$type, "array")
# })

# test_that(".parse_stan_data_block handles row_vector declarations", {
#   code <- "
#   data {
#     int<lower=1> N;
#     row_vector[N] z;
#   }
#   "
#   parsed <- stanworkflow:::.parse_stan_data_block(code)
#   expect_named(parsed, c("N", "z"))
# })

# test_that(".parse_stan_data_block handles cov_matrix and corr_matrix", {
#   code <- "
#   data {
#     int<lower=1> K;
#     cov_matrix[K]  Sigma;
#     corr_matrix[K] Omega;
#   }
#   "
#   parsed <- stanworkflow:::.parse_stan_data_block(code)
#   expect_named(parsed, c("K", "Sigma", "Omega"))
# })

# test_that(".parse_stan_data_block handles comments on same line", {
#   code <- "
#   data {
#     int<lower=1> N;        // number of observations
#     real<lower=0> sigma;   // scale parameter
#   }
#   "
#   parsed <- stanworkflow:::.parse_stan_data_block(code)
#   expect_named(parsed, c("N", "sigma"))
# })

# test_that(".parse_stan_data_block returns empty list for empty data block", {
#   code <- "
#   data { }
#   parameters { real mu; }
#   "
#   parsed <- stanworkflow:::.parse_stan_data_block(code)
#   expect_length(parsed, 0)
# })

# test_that(".parse_stan_data_block errors when no data block found", {
#   code <- "
#   parameters { real mu; }
#   model { mu ~ normal(0, 1); }
#   "
#   expect_error(
#     stanworkflow:::.parse_stan_data_block(code),
#     "No data block found"
#   )
# })


# test_that("get_stan_data() includes skip_likelihood as NA", {
#   code <- "
#   data {
#     int<lower=0> N;
#     vector[N] y;
#     real prior_mu;
#     int<lower=0, upper=1> skip_likelihood;
#   }
#   parameters { real mu; }
#   model { mu ~ normal(prior_mu, 1); y ~ normal(mu, 1); }
#   "
#   m     <- stanmodel() |> compile_from_string(code)
#   sdata <- get_stan_data(m, df = data.frame(y = rnorm(10)), outcome = list(y = "y"))

#   expect_true("skip_likelihood" %in% names(sdata))
#   expect_true(is.na(sdata$skip_likelihood))
# })

# test_that("fit() runs when skip_likelihood is set to 0L", {
#   skip_if_not(cmdstanr::cmdstan_version() > "0", "CmdStan not available")

#   code <- "
#   data {
#     int<lower=0> N;
#     vector[N] y;
#     real prior_mu;
#     int<lower=0, upper=1> skip_likelihood;
#   }
#   parameters { real mu; }
#   model { mu ~ normal(prior_mu, 1); if (!skip_likelihood) { y ~ normal(mu, 1); } }
#   "
#   m     <- stanmodel() |> compile_from_string(code)
#   sdata <- get_stan_data(m, df = data.frame(y = rnorm(10)), outcome = list(y = "y"))
#   sdata$skip_likelihood <- 0L
#   m$priors$prior_mu     <- 0

#   m <- fit(m, sdata)
#   expect_false(is.null(m$fitted_model))
# })



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
    group = rep(c("A", "B"), each = n / 2),
    value = rnorm(n)
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
  fitted_mock <- list(draws = function(...) make_mock_draws(n))
  list(sample = function(...) fitted_mock)
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

# ---- .parse_stan_data_block -------------------------------------------------

test_that(".parse_stan_data_block extracts correct variable names", {
  code <- "
  data {
    int<lower=1> N;
    int<lower=1> n_group;
    array[N] int<lower=1, upper=n_group> group_id;
    vector[N] y;
    real          prior_mu_mean;
    real<lower=0> prior_mu_sd;
    real<lower=0> prior_sigma_rate;
    int<lower=0, upper=1> skip_likelihood;
  }
  "
  parsed <- stanworkflow:::.parse_stan_data_block(code)
  expect_named(
    parsed,
    c("N", "n_group", "group_id", "y",
      "prior_mu_mean", "prior_mu_sd", "prior_sigma_rate", "skip_likelihood")
  )
})

test_that(".parse_stan_data_block handles matrix declarations", {
  code <- "
  data {
    int<lower=1> N;
    int<lower=1> P;
    matrix[N, P] X;
  }
  "
  parsed <- stanworkflow:::.parse_stan_data_block(code)
  expect_named(parsed, c("N", "P", "X"))
  expect_equal(parsed$X$type, "matrix")
})

test_that(".parse_stan_data_block handles 2D array declarations", {
  code <- "
  data {
    int<lower=1> N;
    int<lower=1> P;
    int<lower=1> Q;
    array[N, P] real X;
    array[N, P, Q] int Y;
  }
  "
  parsed <- stanworkflow:::.parse_stan_data_block(code)
  expect_named(parsed, c("N", "P", "Q", "X", "Y"))
  expect_equal(parsed$X$type, "array")
  expect_equal(parsed$Y$type, "array")
})

test_that(".parse_stan_data_block handles row_vector declarations", {
  code <- "
  data {
    int<lower=1> N;
    row_vector[N] z;
  }
  "
  parsed <- stanworkflow:::.parse_stan_data_block(code)
  expect_named(parsed, c("N", "z"))
})

test_that(".parse_stan_data_block handles cov_matrix and corr_matrix", {
  code <- "
  data {
    int<lower=1> K;
    cov_matrix[K]  Sigma;
    corr_matrix[K] Omega;
  }
  "
  parsed <- stanworkflow:::.parse_stan_data_block(code)
  expect_named(parsed, c("K", "Sigma", "Omega"))
})

test_that(".parse_stan_data_block handles comments on same line", {
  code <- "
  data {
    int<lower=1> N;        // number of observations
    real<lower=0> sigma;   // scale parameter
  }
  "
  parsed <- stanworkflow:::.parse_stan_data_block(code)
  expect_named(parsed, c("N", "sigma"))
})

test_that(".parse_stan_data_block returns empty list for empty data block", {
  code <- "
  data { }
  parameters { real mu; }
  "
  parsed <- stanworkflow:::.parse_stan_data_block(code)
  expect_length(parsed, 0)
})

test_that(".parse_stan_data_block errors when no data block found", {
  code <- "
  parameters { real mu; }
  model { mu ~ normal(0, 1); }
  "
  expect_error(
    stanworkflow:::.parse_stan_data_block(code),
    "No data block found"
  )
})

# ---- get_stan_data ----------------------------------------------------------

test_that("get_stan_data() auto-fills N, n_group, group_id, and outcome", {
  m            <- stanmodel()
  m$model_code <- make_minimal_model_code()
  df           <- make_minimal_df(10)

  sd <- suppressMessages(
    get_stan_data(m, df, outcome = list(y = "value"))
  )

  expect_equal(sd$N, nrow(df))
  expect_equal(sd$n_group, dplyr::n_distinct(df$group))
  expect_equal(sd$group_id, as.integer(factor(df$group)))
  expect_equal(sd$y, df$value)
})

test_that("get_stan_data() messages on unfilled variables", {
  m            <- stanmodel()
  m$model_code <- make_minimal_model_code()
  df           <- make_minimal_df(10)

  expect_message(get_stan_data(m, df), regexp = "y")
})

test_that("get_stan_data() includes skip_likelihood as NA", {
  code <- "
  data {
    int<lower=0> N;
    vector[N] y;
    real prior_mu;
    int<lower=0, upper=1> skip_likelihood;
  }
  parameters { real mu; }
  model { mu ~ normal(prior_mu, 1); y ~ normal(mu, 1); }
  "
  skip_if_not(cmdstanr::cmdstan_version() > "0", "CmdStan not available")
  m     <- stanmodel() |> compile_from_string(code)
  sdata <- get_stan_data(m, df = data.frame(y = rnorm(10)), outcome = list(y = "y"))

  expect_true("skip_likelihood" %in% names(sdata))
  expect_true(is.na(sdata$skip_likelihood))
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
             group_id = as.integer(factor(df$group)), y = df$value)

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
    y = df$value
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
             group_id = as.integer(factor(df$group)), y = df$value)

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
             y = df$value, skip_likelihood = 0L)

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
             y = df$value)

  expect_warning(
    fit(m, data = df, stan_data = sd,
        prior_predictive = TRUE, sampler_params = fast_sampler),
    regexp = "skip_likelihood"
  )
})

test_that("fit() runs with skip_likelihood = 0L (CmdStan integration)", {
  skip_if_not(cmdstanr::cmdstan_version() > "0", "CmdStan not available")

  code <- "
  data {
    int<lower=0> N;
    vector[N] y;
    real prior_mu;
    int<lower=0, upper=1> skip_likelihood;
  }
  parameters { real mu; }
  model { mu ~ normal(prior_mu, 1); if (!skip_likelihood) { y ~ normal(mu, 1); } }
  "
  df    <- data.frame(y = rnorm(10))
  m     <- stanmodel() |> compile_from_string(code)
  sdata <- get_stan_data(m, df = df, outcome = list(y = "y"))
  sdata$skip_likelihood <- 0L
  m$priors$prior_mu     <- 0

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
                  y = df$value, prior_mu = NA_real_)
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
             group_id = as.integer(factor(df$group)), y = df$value)

  result <- fit(m, data = df, stan_data = sd, sampler_params = fast_sampler)

  expect_false(is.null(result$draws))
  expect_true(is.character(result$variables))
  expect_false("lp__" %in% result$variables)
})

# ---- gq (CmdStan integration) ------------------------------------------------

.fit_gq_test_model <- function() {
  code <- "
  data {
    int<lower=0> N;
    vector[N] y;
    real prior_mu;
    int<lower=0, upper=1> skip_likelihood;
  }
  parameters { real mu; real<lower=0> sigma; }
  model {
    mu ~ normal(prior_mu, 1);
    if (!skip_likelihood) y ~ normal(mu, sigma);
  }
  generated quantities {
    vector[N] y_rep;
    for (i in 1:N) y_rep[i] = normal_rng(mu, sigma);
  }
  "
  df    <- data.frame(y = rnorm(10))
  m     <- stanmodel() |> compile_from_string(code)
  sdata <- get_stan_data(m, df = df)
  sdata$skip_likelihood <- 0L
  m$priors$prior_mu     <- 0

  fit(m, data = df, stan_data = sdata, sampler_params = fast_sampler)
}

test_that("gq() with n_draws = NULL returns all available draws", {
  skip_if_not(cmdstanr::cmdstan_version() > "0", "CmdStan not available")

  m <- .fit_gq_test_model()
  n_available <- posterior::ndraws(m$fitted_model$draws())

  out <- gq(m, new_data = m$stan_data)

  expect_equal(posterior::ndraws(out), n_available)
})

test_that("gq() with n_draws = k < available subsamples to exactly k draws", {
  skip_if_not(cmdstanr::cmdstan_version() > "0", "CmdStan not available")

  m <- .fit_gq_test_model()

  out <- gq(m, new_data = m$stan_data, n_draws = 20L)

  expect_equal(posterior::ndraws(out), 20L)
})

test_that("gq() cleans up its own temp output_dir after the call", {
  skip_if_not(cmdstanr::cmdstan_version() > "0", "CmdStan not available")

  m <- .fit_gq_test_model()

  before <- list.files(tempdir(), pattern = "^gq_")
  gq(m, new_data = m$stan_data, n_draws = 20L)
  after  <- list.files(tempdir(), pattern = "^gq_")

  expect_equal(after, before)
})

test_that("gq() with cache_compiled = TRUE returns correct draws across repeated calls", {
  skip_if_not(cmdstanr::cmdstan_version() > "0", "CmdStan not available")

  m <- .fit_gq_test_model()

  out1 <- gq(m, new_data = m$stan_data, n_draws = 10L, cache_compiled = TRUE)
  expect_equal(posterior::ndraws(out1), 10L)

  # second call should reuse the memoized compiled model and still work
  out2 <- gq(m, new_data = m$stan_data, n_draws = 10L, cache_compiled = TRUE)
  expect_equal(posterior::ndraws(out2), 10L)
})

test_that("gq() with cache_compiled = FALSE reproduces original per-call recompile behavior", {
  skip_if_not(cmdstanr::cmdstan_version() > "0", "CmdStan not available")

  m <- .fit_gq_test_model()

  out <- gq(m, new_data = m$stan_data, cache_compiled = FALSE)

  expect_equal(posterior::ndraws(out), posterior::ndraws(m$fitted_model$draws()))
})