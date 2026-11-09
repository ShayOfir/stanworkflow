# Counterfactual contrasts on a tiny linear model: y = a[dir] + b[dir] * x.
# mu_rep is linear in x, so the contrast between the two grid points is
# b[dir] * (max(x) - min(x)) in every draw; mu0_rep does not depend on x, so
# its contrast is exactly 0 whenever baseline and grid point use the same draw.

cf_code <- "
data {
  int<lower=1> N;
  int<lower=1> n_direction;
  array[N] int<lower=1, upper=n_direction> direction_id;
  vector[N] x;
  vector[N] y;
  int<lower=0, upper=1> skip_likelihood;
}
parameters { vector[n_direction] a; vector[n_direction] b; real<lower=0> sigma; }
model {
  a ~ normal(0, 1); b ~ normal(0, 1); sigma ~ exponential(1);
  if (!skip_likelihood) y ~ normal(a[direction_id] + b[direction_id] .* x, sigma);
}
generated quantities {
  vector[N] mu_rep  = a[direction_id] + b[direction_id] .* x;
  vector[N] mu0_rep = a[direction_id];
  vector[N] y_rep;
  for (n in 1:N) y_rep[n] = normal_rng(mu_rep[n], sigma);
}
"

cf_fit <- function() {
  set.seed(1)
  df <- data.frame(direction = factor(rep(c("A", "B"), each = 20)), x = runif(40))
  df$y <- ifelse(df$direction == "A", 1 + 2 * df$x, -1 * df$x) + rnorm(40, sd = 0.3)
  sp <- sampler_params(modifyList(.default_sampler_params, list(
    chains = 1L, parallel_chains = 1L, iter_warmup = 200L, iter_sampling = 100L,
    refresh = 0L, show_messages = FALSE)))
  stanmodel() |> compile_from_string(cf_code) |> fit(data = df, sampler_params = sp)
}

cf_contrast <- function(gq_var) {
  contrast(contrast_var = "x", continuous = TRUE, two_way = "direction",
           two_way_levels = c("A", "B"), outcome = outcome_spec(gq_var = gq_var))
}

cf_run <- function(m, gq_var, ...) {
  suppressMessages(capture.output(
    res <- infer_contrasts(m, cf_contrast(gq_var), draws_method = "counterfactual",
                           gq_steps = 2, ...)))
  res$scenario_draws
}

test_that("counterfactual contrasts pair the same draws when gq_n_draws is set", {
  skip_if_not(cmdstanr::cmdstan_version() > "0", "CmdStan not available")
  m <- cf_fit()

  sd0 <- cf_run(m, "mu0_rep", gq_n_draws = 50)
  expect_equal(nrow(sd0), 50L)
  expect_true(all(sd0$contrast_A == 0) && all(sd0$contrast_B == 0))

  # identity check (from the finger paper's ii_tract_counterfactual.R)
  sd1 <- cf_run(m, "mu_rep") |> dplyr::arrange(.draw)
  b  <- posterior::as_draws_df(posterior::subset_draws(m$draws, "b"))
  dx <- diff(range(m$data$x))
  expect_lt(max(abs(sd1$contrast_A - b[["b[1]"]][sd1$.draw] * dx),
                abs(sd1$contrast_B - b[["b[2]"]][sd1$.draw] * dx)), 1e-4)

  expect_warning(cf_run(m, "y_rep", gq_n_draws = 10), "mu_rep")
})

test_that("observation index is read from the brackets, not the variable name", {
  d <- posterior::draws_df(`mu0_rep[1]` = 1:2, `mu0_rep[2]` = 3:4, `mu0_rep[3]` = 5:6)
  m <- list(draws = d, variables = "mu0_rep")
  out <- stanworkflow:::.extract_spec_draws(m, data.frame(g = 1:3),
                                            outcome_spec(gq_var = "mu0_rep"), "g")
  expect_equal(out$.obs_idx, rep(1:3, 2))
  expect_equal(out$g, out$.obs_idx)
})
