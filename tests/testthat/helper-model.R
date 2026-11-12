# One tiny model shared by the CmdStan tests. It is compiled once per test run
# (lazily, so tests that skip without CmdStan never compile) and fitted once.
#
#   y = a[direction] + g[group] + b[direction] * x + e
#
# GQ: mu_rep (expectation), mu0_rep (expectation without the x term),
# y_rep (predictive draw), log_lik; and bg[direction, group] = b[direction] + g[group]
# with mu_bg_rep = mu0_rep + bg[direction, group] * x, a deterministic two-index
# "slope" for averaging tests (no RNG, so the other draws do not change).
# Everything is linear, so counterfactual
# contrasts have closed forms: on x, b[d] * dx; on group, g[k] - g[j].

skip_if_no_cmdstan <- function() {
  skip_if(is.null(suppressWarnings(cmdstanr::cmdstan_version(error_on_NA = FALSE))),
          "CmdStan not available")
}

tm_code <- "
data {
  int<lower=1> N;
  int<lower=1> n_direction;
  array[N] int<lower=1, upper=n_direction> direction_id;
  int<lower=1> n_group;
  array[N] int<lower=1, upper=n_group> group_id;
  vector[N] x;
  vector[N] y;
  real<lower=0> prior_b_sd;
  int<lower=0, upper=1> skip_likelihood;
}
parameters {
  vector[n_direction] a;
  vector[n_group] g;
  vector[n_direction] b;
  real<lower=0> sigma;
}
model {
  a ~ normal(0, 1); g ~ normal(0, 1); b ~ normal(0, prior_b_sd); sigma ~ exponential(1);
  if (!skip_likelihood) y ~ normal(a[direction_id] + g[group_id] + b[direction_id] .* x, sigma);
}
generated quantities {
  vector[N] mu0_rep = a[direction_id] + g[group_id];
  vector[N] mu_rep  = mu0_rep + b[direction_id] .* x;
  matrix[n_direction, n_group] bg = rep_matrix(b, n_group) + rep_matrix(g', n_direction);
  vector[N] mu_bg_rep;
  for (n in 1:N) mu_bg_rep[n] = mu0_rep[n] + bg[direction_id[n], group_id[n]] * x[n];
  vector[N] y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    y_rep[n]   = normal_rng(mu_rep[n], sigma);
    log_lik[n] = normal_lpdf(y[n] | mu_rep[n], sigma);
  }
}
"

# Model folder in the layout regimen()/wflow() expect: <models>/tm/{model_code.stan, priors.csv}
tm_models <- file.path(tempdir(), "swf_models")
dir.create(file.path(tm_models, "tm"), recursive = TRUE, showWarnings = FALSE)
writeLines(tm_code, file.path(tm_models, "tm", "model_code.stan"))
write.csv(data.frame(name = "prior_b_sd", value = 2), file.path(tm_models, "tm", "priors.csv"),
          row.names = FALSE)

tm_data <- function() {
  set.seed(1)
  d <- expand.grid(i = 1:8, group = c("c", "s1", "s2"), direction = c("FLEX", "EXT"),
                   stringsAsFactors = FALSE)
  d$direction <- factor(d$direction, levels = c("FLEX", "EXT"))
  d$x <- round(runif(nrow(d)), 3)
  d$y <- ifelse(d$direction == "FLEX", 1 + 2 * d$x, -d$x) +
    c(c = 0, s1 = -0.5, s2 = 0.5)[d$group] + rnorm(nrow(d), sd = 0.3)
  d$group <- factor(d$group, levels = c("c", "s1", "s2"))   # grouping columns as factors
  d$i <- NULL
  d
}

# Short chains for every fit in the suite (a list, as regimen() takes it).
tm_sampler <- list(chains = 2L, parallel_chains = 2L, iter_warmup = 200L, iter_sampling = 200L,
                   refresh = 0L, show_messages = FALSE, show_exceptions = FALSE, seed = 1L)

.tm <- new.env()

# Compiled stanmodel. compile_from_file() on the same .stan path reuses the
# executable, so wflow()/regimen() do not recompile. The GQ cache is seeded with
# the same executable so gq() does not compile a second copy.
tm_compiled <- function() {
  if (is.null(.tm$compiled)) {
    m <- stanmodel() |> compile_from_file(file.path(tm_models, "tm", "model_code.stan"))
    m$priors$prior_b_sd <- 2
    m$sampler_params <- sampler_params(modifyList(sampler_params(), tm_sampler))
    assign(m$model_code, m$compiled_model, envir = stanworkflow:::.gq_model_cache)
    .tm$compiled <- m
  }
  .tm$compiled
}

# Fitted stanmodel (posterior), fitted once per run.
tm_fit <- function() {
  if (is.null(.tm$fit)) .tm$fit <- tm_compiled() |> fit(data = tm_data())
  .tm$fit
}

tm_draws <- function(m, var) posterior::as_draws_df(posterior::subset_draws(m$draws, var))

# infer_contrasts() prints progress; keep test output clean.
quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(expr))
  out
}
