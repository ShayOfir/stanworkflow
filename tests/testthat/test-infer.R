# Descriptive contrasts (draws_method = "descriptive") are summaries of the
# fitted GQ draws, so every number can be recomputed by hand from the
# S x N mu_rep matrix: contrast = mean over "hi" observations minus mean over
# "lo" observations, draw by draw.

mu_matrix <- function(m) unclass(posterior::as_draws_matrix(posterior::subset_draws(m$draws, "mu_rep")))

# per-draw mean of mu_rep over hi rows minus over lo rows (logical vectors over observations)
by_hand <- function(mu, hi, lo) unname(rowMeans(mu[, hi, drop = FALSE]) - rowMeans(mu[, lo, drop = FALSE]))

# coda::HPDinterval written out: the narrowest window holding round(n * prob) gaps
hpdi_ref <- function(v, prob) {
  v <- sort(v); n <- length(v); gap <- max(1, min(n - 1, round(n * prob)))
  i <- which.min(v[(gap + 1):n] - v[1:(n - gap)])
  c(v[i], v[i + gap])
}

mu_os <- outcome_spec(gq_var = "mu_rep")

# ---- contrast() / outcome_spec() constructors (no Stan) ----------------------

test_that("contrast() validates its arguments", {
  expect_error(contrast("x", absolute_thresholds = TRUE), "continuous")
  expect_error(contrast("x", continuous = TRUE, contrast_levels = c(1.5, 0.2)))
  expect_error(contrast("g", contrast_levels = 1:2))
  expect_error(contrast("g", outcome = list(gq_var = "mu_rep")))
  expect_error(outcome_spec(gq_var = c("a", "b")))
  co <- contrast("g", contrast_levels = c("c", "s1"), grouping_vars = "d")
  expect_s3_class(co, "contrast")
  expect_equal(get_label(co), "s1 - c")
  expect_output(print(co), "categorical")
  expect_output(print(co), "s1 - c")
  expect_equal(get_label(contrast("x", continuous = TRUE)), "75% - 25% (default quantiles)")
  expect_equal(get_label(contrast("x", continuous = TRUE, absolute_thresholds = TRUE,
                                  contrast_levels = c(0.7, 0.3))), "(>=0.7) - (<0.3)")
  expect_match(get_label(cf <- contrast("x", two_way = "d", two_way_levels = c("A", "B"))), "2-way")
  expect_output(print(cf), "2-way")
  expect_output(print(os_hurdle()), "marginal")
  expect_output(print(mu_os), "mu_rep")
})

# ---- outcome_spec fan-out on fake draws (no Stan) ----------------------------

fake_model <- function() {
  # 4 observations, 3 draws; values chosen so contrasts are easy to read
  cols <- c(lapply(1:4, function(j) c(1, 2, 3) * j),        # mu_marginal[j] = draw * j
            lapply(1:4, function(j) rep(j / 10, 3)),        # p_respond[j]   = j / 10
            lapply(1:4, function(j) c(1, 2, 3) + j))        # mu_y[j]        = draw + j
  names(cols) <- c(sprintf("mu_marginal[%d]", 1:4), sprintf("p_respond[%d]", 1:4),
                   sprintf("mu_y[%d]", 1:4))
  list(draws = do.call(posterior::draws_df, cols),
       variables = c("mu_marginal", "p_respond", "mu_y"),
       data = data.frame(g = c("a", "a", "b", "b")))
}

test_that("os_hurdle() fans out into the three margins, each a descriptive contrast", {
  m <- fake_model()
  r <- suppressWarnings(quiet(infer_contrasts(m, contrast("g", contrast_levels = c("a", "b"),
                                                          outcome = os_hurdle()))))
  expect_named(r, c("marginal", "extensive", "intensive", "contrast"))
  for (k in c("marginal", "extensive", "intensive")) expect_s3_class(r[[k]], "inference")
  # b - a: marginal = draw * (3.5 - 1.5), extensive = 0.35 - 0.15, intensive = 2
  expect_equal(r$marginal$posterior$contrast, c(2, 4, 6))
  expect_equal(r$extensive$posterior$contrast, rep(0.2, 3))
  expect_equal(r$intensive$posterior$contrast, rep(2, 3))
})

test_that("outcome_spec transform_fn and filter_fn are applied before contrasting", {
  m <- fake_model()
  os <- outcome_spec(gq_var = "mu_y", transform_fn = function(y) 10 * y,
                     filter_fn = function(df) dplyr::filter(df, .obs_idx != 4))
  r <- suppressWarnings(quiet(infer_contrasts(m, contrast("g", contrast_levels = c("a", "b"),
                                                          outcome = os))))
  # b keeps only obs 3: 10 * (draw + 3) - 10 * (draw + 1.5) = 15
  expect_equal(r$posterior$contrast, rep(15, 3))
  expect_error(suppressWarnings(quiet(infer_contrasts(m, contrast("g", contrast_levels = c("a", "b"),
                                                    outcome = outcome_spec("nope"))))), "nope")
})

test_that("several contrasts return a list of inference objects in input order", {
  m <- fake_model()
  cs <- list(contrast("g", contrast_levels = c("a", "b"), outcome = os_marginal()),
             contrast("g", contrast_levels = c("b", "a"), outcome = os_marginal()))
  r <- suppressWarnings(quiet(infer_contrasts(m, cs)))
  expect_length(r, 2L)
  expect_equal(r[[1]]$posterior$contrast, -r[[2]]$posterior$contrast)
})

# ---- descriptive contrasts on the real fit -----------------------------------

test_that("categorical one-way descriptive contrast: second level minus first, per grouping cell, HPDI via coda", {
  skip_if_no_cmdstan()
  m <- tm_fit(); d <- m$data; mu <- mu_matrix(m)
  co <- contrast("group", contrast_levels = c("c", "s1"), grouping_vars = "direction", outcome = mu_os)
  r <- suppressWarnings(quiet(infer_contrasts(m, co, prob = 0.9)))
  expect_s3_class(r, "inference")
  expect_equal(nrow(r$summary), 2L)
  expect_equal(unique(r$summary$contrast_label), "s1 - c")
  for (dir in c("FLEX", "EXT")) {
    p <- dplyr::filter(r$posterior, direction == dir)
    expect_equal(p$contrast,
                 by_hand(mu, d$group == "s1" & d$direction == dir, d$group == "c" & d$direction == dir)[p$.draw],
                 tolerance = 1e-10)
    s <- dplyr::filter(r$summary, direction == dir)
    expect_equal(c(s$hpdi_lo, s$hpdi_hi), hpdi_ref(p$contrast, 0.9))
    expect_equal(s$P_positive, mean(p$contrast > 0))
    expect_equal(s$contrast_mean, mean(p$contrast))
    expect_equal(s$prob, 0.9)
  }
  expect_output(print(r), "descriptive")
})

test_that("HPDI matches rethinking::HPDI where rethinking is installed", {
  skip_if_not_installed("rethinking")
  set.seed(2); v <- rnorm(401)
  expect_equal(unname(hpdi_ref(v, 0.89)), unname(rethinking::HPDI(v, 0.89)))
})

test_that("categorical contrast without contrast_levels warns and uses the first two levels", {
  skip_if_no_cmdstan()
  m <- tm_fit()
  expect_warning(r <- quiet(infer_contrasts(m, contrast("group", outcome = mu_os))),
                 "contrast_levels not specified")
  expect_equal(unique(r$posterior$contrast_label), "s1 - c")
})

test_that("continuous descriptive contrast compares the top and bottom quantile groups, in either order", {
  skip_if_no_cmdstan()
  m <- tm_fit(); d <- m$data; mu <- mu_matrix(m)
  q <- quantile(d$x, c(0.25, 0.75))
  expected <- by_hand(mu, d$x >= q[2], d$x < q[1])
  for (lv in list(NULL, c(0.75, 0.25), c(0.25, 0.75))) {
    co <- contrast("x", continuous = TRUE, contrast_levels = lv, outcome = mu_os)
    r <- suppressWarnings(quiet(infer_contrasts(m, co)))
    expect_equal(r$posterior$contrast, expected[r$posterior$.draw], tolerance = 1e-10)
  }
})

test_that("absolute thresholds and filter_draws select the expected observations", {
  skip_if_no_cmdstan()
  m <- tm_fit(); d <- m$data; mu <- mu_matrix(m)
  co <- contrast("x", continuous = TRUE, absolute_thresholds = TRUE, contrast_levels = c(0.7, 0.3),
                 grouping_vars = "direction", filter_draws = quote(group == "s1"), outcome = mu_os)
  r <- suppressWarnings(quiet(infer_contrasts(m, co)))
  p <- dplyr::filter(r$posterior, direction == "EXT")
  s1 <- d$group == "s1" & d$direction == "EXT"
  expect_equal(p$contrast, by_hand(mu, s1 & d$x >= 0.7, s1 & d$x < 0.3)[p$.draw], tolerance = 1e-10)
  # the slope in EXT is negative, so high x has lower mu
  expect_lt(mean(p$contrast), 0)
  expect_error(suppressWarnings(quiet(infer_contrasts(
    m, contrast("x", continuous = TRUE, filter_draws = quote(group == "zzz"), outcome = mu_os)))),
    "filter_draws")
})

test_that("two-way descriptive contrast: per-direction contrasts, scenarios, probabilities sum to 1", {
  skip_if_no_cmdstan()
  m <- tm_fit(); d <- m$data; mu <- mu_matrix(m)
  co <- contrast("group", contrast_levels = c("c", "s2"), two_way = "direction",
                 two_way_levels = c("FLEX", "EXT"), outcome = mu_os)
  r <- suppressWarnings(quiet(infer_contrasts(m, co)))
  sd <- r$scenario_draws
  expect_equal(nrow(sd), posterior::ndraws(m$draws))
  for (dir in c("FLEX", "EXT"))
    expect_equal(sd[[paste0("contrast_", dir)]],
                 by_hand(mu, d$group == "s2" & d$direction == dir,
                         d$group == "c" & d$direction == dir)[sd$.draw], tolerance = 1e-10)
  expect_identical(sd$scenario == "Affects FLEX only", sd$contrast_FLEX > 0 & sd$contrast_EXT <= 0)
  expect_equal(sum(r$scenario_probs$probability), 1)
  expect_equal(c(r$lv1, r$lv2, r$hi_lbl, r$lo_lbl), c("FLEX", "EXT", "s2", "c"))
  expect_identical(r$contrast, co)
  expect_output(print(r), "Affects")
})

test_that("two-way descriptive contrast on a continuous variable uses the quantile split", {
  skip_if_no_cmdstan()
  m <- tm_fit(); d <- m$data; mu <- mu_matrix(m)
  co <- contrast("x", continuous = TRUE, two_way = "direction", grouping_vars = "group",
                 outcome = mu_os)
  r <- suppressWarnings(quiet(infer_contrasts(m, co)))
  sd <- dplyr::filter(r$scenario_draws, group == "s1")
  q <- quantile(d$x, c(0.25, 0.75))   # split is over all rows of the two directions
  s1 <- d$group == "s1" & d$direction == "FLEX"
  expect_equal(sd$contrast_FLEX, by_hand(mu, s1 & d$x >= q[2], s1 & d$x < q[1])[sd$.draw],
               tolerance = 1e-10)
  p <- dplyr::summarise(r$scenario_probs, p = sum(probability), .by = group)
  expect_equal(p$p, rep(1, 3))
})
