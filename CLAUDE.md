# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Two things in one git repo:

1. **`stanworkflow`** — an R package (root, `R/`) that wraps `cmdstanr` in a declarative
   compile → prior-predictive → fit → diagnose → infer pipeline.
2. **`dev/`** — the research sandbox that consumes the package (excluded via `.Rbuildignore`).
   It holds the actual study: finger-individuation index (II) in stroke patients vs. controls,
   flexion vs. extension, regressed on lesion load in descending-tract and grey-matter ROIs.

Most day-to-day work happens in `dev/`. The package exists to serve it, so package changes are
usually driven by a gap hit while running an analysis.

## Environment and commands

Windows + Rtools45. The package repo has no renv and no `.Rprofile` (removed 2026-11);
DESCRIPTION is the dependency spec. If Stan compilation cannot find gcc/g++, set `CC`/`CXX`
in your user `~/.Rprofile` / `~/.Renviron` or run `cmdstanr::check_cmdstan_toolchain(fix = TRUE)`.
Development happens in the `stanworkflow-dev` worktree (branch `main`); the original
`stanworkflow` folder is the frozen finger-paper checkout (branch `claude/analysis-fixes`, tag
`finger-paper-2026-10`) and holds `dev/` (~149 GB) — never run `R CMD build/check` there.

```r
devtools::load_all(".")         # load package for interactive work
devtools::document()            # regenerate NAMESPACE + man/ (roxygen2; never hand-edit either)
devtools::check()
cmdstanr::cmdstan_version()     # verify the Stan toolchain before debugging a fit
```

`dev/` scripts load the package by absolute path:
`devtools::load_all("c:/Scripts/Brain-Lab/bayesworkflow/stanworkflow")`, then `setwd()` into
`dev/` because every model/results path in those scripts is relative to it.

Tests:

```r
devtools::test()                                  # whole suite
testthat::test_file("tests/testthat/test-stanmodels.R")
```

Tests: `helper-model.R` holds one tiny linear model (`y ~ a[direction] + g[group] +
b[direction]*x`, GQ `mu_rep`, `mu0_rep`, `y_rep`, `log_lik`) compiled once per run
(`tm_compiled()`, `tm_fit()`; also a `models/tm` folder for `regimen()`/`wflow()`).
Because it is linear, contrasts have closed forms (b·Δx, g[k] − g[j]) and the tests check
them draw by draw. `test-counterfactual.R` (pairing with `gq_n_draws`, identity, y_rep
warning, categorical, the finger-paper two-way pattern), `test-infer.R` (descriptive
contrasts recomputed by hand, HPDI, `outcome_spec`/`os_hurdle` on fake draws),
`test-regimen.R` (lifecycle, checkpoints, `resume`, `quick_fit`, save/load),
`test-parameter.R` (parameter route, `ppc_by_parameter`/`ppc_scenarios` and aliases),
`test-wflow.R` (wflow steps, write/load incl. the JWileymisc-format fixture in
`fixtures/`, ppc/diagnose), `test-loo.R`, `test-stanmodels.R`. CmdStan tests call
`skip_if_no_cmdstan()`. Do not compile further models in tests; extend the helper model.

## Package architecture

Three S3 layers, each wrapping the one below. Read them in this order.

**`stanmodel` (`R/stanmodels.R`)** — thin `cmdstanr` wrapper.
`compile_from_file/string()` → `get_stan_data()` → `fit()` → `gq()` → `diagnose()` / `ppc()`.
After `fit()` the object carries `$draws`, `$variables`, `$data`, `$extra`, `$cond_df`,
`$stan_data`, `$fitted_model`. Also defines the `priors` and `sampler_params` classes
(defaults: 4 chains, 1000/1000, `adapt_delta = 0.8`).

**`wflow` (`R/wflow.R`)** — lifecycle for one model:
`attach_data()` → `tune()` (prior predictive) → `fit_posterior()` → `run_diagnostics()` →
`run_inference()`. Slots `$compiled`, `$prior_fit`, `$post_fit`, `$data`, `$diagnostics`,
`$inference`, `$tuned`. `run.wflow()` chains all of it. Note `tune()` does **not** optimise
anything — it fits with `skip_likelihood = 1` so you can eyeball the prior predictive and
hand-edit `priors.csv`.

**`regimen` (`R/regimen.R`)** — the declarative front end, and what `dev/` scripts use. You
declare `model_name`, `data_list`, `contrasts`, `results_path`, sampler params, then:

```r
reg <- regimen(...) |> run(steps = c("prepare", "tune", "fit", "infer"))
```

- Steps are `prepare` / `tune` / `fit` / `quick_fit` / `infer`. `fit` and `quick_fit` are
  mutually exclusive.
- `quick_fit` = posterior sampling only, skipping prior predictive, PPC, and diagnostics. Use
  it for mass-univariate loops where per-ROI PPC output would be gigabytes.
- Each step checkpoints to `<results_path>/<model_name>_{prior,post,inference}.rds`, and later
  steps reload from disk by default (`load_prior` / `load_post`), so a crashed run resumes
  cheaply. `resume()` inspects which files exist and re-runs only what is missing.

### The naming convention that drives everything

`get_stan_data()` (`R/stanmodels.R:176`) parses the Stan `data` block and fills it **by
variable name** from the attached data frame. This is the single most important thing to
understand before writing a new model:

| Stan declaration | Filled with |
|---|---|
| `N` | `nrow(df)` |
| `n_<col>` | `dplyr::n_distinct(df[[col]])` |
| `<col>_id` | `as.integer(factor(df[[col]]))` |
| `prior_*` | value from `models/<model>/priors.csv` |
| `<col>` | `df[[col]]` verbatim |
| `N_cond_array`, `<col>_cond_array`, `<col>_id_cond_array` | from `data_list$cond_df` |
| anything else | must be supplied in `data_list$extra` |

Consequences worth internalising:

- **`<col>_id` index order is the factor level order of that column.** Always set levels
  explicitly when building data (`dev/helper_functions.R` does this for `finger`, `direction`,
  `group`) — otherwise Stan indices silently reorder and, e.g., `group_id == 1` stops meaning
  "stroke". When mapping indices back to labels, recover them with
  `levels(droplevels(factor(df$col)))` rather than hard-coding.
- Prior data variables **must** be named `prior_*`; `.extract_prior_names()` greps
  `\bprior_[A-Za-z0-9_]+` out of the Stan source to build the `priors` object.
- Every model must declare `int<lower=0,upper=1> skip_likelihood` and gate its likelihood on
  it, or prior predictive checks cannot run.
- `.parse_stan_data_block()` uses a regex (`data\s*\{[^}]*\}`) — a `}` inside the data block
  (including in a comment) will truncate parsing.

### Model folder convention

`models/<name>/` containing `model_code.stan`, `priors.csv` (columns `prior,value`), and
optionally `diagnostics_params.csv` (one column, `parameters`, naming what `diagnose()` should
summarise). `dev/models/` holds `m1` (group × finger × direction), `m2` (adds a single
univariate lesion slope), `m3` (MVNorm over many ROI slopes with an LKJ correlation).

### The inference layer

`R/infer.R` is ~2000 lines and by far the most intricate part.
**`reports/infer-contrasts-explanation.md` is a detailed walkthrough — read it before changing
anything here.** The two orthogonal abstractions:

- **`outcome_spec`** — *what* to extract from draws (`gq_var`, default `"y_rep"`), plus
  `filter_fn` / `transform_fn`. `os_marginal()`, `os_extensive()`, `os_intensive()`,
  `os_hurdle()` exist for hurdle models, where there is no single "the outcome".
- **`contrast`** — *how* to compare: `contrast_var`, `grouping_vars` (facets),
  `continuous` + `contrast_levels` (quantile probabilities, or absolute cutoffs with
  `absolute_thresholds = TRUE`), `effect_positive` (sign convention), `filter_draws` (an
  unevaluated quosure, applied after extraction), and `two_way` + `two_way_levels`.

`draws_method` selects between fundamentally different routes:

- `"descriptive"` — reuse the stored `y_rep` draws. Cheap.
- `"counterfactual"` — call `.simulate_at()` → `gq()` once per grid point (`gq_steps` points
  across the observed range of `contrast_var`), then difference each point against the
  baseline. Very expensive in both time and disk: it re-runs `generate_quantities` per grid
  point and the resulting objects run to hundreds of megabytes per ROI.
- `"parameter"` — summarise a parameter's draws per factor cell (optionally × `delta`), level
  differences, two-way scenarios. No `gq()`; exact for linear, identity-link models.

A `two_way` contrast classifies **each posterior draw** into `"Affects both"` /
`"Affects <lv1> only"` / `"Affects <lv2> only"` / `"Affects neither"` (`.define_scenarios()`),
and reports the proportion of draws in each — a joint sign probability, not four separate
tests.

Virtual contrast variables let a *column of a matrix in `$extra`* act as a contrast variable:
`contrast_var = "roi"` is expanded to `"roi[<colname>]"` per column by
`.expand_numeric_array_contrast()`, and `.simulate_at()` then sets that column to the
counterfactual value **and zeroes all the others**.

### Serialization

`write.wflow()` / `load_wflow()` / `save_regimen()` / `load_regimen()` go through
internal `fst` helpers `.save_rds_fst()` / `.read_rds_fst()` (`R/wflow.R`; ZSTD, level 100,
byte-identical to the former `JWileymisc::saveRDSfst()`); plain RDS files also load.
`save_regimen()` and the inference checkpoint use plain `saveRDS()`. Fitted
objects are large — a single `m2` posterior is ~140 MB, an `m3` counterfactual inference object
several hundred MB. When you only need a few parameters, `posterior::subset_draws()` the
variable you want and `rm()` + `gc()` the wflow rather than holding several in memory.


## Known pitfalls

These cost real debugging time and are not visible from any single file.

**Counterfactual contrasts are algebraically redundant for models that are linear in the
contrast variable.** For `m2`/`m3`, `E[y | lesion = a] − E[y | lesion = 0] = −β·a` exactly, so
the sign-based scenario probabilities *are* the joint sign probabilities of `β` and are
invariant to `a`. Reading `beta_lesion` / `beta_roi_3d` from the draws is exact, and three
orders of magnitude cheaper than the `gq()` sweep. `draws_method = "parameter"`
(`.compute_parameter_contrast()` in `R/infer.R`; inputs `parameter`/`factors`/`delta` on the
`contrast`) implements that route; `test-parameter.R` checks the identity against the
counterfactual `mu_rep` contrast. Counterfactual sweeps only earn their cost for non-linear
models or non-identity links.

**`gq_n_draws` in counterfactual contrasts (fixed on main, 2026-11).** Before the fix,
`gq.stanmodel()` took a fresh `sample.int()` subsample at every grid point and baseline and grid
point were joined on renumbered `.draw`, pairing *different* posterior draws.
`.get_yrep_draws()` now samples once and passes `draw_ids` to `gq()`. Results made with
`gq_n_draws` on older versions (e.g. `dev/results_3.9.26/scenarios_probs_m3.csv`) are not
trustworthy.

**Counterfactual contrasts difference `y_rep`, not the linear predictor.** Since `y_rep` is
`normal_rng(mu, sigma)`, every grid point injects fresh residual noise whose SD is constant in
the counterfactual level while the signal grows with it. Sign-based probabilities are therefore
biased toward 0.25, and the same contrast evaluated at 50% vs. 100% lesion gives different
answers when it mathematically should not. `infer_contrasts()` now warns in this case; use an
expectation GQ (`mu_rep`) via `outcome_spec(gq_var = "mu_rep")`.

**`load_kinematics()` ignores its `scale_to_intact` argument** (`dev/helper_functions.R:17`) —
`y` is always divided by the intact-group 95th percentile per finger × direction × hand,
regardless of what the call site passes.

**Prior structure encodes the scientific hypothesis in `m2`/`m3`.** The hierarchical prior makes
the flexion and extension slopes exchangeable draws around a common mean, which is itself the
"shared substrate" hypothesis, with prior `P(affects both) ≈ 0.38–0.40`. Any scenario-probability
table must be reported against that prior baseline or it is uninterpretable. See
`dev/claude/review_by_opus_12.9.26.qmd` for the full argument, the corrected prior computation
(the pair is a *scale mixture*, so the bivariate-normal orthant formula does not apply), and the
graded plan for the modelling work that follows.

## dev/ orientation

`dev/` accumulates dated scratch scripts; most are superseded. The current entry points:

- `dev/helper_functions.R` — `load_kinematics()`, `load_lesion_data()`,
  `prepare_anatomy_data()` (returns a per-hemisphere list; `extra_option = TRUE` packs ROIs
  into an `extra$roi` matrix for `m3`, `FALSE` keeps one column per ROI for `m2`),
  `filter_rois()` (ROI presets 1/2/3), `scale_na_df()` (q95 scaling, **per hemisphere** — so
  `β` magnitudes are not comparable across hemispheres; use scale-free summaries for any L-vs-R
  claim).
- `dev/analysis_3.9.26.R` — the live driver: `perform_step2()` / `infer_step2()` for the
  mass-univariate `m2` loop over ROI × hemisphere, `perform_step3()` / `infer_step3()` for `m3`.
  Note line 16 (`contrasts <-`) is an incomplete assignment that errors on a clean `source()`;
  the file is meant to be run interactively, region by region.
- `dev/results_3.9.26/` — current results, including
  `MVNorm_and_mass_univariate_results.xlsx` (sheets `multiv`, `mass_univ`).
- `dev/simulate.R` — simulation helpers, the starting point for any design/power or SBC work.

## Dev note: Windows Defender and CmdStan

Defender can block a freshly compiled `model_code.exe` ("system error 225, file contains a virus")
in `tempdir()`; the chain then fails to start. `fit()`/`gq()` now stop on a failed chain instead
of using partial draws. Workaround: add a Defender exclusion for the CmdStan folder and R's
temp dir (`%TEMP%`), then re-run.
