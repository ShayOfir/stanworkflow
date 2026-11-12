# stanworkflow (development version)

* Stan data-block names are read from stanc (`CmdStanModel$variables()`)
  instead of a regex parser that dropped variables after a comment with
  braces, `int a, b;` or a declaration split over lines.
* `<col>_id` / `n_<col>`: a non-factor column is coded with a
  locale-independent (radix) sort, the same on every machine, and a warning
  lists the level order; empty factor levels are dropped with a warning
  (pass `n_<col>` and `<col>_id` via `extra` to keep them). Factors keep
  their level order. Counterfactual and parameter contrasts use the same
  coding. Mixed-case character levels (e.g. "a" vs "B") may now be numbered
  differently than before on Windows.
* Unfilled Stan data variables are reported once, as a message.

# stanworkflow 0.1.0

## Documentation

* README with installation (CmdStan from the Stan r-universe, then
  `remotes::install_github("ShayOfir/stanworkflow")`), an example and a note
  on Windows Defender exclusions.
* Three vignettes: `vignette("stanworkflow")` (getting started),
  `vignette("contrasts")` (descriptive vs counterfactual vs parameter
  contrasts, `mu_rep` vs `y_rep`, the linear identity, `delta`, averaging with
  `weights`) and `vignette("hurdle")` (hurdle outcomes with `outcome_spec()`).
  They are precomputed with CmdStan (sources `vignettes/*.Rmd.orig`, re-run
  with `data-raw/vignettes.R`), so building them needs no Stan.
* Citation metadata: `citation("stanworkflow")`, `CITATION.cff`, `.zenodo.json`.

## Bug fixes

* `fit()` could silently lose chains with all return codes 0. cmdstanr marks a
  chain successful only after parsing its console output; on Windows a fast
  chain can exit before that output is read, and `draws()` then skips its
  complete CSV. `fit()` now rebuilds the fit from all CSVs when that happens
  (and stops if a CSV is really missing). Separately, cmdstanr's CSV names take
  their random part from R's RNG, so after `set.seed()` fits started in the same
  minute shared file names; `fit()` now passes a seed-independent
  `output_basename`.
* `write.wflow()` (used for every checkpoint) now defaults to `compress = 50`
  instead of 100. ZSTD at level 100 needs about 9 times the object's size in
  extra RAM, and saving a large posterior could kill R silently. The file
  format is unchanged and older files load as before.
* `fit.regimen()` no longer warns "not marked as tuned" after `tune()`: the
  prior wflow reloaded from disk is marked tuned (`quick_fit` still warns, since
  it skips the prior checks).
* The grouped PPC `$summary` with a single `.by` variable named its first
  column `...1`; it now carries the variable name.

* Counterfactual contrasts ignored `infer_contrasts(data = )` when running
  `generate_quantities()` (it always used the fitted data; `data` only set the
  grid range and the feature columns). The supplied data are now used, with the
  factor coding of the fitted data, so other covariates can be held at chosen
  values. Calls whose `data` differs from the fitted data only in the contrast
  variable (which the counterfactual overwrites) give the same results as before.

* `fit()` and `gq()` now stop if any CmdStan chain returns a non-zero code.
  Previously cmdstanr only warned and the draws of the surviving chains were
  used silently (e.g. 200 instead of 400 draws when one chain's process failed
  to launch).

## New features

* `draws_method = "parameter"` in `infer_contrasts()`, `regimen()` and
  `run()`: inference straight from the draws of a model parameter, with no
  `generate_quantities()` calls. The contrast carries the inputs
  (`contrast(..., parameter = "b", factors = "direction", delta = "q10q90")`).
  Results hold per-level summaries (median, HPDI, probability of direction),
  differences between levels of a factor (one pair or all pairs), and two-way
  scenarios, in the usual `inference` structure, each labelled with its
  `scale`. `delta` (NULL, a number, `"iqr"`, `"q10q90"`, `"2sd"`, `"range"`)
  multiplies the parameter by a change in the predictor. For a model linear in
  the predictor with an identity link, parameter x `delta = "range"` equals the
  counterfactual contrast on an expectation GQ draw by draw (tested).
  `contrast(link = )` declares a non-identity link, which gives a one-time
  message pointing to the counterfactual route for outcome-scale effects.
* Parameter route: `grouping_vars` now selects the factors kept as cells; the
  other factors (except a categorical `contrast_var` and `two_way`) are averaged
  per draw before summaries, differences and scenarios. New
  `contrast(weights = )`: `"equal"` (default, each level the same) or `"data"`
  (rows per level in the data, joint counts within each kept cell; reproduces
  the counterfactual average over rows for a linear model). The `scale` label
  names what was averaged. `grouping_vars = NULL` (default) is unchanged: every
  factor is a cell.
* `parameter_draws()` and `parameter_scenarios()`: aliases of
  `ppc_by_parameter()` and `ppc_scenarios()` (old names unchanged).

## Bug fixes (results can change)

* `ppc_scenarios(scenario_levels = )` stopped with "object 'col1' not found";
  the two given levels are now used as the outcome columns.
* Counterfactual contrasts with `gq_n_draws` set paired different posterior
  draws: each grid point took a fresh random subsample and baseline and grid
  point were then joined on renumbered `.draw` ids. The subsample is now drawn
  once per contrast group and reused at every grid point (new
  `gq(draw_ids = )` argument). Counterfactual results computed with
  `gq_n_draws` before this fix are not trustworthy; results without
  `gq_n_draws` are unchanged.
* One-way continuous counterfactual contrasts (`continuous = TRUE`, no
  `two_way`) reported the GQ *level* at each grid point instead of its
  difference from the baseline, including a baseline-vs-itself row.
  `grid_summary` / `grid_posterior` now hold grid point minus baseline and
  omit the baseline row. Two-way and categorical contrasts are unchanged.
* The observation index of a generated quantity was taken from the first
  digits in its name, so variables with a digit in the name (e.g. `mu0_rep`,
  `y2_rep`) got index 0 for every observation and observations were mixed
  up in contrasts and grouped PPCs. It is now read from the brackets.
  `y_rep`, `mu_rep` and other digit-free names are unaffected.
* `fit.stanmodel()` now uses its documented `stan_data` argument; it was
  silently ignored and the data always rebuilt from `data`.
* Descriptive contrasts on a continuous variable (`continuous = TRUE`, one-way
  or two-way) split the data at the wrong cut. With the default
  `contrast_levels = c(0.75, 0.25)` "hi" was every value at or above the 25th
  percentile and "lo" the values below it, not the top quartile vs the bottom
  quartile; absolute thresholds given as `c(hi, lo)` were swapped the same
  way. "hi" is now values at or above the larger cut and "lo" values below the
  smaller cut, in whichever order the two levels are given. Calls that passed
  the lower level first (e.g. `c(0.25, 0.75)`) give the same numbers as
  before. Counterfactual contrasts never used this split and are unchanged.
  `get_label()` / `print()` of a categorical contrast now show the levels in
  the order of the computed difference (second minus first).
* Several categorical counterfactual contrasts on the same variable in one
  `infer_contrasts()` call (e.g. `s1 - c` and `s2 - c`) each received the
  draws of all of them: the `contrast_name` column used to pick each
  contrast's rows was dropped before the filter. Each result now holds only
  its own contrast. A single categorical contrast per call was unaffected.
* `print()` of an `inference` object failed for descriptive results (they
  carried no `$contrast`, so the plot code indexed an `NA` column) and for any
  one-way result without `grouping_vars`; categorical counterfactual results
  were labelled "descriptive". Descriptive results now carry `$contrast`
  (added, nothing removed) and print with their plot. Numbers are unchanged.

## Behaviour changes

* `quick_fit()` / `run(steps = "quick_fit")` on an unprepared regimen now
  prepares it, as documented; it stopped with "Regimen not prepared" before.
  `run(steps = c("prepare", "quick_fit"))` is unchanged.
* `filter_draws` given as `rlang::quo(...)` no longer triggers rlang's
  "Subsetting quosures with `[` is deprecated" warning (results unchanged).

* `infer_contrasts(draws_method = "counterfactual")` warns when the outcome is
  `y_rep`: differencing posterior predictive draws adds fresh residual noise at
  each grid point, which biases sign probabilities toward 0.5 (two-way
  scenarios toward 0.25). Add an expectation (e.g. `mu_rep`) to generated
  quantities and use `outcome_spec(gq_var = "mu_rep")`. A warning, not an
  error, so existing scripts still run with unchanged numbers.

## Dependencies and packaging

* `rethinking` removed: HPDI bounds now come from `coda::HPDinterval()`, which
  `rethinking::HPDI()` wraps (identical numbers).
* `JWileymisc` removed: `write.wflow()` / `load_wflow()` use `fst` directly and
  write the same file format, so files saved by earlier versions still load.
* Declared imports that were used but missing (`loo`, `patchwork`, `scales`,
  `stats`, `tibble`); `R (>= 4.2.0)`; licence GPL (>= 3).
* `renv` and the machine-specific `.Rprofile` removed from the repository.
