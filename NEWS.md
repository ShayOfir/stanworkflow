# stanworkflow 0.0.0.9000 (development version)

## Bug fixes (results can change)

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
