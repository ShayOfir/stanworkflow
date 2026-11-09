# stanworkflow 0.0.0.9000 (development version)

## Bug fixes (results can change)

* Counterfactual contrasts with `gq_n_draws` set paired different posterior
  draws: each grid point took a fresh random subsample and baseline and grid
  point were then joined on renumbered `.draw` ids. The subsample is now drawn
  once per contrast group and reused at every grid point (new
  `gq(draw_ids = )` argument). Counterfactual results computed with
  `gq_n_draws` before this fix are not trustworthy; results without
  `gq_n_draws` are unchanged.
* The observation index of a generated quantity was taken from the first
  digits in its name, so variables with a digit in the name (e.g. `mu0_rep`,
  `y2_rep`) got index 0 for every observation and observations were mixed
  up in contrasts and grouped PPCs. It is now read from the brackets.
  `y_rep`, `mu_rep` and other digit-free names are unaffected.
* `fit.stanmodel()` now uses its documented `stan_data` argument; it was
  silently ignored and the data always rebuilt from `data`.

## Behaviour changes

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
