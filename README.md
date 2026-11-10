# stanworkflow
[![R-CMD-check](https://github.com/ShayOfir/stanworkflow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ShayOfir/stanworkflow/actions/workflows/R-CMD-check.yaml)

stanworkflow gives a Stan model you wrote yourself the Bayesian workflow that
formula interfaces such as brms give standard models: **compile → prior
predictive check → fit → diagnose → infer**, declared once in a `regimen` and
run with `cmdstanr`.

- The Stan `data` block is filled from a data frame **by naming convention**
  (`N`, `n_<col>`, `<col>_id`, `<col>`, `prior_*`), so the same model runs on new
  data without glue code.
- Every step **checkpoints** to disk; an interrupted analysis resumes where it
  stopped (`resume()`), and new contrasts reuse the saved fit.
- **Posterior contrasts** draw by draw, by three routes: *descriptive* (from
  the stored predictions), *counterfactual* (set a variable for every row,
  average the change in the expected outcome) and *parameter* (read a
  coefficient, optionally scaled by a change in the predictor); one-way, or
  two-way as joint sign-scenario probabilities. Hurdle outcomes via
  `outcome_spec()`.

## Installation

stanworkflow needs [CmdStan](https://mc-stan.org/cmdstanr/) and a C++
toolchain (Rtools on Windows).

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
remotes::install_github("ShayOfir/stanworkflow")
```

**Windows:** Microsoft Defender can block a freshly compiled model executable
("system error 225"), so the chain fails to start; stanworkflow then stops
rather than continue with fewer chains. Add Defender exclusions for the CmdStan
folder (`cmdstanr::cmdstan_path()`) and for R's temporary folder (`%TEMP%`).

## Example

```r
library(stanworkflow)

# models/lin/model_code.stan declares N, n_group, group_id, x, y, prior_b_sd and
# skip_likelihood, and puts mu_rep (expectation) and y_rep in generated quantities;
# models/lin/priors.csv holds prior_b_sd. See vignette("stanworkflow").
reg <- regimen(
  model_name   = "lin",
  data_list    = list(data = my_data),          # columns group, x, y
  models_folder = "models",
  results_path = "results",
  ppc_by       = "group",
  contrasts    = list(
    group = contrast("group", contrast_levels = c("control", "patient"),
                     outcome = outcome_spec(gq_var = "mu_rep")),
    x     = contrast("x", continuous = TRUE, grouping_vars = "group",
                     outcome = outcome_spec(gq_var = "mu_rep"))
  )
)
reg <- run(reg, steps = c("prepare", "tune", "fit", "infer"))
get_inference(reg, "x")
```

## Documentation

- `vignette("stanworkflow")` — getting started: model folder, naming
  convention, `regimen`, checkpoints.
- `vignette("contrasts")` — descriptive, counterfactual and parameter
  contrasts; why to contrast expectations (`mu_rep`) rather than replicates
  (`y_rep`); the linear identity; `delta`; averaging with `weights`.
- `vignette("mep")` — a non-linear showcase: hierarchical sigmoid MEP
  input-output curves (dataset `mep_io`).
- `vignette("hurdle")` — hurdle outcomes with `outcome_spec()`.

The vignettes are precomputed (building them does not need CmdStan); their
runnable sources are `vignettes/*.Rmd.orig`.

## Citation

```r
citation("stanworkflow")
```

See also [`CITATION.cff`](CITATION.cff).

## Licence

GPL (>= 3). The `mep_io` dataset is from
[ShayOfir/MEP_IOcurve](https://github.com/ShayOfir/MEP_IOcurve)
(Wonderman Bar Sela et al. 2026, *PLOS ONE* 21(2): e0343725).
