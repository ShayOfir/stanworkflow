# `infer_contrasts()` and its helpers in `R/infer.R`

## Context

This walks through `infer_contrasts()`: the stages of execution, the different "contrast"
options available, and the rationale behind the design choices. This is a Bayesian
post-processing layer that sits on top of a fitted Stan model — it turns raw posterior draws
into interpretable group contrasts (e.g. "is the treatment group's predicted outcome higher
than control?") with credible intervals and probability-of-direction statistics.

## 1. The two building-block S3 classes

### `outcome_spec` (what to extract from draws)

Decouples **what quantity to pull out of the posterior draws** from **how to compare groups**.
Rationale: for hurdle/zero-inflated models there isn't one "outcome" — there's `y_rep` (raw
mixture draws), `mu_marginal` (E[y], the usually-correct inferential target per Bürkner 2017 /
Woody et al. 2021), `p_respond` (extensive margin), and `mu_y` (intensive margin, E[y|respond]).
Rather than hard-coding `y_rep` everywhere, `outcome_spec()` (lines 24–41) bundles:
- `gq_var`: which generated-quantities variable to pull from `model$draws`
- `filter_fn`: post-extraction row filter (e.g. drop sentinel/placeholder rows)
- `transform_fn`: element-wise transform (e.g. `exp()` for log-scale models)
- `label`: for plots

Constructors: `os_yrep()` (default, backwards-compatible), `os_marginal()`, `os_extensive()`,
`os_intensive()`, and the compound `os_hurdle()` which bundles all three margins and is
special-cased inside `infer_contrasts()` to fan out into three separate contrast computations
in one call (lines 88–109).

### `contrast` (how to compare groups)

`contrast()` (lines 142–199) specifies the comparison itself, independent of outcome. Key
options and their rationale:
- `contrast_var` / `grouping_vars`: what varies (compared) vs. what's held out as facets (e.g.
  compare `treatment` within each `site`).
- `continuous`: categorical variables split by two named levels; continuous variables split by
  quantile (default 75th/25th pct) or by `absolute_thresholds` (explicit cutoffs) — this
  matters because a continuous contrast_var has no natural "hi"/"lo" groups, so the package
  offers both a robust default (quantile split, adapts to the data) and a precise option
  (fixed thresholds, for a specific clinically-meaningful cutoff).
- `effect_positive`: sign convention for what counts as "the effect direction," used later by
  `.define_scenarios()`.
- `filter_draws`: a quoted dplyr expression to subset draws before contrasting (e.g. restrict
  to one diagnostic group) — kept as a quosure/expression rather than evaluated eagerly so it
  can be applied after the (potentially expensive) draw-extraction step.
- `two_way` / `two_way_levels`: turns a single hi/lo contrast into a **joint 4-scenario
  classification** across two outcome levels of a second variable (e.g. "does the treatment
  help direction A, direction B, both, or neither?"). This is the more elaborate contrast mode,
  explained in section 4.
- `outcome`: an `outcome_spec` object; defaults to `os_yrep()` if unset.

Recent history visible in the comments: `two_way` used to be blocked when `continuous = TRUE`;
that guard was deliberately removed once `.apply_contrast_split()` was generalized to handle
continuous two-way splits uniformly (lines 172–174).

## 2. Draw extraction — two fundamentally different strategies

`infer_contrasts()` supports `draws_method = "descriptive"` (default) or `"counterfactual"`,
implemented by `.get_yrep_draws()` (lines 775–984):

- **Descriptive**: uses the posterior draws already implied by the *observed* data — i.e.
  extracts `y_rep` (or whichever `gq_var`) for the fitted model as-is, groups rows by the
  observed levels/values of `contrast_var`, and compares. Cheap (one draw extraction), but
  the "hi" and "lo" groups differ on every other covariate too, so this is descriptive/
  correlational, not causal.
- **Counterfactual**: holds every other covariate fixed and re-runs `generate_quantities()`
  (via `.simulate_at()`, lines 604–673) at each point of a grid over `contrast_var`, then
  differences each grid point against a baseline (the minimum grid value). This isolates the
  causal effect of `contrast_var` alone, at the cost of running Stan's GQ block once per grid
  point (`gq_steps`, default 5). `.simulate_at()` rebuilds `stan_data` for a synthetic dataset
  with `contrast_var` overwritten, then calls `.restore_training_dimensions()` (lines 289–326)
  to re-map any `*_id` / `*_id_cond_array` index columns back onto the *original* factor
  levels used at fit time — otherwise a synthetic re-encoding could silently shift factor
  codes (e.g. "coil_id = 1" no longer meaning "fig8").

Within counterfactual mode, `.get_yrep_draws()` branches again:
- Continuous `contrast_var`: builds a numeric grid via `seq()` across the observed range and
  keeps **every** grid point (needed downstream for the "full sweep" plots/summaries).
  Categorical `contrast_var`: only simulates at the levels requested in `contrast_levels`,
  then explicitly differences index-level vs. baseline-level draws via a `left_join` keyed on
  `.chain/.iteration/.draw/.obs_idx/feature cols` (lines 940–974) — this differencing must
  happen at the *individual draw* level (not on summaries) to correctly propagate posterior
  uncertainty into the contrast.

`.extract_spec_draws()` (lines 685–748) is the shared low-level worker that pivots a wide
`draws_df` (one column per observation) into a long tibble (`.chain, .iteration, .draw,
.obs_idx, <features>, y`), applies `transform_fn`, and applies `filter_fn` — used by the
descriptive path directly, and conceptually mirrored by `.to_long_draws()` inside the
counterfactual path.

## 3. `infer_contrasts()` orchestration (lines 1157–1355)

Stages, in order:

1. **Normalize input**: wrap a single `contrast` into a list; default `data`/`extra` to the
   model's own stored training data.
2. **Register `numeric_array`s**: lets a contrast reference a column of a matrix-valued
   predictor (e.g. an ROI matrix) that isn't a literal column in `data`.
3. **Expand numeric-array contrasts**: `.expand_numeric_array_contrast()` (lines 1072–1090)
   turns one contrast referencing an array into N contrasts (one per column), injecting a
   virtual column like `"roi[CST_L]"` into `data` on the fly so the rest of the pipeline can
   treat it like any other variable.
4. **Expand hurdle contrasts**: any contrast using `os_hurdle()` is fanned out into 3 sub-
   contrasts (marginal/extensive/intensive), tracked via `hurdle_map` so results can be
   reassembled under one key later.
5. **Collect needed variables**: unions `grouping_vars`, `contrast_var`, `two_way`, and any
   variables referenced inside `filter_draws` expressions (via `all_vars_in_expr()`, a small
   recursive walker over the quoted call) — this is what gets passed to `.get_yrep_draws()` as
   `vars` so the long tibble carries exactly the feature columns needed for grouping/filtering.
6. **Group by `(gq_var, contrast_var)`**: contrasts sharing the same outcome variable *and*
   the same swept variable can share a single, potentially-expensive `.get_yrep_draws()` call
   (this batches numeric-array-expanded contrasts and repeated-variable contrasts without
   re-running GQ simulation redundantly).
7. **Per group**: call `.get_yrep_draws()` once, then dispatch each contrast in the group to
   one of three computation functions, based on `draws_method` and `co$two_way`:
   - `draws_method == "counterfactual"` → `.compute_counterfactual_contrast()`
   - `!is.null(co$two_way)` → `.compute_two_way_contrast()`
   - else → `.compute_one_contrast()`
8. **Reassemble**: flatten results back into the user's original contrast ordering, collapsing
   hurdle sub-results into a named list (`marginal`/`extensive`/`intensive`) under one output
   slot.
9. **Wrap**: single contrast → one `inference` S3 object (`as_inference()`); multiple → a list
   of `inference` objects (or list-of-margins for hurdle contrasts).

## 4. The three contrast computation functions

All three share this shape: apply `filter_draws` → summarise `y` per draw per cell (via
`summary_fn`, default `mean`) → pivot wide → difference → summarize with credible intervals.

- **`.compute_one_contrast()`** (lines 1729–1822): simplest case — one hi/lo split (categorical
  levels via `.resolve_contrast_levels()`, or continuous via quantile/absolute threshold),
  pivot to wide `hi`/`lo` columns, `contrast = hi - lo` per draw, then summarize
  (`mean`, `sd`, `rethinking::HPDI`, `P_positive`, `P_negative`) per `grouping_vars` cell.
  Returns `list(summary, posterior)`.

- **`.compute_two_way_contrast()`** (lines 1828–1934) and its continuous counterfactual sibling
  inside `.compute_counterfactual_contrast()`: answers a richer question — "does the hi/lo
  contrast move *two different outcomes* (levels of `two_way`) in the effect-positive
  direction, and if so which one(s)?" It:
  1. Resolves the two `two_way` levels (`lv1`, `lv2`).
  2. Applies the hi/lo split via the shared `.apply_contrast_split()` (lines 1016–1064),
     which is the same continuous-quantile/absolute-threshold-or-categorical logic factored
     out of `.compute_one_contrast()` so both paths stay consistent.
  3. Pivots to one contrast column per outcome level (`contrast_<lv1>`, `contrast_<lv2>`).
  4. Classifies each posterior draw into one of 4 mutually-exclusive scenarios via
     `.define_scenarios()` (lines 1358–1381): "Affects both", "Affects `lv1` only",
     "Affects `lv2` only", "Affects neither" — sign-per-column check driven by
     `effect_positive`.
  5. Aggregates draw counts into posterior scenario probabilities, optionally faceted by
     `grouping_vars` (minus `two_way` itself, since that's now the thing being crossed).
  Rationale for this design: rather than reporting two separate marginal contrasts (which
  can't tell you about joint/conditional effects), this gives a genuine joint posterior
  answer — e.g. "70% probability the treatment helps direction A only, not B."

- **`.compute_counterfactual_contrast()`** (lines 1410–1655): the counterfactual-draws
  analogue of the above two, but keyed on the full `contrast_var_value` grid rather than a
  single hi/lo pair (since `.get_yrep_draws()` already returns pre-differenced-vs-baseline `y`
  for continuous sweeps). Branches on `co$two_way` (scenario grid across the sweep, reusing
  `.define_scenarios()`) vs. one-way, and again on `co$continuous` (full grid summary,
  `grid_summary`/`grid_posterior`) vs. categorical (pre-differenced pair, `posterior`/
  `summary` — same shape as `.compute_one_contrast()`'s output, for interface consistency).

## 5. Result / plotting layer

`as_inference()` just stamps an `"inference"` class onto the result list for dispatch.
`print.inference()` (lines 347–405) inspects which fields are present (`two_way` set?
`grid_summary`/`scenario_grid` present?) to decide which of the three plot helpers to call —
`.plot_inference_1way()` (density of contrast, faceted), `.plot_inference_2way()` (scatter of
the two per-outcome contrast columns colored by scenario), or `.plot_inference_grid()` (line
plot across the `contrast_var` grid, either scenario probabilities or contrast mean+ribbon).
This mirrors the branching logic inside `infer_contrasts()`/`.compute_*` almost exactly, so the
print/plot dispatch and the computation dispatch stay in lockstep by construction.
