# stanAnimal 0.4.0

## New features

* `ratones_morpho` and `ratones_pedigree` ship with the package: 35 cranial
  distances for 329 mice from the selection experiment of Penna et al. (2017),
  and the pedigree of the 525 individuals of the experiment. `data-raw/ratones.R`
  rebuilds both from the `diogro/ratones` package.
* `plot_G_ellipses()` gains an `X` argument. With `X`, the traits are replaced by
  their residuals after the fixed effects, so the scatter and the dashed
  reference ellipse have roughly the covariance `P = G + E` instead of including
  sex, age or group differences.
* `to_draws_list()` converts the arrays returned by `back_transform()` into one
  `posterior::draws_matrix` per quantity, for the `posterior` and `bayesplot`
  tools.
* `makeA()` builds the additive relationship matrix from a pedigree with the
  tabular method, accounting for inbreeding and accepting the pedigree in any
  order. `animal_design()` uses it when it is given `pedigree`.

## Changes

* `animal_draws()` is now `rename_beta()`, which says what it does: it takes the
  draws from a fit and renames the `beta` variables after the columns of `X`.
  The old name sounded like the way to get draws at all, which is `fit$draws()`.
  This is a rename, not a deprecation: `animal_draws()` is gone.

* `cmdstanr` and `rstan` are suggested rather than imported, so the package
  installs without a Stan engine. A fit uses `cmdstanr` when it is available,
  falls back to `rstan` with a warning when it is not, and stops with the
  installation command when neither is present. The fallback warning used to
  suggest a malformed `install_github()` call.

* The README and `?lmm_animal` now state the cost that grows with the number of
  individuals, which the trait-oriented cost notes did not mention: `A` is dense,
  so each chain decomposes it before sampling. The package targets many traits on
  hundreds to a few thousand individuals, and cohorts of tens of thousands are
  pointed to engines that use the sparse inverse of `A`, such as ASReml.

* The package no longer needs `nadiv`, which was only used to build `A` from a
  pedigree: `makeA()` replaces it, so nothing has to be installed for the
  pedigree path and it works on the same R versions as the rest of the package.
  `nadiv` stays in Suggests for one test, which checks `makeA()` against it.
* `back_transform()` now returns a list of class `size_shape_draws` instead of a
  single `draws_matrix`: `G`, `E`, `P`, `corrG` and `corrP` are
  trait by trait by draw arrays and `h2` is a trait by draw matrix, the same
  layout as `lmm_animal_gibbs()`, with named trait dimensions. Correlations are
  returned for the whole matrix rather than the upper triangle.
* `back_transform()` and `lmm_animal_gibbs()` are documented as sharing that
  layout.
* The `?size_shape` documentation explains what the reparametrization does to
  the LKJ prior, and that a larger `lkj_prior` is often reasonable for shape
  axes.
* `animal_prior_draws()`, `plot_prior_posterior()` and `plot_forest()` no
  longer take `unstructured_E`, which they accepted but never used. The
  residual correlation structure cannot change these draws: heritability uses
  only the diagonal of the residual covariance and the genetic correlations
  come from `G`. `?animal_prior_draws` says so.

## Bug fixes

* `plot_forest()` now checks `prob` and `plot_G_ellipses()` checks `prob`,
  `level` and `n_show`, so a level of 2 or a negative number of draws is a
  message naming the argument instead of a failure inside the interval or
  ellipse calculation. `n_show` must also be a whole number: it used to be
  truncated, so `n_show = 1.5` silently drew one ellipse.

* `animal_prior_draws()` now validates `lkj_prior`, `n_traits` and `n`. A shape
  of 0 returned correlations of exactly +/-1 (and a negative or non-finite
  shape returned `NaN`s with a warning), while a fractional `n_traits` or `n`
  was silently truncated and `0` failed with an unrelated error inside
  `matrix()`. `lmm_animal()` and `lmm_multiCov()` reject an invalid `lkj_prior`
  too, instead of losing every chain inside Stan.
* A count above `.Machine$integer.max` is rejected rather than coerced to `NA`
  with a warning, which used to reach the callers as an invalid matrix size.

* The `control` list now reaches both engines. With `cmdstanr`, only
  `adapt_delta` and `max_treedepth` were forwarded and every other entry was
  silently dropped; all documented sampler controls are now passed on, with
  `step_size` translated to `rstan`'s `stepsize`. Names that neither engine
  supports are an error instead of being ignored.

* The mean of the response is now restored for any design whose columns can
  represent a constant, not only when there is a literal intercept column. A
  formula such as `y ~ 0 + factor(g)` used to return coefficients on the
  centred response, missing the overall mean; they are now the group means.
  A design that cannot represent a constant (`y ~ 0 + x`) is rejected with a
  message, instead of silently returning coefficients that no longer match the
  data.

* A response trait that does not vary is now rejected with a clear message
  before the model is fitted. The univariate and multivariate Stan models
  standardize the responses by their standard deviation, so a constant trait
  used to fail inside Stan, after compiling and starting the chains.

* The size and shape example in the README built its response with a fixed
  number of shape axes, so `back_transform()` failed whenever the number of
  traits was not four.
* `back_transform()` called `stats::cov2cor()` through an unqualified name,
  which `R CMD check` reported.
* The package declared `R (>= 3.5.0)` but needs 4.0.0 (`tools::R_user_dir()`,
  and `cmdstanr`).

## Testing

* New tests cover the bundled datasets, the `X` option of `plot_G_ellipses()`
  and `to_draws_list()`, the sampler `control` pass-through of both engines, and
  the choice between `cmdstanr` and `rstan`.
* The Stan-specific branches of `animal_parameter_draws()`,
  `plot_prior_posterior()`, `plot_forest()` and `plot_G_ellipses()` are exercised
  on a real `CmdStanMCMC` fit, not only on synthetic draws, and the cache
  directory is tested through `STANANIMAL_CACHE_DIR` as well as the option.
* The code blocks of the README are now evaluated as a test, so documentation
  cannot drift from the API unnoticed (skipped outside a checkout and on CRAN).
