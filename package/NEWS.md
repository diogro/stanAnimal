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

## Bug fixes

* The size and shape example in the README built its response with a fixed
  number of shape axes, so `back_transform()` failed whenever the number of
  traits was not four.
* `back_transform()` called `stats::cov2cor()` through an unqualified name,
  which `R CMD check` reported.
* The package declared `R (>= 3.5.0)` but needs 4.0.0 (`tools::R_user_dir()`,
  and `cmdstanr`).

## Testing

* New tests cover the bundled datasets, the `X` option of `plot_G_ellipses()`
  and `to_draws_list()`.
* The code blocks of the README are now evaluated as a test, so documentation
  cannot drift from the API unnoticed (skipped outside a checkout and on CRAN).
