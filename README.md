# stanAnimal

Bayesian animal models in Stan, with a focus on estimating the **additive genetic covariance matrix (G)** of one or more traits from a known relationship matrix.

Fits use [CmdStan](https://mc-stan.org/cmdstanr/) through `cmdstanr` (default) or `rstan`.

# Install

```r
if (!requireNamespace("pak", quietly = TRUE))
    install.packages("pak")
pak::pkg_install("stan-dev/cmdstanr")
cmdstanr::install_cmdstan()
install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
pak::pkg_install("diogro/stanAnimal/package")
```

# Usage

```r
library(stanAnimal)

# Y: trait matrix (individuals in rows), one record per individual
# X: fixed effects model matrix (full column rank, include an intercept column)
# A: relationship matrix, same individual order as Y and X
A <- as.matrix(nadiv::makeA(pedigree))
fit <- lmm_animal(Y, X, A)
```

`Y` can be a vector (univariate) or a matrix (multivariate). `A` must be symmetric and positive semidefinite. When `A` and the rows of `Y` are both named, `A` is put in the row order of `Y`; names that disagree are an error. Without names, the rows are assumed to be in the same order and only the dimensions are checked.

## Formula interface

`animal_design()` builds `Y` and `X` from a formula and a data frame, for when you would rather write the model the usual way than assemble the matrices by hand:

```r
design <- animal_design(cbind(t1, t2) ~ age + sex, data = dat)
fit <- lmm_animal(design$Y, design$X, A)
```

- the left-hand side is a single response or a `cbind()` of several traits, and the columns of `Y` are named after them;
- `X` is the fixed-effects design matrix, so its columns match the `beta` draws reported by the fit; factors are dummy coded and the intercept is kept unless the formula removes it;
- random-effect terms such as `(1 | animal)` are ignored, so `lme4`-style formulas can be reused unchanged — relatedness is supplied through `A`;
- rows with missing values are dropped (`na.action = stats::na.omit` by default), and the individuals kept, in row order, are returned in `attr(design, "individuals")`. Because `Y` keeps them as row names, a named `A` is reordered to match automatically:

```r
rownames(dat) <- dat$id
design <- animal_design(cbind(t1, t2) ~ age + sex, data = dat)
fit <- lmm_animal(design$Y, design$X, A)   # A is reordered to the rows of Y
```

With many traits, the `cbind()` response is easier to assemble from the column names:

```r
traits <- c("t1", "t2", "t3")
f <- reformulate(c("age", "sex"),
                 response = sprintf("cbind(%s)", paste(traits, collapse = ", ")))
design <- animal_design(f, data = dat)
```

The result feeds `lmm_animal()`, `lmm_multiCov()` and `lmm_animal_gibbs()` alike.

`animal_design()` can also build `A`, either from a pedigree (needs `nadiv`) or from a genotype matrix (through `genomic_relationship()`), as long as `data` has the individual IDs as row names. `A` is restricted to the rows kept in `Y`, so ancestors and individuals without records are handled:

```r
rownames(dat) <- dat$id

# from a pedigree (first three columns: id, sire, dam); needs nadiv
design <- animal_design(cbind(t1, t2) ~ age + sex, dat, pedigree = ped)
fit <- lmm_animal(design$Y, design$X, design$A)

# or from genotypes (individuals in rows, markers in columns)
design <- animal_design(cbind(t1, t2) ~ age + sex, dat, genotypes = M, blend = 0.01)
fit <- lmm_animal(design$Y, design$X, design$A)
```

## Genomic relationship matrix

If you have marker genotypes instead of a pedigree, `genomic_relationship()` builds a VanRaden (2008) GRM to use as `A`:

```r
# M: individuals in rows, markers in columns, coded 0/1/2 (NA allowed)
A <- genomic_relationship(M, blend = 0.01)
fit <- lmm_animal(Y, X, A)
```

Markers with low minor allele frequency (`maf`), too much missing data (`max_missing`) or no variation are dropped, and missing genotypes are set to their expected value. Known base-population allele frequencies can be passed with `p`. If there are fewer markers than individuals, or duplicated genotypes, the GRM is singular; a small `blend` mixes in the identity to make it positive definite.

## Output

With the default `cmd_stan = TRUE`, `lmm_animal()` returns a `CmdStanMCMC` object, so the usual cmdstanr, posterior and bayesplot tools work directly:

```r
fit$summary(c("G", "E", "h2"))
fit$draws("G", format = "draws_matrix")

bayesplot::mcmc_trace(fit$draws("h2"))
bayesplot::mcmc_intervals(fit$draws("corrG"))
bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(fit))
```

Variables are reported on the **original scale** of the data:

| Model | Variables |
|---|---|
| Multivariate | `G`, `E`, `P = G + E`, `corrG`, `h2`, `sigma_G`, `sigma_E`, `beta[trait, j]` |
| Univariate | `h2`, `sigma_G`, `sigma_E` (standard deviations), `beta[j]` |

Fixed effects are `beta[j]` (univariate) or `beta[trait, j]` for column `j` of `X`. To get draws named after the columns of `X` (for example `beta[2,age]`), use:

```r
animal_draws(fit, "beta")
```

With `cmd_stan = FALSE` the fit is a `stanfit` object from `rstan`, with `beta` already named after the columns of `X`.

## Multivariate model

The multivariate model estimates `G` and the residual covariance `E` directly. Each is a vector of scales times an LKJ-distributed correlation matrix, with half-normal priors on the scales (of the standardized responses). There is no prior on h2; `h2 = diag(G) / diag(P)` is derived.

- `lkj_prior` sets the LKJ shape for both correlation matrices (default 2; 1 is uniform, larger shrinks toward zero correlation). A shape above 1 is deliberately recommended for many traits: at high dimensionality, the uniform prior puts most of its mass near ±1 correlations and overestimates their magnitude, which is a well-known problem when estimating genetic correlations in quantitative genetics.
- `residual_structure = "unstructured"` (default) fits a full residual covariance. `"diagonal"` assumes independent residuals, which pushes environmental covariance between traits into `G` and inflates heritability, so it is not recommended when estimating `G`.
- Breeding values are marginalized using the eigendecomposition of `A`, so the model has no per-individual random effects. This is much faster and better mixing than sampling breeding values, and it only requires `A` to be positive semidefinite.
- The per-individual likelihoods are also whitened by the residual Cholesky factor and rotated into a single eigendecomposition of $E^{-1/2}GE^{-1/2}$, so an evaluation costs $O(NK^2 + K^3)$ rather than $O(NK^3)$. The model scales to dozens of traits (e.g. 40) without per-individual Cholesky factorizations. When the eigenvalues of $E^{-1/2}GE^{-1/2}$ are nearly coincident (at initialization, or when `G` is roughly proportional to `E`), the model automatically falls back to a per-individual Cholesky evaluation, whose gradient is well behaved.
- Fitting many traits is still computationally heavy even though the likelihood is cheap: an unstructured `G` at `K = 40` has roughly 1,700 parameters, and NUTS needs many leapfrog steps, so one fit can take an hour or more. Fewer `iter`, a larger `lkj_prior`, or `residual_structure = "diagonal"` make fits faster.

## Size and shape reparametrization

Morphometric traits are usually dominated by a single size gradient, so their covariance is a long, thin ridge that is slow to sample. `size_shape()` rewrites them as log size plus orthogonal shape axes, and `back_transform()` maps the posterior draws back to the trait scale:

```r
ss <- size_shape(log_traits)                  # individuals in rows, traits in columns
dat <- cbind(morphology, ss$scores)           # size, shape1, shape2, ...
design <- animal_design(cbind(size, shape1, shape2, shape3) ~ sex + site, dat)
fit <- lmm_animal(design$Y, design$X, A)

post <- back_transform(fit, ss)               # G, E, P, corrG and h2 on the log scale
posterior::summarise_draws(
  posterior::subset_draws(post, variable = c("G[1,1]", "corrG[1,2]", "h2[1]")))
```

Size is `rowMeans(Y)` and the shape axes are the principal components of `Y - size`, so `Y = scores %*% t(loadings)` up to the means and covariance matrices map exactly as $B M B^\top$. Keeping fewer axes (`n_shape`) gives a lower-rank, singular covariance in the original trait space. Heritabilities and genetic correlations are recomputed after the transformation, since neither is invariant under a change of traits; `beta` stays on the model scale. `back_transform()` also accepts a `lmm_multiCov()` fit, in which case it returns `P`, `corrP` and `sigma_P`.

## Phenotypic covariance only (P)

With no pedigree or genomic data, `lmm_multiCov()` estimates just the phenotypic covariance matrix `P`. It takes the same arguments and priors as `lmm_animal()` (without `A`), and returns the same kind of object, with `P`, `corrP`, `sigma_P` and `beta`. The likelihood is a single whitened triangular solve, $O(NK^2 + K^3)$, so it also scales to dozens of traits:

```r
fit <- lmm_multiCov(Y, X, lkj_prior = 2)
fit$summary("P")
```

## Univariate models

Two parametrizations are available through `parametrization`:

- `"variance"` (default): estimates the genetic and residual standard deviations with half-normal priors and derives h2.
- `"h2"`: estimates a total standard deviation and h2, with a `beta(beta_prior)` prior on h2. The default `beta_prior = c(1, 1)` is flat on h2.

If the `"variance"` version has chains with low E-BFMI, `lmm_animal()` warns and suggests switching to `parametrization = "h2"`. This tends to happen when the genetic variance is weakly identified, for example with few families.

```r
fit <- lmm_animal(y, X, A, parametrization = "h2")
```

## Scaling

All models standardize the response and divide continuous predictors by their standard deviation (constant and 0/1 columns are not scaled) for better sampling. Coefficients are back-transformed, so `beta` is in the original units of `X` and `Y`; the intercept is recovered when `X` has a constant column.

## Compilation cache

Models fit with cmdstanr are compiled once and the executable is reused in later R sessions. The cache is in `tools::R_user_dir("stanAnimal", "cache")` by default; change it with `options(stanAnimal.cache_dir = "path")` or the `STANANIMAL_CACHE_DIR` environment variable. Use `stan_cache_dir()` to see the location, and run `stanAnimal_clear_cache()` after upgrading CmdStan.

## Other options

- `cmd_stan = FALSE` uses `rstan` instead of `cmdstanr`.
- `blend_A = TRUE` replaces `A` with `0.99 * A + 0.01 * I`.
- `model_file` runs a custom Stan file in place of the default model.
- Further arguments (`...`) are passed to the sampler, for example `refresh = 0`.

## Gibbs sampler

`lmm_animal_gibbs()` is a quick single-chain Gibbs sampler for simulations and sanity checks. Its priors (inverse-Wishart) differ from those in `lmm_animal()`, so the results are not directly interchangeable.

# Example

```r
set.seed(1)
n_fam <- 40; fam_size <- 10; n <- n_fam * fam_size
family <- rep(seq_len(n_fam), each = fam_size)
A <- outer(family, family, "==") * 0.5   # half-sib families
diag(A) <- 1

G <- matrix(c(1.0, 0.5, 0.5, 1.5), 2)
E <- matrix(c(1.0, 0.3, 0.3, 1.0), 2)
X <- cbind(Intercept = 1, age = rnorm(n, 50, 10), sex = rbinom(n, 1, 0.5))
B <- cbind(c(10, 0.1, 1), c(-5, 0.05, -0.5))
Y <- X %*% B +
  t(chol(A)) %*% matrix(rnorm(n * 2), n) %*% chol(G) +
  matrix(rnorm(n * 2), n) %*% chol(E)
colnames(Y) <- c("t1", "t2")

fit <- lmm_animal(Y, X, A)
fit$summary(c("G", "corrG", "h2"))
animal_draws(fit, "beta")
```

# Notes

- The models assume one record per individual.
- Under development: the API may still change.
