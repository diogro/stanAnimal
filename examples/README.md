# Examples

Runnable scripts for the models in `stanimal`. Run them from the project root,
because the figures are written to `examples/figures/`:

```sh
Rscript examples/univariate.R
```

| Script | Shows |
|---|---|
| `univariate.R` | `lmm_animal()` with a single trait, both parametrisations (`"variance"` and `"h2"`), compared with a simulated truth |
| `multivariate.R` | `lmm_animal()` with five correlated traits: `G`, `E`, `P`, `corrG`, `h2` and the fixed effects |
| `ratones.R` | `lmm_animal()` on real data: four log cranial distances of 329 mice, with `A` from the bundled pedigree (`ratones_morpho` and `ratones_pedigree`) |
| `phenotypic_covariance.R` | `lm_phen_cov()` with no pedigree, plus `size_shape()` and `back_transform()`, on the same mouse data |
| `gibbs.R` | `lmm_animal_gibbs()`, the single-chain Gibbs sampler; its inverse-Wishart priors differ from `lmm_animal()` |
| `multivariate_power_simulation.R` | Parallel simulation study of genetic-correlation recovery across sample sizes, writing summaries to `results/` |
| `package_comparison.R` | The same bivariate animal model fitted with `stanimal`, `MCMCglmm`, `sommer` and ASReml-R, estimated `G` and `E` compared with the simulated truth |

In the interval plots (`bayesplot::mcmc_recover_intervals()`), the filled circle
is the posterior mean, the thick and thin segments are the 80% and 99% credible
intervals, and the simulated value is the hollow triangle. The true values are
matched to the draws by name, not by position.

All scripts need `stanimal` installed and CmdStan available
(`cmdstanr::install_cmdstan()`), except `gibbs.R`, which needs neither. The
mouse data used by `ratones.R` and `phenotypic_covariance.R` ships with the
package, so nothing is downloaded. `package_comparison.R` also needs
`MCMCglmm`, `sommer` and a licensed ASReml-R.
