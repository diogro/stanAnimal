# Animal model on real data: mice from a long-term selection experiment
#
# ratones_morpho ships with the package (Penna et al. 2017): 35 cranial
# distances for 329 mice from five artificial selection lines, with sex, age,
# generation and line, plus ratones_pedigree, the pedigree of the 525
# individuals in the experiment. Four distances are fitted here, with the
# relationship matrix built from the pedigree by animal_design().
#
# Run from the project root:
#   Rscript examples/ratones.R

# Install (if needed) and load all dependencies in one go
pkgs <- c("rstan", "bayesplot", "posterior", "ggplot2", "patchwork", "evolqg")
if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install(pkgs, upgrade = FALSE)
invisible(lapply(pkgs, library, character.only = TRUE))

# cmdstanr is not on CRAN; install it only if missing
if(!requireNamespace("cmdstanr", quietly = TRUE)) pak::pkg_install("stan-dev/cmdstanr")
library(cmdstanr)

options(mc.cores = 4)

library(stanAnimal)

## Data
data(ratones_morpho)
data(ratones_pedigree)
dat <- ratones_morpho
cat(nrow(dat), "mice,", ncol(dat) - 11, "cranial distances\n")
table(dat$selection)

# Four distances spanning the skull, on the log scale as usual for
# morphometrics. The 35 distances are moderately integrated (one principal
# component explains about 40% of the variance, mean r = 0.35), so a handful of
# them already carry a correlated genetic signal.
traits <- c("IS_PM", "NSL_NA", "PT_BA", "BA_OPI")
dat[traits] <- log(dat[traits])
round(cor(dat[traits]), 2)

# animal_design() builds Y, X and A in one call: the formula gives the fixed
# effects (sex, age at death and selection line -- the lines differ in body weight
# by construction, so their means must be modelled) and pedigree =
# ratones_pedigree adds the relationship matrix. The rows of ratones_morpho are
# named with the individual IDs, so the pedigree -- which also contains the
# unmeasured ancestors -- is matched to the rows of Y automatically.
design <- animal_design(cbind(IS_PM, NSL_NA, PT_BA, BA_OPI) ~ SEX + AGE + LIN,
                        dat, pedigree = ratones_pedigree)
cat(nrow(design$Y), "mice with records,", ncol(design$Y), "traits,",
    ncol(design$X), "fixed effects\n")
cat("pedigree:", nrow(ratones_pedigree), "individuals,",
    nrow(ratones_pedigree) - nrow(design$Y), "of them ancestors without records\n")
dim(design$A)
range(design$A)

## Model fit
time_traits <- system.time(
  fit <- lmm_animal(design$Y, design$X, design$A, lkj_prior = 2,
                    iter = 3000, warmup = 1500, chains = 4, cores = 4))

# Heritability, genetic and residual standard deviations, and genetic
# correlations, with sampler diagnostics
fit$summary(c("h2", "sigma_G", "sigma_E", "corrG"))
fit$diagnostic_summary()

# Fixed effects, named after the columns of X: beta[trait, term]
beta_draws <- rename_beta(fit, "beta")
posterior::summarise_draws(beta_draws)

# Posterior means of G and E. The traits are log distances, so variances are
# small numbers: h2 and the correlations are scale free and easier to read.
post_mean <- function(name, k = ncol(design$Y)) {
  s <- fit$summary(name, "mean")
  matrix(s$mean[match(sprintf("%s[%d,%d]", name, rep(seq_len(k), k),
                              rep(seq_len(k), each = k)), s$variable)], k,
         dimnames = list(traits, traits))
}
signif(post_mean("G"), 2)
signif(post_mean("E"), 2)
round(post_mean("corrG"), 2)
# P = G + E was estimated after the fixed effects, so it is smaller than the raw
# covariance of the traits, which also contains the sex, age and line effects
head(signif(post_mean("P"), 2))
head(signif(cov(residuals(lm(cbind(IS_PM, NSL_NA, PT_BA, BA_OPI) ~ SEX + AGE + LIN,
                        dat))), 2))

## Plots
dir.create("examples/figures", recursive = TRUE, showWarnings = FALSE)
pdf("examples/figures/ratones_plots.pdf", width = 10, height = 7)
bayesplot::color_scheme_set("blue")

print(bayesplot::mcmc_trace(fit$draws(c("h2", "corrG"))))
print(bayesplot::mcmc_intervals(beta_draws, prob = 0.8, prob_outer = 0.99))

# Where the prior and the posterior overlap for h2 and corrG, the data add
# little information
print(plot_prior_posterior(fit, trait_names = traits, ncol = 4))
print(plot_forest(fit, trait_names = traits))

# Residuals after the fixed effects, with posterior draws of the G ellipses and
# the genetic correlations. X = design$X removes the sex, age and line effects,
# so the cloud and the dashed ellipse have roughly the covariance P = G + E.
print(plot_G_ellipses(fit, design$Y, trait_names = traits, n_show = 100,
                      X = design$X))

# Estimated genetic correlations next to the observed phenotypic correlations
print(plot_corr(post_mean("corrG"), "Posterior mean corrG") +
        plot_corr(cor(dat[traits]), "Observed P correlations"))

dev.off()

## Extra: the same model in size and shape coordinates
# Size plus the three shape axes span the same space as the four traits, so this
# is the same model in another parametrisation: same number of traits, same X
# and A. The prior is not invariant, though, and here that helps: the LKJ prior
# acts on correlations between size and shape rather than between the raw
# traits, and the shape axes are principal components, so their correlations are
# expected to be small. lkj_prior stays at 2 to keep the two fits comparable;
# for a size and shape fit on its own, a larger value is worth trying
# (see ?size_shape).
ss <- size_shape(dat[traits], n_shape = length(traits) - 1L)
dat <- cbind(dat, ss$scores)
design_ss <- animal_design(cbind(size, shape1, shape2, shape3) ~ SEX + AGE + LIN,
                           dat, pedigree = ratones_pedigree)

time_size_shape <- system.time(
  fit_ss <- lmm_animal(design_ss$Y, design_ss$X, design_ss$A, lkj_prior = 2,
                       iter = 3000, warmup = 1500, chains = 4, cores = 4))

# Running time. Both fits use the same compiled model, so the sampling column is
# the direct comparison; elapsed also covers R start-up and, on a cold cache,
# compilation. A single run is noisy, so repeat it for a reliable comparison.
rbind(traits = c(elapsed = time_traits[["elapsed"]],
                 sampling = fit$time()$total),
      size_shape = c(elapsed = time_size_shape[["elapsed"]],
                     sampling = fit_ss$time()$total))

# G on the measured trait scale, from both fits. back_transform() returns G as a
# traits by traits by draws array, already named after the traits, the same
# layout lmm_animal_gibbs() uses
post_ss <- back_transform(fit_ss, ss)
G_traits <- post_mean("G")
G_size_shape <- apply(post_ss$G, c(1, 2), mean)
round(G_traits, 4)
round(G_size_shape, 4)
# Largest difference between the two estimates, relative to the largest element
# of G: Monte Carlo error, plus the prior difference described above
round(max(abs(G_traits - G_size_shape)) / max(abs(G_traits)), 3)

evolqg::MatrixCompare(G_traits, G_size_shape)
