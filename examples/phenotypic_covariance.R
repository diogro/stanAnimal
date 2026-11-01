# Phenotypic covariance example (no pedigree)
#
# Real mouse morphometric data (35 log-transformed traits), bundled with the
# package as ratones_morpho. lm_phen_cov() fits the phenotypic covariance
# matrix P after the fixed effects, so there is no relatedness and no G. Because
# the distances are strongly integrated, the traits are reparametrized
# with size_shape() and every draw is mapped back to the log-trait scale with
# back_transform(). Run from the project root:
#   Rscript examples/phenotypic_covariance.R

list_pkgs <- c("dplyr", "bayesplot", "posterior", "cmdstanr", "patchwork")
new_pkgs <- list_pkgs[!(list_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs) > 0){ install.packages(new_pkgs) }
invisible(lapply(sub("^.*/", "", list_pkgs), library, character.only = TRUE))


library(stanimal)

# ratones_morpho ships with the package: 329 mice, 35 cranial distances and the
# individual information (sex, age, generation, line). No network access needed.
ratonesdf <- ratones_morpho

# lm_phen_cov estimates the phenotypic covariance matrix P after fixed effects;
# there is no relatedness, so no G or breeding values.
all_traits = names(ratonesdf)[12:46]
traits = all_traits[1:35]  # subset to keep the example small (50s on full matrix)

# Size and shape reparametrization: the measurements share a single size
# gradient, so the model is fit to log size and orthogonal shape axes, and every
# draw is transformed back to the log-trait scale below. Keeping all k - 1 shape
# axes makes the reparametrization exact; a smaller n_shape gives a lower-rank,
# and much cheaper, model.
log_traits = log(as.matrix(ratonesdf[, traits]))
ss = size_shape(log_traits)
ratonesdf = cbind(ratonesdf, ss$scores)
cat(length(traits), "measured traits ->", ncol(ss$scores), "model traits\n")

# animal_design() builds Y (one column per model trait) and X from the data
# frame, so the model can be written as a formula. With many traits the cbind()
# response is assembled from the column names.
ratones_formula = stats::reformulate(c("SEX", "AGE", "LIN"),
                                     response = sprintf("cbind(%s)",
                                                        paste(ss$traits, collapse = ", ")))
ratones = animal_design(ratones_formula, data = ratonesdf)
Y = ratones$Y
X = ratones$X
cat(nrow(Y), "individuals,", ncol(Y), "model traits,", ncol(X), "fixed effects\n")

cov_fit = lm_phen_cov(Y, X, lkj_prior = 6, iter = 800,
                       warmup = 400, cores = 4,  chains = 4)

# Summaries and sampler diagnostics (model scale: size and shape axes)
cov_fit$summary(c("sigma_P[1]", "P[1,1]", "P[2,2]", "corrP[1,2]"))
cov_fit$diagnostic_summary()

# Fixed effects, named after the columns of X (model scale: one coefficient per
# size or shape trait)
beta_draws = rename_beta(cov_fit, "beta")
posterior::summarise_draws(beta_draws)

# Back-transform every draw to the log-trait scale. The result is a list of
# arrays indexed by trait first and draw last, as lmm_animal_gibbs() returns
# them: P and corrP are traits x traits x draws and sigma_P is traits x draws,
# all named after the traits
post = back_transform(cov_fit, ss)
post
k = length(traits)

round(apply(post$sigma_P, 1, mean), 3)
round(apply(post$corrP, c(1, 2), mean)[1:3, 1:3], 3)

# Posterior mean P on the log-trait scale, compared with the residual covariance
P_hat = apply(post$P, c(1, 2), mean)
round(P_hat, 3)
round(cov(residuals(lm(log_traits ~ X - 1))), 3)

library(evolqg)
MatrixCompare(P_hat, cov(residuals(lm(log_traits ~ X - 1))))


# Plots, all written to a single PDF
dir.create("examples/figures", recursive = TRUE, showWarnings = FALSE)
pdf("examples/figures/phenotypic_covariance_plots.pdf", width = 10, height = 7)
bayesplot::color_scheme_set("blue")
print(bayesplot::mcmc_trace(cov_fit$draws(c("P[1,1]", "P[2,2]", "P[3,3]"))))
print(bayesplot::mcmc_intervals(beta_draws))
print(bayesplot::mcmc_intervals(
  t(post$sigma_P[1:3, ]), prob = 0.8, prob_outer = 0.99, point_est = "mean"))
print(bayesplot::mcmc_intervals(
  cbind(`corrP[1,2]` = post$corrP[1, 2, ], `corrP[1,3]` = post$corrP[1, 3, ],
        `corrP[2,3]` = post$corrP[2, 3, ]),
  prob = 0.8, prob_outer = 0.99, point_est = "mean"))
print(plot_corr(cov2cor(P_hat), "Posterior mean phenotypic correlations") +
plot_corr(cov2cor(cov(residuals(lm(log_traits ~ X - 1)))), "Residual correlations"))
dev.off()
