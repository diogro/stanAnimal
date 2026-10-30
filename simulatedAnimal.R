if(!require(MCMCglmm)){install.packages("MCMCglmm"); library(MCMCglmm)}
if(!require(mvtnorm)){install.packages("mvtnorm"); library(mvtnorm)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(mvtnorm)){install.packages("mvtnorm"); library(mvtnorm)}
if(!require(nadiv)){install.packages("nadiv"); library(nadiv)}
if(!require(tictoc)){install.packages("tictoc"); library(tictoc)}
if(!require(pedtools)){pak::pkg_install("cran/pedtools"); library(pedtools)}

pak::pkg_install(c("tidyverse", "devtools",  "roxygen2"))
library(devtools)
rstan_options(auto_write = TRUE)
options(mc.cores = 12)

library(stanAnimal)
"cmdstanr" %in% installed.packages()[,"Package"]
#cmdstanr::install_cmdstan(cores = 2)
library(AtchleyMice)

#ped <- read.table("volesPED.txt",header=T)
ped = randomPed(100, 10, seed = 2)
plot(ped)
dev.off()
ped = data.frame(ID = ped$ID, sire = ped$FIDX, dam = ped$MIDX)
ped = prepPed(ped)
A <- as.matrix(nadiv::makeA(ped))
#ped <- mice_pedigree
#F6_ID = mice_info$F6$ID
#A <- as.matrix(nadiv::makeA(ped))[F6_ID, F6_ID]
L_A = chol(A)

corrG <- matrix(c(1, 0.7, 0.2,
                  0.7, 1, 0.0,
                  0.2, 0,   1), 3, 3, byrow = T)
corrE <- matrix(c(  1, 0.2, 0.0,
                  0.2,   1, 0.0,
                  0.0, 0.0,  1), 3, 3, byrow = T)
varG = rep(1, 3)
varE = 2*varG
G = sqrt(varG) %*% t(sqrt(varG)) * corrG
E = sqrt(varE) %*% t(sqrt(varE)) * corrE
varG / (varG + varE)
n = nrow(A)
p = nrow(G)
a = t(L_A) %*% matrix(rnorm(n*p), n, p) %*% chol(G)


beta = matrix(c(1, 2, 3,
                0.1, 0.2, 0.5,
                0.05, 0.1, 0.3), 3, 3, byrow = TRUE)
colnames(beta) = c("x", "y", "z")
rownames(beta) = c("Intercept", "sex", "Z")


#sex = as.numeric(as.factor(mice_info$F6$Sex))-1
sex = sample(c(0,1), nrow(a), replace = TRUE)
sex = sex[1:nrow(a)]
Z = rnorm(nrow(a))
Intercept = rep(1, nrow(a))
X = cbind(Intercept, sex, Z)
rownames(X) = rownames(a)
e = rmvnorm(nrow(a), sigma = E)

Y = X %*% beta + a + e

colnames(Y) = c("x", "y", "z")

sim_data = data.frame(Y, X, animal = rownames(a))

# prior_bi <- list(G = list(G1 = list(V = diag(3), nu = 1.002)),
#                           R = list(V = diag(3), nu = 1.002))
# tic()
# model_bi <- MCMCglmm(cbind(x, y, z) ~ trait + trait:sex + trait:Z,
#                      random = ~us(trait):animal,
#                      rcov = ~us(trait):units, family = c("gaussian", "gaussian", "gaussian"),
#                      pedigree = ped, data = sim_data, prior = prior_bi,
#                      nitt = 130000, thin = 100, burnin = 30000, verbose = TRUE)
# toc()
# summary(model_bi)
# colMeans(model_bi$VCV[,c("traitx:traitx.animal", "traity:traity.animal", "traitz:traitz.animal")])
# G_mcmc = matrix(colMeans(model_bi$VCV[, grep("animal", colnames(model_bi$VCV))]), 3, 3)
# corrG_mcmc = cov2cor(G_mcmc)


stan_model = lmm_animal(Y, X, A, chains = 4, iter = 2000, warmup = 1000,
                        cores = 4, lkj_prior = 2)

# The fit is a CmdStanMCMC object: summaries and diagnostics use its methods,
# and draws work with posterior and bayesplot.

## Summaries
stan_model$summary(c("h2", "G", "E", "corrG"))
stan_model$summary(c("h2", "G", "E", "corrG"),
                   posterior::default_convergence_measures())

## Sampler diagnostics (divergences, treedepth, E-BFMI)
stan_model$diagnostic_summary()

# Fixed effects are named after the columns of X: beta[trait, column]
beta_draws = animal_draws(stan_model, "beta")
posterior::summarise_draws(beta_draws)

# Residual correlations are not a model output; derive them from the draws
draws_df = posterior::as_draws_df(stan_model$draws(c("G", "E", "corrG", "h2")))
draws_df = posterior::mutate_variables(
  draws_df,
  `corrE[1,2]` = `E[1,2]` / sqrt(`E[1,1]` * `E[2,2]`),
  `corrE[1,3]` = `E[1,3]` / sqrt(`E[1,1]` * `E[3,3]`),
  `corrE[2,3]` = `E[2,3]` / sqrt(`E[2,2]` * `E[3,3]`))
posterior::summarise_draws(posterior::subset_draws(draws_df, variable = "corrE"))

# Posterior means as matrices, compared with the simulated truth
post_mean = function(name, k = 3) {
  s = stan_model$summary(name, "mean")
  matrix(s$mean[match(sprintf("%s[%d,%d]", name, rep(1:k, k), rep(1:k, each = k)),
                      s$variable)], k)
}
G_hat = post_mean("G")
E_hat = post_mean("E")
round(G_hat, 2); round(G, 2)
round(E_hat, 2); round(E, 2)
round(cov2cor(post_mean("P")), 2)

# Fixed effects, original scale: rows are traits, as in t(beta)
beta_hat = stan_model$summary("beta", "mean")$mean
round(matrix(beta_hat, 3), 2)
t(beta)

## Visualization
bayesplot::color_scheme_set("blue")

# Traces and densities
bayesplot::mcmc_trace(stan_model$draws("G"))
bayesplot::mcmc_trace(stan_model$draws("E"))
bayesplot::mcmc_dens_overlay(stan_model$draws("h2"))

# Convergence and sampler checks
all_draws = stan_model$draws()
bayesplot::mcmc_rhat(posterior::summarise_draws(all_draws, "rhat")$rhat)
bayesplot::mcmc_neff(posterior::summarise_draws(all_draws, "ess_bulk")$ess_bulk /
                       posterior::ndraws(all_draws))
bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(stan_model))

# Fixed effects
bayesplot::mcmc_intervals(beta_draws)

# Intervals with the simulated truth as a reference line
interval_plot = function(pars, truth) {
  bayesplot::mcmc_intervals(
    stan_model$draws(), pars = pars,
    prob = 0.8, prob_outer = 0.99, point_est = "mean") +
    ggplot2::geom_vline(xintercept = truth)
}
interval_plot(c("h2[1]", "h2[2]", "h2[3]"), diag(G) / (diag(G) + diag(E)))
interval_plot(c("G[1,1]", "G[2,2]", "G[3,3]"), diag(G))
interval_plot(c("E[1,1]", "E[2,2]", "E[3,3]"), diag(E))
interval_plot(c("G[1,2]", "G[1,3]", "G[2,3]"), G[lower.tri(G)])
interval_plot(c("corrG[1,2]", "corrG[1,3]", "corrG[2,3]"), corrG[lower.tri(corrG)])

# Estimated vs simulated genetic correlations
library(corrplot)
par(mfrow = c(1, 3))
corrplot.mixed(cor(a), upper = "ellipse")      # realized breeding values
corrplot.mixed(post_mean("corrG"), upper = "ellipse")
corrplot.mixed(corrG, upper = "ellipse")       # truth
par(mfrow = c(1, 1))

