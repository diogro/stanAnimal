# Gibbs sampler example
#
# lmm_animal_gibbs() is a quick single-chain Gibbs sampler, intended for
# simulations and sanity checks rather than for inference. Its priors are
# inverse-Wishart on G and R, which differ from the priors used by
# lmm_animal(), so the two posteriors are close but not interchangeable.
#
# This script simulates a bivariate trait on half-sib families, fits the Gibbs
# sampler, and compares the posterior with the simulated truth. Run from the
# project root:
#   Rscript examples/gibbs.R

# Install (if needed) and load all dependencies in one go
pkgs <- c("posterior", "bayesplot")
if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install(pkgs, upgrade = FALSE)
invisible(lapply(pkgs, library, character.only = TRUE))

library(stanAnimal)

## Simulation
# Half-sib families give a block relationship matrix without needing a package.
set.seed(1)
n_fam <- 40
fam_size <- 8
n <- n_fam * fam_size
family <- rep(seq_len(n_fam), each = fam_size)

A <- outer(family, family, "==") * 0.5
diag(A) <- 1

G <- matrix(c(0.60, 0.30,
              0.30, 0.50), 2, byrow = TRUE)
E <- matrix(c(0.80, 0.20,
              0.20, 0.70), 2, byrow = TRUE)

X <- cbind(Intercept = 1, sex = rbinom(n, 1, 0.5))
B <- cbind(c(1.0, 0.5),     # intercepts
           c(-1.0, -0.2))   # sex effects

Y <- X %*% B +
  t(chol(A)) %*% matrix(rnorm(n * 2), n) %*% chol(G) +
  matrix(rnorm(n * 2), n) %*% chol(E)
colnames(Y) <- c("t1", "t2")

h2_true <- diag(G) / diag(G + E)
round(G, 2)
round(E, 2)
round(h2_true, 2)

## Model fit
# iter/warmup/thin are counts of Gibbs iterations, not Hamiltonian trajectories,
# so a few thousand retained draws already give a usable sample. Successive
# draws are autocorrelated, so the effective sample size is much smaller than
# the number of draws; the example keeps 5000 of 20000 iterations and runs in
# a few seconds.
gibbs <- lmm_animal_gibbs(Y, X, A, iter = 20000, warmup = 5000, thin = 3)

# Draws are plain arrays: G and R are k x k x draws, beta is
# traits x fixed effects x draws, h2 is traits x draws.
str(gibbs, max.level = 1)
length(gibbs$iterations)
dim(gibbs$G)
dim(gibbs$beta)
dim(gibbs$h2)

## Posterior means compared with the truth
G_hat <- apply(gibbs$G, 1:2, mean)
R_hat <- apply(gibbs$R, 1:2, mean)   # residual covariance, the analogue of E
round(G_hat, 2)
round(G, 2)
round(R_hat, 2)
round(E, 2)

# Fixed effects, original scale: rows are traits, columns the columns of X
beta_hat <- apply(gibbs$beta, 1:2, mean)
round(beta_hat, 2)
t(B)

round(rowMeans(gibbs$h2), 2)
round(h2_true, 2)

## Summaries as posterior draws
# Flatten the arrays into a draws matrix so posterior and bayesplot apply.
gibbs_draws <- posterior::as_draws_matrix(cbind(
  G11 = gibbs$G[1, 1, ], G12 = gibbs$G[1, 2, ], G22 = gibbs$G[2, 2, ],
  R11 = gibbs$R[1, 1, ], R12 = gibbs$R[1, 2, ], R22 = gibbs$R[2, 2, ],
  `h2[t1]` = gibbs$h2[1, ], `h2[t2]` = gibbs$h2[2, ]))

posterior::summarise_draws(gibbs_draws)

## Plots
dir.create("examples/figures", recursive = TRUE, showWarnings = FALSE)
pdf("examples/figures/gibbs_plots.pdf", width = 10, height = 7)
bayesplot::color_scheme_set("blue")

# Single chain, so mcmc_trace shows one line per parameter and no Rhat
print(bayesplot::mcmc_trace(gibbs_draws,
                            pars = c("G11", "G12", "G22", "R11", "R22")))

# Posterior intervals with the simulated truth as a reference line
truth <- c(
  G11 = G[1, 1], G12 = G[1, 2], G22 = G[2, 2],
  R11 = E[1, 1], R12 = E[1, 2], R22 = E[2, 2],
  `h2[t1]` = h2_true[1], `h2[t2]` = h2_true[2])
print(bayesplot::mcmc_intervals(gibbs_draws, prob = 0.8, prob_outer = 0.99,
                                point_est = "mean") +
        ggplot2::geom_vline(xintercept = truth[posterior::variables(gibbs_draws)]))
dev.off()
