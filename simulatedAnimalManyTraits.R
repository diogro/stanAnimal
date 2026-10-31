# Install (if needed) and load all dependencies in one go
pkgs <- c("mvtnorm", "rstan", "nadiv", "evolqg", "corrplot", "patchwork",
          "cran/pedtools")  # pedtools is no longer on CRAN -> install from GitHub
if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install(pkgs, upgrade = FALSE)
invisible(lapply(sub("^.*/", "", pkgs), library, character.only = TRUE))

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

library(stanAnimal)

lt = function(x, diag = T) x[lower.tri(x, diag = diag)]

ped = randomPed(1500, 50, seed = 2)
ped = data.frame(ID = ped$ID, sire = ped$FIDX, dam = ped$MIDX)
ped = prepPed(ped)
A <- as.matrix(nadiv::makeA(ped))
nrow(A)
L_A = chol(A)

n_traits = 5

set.seed(41)
(G <- RandomMatrix(n_traits, LKJ = FALSE))
E = diag(n_traits)

n = nrow(A)
p = nrow(G)
a = t(L_A) %*% matrix(rnorm(n*p), n, p) %*% chol(G)
cor(a)

beta = matrix(c(rep(1, n_traits),
                rnorm(n_traits, 0, 1),
                rnorm(n_traits, 1, 1)), 3, n_traits, byrow = TRUE)
rownames(beta) = c("Intercept", "sex", "Z")

sex = sample(c(0, 1), nrow(a), replace = TRUE)
Z = rnorm(nrow(a))
Intercept = rep(1, nrow(a))
X = cbind(Intercept, sex, Z)
rownames(X) = rownames(a)
e = rmvnorm(nrow(a), sigma = E)

Y = X %*% beta + a + e
   
stan_model = lmm_animal(Y, X, A, lkj_prior = 2, 
                        iter = 400, warmup = 200,
                        cores = 4, chains = 4, residual_structure = "diagonal")

stan_model$summary(c("h2", "G", "E", "corrG"))
stan_model$summary(c("h2", "G", "E", "corrG"),
                   posterior::default_convergence_measures())

## Sampler diagnostics (divergences, treedepth, E-BFMI)
stan_model$diagnostic_summary()

# Posterior means as matrices, compared with the simulated truth
post_mean = function(name, k = n_traits) {
  s = stan_model$summary(name, "mean")
  matrix(s$mean[match(sprintf("%s[%d,%d]", name, rep(1:k, k), rep(1:k, each = k)),
                      s$variable)], k)
}
G_hat = post_mean("G")
E_hat = post_mean("E")
round(G_hat, 2); round(G, 2)
round(E_hat, 2); round(E, 2)
round(cov2cor(post_mean("P")), 2); round(cov2cor(G + E), 2)

# Data with posterior draws of the G ellipses and genetic correlations
print(plot_G_ellipses(stan_model, Y, n_show = 100))

par(mfrow = c(1, 2))
print(corrplot.mixed(post_mean("corrG"), upper = "ellipse"), main = "Estimated Genetic Correlations")
print(corrplot.mixed(cov2cor(G), upper = "ellipse"), main = "True Genetic Correlations")
par(mfrow = c(1, 1))

plot_corr(post_mean("corrG"), "Estimated") +
plot_corr(cov2cor(G), "True")

# Posterior distribution of the eigenvalue variance of G
G_draws <- as.array(stan_model$draws("G"))  # iterations x chains x traits x traits
eig_var <- apply(G_draws, c(1, 2), function(d) CalcEigenVar(matrix(d, n_traits)))
eig_var_true <- CalcEigenVar(G)

ggplot2::ggplot(data.frame(value = as.vector(eig_var)), ggplot2::aes(value)) +
  ggplot2::geom_density(fill = "grey70", colour = NA) +
  ggplot2::geom_vline(xintercept = eig_var_true, colour = "red", linetype = "dashed") +
  ggplot2::labs(x = "Eigenvalue variance of G", y = "Density",
                title = "Posterior eigenvalue variance of G (dashed line = true value)")
