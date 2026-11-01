# One simulated data set, four engines for a bivariate animal model
#
# Simulates two genetically correlated traits with a known G and E on a
# five-generation pedigree, then fits the same bivariate animal model with
# stanimal, MCMCglmm, sommer and ASReml-R, and compares the estimates with the
# simulated truth and the time each fit took. Run from the project root:
#   Rscript examples/package_comparison.R
#
# stanimal and MCMCglmm sample the posterior; sommer and ASReml-R maximize the
# REML likelihood. All four need to be installed (ASReml-R is licensed).

library(stanimal)

## Simulation ----------------------------------------------------------------

set.seed(1234)

n_gen <- 5; n_breed <- 100
ped <- data.frame(id = seq_len(n_gen * n_breed),
                  sire = NA_integer_, dam = NA_integer_)
for (g in 2:n_gen) {
  kids <- ((g - 1) * n_breed + 1):(g * n_breed)
  pool <- ped[ped$id <= (g - 1) * n_breed, ]
  ped$sire[kids] <- sample(pool$id[pool$id %% 2 == 1], n_breed, TRUE)
  ped$dam[kids]  <- sample(pool$id[pool$id %% 2 == 0], n_breed, TRUE)
}
A <- makeA(ped)

G <- matrix(c(0.8, 0.4, 0.4, 1.2), 2, dimnames = list(c("t1", "t2"), c("t1", "t2")))
E <- matrix(c(1.0, 0.3, 0.3, 0.9), 2, dimnames = list(c("t1", "t2"), c("t1", "t2")))

n <- nrow(A)
X <- cbind(Intercept = 1, sex = rbinom(n, 1, 0.5), age = round(rnorm(n, 3, 0.5), 2))
B <- rbind(c(10, -5), c(1.5, 0.8), c(0.4, -0.2))

Y <- X %*% B +
  t(chol(A)) %*% matrix(rnorm(n * 2), n) %*% chol(G) +   # breeding values
  matrix(rnorm(n * 2), n) %*% chol(E)                    # residuals
colnames(Y) <- colnames(G)

dat <- data.frame(id = factor(ped$id), sex = X[, "sex"], age = X[, "age"], Y)

## Comparison helpers --------------------------------------------------------

# One row per engine: the three distinct entries of G and of E, heritabilities,
# the genetic correlation and the time the fit took in seconds.
compare <- function(engine, G_hat, E_hat, elapsed = NA_real_) {
  data.frame(engine = engine,
             G11 = G_hat[1, 1], G12 = G_hat[1, 2], G22 = G_hat[2, 2],
             E11 = E_hat[1, 1], E12 = E_hat[1, 2], E22 = E_hat[2, 2],
             h2_1 = G_hat[1, 1] / (G_hat[1, 1] + E_hat[1, 1]),
             h2_2 = G_hat[2, 2] / (G_hat[2, 2] + E_hat[2, 2]),
             rG = cov2cor(G_hat)[1, 2],
             time_s = elapsed)
}

# lower triangle, column order -> symmetric matrix
lower2mat <- function(x, k = 2) {
  m <- matrix(0, k, k)
  m[lower.tri(m, diag = TRUE)] <- x
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

results <- list(compare("truth", G, E))

## stanimal: Bayesian, breeding values marginalized out of the likelihood ----

design <- animal_design(cbind(t1, t2) ~ sex + age, dat, pedigree = ped)
time_stan <- system.time(
  fit_stan <- lmm_animal(design$Y, design$X, design$A,
                         iter = 2000, warmup = 1000, chains = 4, cores = 4,
                         refresh = 0))

# G[i,j] and E[i,j] draws, summarised as trait by trait posterior means
results$stanimal <- compare("stanimal", draws_mean(fit_stan, "G"),
                            draws_mean(fit_stan, "E"),
                            elapsed = time_stan[["elapsed"]])

## MCMCglmm: Bayesian, pedigree supplied as the inverse of A ----------------

library(MCMCglmm)
prior <- list(G = list(G = list(V = diag(2) * 0.5, nu = 2)),
              R = list(R = list(V = diag(2) * 0.5, nu = 2)))
time_mcmc <- system.time(
  fit_mcmc <- MCMCglmm(cbind(t1, t2) ~ trait + trait:sex + trait:age - 1,
                       random   = ~ us(trait):id,
                       rcov     = ~ us(trait):units,
                       ginverse = list(id = MCMCglmm::inverseA(ped)$Ainv),
                       family   = c("gaussian", "gaussian"),
                       prior = prior, data = dat,
                       nitt = 13000, burnin = 3000, thin = 10, verbose = FALSE))
# VCV columns are G[t1,t1], G[t2,t1], G[t1,t2], G[t2,t2], then the same for units
results$MCMCglmm <- compare("MCMCglmm",
                            matrix(colMeans(fit_mcmc$VCV[, 1:4]), 2),
                            matrix(colMeans(fit_mcmc$VCV[, 5:8]), 2),
                            elapsed = time_mcmc[["elapsed"]])

## sommer: REML, relationship matrix passed directly to vsr() ----------------

library(sommer)
# mmer() still works; sommer >= 4.5 also points to mmes() with the same syntax
time_sommer <- system.time(
  fit_sommer <- mmer(cbind(t1, t2) ~ sex + age,
                     random = ~ vsr(id, Gu = A),
                     rcov   = ~ vsr(units, Gtc = unsm(2)),
                     data = dat, verbose = FALSE))
results$sommer <- compare("sommer", fit_sommer$sigma$`u:id`,
                          fit_sommer$sigma$`u:units`,
                          elapsed = time_sommer[["elapsed"]])

## ASReml-R: REML, inverse of A passed through vm() --------------------------

library(asreml)
ped0 <- ped; ped0[is.na(ped0)] <- 0                 # ASReml codes founders as 0
# asreml::ainverse() is not MCMCglmm::inverseA(): same pedigree, two packages
ainv <- try(asreml::ainverse(ped0), silent = TRUE)  # ASReml's sparse inverse of A
if (inherits(ainv, "try-error")) {
  message("asreml::ainverse() failed to form A^-1, skipping the ASReml-R fit")
} else {
  # Time both passes: bivariate ASReml-R fits are estimated in two steps
  time_asr <- system.time({
    fit_asr <- asreml(fixed  = cbind(t1, t2) ~ trait + trait:sex + trait:age,
                      random = ~ vm(id, ainv):us(trait, init = c(1, 0.2, 1)),
                      rcov   = ~ units:us(trait, init = c(1, 0.2, 1)),
                      data = dat)
    fit_asr <- update(fit_asr)
  })
  vc <- summary(fit_asr)$varcomp
  results$asreml <- compare("asreml",
                            lower2mat(vc[grep("^vm\\(", rownames(vc)), 1]),
                            lower2mat(vc[grep("^units.*trait!trait", rownames(vc)), 1]),
                            elapsed = time_asr[["elapsed"]])
}

## The comparison -----------------------------------------------------------

results <- do.call(rbind, results)
rownames(results) <- results$engine
out <- round(results[, -1], 3)
out[, "time_s"] <- round(results$time_s, 1)
print(out)

# Timings are a single run each, in seconds, and cover the fit call only.
# Setup outside the call is not timed: animal_design() builds A for stanimal and
# ainverse() is run before the ASReml-R fit, while the inverse of A that
# MCMCglmm needs is evaluated inside its call and so is included. stanimal's
# time also includes compiling the model when the cache is cold; after the first
# fit the compiled executable is reused and only sampling is timed. These are
# not like-for-like posteriors either: MCMCglmm and stanimal sample and stop at
# different effective sample sizes, while sommer and ASReml-R maximize the REML
# likelihood.
