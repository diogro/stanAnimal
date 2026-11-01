# One simulated data set, six engines for a bivariate animal model
#
# Simulates two genetically correlated traits with a known G and E on a
# five-generation pedigree, then fits the same bivariate animal model with
# stanimal, MCMCglmm, BGLR, sommer, lme4breeding and ASReml-R, and compares the
# estimates with the simulated truth and the time each fit took. Run from the
# project root:
#   Rscript examples/package_comparison.R
#
# stanimal, MCMCglmm and BGLR sample the posterior; sommer, lme4breeding and
# ASReml-R maximize the REML likelihood. The open-source engines are installed
# below if they are missing; ASReml-R is licensed and has to be installed
# separately.

# Install the engines that are missing (stanimal is assumed installed)
pkgs <- c("MCMCglmm", "BGLR", "sommer", "lme4breeding")
new <- setdiff(pkgs, rownames(utils::installed.packages()))
if (length(new))
  install.packages(new, repos = "https://cloud.r-project.org")

# ASReml-R is licensed and not on CRAN: only check whether it is available,
# never try to install it
has_asreml <- requireNamespace("asreml", quietly = TRUE)

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
             rG = if (anyNA(G_hat)) NA_real_ else cov2cor(G_hat)[1, 2],
             time_s = elapsed)
}

# ASReml reports us() parameters as the lower triangle, row by row
lower2mat <- function(x, k = 2) {
  m <- matrix(0, k, k)
  rows <- rep(seq_len(k), seq_len(k))
  cols <- sequence(seq_len(k))
  m[cbind(rows, cols)] <- x
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

results <- list(compare("truth", G, E))

# One row per engine, filled in as each fit runs, so an engine that is not
# installed (or whose fit is skipped) still appears in the table, with NA
# estimates and no time.
engines <- c("stanimal", "MCMCglmm", "BGLR", "sommer", "lme4breeding",
             "asreml")
na_mat <- matrix(NA_real_, 2, 2)
for (engine in engines)
  results[[engine]] <- compare(engine, na_mat, na_mat)

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

## BGLR: Bayesian multitrait, pedigree relationship as an RKHS kernel -------
# Multitrait gives each trait its own fixed effects, and a K kernel supplies  
# the between-individual covariance: K = A makes the random effect u ~ N(0, A
# (x) G), the animal model. The residual covariance is fitted unstructured.
if (requireNamespace("BGLR", quietly = TRUE)) {
  X_bglr <- model.matrix(~ sex + age, dat)
  A_dat <- A[match(as.character(dat$id), rownames(A)),
             match(as.character(dat$id), rownames(A))]
  # S0 and df0 are set explicitly so the priors are in the same spirit as the
  # other engines and Multitrait does not print the values it would pick
  prior_cov <- list(type = "UN", df0 = 5, S0 = 2 * diag(2))
  time_bglr <- system.time(
    fit_bglr <- BGLR::Multitrait(
      y = cbind(dat$t1, dat$t2),
      ETA = list(list(X = X_bglr, model = "FIXED"),
                 list(K = A_dat, model = "RKHS", Cov = prior_cov)),
      intercept = FALSE, resCov = prior_cov,
      nIter = 13000, burnIn = 3000, thin = 10,
      saveAt = file.path(tempdir(), "bglr_"), verbose = FALSE))
  # Omega is the trait by trait genetic covariance, R the residual one
  results$BGLR <- compare("BGLR", fit_bglr$ETA[[2]]$Cov$Omega,
                          fit_bglr$resCov$R,
                          elapsed = time_bglr[["elapsed"]])
} else {
  message("BGLR is not installed; install.packages(\"BGLR\") to add it")
}

## sommer: REML through mmes(), in long format -------------------------------

library(sommer)
# mmes() takes the inverse relationship matrix, marked as such, and a
# long-format response: usm(trait) estimates the trait covariance and ism(id)
# applies A^-1 to the animal effect.
Ainv <- solve(A)
attr(Ainv, "inverse") <- TRUE
long_sommer <- data.frame(
  id = factor(rep(dat$id, 2)),
  trait = factor(rep(c("t1", "t2"), each = nrow(dat))),
  value = c(dat$t1, dat$t2),
  sex = rep(dat$sex, 2),
  age = rep(dat$age, 2),
  record = rep(seq_len(nrow(dat)), 2))
time_sommer <- system.time(
  fit_sommer <- mmes(
    value ~ 0 + trait + trait:sex + trait:age,
    random = ~ vsm(usm(trait), ism(id), Gu = Ainv),
    rcov   = ~ vsm(usm(trait), ism(record)),
    data = long_sommer,
    control = mmes_control(diagnostics = list(verbose = FALSE))))
# covmatrix_mmes() returns the trait covariance of each term: 1 is the animal
# effect (G) and 2 the residual (E)
results$sommer <- compare("sommer",
                          covmatrix_mmes(fit_sommer, 1)$covariance,
                          covmatrix_mmes(fit_sommer, 2)$covariance,
                          elapsed = time_sommer[["elapsed"]])

## lme4breeding: REML through lme4, relationship matrix in relmat() ----------
# The multivariate model is written in long format: one row per trait and
# individual, with trait-specific intercepts and slopes. (0 + trait | id) is
# the genetic term, with relmat supplying A, and (0 + trait | units) is the
# residual covariance. lme4 keeps a scalar residual variance beside that term.
if (requireNamespace("lme4breeding", quietly = TRUE)) {
  # lmeb() calls lme4 helpers from the search path, so attach the package
  library(lme4breeding)
  long <- data.frame(
    id = factor(rep(dat$id, 2)),
    trait = factor(rep(c("t1", "t2"), each = nrow(dat))),
    value = c(dat$t1, dat$t2),
    sex = rep(dat$sex, 2),
    age = rep(dat$age, 2))
  long$units <- long$id
  time_l4b <- system.time(
    fit_l4b <- lmeb(
      value ~ 0 + trait + trait:sex + trait:age +
        (0 + trait | id) + (0 + trait | units),
      relmat = list(id = A), data = long,
      REML = TRUE, rotation = FALSE, verbose = 0, trace = 0))
  vc_l4b <- VarCorr(fit_l4b)
  E_l4b <- as.matrix(vc_l4b$units)
  diag(E_l4b) <- diag(E_l4b) + attr(vc_l4b, "sc")^2
  results$lme4breeding <- compare("lme4breeding", as.matrix(vc_l4b$id),
                                  E_l4b, elapsed = time_l4b[["elapsed"]])
} else {
  message("lme4breeding is not installed; install.packages(\"lme4breeding\") ",
          "to add it")
}

## ASReml-R: REML, inverse of A passed through vm() --------------------------

if (!has_asreml) {
  message("asreml is not installed; it is licensed and not on CRAN, so the ",
          "ASReml-R fit is skipped")
} else {
  library(asreml)
  ped0 <- ped; ped0[is.na(ped0)] <- 0               # ASReml codes founders as 0
  # asreml::ainverse() is not MCMCglmm::inverseA(): same pedigree, two packages
  ainv <- try(asreml::ainverse(ped0), silent = TRUE)  # ASReml's sparse inverse
  if (inherits(ainv, "try-error")) {
    message("asreml::ainverse() failed to form A^-1 (",
            conditionMessage(attr(ainv, "condition")),
            "), skipping the ASReml-R fit")
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
                              lower2mat(vc[grep("^units.*trait!trait",
                                                rownames(vc)), 1]),
                              elapsed = time_asr[["elapsed"]])
  }
}

## The comparison -----------------------------------------------------------

results <- do.call(rbind, results)
rownames(results) <- results$engine
out <- round(results[, -1], 3)
out[, "time_s"] <- round(results$time_s, 1)
print(out)

# Timings are a single run each, in seconds, and cover the fit call only.
# Setup outside the call is not timed: animal_design() builds A for stanimal,
# solve(A) inverts it for sommer, and ainverse() is run before the ASReml-R
# fit, while the inverse of A that MCMCglmm needs is evaluated inside its call
# and so is included. stanimal's time also includes compiling the model when the
# cache is cold; after the first fit the compiled executable is reused and only
# sampling is timed. These are not like-for-like posteriors either: MCMCglmm,
# BGLR and stanimal sample and stop at different effective sample sizes, while
# sommer, lme4breeding and ASReml-R maximize the REML likelihood.
