# The same animal model fitted with several engines, for a chosen number of
# traits
#
# Simulates K genetically correlated traits with a known G and E on a
# five-generation pedigree, fits the same multivariate animal model with
# stanimal, MCMCglmm, BGLR, sommer, lme4breeding and ASReml-R, and compares the
# estimates with the simulated truth and the time each fit took.
#
# K is the first command-line argument (default 5):
#   Rscript examples/trait_comparison.R       # five traits
#   Rscript examples/trait_comparison.R 3     # three traits
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

## Number of traits ----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
n_traits <- if (length(args)) as.integer(args[1L]) else 5L
if (is.na(n_traits) || n_traits < 2L)
  stop("the number of traits must be a whole number of at least 2.",
       call. = FALSE)
traits <- paste0("t", seq_len(n_traits))

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

# Simulated truth from the package's LKJ construction: correlations are shrunk
# toward zero, so G and E are well conditioned and REML does not sit on the
# boundary.
set.seed(4321)
G <- random_covariance(n_traits, min_var = 0.7, max_var = 1.3,
                       trait_names = traits)
E <- random_covariance(n_traits, min_var = 0.7, max_var = 1.3,
                       trait_names = traits)

n <- nrow(A)
X <- cbind(Intercept = 1, sex = rbinom(n, 1, 0.5), age = round(rnorm(n, 3, 0.5), 2))
B <- rbind(rep(10, n_traits), rep(1.5, n_traits), rep(0.4, n_traits))

Y <- X %*% B +
  t(chol(A)) %*% matrix(rnorm(n * n_traits), n) %*% chol(G) +
  matrix(rnorm(n * n_traits), n) %*% chol(E)
colnames(Y) <- traits

dat <- data.frame(id = factor(ped$id), sex = X[, "sex"], age = X[, "age"], Y)

## Comparison helpers --------------------------------------------------------

# ASReml reports us() parameters as the lower triangle, row by row
lower2mat <- function(x, k) {
  m <- matrix(0, k, k)
  rows <- rep(seq_len(k), seq_len(k))
  cols <- sequence(seq_len(k))
  m[cbind(rows, cols)] <- x
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

corr_of <- function(S) {
  if (nrow(S) < 2L) return(numeric(0))
  cov2cor(S)[upper.tri(S)]
}

h2_true <- diag(G) / diag(G + E)
rG_true <- corr_of(G)

# One row per engine: mean squared error of G, E, the heritabilities and the
# genetic correlations, and the time the fit took.
compare <- function(engine, G_hat, E_hat, elapsed = NA_real_) {
  ok <- !is.null(G_hat) && !is.null(E_hat) &&
    !anyNA(G_hat) && !anyNA(E_hat)
  if (!ok)
    return(data.frame(engine = engine, G_mse = NA_real_, E_mse = NA_real_,
                      h2_mse = NA_real_, corrG_mse = NA_real_,
                      time_s = elapsed))
  data.frame(engine = engine,
             G_mse = mean((G_hat - G)^2),
             E_mse = mean((E_hat - E)^2),
             h2_mse = mean((diag(G_hat) / diag(G_hat + E_hat) - h2_true)^2),
             corrG_mse = mean((corr_of(G_hat) - rG_true)^2),
             time_s = elapsed)
}

results <- list(compare("truth", G, E))

# One row per engine, filled in as each fit runs, so an engine that is not
# installed (or whose fit is skipped) still appears in the table, with NA
# estimates and no time.
engines <- c("stanimal", "MCMCglmm", "BGLR", "sommer", "lme4breeding",
             "asreml")
na_mat <- matrix(NA_real_, n_traits, n_traits)
for (engine in engines)
  results[[engine]] <- compare(engine, na_mat, na_mat)

## stanimal: Bayesian, breeding values marginalized out of the likelihood ----

response <- sprintf("cbind(%s)", paste(traits, collapse = ", "))
design <- animal_design(reformulate(c("sex", "age"), response = response),
                        dat, pedigree = ped)
time_stan <- system.time(
  fit_stan <- lmm_animal(design$Y, design$X, design$A,
                         iter = 2000, warmup = 1000, chains = 4, cores = 4,
                         refresh = 0))
results$stanimal <- compare("stanimal", draws_mean(fit_stan, "G"),
                            draws_mean(fit_stan, "E"),
                            elapsed = time_stan[["elapsed"]])

## MCMCglmm: Bayesian, pedigree supplied as the inverse of A ----------------

library(MCMCglmm)
prior <- list(G = list(G = list(V = diag(n_traits) * 0.5, nu = n_traits)),
              R = list(R = list(V = diag(n_traits) * 0.5, nu = n_traits)))
time_mcmc <- system.time(
  fit_mcmc <- MCMCglmm(as.formula(paste(response,
                                        "~ trait + trait:sex + trait:age - 1")),
                       random   = ~ us(trait):id,
                       rcov     = ~ us(trait):units,
                       ginverse = list(id = MCMCglmm::inverseA(ped)$Ainv),
                       family   = rep("gaussian", n_traits),
                       prior = prior, data = dat,
                       nitt = 13000, burnin = 3000, thin = 10, verbose = FALSE))
# VCV holds the K x K G matrix in column order, then the K x K residual matrix
k2 <- n_traits^2
results$MCMCglmm <- compare(
  "MCMCglmm",
  matrix(colMeans(fit_mcmc$VCV[, seq_len(k2)]), n_traits),
  matrix(colMeans(fit_mcmc$VCV[, k2 + seq_len(k2)]), n_traits),
  elapsed = time_mcmc[["elapsed"]])

## BGLR: Bayesian multitrait, pedigree relationship as an RKHS kernel -------

if (requireNamespace("BGLR", quietly = TRUE)) {
  X_bglr <- model.matrix(~ sex + age, dat)
  A_dat <- A[match(as.character(dat$id), rownames(A)),
             match(as.character(dat$id), rownames(A))]
  # Weak inverse-Wishart priors with prior mean I, for G and for the residual
  prior_cov <- list(type = "UN", df0 = n_traits + 3,
                    S0 = 2 * diag(n_traits))
  time_bglr <- system.time(
    fit_bglr <- BGLR::Multitrait(
      y = as.matrix(dat[traits]),
      ETA = list(list(X = X_bglr, model = "FIXED"),
                 list(K = A_dat, model = "RKHS", Cov = prior_cov)),
      intercept = FALSE, resCov = prior_cov,
      nIter = 13000, burnIn = 3000, thin = 10,
      saveAt = file.path(tempdir(), "bglr_"), verbose = FALSE))
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
  id = factor(rep(dat$id, n_traits)),
  trait = factor(rep(traits, each = nrow(dat)), levels = traits),
  value = unlist(dat[traits], use.names = FALSE),
  sex = rep(dat$sex, n_traits),
  age = rep(dat$age, n_traits),
  record = rep(seq_len(nrow(dat)), n_traits))
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

if (requireNamespace("lme4breeding", quietly = TRUE)) {
  # lmeb() calls lme4 helpers from the search path, so attach the package
  library(lme4breeding)
  long <- data.frame(
    id = factor(rep(dat$id, n_traits)),
    trait = factor(rep(traits, each = nrow(dat)), levels = traits),
    value = unlist(dat[traits], use.names = FALSE),
    sex = rep(dat$sex, n_traits),
    age = rep(dat$age, n_traits))
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
    # One initial variance per trait and one covariance per pair, in the
    # row-wise lower-triangle order us() expects
    init_cov <- function(k) {
      m <- matrix(0.2, k, k)
      diag(m) <- 1
      m[cbind(rep(seq_len(k), seq_len(k)), sequence(seq_len(k)))]
    }
    time_asr <- system.time({
      fit_asr <- asreml(fixed  = as.formula(paste(response,
                                                  "~ trait + trait:sex + trait:age")),
                        random = ~ vm(id, ainv):us(trait, init = init_cov(n_traits)),
                        rcov   = ~ units:us(trait, init = init_cov(n_traits)),
                        data = dat)
      fit_asr <- update(fit_asr)
    })
    vc <- summary(fit_asr)$varcomp
    results$asreml <- compare(
      "asreml",
      lower2mat(vc[grep("^vm\\(", rownames(vc)), 1], n_traits),
      lower2mat(vc[grep("^units.*trait!trait", rownames(vc)), 1], n_traits),
      elapsed = time_asr[["elapsed"]])
  }
}

## The comparison -----------------------------------------------------------

results <- do.call(rbind, results)
rownames(results) <- results$engine
out <- round(results[, -1], 4)
out[, "time_s"] <- round(results$time_s, 1)
cat("\n", n_traits, "traits ,", n, "individuals \n", sep = "")
print(out)

# G_mse and E_mse are the elementwise mean squared errors against the simulated
# matrices; h2_mse and corrG_mse compare the trait-level summaries (corrG uses
# the upper triangle, since the diagonal is always 1). Timings are a single run
# each, in seconds, and cover the fit call only. Setup outside the call is not
# timed: animal_design() builds A for stanimal, solve(A) inverts it for sommer,
# and ainverse() is run before the ASReml-R fit, while the inverse of A that
# MCMCglmm needs is evaluated inside its call and so is included. stanimal's
# time also includes compiling the model when the cache is cold; after the first
# fit the compiled executable is reused and only sampling is timed. These are
# not like-for-like posteriors either: MCMCglmm, BGLR and stanimal sample and
# stop at different effective sample sizes, while sommer, lme4breeding and
# ASReml-R maximize the REML likelihood.
