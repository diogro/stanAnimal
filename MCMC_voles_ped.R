#1. read in the main data file
voles<-read.table("voles2.txt", header=T)
voles$ID<-as.factor(voles$ID)
voles$sex<-as.factor(voles$sex)
voles$SIRE<-as.factor(voles$SIRE)
voles$DAM<-as.factor(voles$DAM)
voles$animal<-as.factor(voles$ID)

#read in the pedigree file
volesPED<-read.table("volesPED.txt",header=T)

# scale data
voles_scaled = voles
voles_scaled[, c("aggression", "size")] = scale(voles[,c("aggression", "size")])
narrow_voles = tidyr::gather(voles_scaled, trait, value, aggression:size)

#2. Automatically installs packages in this list that user does not already have

list_pkgs <- c("lme4","lmerTest","MCMCglmm","pedantics","pedigreemm", "brms", "tidyr")
new_pkgs <- list_pkgs[!(list_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs) > 0){ install.packages(new_pkgs) }

library(pedantics)
library(MCMCglmm)
library(lme4)
library(lmerTest)
library(pedigreemm)
library(plyr)
library(tidyr)
library(brms)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 10)

prior_bi <- list(G = list(G1 = list(V = diag(2), nu = 1.002)),
                          R = list(V = diag(2), nu = 1.002))
model_bi <- MCMCglmm(cbind(aggression, size) ~ trait + trait:sex + trait:forage,
                    random = ~us(trait):animal,
                    rcov = ~us(trait):units, family = c("gaussian", "gaussian"),
                    pedigree = volesPED, data = voles_scaled, prior = prior_bi,
                    nitt = 65000, thin = 50, burnin = 15000, verbose = TRUE)
summary(model_bi)

ntraits = 2
inv.phylo <- MCMCglmm::inverseA(volesPED, scale = TRUE)
A <- solve(inv.phylo$Ainv)
A = (A + t(A))/2
rownames(A) <- rownames(inv.phylo$Ainv)
pos = aaply(as.numeric(voles_scaled$ID), 1, function(x) which(x == rownames(A)))

lm_model = lm(cbind(aggression, size) ~ sex + forage, data = voles_scaled)
X = model.matrix(lm_model)
Y = lm_model$model$`cbind(aggression, size)`

library(stanAnimal)
cov_fit = lmm_multiCov(Y, X)
cov_fit$summary("P", "mean") 

animal_fit = lmm_animal(Y, X, A)
animal_fit$summary("G")

post_mean = function(name, k = 2) {
  s = animal_fit$summary(name, "mean")
  matrix(s$mean[match(sprintf("%s[%d,%d]", name, rep(1:k, k), rep(1:k, each = k)),
                      s$variable)], k)
}
G_hat = post_mean("G")
plot_G_ellipses(animal_fit, Y, n_show = 100)
dev.off()

print(plot_forest(animal_fit, trait_names = colnames(Y), lkj_prior = 2))

model = rstan::extract(animal_fit)
colMeans(model$G)
print(animal_fit, pars = c("G","E", "beta"))

set.seed(1)
n_fam <- 400; fam_size <- 10; n <- n_fam * fam_size
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

fit <- lmm_animal(Y[,1], X, A)
fit$summary(c("h2"))
animal_draws(fit, "beta")
