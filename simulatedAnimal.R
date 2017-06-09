list_pkgs <- c("lme4","lmerTest","MCMCglmm","pedantics","pedigreemm", "brms", "tidyr", "rstan")
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
library(mvtnorm)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 10)

ped <- read.table("volesPED.txt",header=T)
G <- matrix(c(1, 0.7, 0.7, 2), 2, 2)
E <- matrix(c(1.5, 0.2, 0.2, 3), 2, 2)
a = rbv(ped, G)

beta = matrix(c(1, 2,
                0.1, 0.2,
                0.05, 0.1), 3, 2, byrow = TRUE)
colnames(beta) = c("x", "y")
rownames(beta) = c("Intercept", "sex", "Z")

sex = numeric(nrow(a))
sex = sample(c(0, 1), length(sex), replace = TRUE)
sex[rownames(a) %in% ped$SIRE] <- 1
sex[rownames(a) %in% ped$DAM] <- 0
Z = rnorm(nrow(a))
Intercept = rep(1, nrow(a))
X = cbind(Intercept, sex, Z)
rownames(X) = rownames(a)
e = rmvnorm(nrow(a), sigma = E)

Y = X %*% beta + a + e

colnames(Y) = c("x", "y")

sim_data = data.frame(Y, X, animal = rownames(a))

prior_bi <- list(G = list(G1 = list(V = diag(2), n = 1.002)),
                          R = list(V = diag(2), n = 1.002))
model_bi <- MCMCglmm(cbind(x, y) ~ trait + trait:sex + trait:Z,
                     random = ~us(trait):animal,
                     rcov = ~us(trait):units, family = c("gaussian", "gaussian"),
                     pedigree = ped, data = sim_data, prior = prior_bi,
                     nitt = 130000, thin = 100, burnin = 30000, verbose = TRUE)
summary(model_bi)

inv.phylo <- MCMCglmm::inverseA(ped, scale = TRUE)
A <- solve(inv.phylo$Ainv)
A = (A + t(A))/2
rownames(A) <- rownames(inv.phylo$Ainv)

stan_data = list(K = ncol(Y),
                 J = ncol(X),
                 N = nrow(Y),
                 X = X,
                 Y = Y,
                 A = as.matrix(A))
stan_model = stan(file = "./animalModel.stan", data = stan_data, chains = 4, 
                  iter = 2000,
                  control = list(adapt_delta = 0.99))
model = rstan::extract(stan_model)
rstan::summary(stan_model, pars = "G")[[1]]
rstan::plot(stan_model, pars = "beta")
rstan::traceplot(stan_model, pars = "G")
rstan::traceplot(stan_model, pars = "E")
colMeans(model$G)
colMeans(model$aM)
colMeans(model$E)
colMeans(model$beta)
t(beta)
