if(!require(MCMCglmm)){install.packages("MCMCglmm"); library(MCMCglmm)}
if(!require(mvtnorm)){install.packages("mvtnorm"); library(mvtnorm)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(mvtnorm)){install.packages("mvtnorm"); library(mvtnorm)}
if(!require(nadiv)){install.packages("nadiv"); library(nadiv)}
if(!require(tictoc)){install.packages("tictoc"); library(tictoc)}
if(!require(pedtools)){install.packages("pedtools"); library(pedtools)}
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

library(stanAnimal)
library(AtchleyMice)

lt = function(x, diag = T) x[lower.tri(x, diag = diag)]

#ped <- read.table("volesPED.txt",header=T)
ped = randomPed(1500, 50, seed = 2)
#ped = ped$`_comp1`
ped = data.frame(ID = ped$ID, sire = ped$FIDX, dam = ped$MIDX)
ped = prepPed(ped)
A <- as.matrix(nadiv::makeA(ped))
nrow(A)
#ped <- mice_pedigree
#F6_ID = mice_info$F6$ID
#A <- as.matrix(nadiv::makeA(ped))[F6_ID, F6_ID]
L_A = chol(A)

n_traits = 20

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

sex = as.numeric(as.factor(mice_info$F6$Sex))-1
sex = sex[1:nrow(a)]
Z = rnorm(nrow(a))
Intercept = rep(1, nrow(a))
X = cbind(Intercept, sex, Z)
rownames(X) = rownames(a)
e = rmvnorm(nrow(a), sigma = E)

#Y = X %*% beta + a + e
Y = 1 + a + e


apply(Y, 2, mean)
apply(Y, 2, sd)
apply(a, 2, sd)
stan_model = lmm_animal(Y, X[,1, drop = F], A, lkj_prior = 4, chains = 4, iter = 4000, warmup = 2000, 
                        cores = 4, control = list(max_treedepth = 10))

model = rstan::extract(stan_model, permuted = TRUE)
str(model)
rstan::summary(stan_model, pars = c("h2", "L_sigma", "L_sigma_G"))[[1]]
rstan::summary(stan_model, pars = "G")[[1]]
dimnames(model)
colMeans(model[,1,grep("^L_sigma", dimnames(model)[[3]])])
colMeans(model[,2,grep("^L_sigma", dimnames(model)[[3]])])
colMeans(model[,3,grep("^L_sigma", dimnames(model)[[3]])])
colMeans(model[,4,grep("^L_sigma", dimnames(model)[[3]])])
diag(G_stan)
plot(colMeans(model$a), a, pch = 19)

par(mfrow = c(1, 1))
plot(lt(cov(a)), lt(colMeans(model$G)), pch = 19)
abline(0, 1)

colMeans(model$beta)
t(beta)

mcmc_intervals(
  as.array(stan_model),
  pars = c("h2[1]","h2[2]","h2[3]"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean") + geom_vline(xintercept = 0.33)

mcmc_intervals(
  as.array(stan_model),
  pars = c("G[1,2]", "G[1,3]", "G[2,3]"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean") + geom_vline(xintercept = G[lower.tri(G)])

mcmc_intervals(
  as.array(stan_model),
  pars = c("corrG[1,2]", "corrG[1,3]", "corrG[2,3]"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean") + geom_vline(xintercept = cor(a)[lower.tri(G)])

library(corrplot)
par(mfrow=c(1,3))
corrplot.mixed(cor(a), upper = "ellipse")
corrplot.mixed(colMeans(model$corrG), upper = "ellipse")
corrplot.mixed(colMeans(model$corrE), upper = "ellipse")


