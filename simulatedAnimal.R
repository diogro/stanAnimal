if(!require(MCMCglmm)){install.packages("MCMCglmm"); library(MCMCglmm)}
if(!require(mvtnorm)){install.packages("mvtnorm"); library(mvtnorm)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(mvtnorm)){install.packages("mvtnorm"); library(mvtnorm)}
if(!require(nadiv)){install.packages("nadiv"); library(nadiv)}
if(!require(tictoc)){install.packages("tictoc"); library(tictoc)}
if(!require(pedtools)){install.packages("pedtools"); library(pedtools)}

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

library(stanAnimal)
"cmdstanr" %in% installed.packages()[,"Package"]
library(AtchleyMice)

#ped <- read.table("volesPED.txt",header=T)
ped = randomPed(100, 10, seed = 2)
plot(ped)
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
varG = 1:3*10
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

sex = as.numeric(as.factor(mice_info$F6$Sex))-1
sex = sex[1:nrow(a)]
Z = rnorm(nrow(a))
Intercept = rep(1, nrow(a))
X = cbind(Intercept, sex, Z)
rownames(X) = rownames(a)
e = rmvnorm(nrow(a), sigma = E)

Y = X %*% beta + a + e

colnames(Y) = c("x", "y", "z")

sim_data = data.frame(Y, X, animal = rownames(a))

prior_bi <- list(G = list(G1 = list(V = diag(3), n = 1.002)),
                          R = list(V = diag(3), n = 1.002))
tic()
model_bi <- MCMCglmm(cbind(x, y, z) ~ trait + trait:sex + trait:Z,
                     random = ~us(trait):animal,
                     rcov = ~us(trait):units, family = c("gaussian", "gaussian", "gaussian"),
                     pedigree = ped, data = sim_data, prior = prior_bi,
                     nitt = 130000, thin = 100, burnin = 30000, verbose = TRUE)
toc()
summary(model_bi)
colMeans(model_bi$VCV[,c("traitx:traitx.animal", "traity:traity.animal", "traitz:traitz.animal")])
G_mcmc = matrix(colMeans(model_bi$VCV[, grep("animal", colnames(model_bi$VCV))]), 3, 3)
corrG_mcmc = cov2cor(G_mcmc)


stan_model = lmm_animal(Y, X, A, chains = 4, iter = 200, warmup = 100, 
                        cores = 4)



model = rstan::extract(stan_model)
rstan::summary(stan_model, pars = "h2")[[1]]
rstan::summary(stan_model, pars = "G")[[1]]
rstan::plot(stan_model, pars = "beta")
rstan::traceplot(stan_model, pars = "G")
rstan::traceplot(stan_model, pars = "E")
colMeans(model$L_sigma)
colMeans(model$corrG)
plot(colMeans(model$a), a)
colMeans(model$corrE)
cov2cor(colMeans(model$P))

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
  pars = c("G[1,1]", "G[2,2]", "G[3,3]"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean") + geom_vline(xintercept = colMeans(model_bi$VCV[,c("traitx:traitx.animal",
                                                                         "traity:traity.animal",
                                                                         "traitz:traitz.animal")]))
mcmc_intervals(
  as.array(stan_model),
  pars = c("E[1,1]", "E[2,2]", "E[3,3]"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean") + geom_vline(xintercept = colMeans(model_bi$VCV[,c("traitx:traitx.units",
                                                                         "traity:traity.units",
                                                                         "traitz:traitz.units")]))
mcmc_intervals(
  as.array(stan_model),
  pars = c("G[1,2]", "G[1,3]", "G[2,3]"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean") + geom_vline(xintercept = G[lower.tri(corrG)])

mcmc_intervals(
  as.array(stan_model),
  pars = c("corrG[1,2]", "corrG[1,3]", "corrG[2,3]"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean") + geom_vline(xintercept = corrG[lower.tri(corrG)])

library(corrplot)
par(mfrow=c(1,3))
corrplot.mixed(cor(a), upper = "ellipse")
corrplot.mixed(colMeans(model$corrG), upper = "ellipse")
corrplot.mixed(corrG_mcmc, upper = "ellipse")

