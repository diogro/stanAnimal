list_pkgs <- c("lme4","lmerTest","MCMCglmm","pedantics","pedigreemm", "rstan", "mvtnorm", "bayesplot", "shinystan")
new_pkgs <- list_pkgs[!(list_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs) > 0){ install.packages(new_pkgs) }

library(pedantics)
library(MCMCglmm)
library(lme4)
library(lmerTest)
library(pedigreemm)
library(rstan)
library(mvtnorm)
library(bayesplot)
library(shinystan)
rstan_options(auto_write = TRUE)
options(mc.cores = 10)

library(stanAnimal)

ped <- read.table("volesPED.txt",header=T)
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
a = rbv(ped, G)

beta = matrix(c(1, 2, 3, 
                0.1, 0.2, 0.5,
                0.05, 0.1, 0.3), 3, 3, byrow = TRUE)
colnames(beta) = c("x", "y", "z")
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

colnames(Y) = c("x", "y", "z")

sim_data = data.frame(Y, X, animal = rownames(a))

prior_bi <- list(G = list(G1 = list(V = diag(3), n = 1.002)),
                          R = list(V = diag(3), n = 1.002))
model_bi <- MCMCglmm(cbind(x, y, z) ~ trait + trait:sex + trait:Z,
                     random = ~us(trait):animal,
                     rcov = ~us(trait):units, family = c("gaussian", "gaussian", "gaussian"),
                     pedigree = ped, data = sim_data, prior = prior_bi,
                     nitt = 130000, thin = 100, burnin = 30000, verbose = TRUE)
summary(model_bi)
colMeans(model_bi$VCV[,c("traitx:traitx.animal", "traity:traity.animal", "traitz:traitz.animal")])

inv.phylo <- MCMCglmm::inverseA(ped, scale = TRUE)
A <- solve(inv.phylo$Ainv)
A = as.matrix((A + t(A))/2)
rownames(A) <- rownames(inv.phylo$Ainv)

stan_model = lmm_animal(Y, X, A, chains = 4, iter = 2000, warmup = 1000)

model = rstan::extract(stan_model)
rstan::summary(stan_model, pars = "G")[[1]]
rstan::plot(stan_model, pars = "beta")
rstan::traceplot(stan_model, pars = "G")
rstan::traceplot(stan_model, pars = "E")
colMeans(model$G)
plot(colMeans(model$a), a)
colMeans(model$E)
colMeans(model$beta)
t(beta)
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
