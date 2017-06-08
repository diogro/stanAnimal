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

set.seed(3)

# Pedigree
ped <- read.table("https://raw.githubusercontent.com/diogro/QGcourse/master/tutorials/volesPED.txt",header=T)

# Genetic variance
G <- 1.

# Residual variance
E <- 2.

# Create relationship matrix
inv.phylo <- MCMCglmm::inverseA(ped, scale = TRUE)
A <- solve(inv.phylo$Ainv)
A = (A + t(A))/2 # Not always symmetric after inversion
rownames(A) <- rownames(inv.phylo$Ainv)
A[A < 1e-10] = 0

# Simulated breeding values
#a = G * (as.matrix(chol(A)) %*% rnorm(nrow(A)))
#a = rbv(pedigree = ped, G = G)
a = t(rmvnorm(1, sigma = G*as.matrix(A)))

# Fixed effects
beta = matrix(c(1,
                0.3), 2, 1, byrow = TRUE)
rownames(beta) = c("Intercept", "sex")

sex = sample(c(0, 1), nrow(a), replace = TRUE)
sex[rownames(a) %in% ped$SIRE] <- 1
sex[rownames(a) %in% ped$DAM] <- 0
Intercept = rep(1, nrow(a))
X = cbind(Intercept, sex)
rownames(X) = rownames(a)

# Correlated noise
e = rnorm(nrow(a), sd = sqrt(E))

# Simulated data
Y = X %*% beta + a + e
colnames(Y) = c("x")

ped2 <- pedigree(ped$SIRE, ped$DAM, ped$ID)  #restructure ped file
data = data.frame(Y = Y, sex = sex, ID = rownames(A))
mod_animalREML<-pedigreemm(Y ~ sex + (1|ID), pedigree=list(ID=ped2), 
                           data = data, REML=TRUE, 
                           control = lmerControl(check.nobs.vs.nlev="ignore",
                                                 check.nobs.vs.nRE="ignore"))
summary(mod_animalREML)
#launch_shinystan(stan_model)
stan_data = list(J = ncol(X),
                 N = nrow(Y),
                 X = X,
                 Y = as.numeric(Y),
                 A = as.matrix(A))
stan_model = stan(file = "./animalModelUni.stan", data = stan_data, iter = 2000, control = list(adapt_delta = 0.99))
rstan::summary(stan_model, pars = c("sigma_G", "sigma_E"))[[1]]
model = rstan::extract(stan_model, permuted = FALSE, pars = c("sigma_G", "sigma_R"))
model_list = vector("list", 2)
names(model_list) = c("sigma_G", "sigma_R")
model_list$sigma_G = matrix(model[,,"sigma_G"], ncol = 1)
model_list$sigma_R = matrix(model[,,"sigma_R"], ncol = 1)
sampler_params <- get_sampler_params(stan_model, inc_warmup = FALSE)
energy <- sapply(sampler_params, function(x) x[, "energy__"])
model_list$energy__ = matrix(energy, ncol = 1)
pairs(model_list)


mcmc_intervals( 
  as.array(stan_model),  
  pars = c("sigma_G", "sigma_E"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean")
