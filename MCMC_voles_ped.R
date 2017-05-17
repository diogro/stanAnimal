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


prior_bi <- list(G = list(G1 = list(V = diag(2), n = 1.002)),
                          R = list(V = diag(2), n = 1.002))
model_bi <- MCMCglmm(cbind(aggression, size) ~ trait + trait:sex + trait:forage,
                    random = ~us(trait):animal,
                    rcov = ~us(trait):units, family = c("gaussian", "gaussian"),
                    pedigree = volesPED, data = voles_scaled, prior = prior_bi,
                    nitt = 65000, thin = 50, burnin = 15000, verbose = TRUE)
summary(model_bi)

inv.phylo <- MCMCglmm::inverseA(volesPED, scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)
model_simple <- brm(value ~ trait + trait:sex + trait:forage + (trait|animal) + (trait|ID), 
                    data = narrow_voles, 
                    family = gaussian(), cov_ranef = list(animal = A),
                    prior = c(prior(normal(0, 10), "b"),
                              prior(normal(0, 10), "Intercept"),
                              prior(student_t(3, 0, 20), "sd"),
                              prior(student_t(3, 0, 20), "sigma")),
                    iter = 1000, chains = 4, control = list(adapt_delta = 0.95))
summary(model_simple)
stancode(model_simple)

ntraits = 2
inv.phylo <- MCMCglmm::inverseA(volesPED, scale = TRUE)
A <- solve(inv.phylo$Ainv)
A = (A + t(A))/2
rownames(A) <- rownames(inv.phylo$Ainv)
pos = aaply(as.numeric(voles_scaled$ID), 1, function(x) which(x == rownames(A)))

lm_model = lm(cbind(aggression, size) ~ sex + forage, data = voles_scaled)
X = model.matrix(lm_model)

stan_data = list(K = 2,
                 J = ncol(X),
                 N = nrow(voles_scaled),
                 X = X,
                 Y = voles_scaled[,c("aggression", "size")],
                 Z = pos,
                 A = as.matrix(chol(A)))

stan_model = stan(file = "./animalModel.stan", data = stan_data, chains = 4, iter = 1000)
model = rstan::extract(stan_model)
colMeans(model$G)
