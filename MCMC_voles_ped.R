#1. read in the main data file
voles<-read.table("voles2.txt", header=T)
voles$ID<-as.factor(voles$ID)
voles$sex<-as.factor(voles$sex)
voles$SIRE<-as.factor(voles$SIRE)
voles$DAM<-as.factor(voles$DAM)

#read in the pedigree file
volesPED<-read.table("volesPED.txt",header=T)

#2. Automatically installs packages in this list that user does not already have

list_pkgs <- c("lme4","lmerTest","MCMCglmm","pedantics","pedigreemm")
new_pkgs <- list_pkgs[!(list_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs) > 0){ install.packages(new_pkgs) }

library(pedantics)
library(MCMCglmm)
library(lme4)
library(lmerTest)
library(pedigreemm)

voles$animal<-as.factor(voles$ID)

prior_bi <- list(G = list(G1 = list(V = diag(2), n = 1.002)),
                          R = list(V = diag(2), n = 1.002))


model_bi <- MCMCglmm(cbind(size, aggression) ~ trait + trait:sex +trait:forage- 1,
                    random = ~us(trait):animal,
                    rcov = ~us(trait):units, family = c("gaussian", "gaussian"),
                    pedigree = volesPED, data = voles, prior = prior_bi,
                    nitt = 65000, thin = 50, burnin = 15000, verbose = TRUE)

summary(model_bi)


inv.phylo <- MCMCglmm::inverseA(volesPED, scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

library(tidyr)
narrow_voles = tidyr::gather(voles, trait, value, aggression:size)

library(brms)
options(mc.cores = parallel::detectCores())
model_simple <- brm(value ~ trait + trait:sex + trait:forage + (trait|ID), 
                    data = narrow_voles, 
                    family = gaussian(), cov_ranef = list(ID = A),
                    prior = c(prior(normal(0, 10), "b"),
                              prior(normal(0, 50), "Intercept"),
                              prior(student_t(3, 0, 20), "sd"),
                              prior(student_t(3, 0, 20), "sigma")),
                    iter = 200)

