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
library(ratones)
rstan_options(auto_write = TRUE)
options(mc.cores = 10)


formula = paste0("cbind(",
                paste(names(dplyr::select(ratonesdf, IS_PM:IS_PNS)),
                      collapse = ", "), ") ~ SEX + AGE + LIN")
stan_data = generateAnimalModelInput(formula, ratonesdf,
                                     ratones_ped$A,
                                     line = "all")
stan_model = stan(file = "./animalModel.stan", data = stan_data, chains = 1, 
                  iter = 800,
                  control = list(adapt_delta = 0.99))
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
  point_est = "mean")
mcmc_intervals( 
  as.array(stan_model),  
  pars = c("corrG[1,2]", "corrG[1,3]", "corrG[2,3]"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean")