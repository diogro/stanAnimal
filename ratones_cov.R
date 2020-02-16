library(rstan)
library(stanAnimal)
rstan_options(auto_write = TRUE)
options(mc.cores = 10)

devtools::install_github("diogro/ratones")
library(ratones)

data("ratones")
trait_1 = which(names(ratonesdf) == "IS_PM")
n_traits = 15 
Y = sqrt(10)*as.matrix(ratonesdf[, trait_1:(n_traits + trait_1 - 1)])
X = model.matrix(lm(IS_PM ~ SEX + LIN, data = ratonesdf))
cov_fit = lmm_multiCov(Y, X, LKJ = 4, warmup = 3000, iter = 6000, 
                       control = list(adapt_delta = 0.99, max_treedepth = 12))              
model = rstan::extract(cov_fit)
colMeans(model$P)
cov(residuals(lm(Y~X +0)))
dim(Y)
