library(stanAnimal)
rstan_options(auto_write = TRUE)
options(mc.cores = 10)

devtools::install_github("diogro/ratones")
library(ratones)

data("ratones")
Y = sqrt(10)*as.matrix(dplyr::select(ratonesdf, IS_PM:NSL_NA))
X = model.matrix(lm(IS_PM ~ SEX + LIN, data = ratonesdf))
cov_fit = lmm_multiCov(Y, X)              
model = rstan::extract(cov_fit)
colMeans(model$P)
cov(residuals(lm(Y~X +0)))