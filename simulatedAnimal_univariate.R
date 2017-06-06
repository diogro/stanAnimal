library(MCMCglmm)
library(rstan)
library(mvtnorm)

rstan_options(auto_write = TRUE)
options(mc.cores = 10)

set.seed(17)

# Pedigree
ped <- read.table("https://raw.githubusercontent.com/diogro/QGcourse/master/tutorials/volesPED.txt",header=T)

# Genetic variance
G <- 0.7

# Residual variance
E <- 0.3

# Simulated breeding values
a = rbv(ped, G)

# Fixed effects
beta = matrix(c(1,
                0.3), 2, 1, byrow = TRUE)
rownames(beta) = c("Intercept", "sex")

sex = numeric(nrow(a))
sex = sample(c(0, 1), length(sex), replace = TRUE)
sex[rownames(a) %in% ped$SIRE] <- 1
sex[rownames(a) %in% ped$DAM] <- 0
Intercept = rep(1, nrow(a))
X = cbind(Intercept, sex)
rownames(X) = rownames(a)

# Correlated noise
e = rnorm(nrow(a), sqrt(E))

# Simulated data
Y = X %*% beta + a + e
colnames(Y) = c("x")

# Create relationship matrix
inv.phylo <- MCMCglmm::inverseA(ped, scale = TRUE)
A <- solve(inv.phylo$Ainv)
A = (A + t(A))/2 # Not always symmetric after inversion
rownames(A) <- rownames(inv.phylo$Ainv)

# Relate individuals to columns of the A matrix
pos = unlist(lapply(rownames(Y), function(x) which(x == rownames(A))))

stan_data = list(J = ncol(X),
                 N = nrow(Y),
                 X = X,
                 Y = as.numeric(Y),
                 A = as.matrix(chol(A)))
stan_model = stan(file = "./animalModelUni.stan", data = stan_data)
model = rstan::extract(stan_model)
rstan::summary(stan_model, pars = "sigma_E")[[1]]
rstan::traceplot(stan_model, pars = "sigma_E")
colMeans(model$G)
