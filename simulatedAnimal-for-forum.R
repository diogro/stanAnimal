library(MCMCglmm)
library(rstan)
library(mvtnorm)

set.seed(17)

# Pedigree
ped <- read.table("https://raw.githubusercontent.com/diogro/QGcourse/master/tutorials/volesPED.txt",header=T)

# G matrix
G <- matrix(c(1, 0.7, 0.7, 2), 2, 2)

# Residual matrix
E <- matrix(c(1, 0.2, 0.2, 2), 2, 2)

# Simulated breeding values
a = rbv(ped, G)

# Fixed effects
beta = matrix(c(1, 2,
                0.1, 0.2,
                0.05, 0.1), 3, 2, byrow = TRUE)
colnames(beta) = c("x", "y")
rownames(beta) = c("Intercept", "sex", "Z")

sex = numeric(nrow(a))
sex = sample(c(0, 1), length(sex), replace = TRUE)
sex[rownames(a) %in% ped$SIRE] <- 1
sex[rownames(a) %in% ped$DAM] <- 0
Z = rnorm(nrow(a))
Intercept = rep(1, nrow(a))
X = cbind(Intercept, sex, Z)
rownames(X) = rownames(a)

# Correlated noise
e = rmvnorm(nrow(a), sigma = E)

# Simulated data
Y = X %*% beta + a + e
colnames(Y) = c("x", "y")

# Create relationship matrix
inv.phylo <- MCMCglmm::inverseA(ped, scale = TRUE)
A <- solve(inv.phylo$Ainv)
A = (A + t(A))/2 # Not always symmetric after inversion
rownames(A) <- rownames(inv.phylo$Ainv)

# Relate individuals to columns of the A matrix
pos = unlist(lapply(rownames(Y), function(x) which(x == rownames(A))))

stan_data = list(K = ncol(Y),
                 J = ncol(X),
                 N = nrow(Y),
                 X = X,
                 Y = Y,
                 Z = pos,
                 A = as.matrix(chol(A)))
stan_model = stan(file = "./animalModel.stan", data = stan_data)
model = rstan::extract(stan_model)
rstan::summary(stan_model, pars = "G")[[1]]
rstan::traceplot(stan_model, pars = "G")
colMeans(model$G)tan_model, pars = "G")
colMeans(model$G)
