library(ape)
library(MCMCglmm)

# uncomment and set your path to the folder where you stored the data
my_path <- "./"
phylo <- ape::read.nexus(paste0(my_path, "phylo.nex"))
data_simple <- read.table(paste0(my_path, "data_simple.txt"), header = TRUE)
head(data_simple)

inv.phylo <- MCMCglmm::inverseA(phylo, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

library(brms)
model_simple <- brm(phen ~ cofactor + (1|phylo), data = data_simple, 
                    family = gaussian(), cov_ranef = list(phylo = A),
                    prior = c(prior(normal(0, 10), "b"),
                              prior(normal(0, 50), "Intercept"),
                              prior(student_t(3, 0, 20), "sd"),
                              prior(student_t(3, 0, 20), "sigma")))
                              
stancode(model_simple)

