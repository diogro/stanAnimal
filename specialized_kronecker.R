if(!require(microbenchmark)){install.packages("microbenchmark"); library(microbenchmark)}

ped = data.frame(ID = 1:14,
                 SIRE = c(NA, NA, NA, NA, NA, 1, 3, 3, NA, NA, 4, 6, 8, 10), 
                 DAM  = c(NA, NA, NA, NA, NA, 2, 2, 2, NA, NA, 5, 5, 9, 9 ))


inv.phylo <- MCMCglmm::inverseA(ped, scale = TRUE)
A <- solve(inv.phylo$Ainv)
A = (A + t(A))/2 # Not always symmetric after inversion
rownames(A) <- rownames(inv.phylo$Ainv)
A[A < 1e-10] = 0
LA = as.matrix(t(chol(A)))
G = matrix(c(1, 0.8, 0.8, 1), 2, 2)
LG = t(chol(G))
a = rnorm(28)

kroneker_product <- function(LA, LG, a){
  nA = nrow(LA)
  nG = nrow(LG)
  new_a = numeric(length(a))
  for(iA in 1:nA){
    for(jA in 1:iA){
      if(LA[iA, jA] > 10e-10){
        for(iG in 1:nG){
          for(jG in 1:iG){
            new_a[(nG*(iA-1))+iG] = new_a[(nG*(iA-1))+iG] + LA[iA, jA] * LG[iG, jG] * a[(nG*(jA-1))+jG]
          }
        }
      }
    }
  }
  return(new_a)
}

asw = matrix((kronecker(LA, LG) %*% a)[,1], 14, 2, byrow = T)
new_a = matrix(kroneker_product(LA, LG, a), 14, 2, byrow = T)
fst_a = matrix(LA %*% matrix(a, 14, 2, byrow = T) %*% t(LG), 14, 2)
microbenchmark(asw = matrix((kronecker(LA, LG) %*% a)[,1], 14, 2, byrow = T),
               new_a = matrix(kroneker_product(LA, LG, a), 14, 2, byrow = T),
               fst_a = LA %*% matrix(a, 14, 2, byrow = T) %*% t(LG))
all.equal(asw, new_a)
all.equal(asw, fst_a)
