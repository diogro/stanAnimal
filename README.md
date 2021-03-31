# Stan animal model

Animal model implementation in Stan.

Somewhat experimental, univariate works better, but multivariate is pretty good too.

# Install

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rhdf5")
install.packages("remotes")
remotes::install_github("stan-dev/cmdstanr")
cmdstanr::install_cmdstan() 
install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
remotes::install_github("diogro/stanAnimal", subdir = "package")
```

# Usage

```r
# Y trait matrix
# X fixed effect model matrix
# A relatedness matrix
A = nadiv::makeA(pedigree)
lmm_animal(Y, X, A)
```