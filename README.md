# Stan animal model

Animal model implementation in Stan.

Somewhat experimental, univariate works better, but multivariate is pretty good too.

# Install

```r
if (!requireNamespace("pak", quietly = TRUE))
    install.packages("pak")
pak::pkg_install(c("rhdf5", "stan-dev/cmdstanr"))
cmdstanr::install_cmdstan() 
install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
pak::pkg_install("diogro/stanAnimal/package")
```

# Usage

```r
# Y trait matrix
# X fixed effect model matrix
# A relatedness matrix
A = nadiv::makeA(pedigree)
lmm_animal(Y, X, A)
```