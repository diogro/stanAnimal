# Stan animal model

Animal model implementation in Stan.

Somewhat experimental, univariate works better, but multivariate is pretty good too.

# Install

```r
devtools::install_github("diogro/stanAnimal", subdir = "package")
```

# Usage

```r
# Y trait matrix
# X fixed effect model matrix
# K relatedness matrix
lmm_animal(Y, X, K)
```