# Stan animal model

Animal model implementation in Stan
Very experimental, univariate works better.

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