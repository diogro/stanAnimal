// Univariate animal model parametrized by the genetic and residual variances.
// Breeding values are marginalized: with A = Q diag(lambda) Q', the rotated
// records are independent, y*_i ~ N(X*_i beta, lambda_i sigma_G^2 + sigma_E^2).
data {
  int<lower=1> J; // number of fixed effects
  int<lower=2> N; // number of individuals
  matrix[N, J] X; // fixed effects design matrix
  vector[N] Y; // response variable
  matrix[N, N] A; // known relationship matrix (symmetric, PSD)
}
transformed data {
  vector[N] lambda;
  matrix[N, N] Q = eigenvectors_sym(A);
  real y_sd = sd(Y);
  real y_mean = mean(Y);
  vector[J] x_scale;
  vector[J] x_scale_inv;
  array[J] int x_is_constant;
  vector[J] x_constant_value;
  vector[N] Ys;
  matrix[N, J] Xs;

  {
    vector[N] ev = eigenvalues_sym(A);
    for (n in 1:N)
      lambda[n] = fmax(ev[n], 0);
  }

  {
    matrix[N, J] X_std;
    for (j in 1:J) {
      int is_indicator = 1;
      int is_constant = 1;
      for (n in 1:N) {
        if (X[n, j] != X[1, j])
          is_constant = 0;
        if (X[n, j] != 0 && X[n, j] != 1)
          is_indicator = 0;
      }
      x_is_constant[j] = is_constant;
      x_constant_value[j] = X[1, j];
      x_scale[j] = 1;
      if (is_indicator == 0 && sd(X[, j]) > 0)
        x_scale[j] = sd(X[, j]);
      x_scale_inv[j] = 1 / x_scale[j];
      X_std[, j] = X[, j] * x_scale_inv[j];
    }
    Ys = Q' * ((Y - y_mean) / y_sd);
    Xs = Q' * X_std;
  }
}
parameters {
  vector[J] beta_scaled; // fixed effects on scaled predictors
  real<lower=0> s_G; // genetic sd, standardized response scale
  real<lower=0> s_E; // residual sd, standardized response scale
}
transformed parameters {
  vector[J] beta; // fixed effects on original predictor and response scales
  for (j in 1:J) {
    beta[j] = y_sd * beta_scaled[j] * x_scale_inv[j];
    if (x_is_constant[j] == 1 && x_constant_value[j] != 0)
      beta[j] += y_mean / x_constant_value[j];
  }
}
model {
  Ys ~ normal(Xs * beta_scaled, sqrt(lambda * square(s_G) + square(s_E)));

  beta_scaled ~ normal(0, 1);
  s_G ~ normal(0, 1);
  s_E ~ normal(0, 1);
}
generated quantities {
  real sigma_G = y_sd * s_G;
  real sigma_E = y_sd * s_E;
  real h2 = square(s_G) / (square(s_G) + square(s_E));
}
