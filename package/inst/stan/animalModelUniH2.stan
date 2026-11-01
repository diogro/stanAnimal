// Univariate animal model parametrized by total variance and h2.
// Breeding values are marginalized: with A = Q diag(lambda) Q', the rotated
// records are independent, y*_i ~ N(X*_i beta, sigma^2 (lambda_i h2 + 1 - h2)).
data {
  int<lower=1> J; // number of fixed effects
  int<lower=2> N; // number of individuals
  matrix[N, J] X; // fixed effects design matrix
  vector[N] Y; // response variable
  matrix[N, N] A; // known relationship matrix (symmetric, PSD)
  vector<lower=0>[2] beta_prior; // beta prior on h2; c(1, 1) is flat
  vector[J] x_mean_weight; // X * x_mean_weight is a column of ones
}
transformed data {
  vector[N] lambda;
  matrix[N, N] Q = eigenvectors_sym(A);
  real y_sd = sd(Y);
  real y_mean = mean(Y);
  vector[J] x_center;
  vector[J] x_scale;
  vector[J] x_scale_inv;
  vector[J] beta_prior_sd;
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
      x_center[j] = 0;
      x_scale[j] = 1;
      if (is_indicator == 0 && is_constant == 0) {
        x_center[j] = mean(X[, j]);
        x_scale[j] = sd(X[, j]);
      }
      x_scale_inv[j] = 1 / x_scale[j];
      beta_prior_sd[j] = 1;
      if (is_constant == 1 && x_constant_value[j] != 0)
        beta_prior_sd[j] = 2.5 / abs(x_constant_value[j]);
      X_std[, j] = (X[, j] - x_center[j]) * x_scale_inv[j];
    }
    Ys = Q' * ((Y - y_mean) / y_sd);
    Xs = Q' * X_std;
  }
}
parameters {
  vector[J] beta_scaled; // fixed effects on scaled predictors
  real<lower=0> s_total; // total sd, standardized response scale
  real<lower=0, upper=1> h2; // heritability
}
transformed parameters {
  vector[J] beta; // fixed effects on original predictor and response scales
  for (j in 1:J) {
    beta[j] = y_sd * beta_scaled[j] * x_scale_inv[j];
  }
  // The response is centred and scaled, so the model predicts deviations from
  // its mean. Put the mean back where the design can hold it: X x_mean_weight is
  // a column of ones, so the shift lands on every prediction, whether the design
  // has an intercept column or the indicators of a factor.
  {
    real mean_shift = y_mean - dot_product(x_center, beta);
    for (j in 1:J)
      beta[j] += mean_shift * x_mean_weight[j];
  }
}
model {
  Ys ~ normal(Xs * beta_scaled, s_total * sqrt(lambda * h2 + (1 - h2)));

  for (j in 1:J)
    beta_scaled[j] ~ normal(0, beta_prior_sd[j]);
  s_total ~ normal(0, 1);
  h2 ~ beta(beta_prior[1], beta_prior[2]);
}
generated quantities {
  real sigma_G = y_sd * s_total * sqrt(h2);
  real sigma_E = y_sd * s_total * sqrt(1 - h2);
}
