// Multivariate animal model parametrized directly by G and E.
// Breeding values are marginalized: with A = Q diag(lambda) Q', the rotated
// records are independent, y*_i ~ MVN(X*_i beta, lambda_i G + E).
data {
  int<lower=1> K; // number of traits
  int<lower=1> J; // number of fixed effects
  int<lower=2> N; // number of individuals
  matrix[N, J] X; // fixed effects design matrix
  matrix[N, K] Y; // response variables
  matrix[N, N] A; // known relationship matrix (symmetric, PSD)
  real<lower=0> lkj_prior; // LKJ shape for the correlation matrices
  int<lower=0, upper=1> unstructured_E; // 1: full residual covariance, 0: diagonal
}
transformed data {
  int E_dim = unstructured_E ? K : 0;
  vector[N] lambda;
  matrix[N, N] Q = eigenvectors_sym(A);
  vector[K] y_sd;
  vector[K] y_mean;
  vector[J] x_center;
  vector[J] x_scale;
  vector[J] x_scale_inv;
  vector[J] beta_prior_sd;
  array[J] int x_is_constant;
  vector[J] x_constant_value;
  matrix[N, K] Ys;
  matrix[N, J] Xs;

  {
    vector[N] ev = eigenvalues_sym(A);
    for (n in 1:N)
      lambda[n] = fmax(ev[n], 0);
  }

  {
    matrix[N, K] Y_std;
    matrix[N, J] X_std;
    for (k in 1:K) {
      y_sd[k] = sd(Y[, k]);
      y_mean[k] = mean(Y[, k]);
      Y_std[, k] = (Y[, k] - y_mean[k]) / y_sd[k];
    }
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
    Ys = Q' * Y_std;
    Xs = Q' * X_std;
  }
}
parameters {
  matrix[K, J] beta_scaled; // fixed effects on scaled predictors
  vector<lower=0>[K] s_G; // genetic sds, standardized response scale
  vector<lower=0>[K] s_E; // residual sds, standardized response scale
  cholesky_factor_corr[K] L_Omega_G;
  cholesky_factor_corr[E_dim] L_Omega_E;
}
transformed parameters {
  matrix[K, J] beta; // fixed effects on original predictor and response scales
  for (k in 1:K) {
    for (j in 1:J) {
      beta[k, j] = y_sd[k] * beta_scaled[k, j] * x_scale_inv[j];
    }
  }
  for (k in 1:K) {
    for (j in 1:J) {
      if (x_is_constant[j] == 1 && x_constant_value[j] != 0) {
        real center_shift = 0;
        for (m in 1:J)
          center_shift += x_center[m] * beta[k, m];
        beta[k, j] += (y_mean[k] - center_shift) / x_constant_value[j];
      }
    }
  }
}
model {
  matrix[K, K] G = multiply_lower_tri_self_transpose(diag_pre_multiply(s_G, L_Omega_G));
  matrix[K, K] E;
  matrix[N, K] mu = Xs * beta_scaled';

  if (unstructured_E)
    E = multiply_lower_tri_self_transpose(diag_pre_multiply(s_E, L_Omega_E));
  else
    E = diag_matrix(square(s_E));

  for (n in 1:N)
    Ys[n] ~ multi_normal_cholesky(mu[n], cholesky_decompose(lambda[n] * G + E));

  for (k in 1:K)
    for (j in 1:J)
      beta_scaled[k, j] ~ normal(0, beta_prior_sd[j]);
  s_G ~ normal(0, 1);
  s_E ~ normal(0, 1);
  L_Omega_G ~ lkj_corr_cholesky(lkj_prior);
  if (unstructured_E)
    L_Omega_E ~ lkj_corr_cholesky(lkj_prior);
}
generated quantities {
  matrix[K, K] G; // additive genetic covariance, original response scale
  matrix[K, K] E; // residual covariance, original response scale
  matrix[K, K] P; // phenotypic covariance
  corr_matrix[K] corrG;
  vector[K] sigma_G = y_sd .* s_G;
  vector[K] sigma_E = y_sd .* s_E;
  vector[K] h2;
  {
    matrix[K, K] G_std = multiply_lower_tri_self_transpose(diag_pre_multiply(s_G, L_Omega_G));
    matrix[K, K] E_std;
    if (unstructured_E)
      E_std = multiply_lower_tri_self_transpose(diag_pre_multiply(s_E, L_Omega_E));
    else
      E_std = diag_matrix(square(s_E));
    G = quad_form_diag(G_std, y_sd);
    E = quad_form_diag(E_std, y_sd);
  }
  P = G + E;
  corrG = multiply_lower_tri_self_transpose(L_Omega_G);
  for (k in 1:K)
    h2[k] = G[k, k] / P[k, k];
}
