// Multivariate linear model for the phenotypic (residual) covariance matrix P,
// for data with no relatedness or genomic information. Parametrized like the
// animal model: scales times an LKJ-distributed correlation matrix.
data {
  int<lower=1> K; // number of traits
  int<lower=1> J; // number of fixed effects
  int<lower=2> N; // number of individuals
  matrix[N, J] X; // fixed effects design matrix
  matrix[N, K] Y; // response variables
  real<lower=0> lkj_prior; // LKJ shape for the correlation matrix
}
transformed data {
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

  for (k in 1:K) {
    y_sd[k] = sd(Y[, k]);
    y_mean[k] = mean(Y[, k]);
    Ys[, k] = (Y[, k] - y_mean[k]) / y_sd[k];
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
    Xs[, j] = (X[, j] - x_center[j]) * x_scale_inv[j];
  }
}
parameters {
  matrix[K, J] beta_scaled; // fixed effects on scaled predictors
  vector<lower=0>[K] s_P; // phenotypic sds, standardized response scale
  cholesky_factor_corr[K] L_Omega_P;
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
  // Cholesky factor of P on the standardized response scale.
  matrix[K, K] L_P = diag_pre_multiply(s_P, L_Omega_P);

  // log det P uses the true diagonal. The diagonal is floored only for the
  // whitening below, so that an extreme proposal (a correlation at the
  // boundary, or a scale underflowing to zero) yields a finite, hugely
  // negative log density instead of NaN/Inf, which would abort the chain.
  vector[K] log_diag_P = log(diagonal(L_P));
  for (k in 1:K)
    L_P[k, k] = fmax(L_P[k, k], 1e-12);

  // Whitened residuals: column n of W is L_P^{-1} r_n, so sum(square(W)) is
  // the quadratic form of the multivariate normal log density. Computing it
  // with one triangular solve costs O(N K^2) and avoids a per-individual loop.
  matrix[N, K] resid = Ys - Xs * beta_scaled';
  matrix[K, N] W = mdivide_left_tri_low(L_P, resid');
  target += -0.5 * (N * K * log(2 * pi())
                    + 2 * N * sum(log_diag_P)
                    + sum(square(W)));

  for (k in 1:K)
    for (j in 1:J)
      beta_scaled[k, j] ~ normal(0, beta_prior_sd[j]);
  s_P ~ normal(0, 1);
  L_Omega_P ~ lkj_corr_cholesky(lkj_prior);
}
generated quantities {
  vector[K] sigma_P = y_sd .* s_P; // phenotypic sds, original response scale
  matrix[K, K] P = quad_form_diag(
    multiply_lower_tri_self_transpose(diag_pre_multiply(s_P, L_Omega_P)), y_sd);
  corr_matrix[K] corrP = multiply_lower_tri_self_transpose(L_Omega_P);
}
