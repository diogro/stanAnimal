// Multivariate animal model parametrized directly by G and E.
// Breeding values are marginalized: with A = Q diag(lambda) Q', the rotated
// records are independent, y*_i ~ MVN(X*_i beta, lambda_i G + E).
//
// The N per-individual covariances are never factorized one at a time. With
// E = L_E L_E' and C = L_E^{-1} L_G (lower triangular), for any lambda
//   lambda G + E = L_E (I + lambda C C') L_E',
// and writing C C' = U diag(d) U' gives, for every individual,
//   log det(lambda G + E) = log det(E) + sum_k log(1 + lambda d_k)
//   r' (lambda G + E)^{-1} r = sum_k (u_k' L_E^{-1} r)^2 / (1 + lambda d_k).
// Whitening the residuals by L_E and rotating them once by U therefore costs
// O(K) per individual, so an evaluation is O(N K^2 + K^3) instead of the
// O(N K^3) of a Cholesky decomposition for each of the N individuals. When the
// eigenvalues of E^{-1/2} G E^{-1/2} are nearly coincident (where eigenvector
// gradients are undefined) the model falls back to the per-individual
// Cholesky route, which has a well-behaved gradient.
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
  // Cholesky factors and covariances of G and E.
  matrix[K, K] L_G = diag_pre_multiply(s_G, L_Omega_G);
  matrix[K, K] L_E;
  if (unstructured_E)
    L_E = diag_pre_multiply(s_E, L_Omega_E);
  else
    L_E = diag_matrix(s_E);

  // log det E uses the true diagonal. The diagonal is floored only for the
  // whitening below, so that an extreme HMC proposal (e.g. a residual SD that
  // underflows to zero) yields a finite, hugely negative log density instead
  // of NaN/Inf, which would abort the chain.
  vector[K] log_diag_E = log(diagonal(L_E));
  for (k in 1:K)
    L_E[k, k] = fmax(L_E[k, k], 1e-12);

  matrix[K, K] G = L_G * L_G';
  matrix[K, K] E = L_E * L_E';

  // C = L_E^{-1} L_G is lower triangular, so C C' = L_E^{-1} G L_E^{-T}.
  matrix[K, K] C = mdivide_left_tri_low(L_E, L_G);
  matrix[K, K] M = C * C';
  M = (M + M') / 2; // symmetrize against rounding for the eigendecomposition
  for (i in 1:K)
    for (j in 1:K)
      M[i, j] = fmin(fmax(M[i, j], -1e150), 1e150); // guard overflow in the tails
  vector[K] d = eigenvalues_sym(M);
  for (k in 1:K)
    d[k] = fmax(d[k], 0); // M is PSD; guard rounding before log()

  // Eigenvector gradients are only well defined when the eigenvalues are well
  // separated. Near-coincident eigenvalues (the identity starting point, or G
  // proportional to E) fall back to the per-individual Cholesky route.
  real gap = positive_infinity();
  for (k in 2:K)
    gap = fmin(gap, d[k] - d[k - 1]);

  if (gap > 1e-4 * fmax(d[K], 1.0)) {
    // Fast path: rotate once into the eigenvectors of E^{-1/2} G E^{-1/2}.
    // Column n of W is L_E^{-1} r_n, and Z[n, k] = u_k' L_E^{-1} r_n.
    matrix[K, K] U = eigenvectors_sym(M);
    matrix[N, K] resid = Ys - Xs * beta_scaled';
    matrix[K, N] W = mdivide_left_tri_low(L_E, resid');
    matrix[N, K] Z = W' * U;

    // denom[n, k] = 1 + lambda_n d_k
    matrix[N, K] denom = rep_matrix(1, N, K) + lambda * d';
    target += -0.5 * (N * K * log(2 * pi())
                      + 2 * N * sum(log_diag_E)
                      + sum(log(denom))
                      + sum(square(Z) ./ denom));
  } else {
    // Robust fallback: per-individual Cholesky of lambda_n G + E.
    matrix[N, K] mu = Xs * beta_scaled';
    for (n in 1:N)
      target += multi_normal_cholesky_lpdf(Ys[n] | mu[n],
                                           cholesky_decompose(lambda[n] * G + E));
  }

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
