data {
  int<lower=1>    K; // number of traits
  int<lower=1>    J; // number of fixed effects
  int<lower=0>    N; // number of individuals
  array[N] vector[J] X;// Fixed effects design matrix
  array[N] vector[K] Y; // response variable
  matrix[N, N]    A; // known covariance matrix
  int<lower=1>  lkj_prior; // LKJ prior, controls shrinkage of the correlations
  vector<lower=0>[2] beta_prior; // h2 beta prior parameters
}
transformed data{
  matrix[N, N] LA;
  vector[K] y_sd;
  vector[K] y_var;
  vector[K] y_mean;
  array[N] vector[K] Y_std; // response variable
  vector[J] x_scale;
  vector[J] x_scale_inv;
  array[J] int x_is_constant;
  vector[J] x_constant_value;
  LA = cholesky_decompose(A);

  for(j in 1:J){
    vector[N] x_j;
    int is_indicator = 1;
    int is_constant = 1;
    for(n in 1:N){
      x_j[n] = X[n][j];
      if(X[n][j] != X[1][j])
        is_constant = 0;
      if(X[n][j] != 0 && X[n][j] != 1)
        is_indicator = 0;
    }
    x_is_constant[j] = is_constant;
    x_constant_value[j] = X[1][j];
    x_scale[j] = 1;
    if(is_indicator == 0 && sd(x_j) > 0)
      x_scale[j] = sd(x_j);
  }
  x_scale_inv = 1 ./ x_scale;

  for(k in 1:K){
    y_sd[k] = sd(Y[,k]);
    y_mean[k] = mean(Y[,k]);
    for(n in 1:N)
      Y_std[n,k] = (Y[n,k] - y_mean[k]) / y_sd[k];
  }
}
parameters {
  matrix[K, J]    beta_scaled; // fixed effects on scaled predictors
  matrix[N, K] a_tilde; // breeding values precursor
  vector<lower=0, upper = 1>[K] h2; //variance partition

// total variances
  vector<lower=0>[K] L_sigma;

// G matrix
  cholesky_factor_corr[K] L_Omega_G;

}
transformed parameters {
  matrix[N, K] a;
  matrix[K, J] beta; // fixed effects on original predictor scales
  vector<lower=0>[K] L_sigma_G;
  vector<lower=0>[K] L_sigma_R;

  for(k in 1:K){
    for(j in 1:J){
      beta[k,j] = y_sd[k] * beta_scaled[k,j] * x_scale_inv[j];
      if(x_is_constant[j] == 1 && x_constant_value[j] != 0)
        beta[k,j] += y_mean[k] / x_constant_value[j];
    }
  }
  for(k in 1:K){
    L_sigma_G[k] = L_sigma[k] * sqrt(h2[k]);
    L_sigma_R[k] = L_sigma[k] * sqrt(1-h2[k]);
  }
  a = (LA * a_tilde) * diag_pre_multiply(L_sigma_G, L_Omega_G)'; // a ~ N(0, A x G)
}
model {
    array[N] vector[K] mu;

    for(n in 1:N){
      mu[n] = beta_scaled * (X[n] ./ x_scale) + to_vector(a[n]);
      Y_std[n] ~ normal(mu[n], L_sigma_R);
    }

    to_vector(beta_scaled) ~ normal(0, 1);
    to_vector(a_tilde) ~ normal(0, 1);
    h2 ~ beta(beta_prior[1], beta_prior[2]);
    L_sigma ~ normal(0, 0.5);
    L_Omega_G ~ lkj_corr_cholesky(lkj_prior);
}
generated quantities {
    vector[K] sigma_P;
    vector[K] sigma_G;
    vector[K] sigma_R;
    cov_matrix[K] P;
    cov_matrix[K] G;
    cov_matrix[K] E;
    corr_matrix[K] corrG;

    sigma_G = y_sd .* L_sigma_G;
    sigma_R = y_sd .* L_sigma_R;

    G = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_G, L_Omega_G));
    E = diag_matrix(sigma_R);
    P = G + E;

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
}
