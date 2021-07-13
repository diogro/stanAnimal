data {
  int<lower=1>    K; // number of traits
  int<lower=1>    J; // number of fixed effects
  int<lower=0>    N; // number of individuals
  vector[J]    X[N]; // Fixed effects design matrix
  vector[K]    Y[N]; // response variable
  matrix[N, N]    A; // known covariance matrix
  int<lower=1>  lkj_prior; // LKJ prior, controls shrinkage of the correlations
  vector<lower=0>[2] beta_prior; // h2 beta prior parameters
}
transformed data{
  matrix[N, N] LA;
  vector[K] y_sd;
  vector[K] y_var;
  vector[K] y_mean;
  vector[K] Y_std[N]; // response variable
  LA = cholesky_decompose(A);
  
  for(k in 1:K){
    y_sd[k] = sd(Y[,k]);
    y_mean[k] = mean(Y[,k]);
    for(n in 1:N)
      Y_std[n,k] = (Y[n,k] - y_mean[k]) / y_sd[k];
  }
}
parameters {
  matrix[K, J]    beta; // fixed effects
  matrix[N, K] a_tilde; // breeding values precursor
  vector<lower=0, upper = 1>[K] h2; //variance partition

// total variances
  vector<lower=0>[K] L_sigma;

// G matrix
  cholesky_factor_corr[K] L_Omega_G;

// R matrix
  cholesky_factor_corr[K] L_Omega_R;

}
transformed parameters {
  matrix[N, K] a;
  vector<lower=0>[K] L_sigma_G;
  vector<lower=0>[K] L_sigma_R;

  for(k in 1:K){
    L_sigma_G[k] = L_sigma[k] *    h2[k];
    L_sigma_R[k] = L_sigma[k] * (1-h2[k]);
  }
  a = diag_pre_multiply(L_sigma_G, L_Omega_G) * a_tilde * LA';  // a ~ N(0, A x G)
  //a = (LA * a_tilde) * diag_pre_multiply(L_sigma_G, L_Omega_G)'; //old
}
model {
    vector[K] mu[N];
    matrix[K, K] L_Sigma_R;

    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    for(n in 1:N)
      mu[n] = beta * X[n] + to_vector(a[n]);

    Y_std ~ multi_normal_cholesky(mu, L_Sigma_R);

    to_vector(beta) ~ normal(0, 1);
    to_vector(a_tilde) ~ normal(0, 1);
    h2 ~ beta(beta_prior[1], beta_prior[2]);
    L_sigma ~ normal(0, 0.5);
    L_Omega_G ~ lkj_corr_cholesky(lkj_prior);
    L_Omega_R ~ lkj_corr_cholesky(lkj_prior);
}
generated quantities {
    vector[K] sigma_P;
    vector[K] sigma_G;
    vector[K] sigma_R;
    cov_matrix[K] P;
    cov_matrix[K] G;
    cov_matrix[K] E;
    corr_matrix[K] corrG;
    corr_matrix[K] corrE;
  
    sigma_P = y_sd .* L_sigma;
    sigma_G = y_sd .* L_sigma_G;
    sigma_R = y_sd .* L_sigma_R;

    G = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_G, L_Omega_G));
    E = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_R, L_Omega_R));
    P = G + E;

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
    corrE = multiply_lower_tri_self_transpose(L_Omega_R);
}
