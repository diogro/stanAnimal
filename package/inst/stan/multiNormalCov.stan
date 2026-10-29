data {
  int<lower=1>    K; // number of traits
  int<lower=1>    J; // number of fixed effects
  int<lower=0>    N; // number of individuals
  array[N] vector[J] X; // Fixed effects design matrix
  array[N] vector[K] Y; // response variable
  int<lower=1>    c; // LKJ parameter
}

transformed data{
  matrix[N, N] LA;
  vector[K] y_sd;
  vector[K] y_var;
  vector[K] y_mean;
  array[N] vector[K] Y_std; // response variable

  for(k in 1:K){
    y_sd[k] = sd(Y[,k])/2;
    y_mean[k] = mean(Y[,k]);
    for(n in 1:N)
      Y_std[n,k] = (Y[n,k] - y_mean[k]) / y_sd[k];
  }
}
parameters {
  matrix[K,J] beta; // fixed effects

// R matrix
  cholesky_factor_corr[K] L_Omega_R;
  vector<lower=0>[K] L_sigma_R;

}

model {
    array[N] vector[K] mu;
    matrix[K, K] L_Sigma_R;

    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    for(n in 1:N)
      mu[n] = beta * X[n];

    Y_std ~ multi_normal_cholesky(mu, L_Sigma_R);

    to_vector(beta) ~ normal(0, 1);

    L_Omega_R ~ lkj_corr_cholesky(c);
    L_sigma_R ~ normal(0, 1);
}
generated quantities {
    vector[K] sigma_R;
    cov_matrix[K] P;
    corr_matrix[K] corrP;
  
    sigma_R = y_sd .* L_sigma_R;

    P = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_R, L_Omega_R));

    corrP = multiply_lower_tri_self_transpose(L_Omega_R);
}
