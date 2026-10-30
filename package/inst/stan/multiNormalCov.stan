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
  vector[J] x_scale;
  vector[J] x_scale_inv;
  array[J] int x_is_constant;
  vector[J] x_constant_value;

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
    y_sd[k] = sd(Y[,k])/2;
    y_mean[k] = mean(Y[,k]);
    for(n in 1:N)
      Y_std[n,k] = (Y[n,k] - y_mean[k]) / y_sd[k];
  }
}
parameters {
  matrix[K,J] beta_scaled; // fixed effects on scaled predictors

// R matrix
  cholesky_factor_corr[K] L_Omega_R;
  vector<lower=0>[K] L_sigma_R;

}
transformed parameters {
  matrix[K,J] beta; // fixed effects on original predictor scales
    for(k in 1:K){
      for(j in 1:J){
        beta[k,j] = y_sd[k] * beta_scaled[k,j] * x_scale_inv[j];
        if(x_is_constant[j] == 1 && x_constant_value[j] != 0)
          beta[k,j] += y_mean[k] / x_constant_value[j];
      }
    }
}

model {
    array[N] vector[K] mu;
    matrix[K, K] L_Sigma_R;

    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    for(n in 1:N)
      mu[n] = beta_scaled * (X[n] ./ x_scale);

    Y_std ~ multi_normal_cholesky(mu, L_Sigma_R);

    to_vector(beta_scaled) ~ normal(0, 1);

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
