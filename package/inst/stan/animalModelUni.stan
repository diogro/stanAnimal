data {
  int<lower=1>    J; // number of fixed effects
  int<lower=0>    N; // number of individuals
  array[N] vector[J] X; // Fixed effects design matrix
  array[N] real Y; // response variable
  cov_matrix[N]   A; // known covariance matrix
}
transformed data{
  matrix[N, N] LA;
  real y_sd;
  real y_mean;
  vector[N] Y_std; // standardized response
  vector[J] x_scale;
  vector[J] x_scale_inv;
  array[J] int x_is_constant;
  vector[J] x_constant_value;
  LA = cholesky_decompose(A);
  y_sd = sd(Y);
  y_mean = mean(Y);
  for(n in 1:N)
    Y_std[n] = (Y[n] - y_mean) / y_sd;

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
}
parameters {
  vector[N]  a_tilde; // breeding values
  row_vector[J] beta_scaled; // fixed effects on scaled predictors
  real<lower=0, upper = 1> h2; // variance partition
  real<lower=0> sigma; // total residual + genetic sd on the standardized scale
}
transformed parameters {
  row_vector[J] beta; // fixed effects on original predictor scales
  for(j in 1:J){
    beta[j] = y_sd * beta_scaled[j] * x_scale_inv[j];
    if(x_is_constant[j] == 1 && x_constant_value[j] != 0)
      beta[j] += y_mean / x_constant_value[j];
  }
}
model {
    vector[N] mu;
    vector[N] a;
    
    a_tilde ~ normal(0, 1);
    a = sigma * sqrt(h2) * (LA * a_tilde);
 
    for(n in 1:N)
      mu[n] = beta_scaled * (X[n] ./ x_scale) + a[n];

    Y_std ~ normal(mu, sigma * sqrt(1 - h2));

    to_vector(beta_scaled) ~ normal(0, 1);
    h2 ~ beta(2, 4);
    sigma ~ normal(0, 0.5);

}
generated quantities{
  real sigma_E;
  real sigma_G;
  sigma_E = square(y_sd * sigma) * (1 - h2);
  sigma_G = square(y_sd * sigma) * h2;
}
