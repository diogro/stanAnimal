data {
  int<lower=1>    J; // number of fixed effects
  int<lower=0>    N; // number of individuals
  vector[J]    X[N]; // Fixed effects design matrix
  real         Y[N]; // response variable
  cov_matrix[N]   A; // known covariance matrix
}
transformed data{
  matrix[N, N] LA;
  LA = cholesky_decompose(A);
}
parameters {
  vector[N]  a_tilde; // breeding values
  row_vector[J] beta; // fixed effects

// Genetic variance
  real<lower=0> sigma_G;

// Residual variance
  real<lower=0> sigma_R;
}
model {
    vector[N] mu;
    vector[N] a;
    
    a_tilde ~ normal(0, 1);
    a = sqrt(sigma_G) * (LA * a_tilde);
 
    for(n in 1:N)
      mu[n] = beta * X[n] + a[n];

    Y ~ normal(mu, sigma_R);

    to_vector(beta) ~ normal(0, 1);
    
    sigma_G ~ normal(0, 1);
    sigma_R ~ normal(0, 1);
}
generated quantities{
  real sigma_E;
  sigma_E = sigma_R * sigma_R;
}
