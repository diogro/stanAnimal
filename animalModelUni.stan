data {
  int<lower=1>    J; // number of fixed effects
  int<lower=0>    N; // number of individuals
  vector[J]    X[N]; // Fixed effects design matrix
  real         Y[N]; // response variable
  matrix[N, N]    A; // cholesky factor of known covariance matrix
}
parameters {
  row_vector[J] beta; // fixed effects
  vector[N]    a; // breeding values

# Genetic variance
  real<lower=0> sigma_G;

# Residual variance
  real<lower=0> sigma_R;

}
transformed parameters {
  vector[N] aM;
  aM = (sqrt(sigma_G) * A) * a;
}
model {
    vector[N] mu;
 
    for(n in 1:N)
      mu[n] = beta * X[n] + aM[n];

    Y ~ normal(mu, sigma_R);

    to_vector(beta) ~ normal(0, 1);
    a ~ normal(0, 1);
    sigma_G ~ normal(0, 1);
    sigma_R ~ normal(0, 1);
}
