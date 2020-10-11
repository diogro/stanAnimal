data {
  int<lower=1>    J; // number of fixed effects
  int<lower=0>    N; // number of individuals
  vector[J]    X[N]; // Fixed effects design matrix
  real         Y[N]; // response variable
  cov_matrix[N]   A; // known covariance matrix
}
transformed data{
  matrix[N, N] LA;
  real<lower=0> sigma;
  LA = cholesky_decompose(A);
  sigma = sd(Y) * sd(Y);
}
parameters {
  vector[N]  a_tilde; // breeding values
  row_vector[J] beta; // fixed effects
  real<lower=0, upper = 1> h2; // variance partition

// Total variance
  

}
model {
    vector[N] mu;
    vector[N] a;
    
    a_tilde ~ normal(0, 1);
    a = sqrt(sigma*h2) * (LA * a_tilde);
 
    for(n in 1:N)
      mu[n] = beta * X[n] + a[n];

    Y ~ normal(mu, sqrt(sigma*(1 - h2)));

    to_vector(beta) ~ normal(0, 1);
    h2 ~ beta(2, 4);

}
generated quantities{
  real sigma_E;
  real sigma_G;
  sigma_E = sigma*(1 - h2);
  sigma_G = sigma*(    h2);
}
