functions {
  matrix as_matrix(vector X, int N, int K) { 
    matrix[N, K] Y; 
    for (i in 1:N) {
      Y[i] = to_row_vector(X[((i - 1) * K + 1):(i * K)]); 
    }
    return Y; 
  }
  vector chol_kronecker_product(matrix LA, matrix LG, vector a) {
    vector[num_elements(a)] new_a;
    new_a = rep_vector(0, num_elements(a));
    for(iA in 1:cols(LA)){
      for(jA in 1:iA){
        if(LA[iA, jA] > 10e-10){
          for(iG in 1:cols(LG)){
            for(jG in 1:iG){
              new_a[(cols(LG)*(iA-1))+iG] = new_a[(cols(LG)*(iA-1))+iG] + LA[iA, jA] * LG[iG, jG] * a[(cols(LG)*(jA-1))+jG];
            }
          }
        }
      }
    }
    return new_a;
  }
}
data {
  int<lower=1>    K; // number of traits
  int<lower=1>    J; // number of fixed effects
  int<lower=0>    N; // number of individuals
  vector[J]    X[N]; // Fixed effects design matrix
  vector[K]    Y[N]; // response variable
  matrix[N, N]    A; // known covariance matrix
}
transformed data{
  matrix[N, N] LA;
  LA = cholesky_decompose(A);
}
parameters {
  matrix[K,J]    beta; // fixed effects
  vector[N*K] a_tilde; // breeding values

# G matrix
  cholesky_factor_corr[K] L_Omega_G;
  vector<lower=0>[K] L_sigma_G;

# R matrix
  cholesky_factor_corr[K] L_Omega_R;
  vector<lower=0>[K] L_sigma_R;

}
transformed parameters {
  matrix[N, K] a;
  a = as_matrix(chol_kronecker_product(LA, diag_pre_multiply(L_sigma_G, L_Omega_G), a_tilde), N, K);
}
model {
    vector[K] mu[N];
    matrix[K, K] L_Sigma_R;

    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    for(n in 1:N)
      mu[n] = beta * X[n] + to_vector(a[n]);

    Y ~ multi_normal_cholesky(mu, L_Sigma_R);

    to_vector(beta) ~ normal(0, 1);
    a_tilde   ~ normal(0, 1);
    L_Omega_G ~ lkj_corr_cholesky(4);
    L_sigma_G ~ cauchy(0, 2.5);
    L_Omega_R ~ lkj_corr_cholesky(4);
    L_sigma_R ~ cauchy(0, 2.5);
}
generated quantities {
    cov_matrix[K] P;
    cov_matrix[K] G;
    cov_matrix[K] E;
    corr_matrix[K] corrG;
    corr_matrix[K] corrE;

    G = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_G, L_Omega_G));
    E = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R, L_Omega_R));
    P = G + E;

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
    corrE = multiply_lower_tri_self_transpose(L_Omega_R);
}
