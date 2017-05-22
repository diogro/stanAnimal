functions {
  matrix as_matrix(vector X, int N, int K) { 
    matrix[N, K] Y; 
    for (i in 1:N) {
      Y[i] = to_row_vector(X[((i - 1) * K + 1):(i * K)]); 
    }
    return Y; 
  }
  matrix kronecker(matrix A, matrix B) {
    matrix[rows(A)*rows(B), cols(A)*cols(B)] kron;
    for (i in 1:cols(A)) {
      for (j in 1:rows(A)) {
        kron[((j-1)*rows(B)+1):(j*rows(B)), ((i-1)*cols(B)+1):(i*cols(B))] = A[j,i] * B;
      }
    }
    return kron;
  }
}
data {
  int<lower=1>    K; // number of traits
  int<lower=1>    J; // number of fixed effects
  int<lower=0>    N; // number of individuals
  vector[J]    X[N]; // Fixed effects design matrix
  vector[K]    Y[N]; // response variable
  int          Z[N]; // Random effect design matrix
  matrix[N, N]    A; // cholesky factor of known covariance matrix
}
parameters {
  matrix[K,J] beta; // fixed effects
  vector[N*K]    a; // breeding values

# G matrix
  cholesky_factor_corr[K] L_Omega_G;
  vector<lower=0>[K] L_sigma_G;

# R matrix
  cholesky_factor_corr[K] L_Omega_R;
  vector<lower=0>[K] L_sigma_R;

}
transformed parameters {
  matrix[N, K] aM;
  aM = as_matrix(kronecker(A, diag_pre_multiply(L_sigma_G, L_Omega_G)) * a, N, K);
}
model {
    vector[K] mu[N];
    matrix[K,K] L_Sigma_R;

    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    for(n in 1:N)
      mu[n] = beta * X[n] + to_vector(aM[Z[n]]);

    Y ~ multi_normal_cholesky(mu, L_Sigma_R);

    to_vector(beta) ~ normal(0, 1);
    a ~ normal(0, 1);
    L_Omega_G ~ lkj_corr_cholesky(4);
    L_sigma_G ~ cauchy(0, 2.5);
    L_Omega_R ~ lkj_corr_cholesky(4);
    L_sigma_R ~ cauchy(0, 2.5);
}
generated quantities {
    matrix[K, K] P;
    matrix[K, K] G;
    matrix[K, K] E;
    corr_matrix[K] corrG;
    corr_matrix[K] corrE;

    G = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_G, L_Omega_G));
    E = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R, L_Omega_R));
    P = G + E;

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
    corrE = multiply_lower_tri_self_transpose(L_Omega_R);
}
