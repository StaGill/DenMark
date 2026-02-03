//////////////////////////////////////////////////////////////
// Rstan program for the LMC  
// With log(y) included, V=AA^T reparametrized as LKJ prior 
// Exact GP Poisson-Poisson
/////////////////////////////////////////////////////////////
functions {
  // Matern 3/2 covariance and its Cholesky factor
  matrix matern32_chol(matrix D, real rho) {
    int N = rows(D);
    matrix[N, N] K;
    real kappa = sqrt(3) / rho;
    K = (1 + kappa * D) .* exp(-kappa * D);
    for (n in 1:N) K[n, n] += 1e-10; // jitter
    return cholesky_decompose(K);
  }
}

data {
  int<lower=0> N;
  array[N] int Y;
  array[N] int M;
  matrix[N, N] dist_matrix;
  real<lower=0> grid_area; 
}

parameters {
  real beta0;
  real beta1;
  //real beta_y;
  vector<lower=0.0001>[2] sigma;
  cholesky_factor_corr[2] L;
  real<lower=0> rho1;
  real<lower=0> rho2;
  vector[N] z1;
  vector[N] z2;
}

transformed parameters {
  matrix[2, 2] A = diag_pre_multiply(sigma, L);
  matrix[2, 2] R = multiply_lower_tri_self_transpose(L);
  real<lower=-1, upper=1> corr = R[2, 1];

  // compute log offset once
  real log_grid_area = log(grid_area);
  
  // Cholesky factors of Matern covariances
  matrix[N, N] L_Sigma1 = matern32_chol(dist_matrix, rho1);
  matrix[N, N] L_Sigma2 = matern32_chol(dist_matrix, rho2);

  // latent GPs
  vector[N] omega1 = L_Sigma1 * z1;
  vector[N] omega2 = L_Sigma2 * z2;

  // linear predictors
  vector[N] eta1 = beta0 + log_grid_area + A[1, 1] * omega1;
  vector[N] eta2 = beta1 + A[2, 1] * omega1 + A[2, 2] * omega2;

  // means
  //vector[N] lambda = exp(eta1);
  //vector[N] mu = exp(eta2 +  log(to_vector(Y) + 0.05));
}

model {
  // priors
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  //beta_y ~ normal(0, 2);
  sigma ~ std_normal();
  //sigma[1] ~ normal(1, 2);
  //sigma[2] ~ normal(1, 2);
  L ~ lkj_corr_cholesky(1.2);
  //rho1 ~ lognormal(log(0.5), 0.4);
  //rho2 ~ lognormal(log(0.3), 0.4);
  rho1 ~ inv_gamma(2, 0.5);
  rho2 ~ inv_gamma(2, 0.3);
  z1 ~ std_normal();
  z2 ~ std_normal();

  // likelihood
  //Y ~ poisson(lambda);
  //M ~ poisson(mu);
  
  Y ~ poisson_log(eta1);
  M ~ poisson_log(eta2 + log(to_vector(Y) + 0.05));
}

generated quantities {
  //vector[N] log_lik;
  //for (s in 1:N) {
  //  log_lik[s] = poisson_lpmf(Y[s] | lambda[s]) +
  //               poisson_lpmf(M[s] | mu[s]);
  //}
}
