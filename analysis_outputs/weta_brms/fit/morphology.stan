// generated with brms 2.21.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_fore;  // number of observations
  vector[N_fore] Y_fore;  // response variable
  int<lower=1> K_fore;  // number of population-level effects
  matrix[N_fore, K_fore] X_fore;  // population-level design matrix
  int<lower=1> N_mid;  // number of observations
  vector[N_mid] Y_mid;  // response variable
  int<lower=1> K_mid;  // number of population-level effects
  matrix[N_mid, K_mid] X_mid;  // population-level design matrix
  int<lower=1> N_hind;  // number of observations
  vector[N_hind] Y_hind;  // response variable
  int<lower=1> K_hind;  // number of population-level effects
  matrix[N_hind, K_hind] X_hind;  // population-level design matrix
  int<lower=1> N_earlin;  // number of observations
  vector[N_earlin] Y_earlin;  // response variable
  int<lower=1> K_earlin;  // number of population-level effects
  matrix[N_earlin, K_earlin] X_earlin;  // population-level design matrix
  int<lower=1> N_eye;  // number of observations
  vector[N_eye] Y_eye;  // response variable
  int<lower=1> K_eye;  // number of population-level effects
  matrix[N_eye, K_eye] X_eye;  // population-level design matrix
  int<lower=1> nresp;  // number of responses
  int nrescor;  // number of residual correlations
  int prior_only;  // should the likelihood be ignored?
  real prior_scale;
}
transformed data {
  array[N] vector[nresp] Y;  // response array
  for (n in 1:N) {
    Y[n] = transpose([Y_fore[n], Y_mid[n], Y_hind[n], Y_earlin[n], Y_eye[n]]);
  }
}
parameters {
  vector[K_fore] b_fore;  // regression coefficients
  real<lower=0> sigma_fore;  // dispersion parameter
  vector[K_mid] b_mid;  // regression coefficients
  real<lower=0> sigma_mid;  // dispersion parameter
  vector[K_hind] b_hind;  // regression coefficients
  real<lower=0> sigma_hind;  // dispersion parameter
  vector[K_earlin] b_earlin;  // regression coefficients
  real<lower=0> sigma_earlin;  // dispersion parameter
  vector[K_eye] b_eye;  // regression coefficients
  real<lower=0> sigma_eye;  // dispersion parameter
  cholesky_factor_corr[nresp] Lrescor;  // parameters for multivariate linear models
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(b_fore[1] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_fore[2] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_fore[3] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_fore[4] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_fore[5] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_fore[6] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_fore[7] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_fore[8] | 0, 1 * prior_scale);
  lprior += student_t_lpdf(sigma_fore | 3, 0, 0.15 * prior_scale)
    - 1 * student_t_lccdf(0 | 3, 0, 0.15 * prior_scale);
  lprior += normal_lpdf(b_mid[1] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_mid[2] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_mid[3] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_mid[4] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_mid[5] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_mid[6] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_mid[7] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_mid[8] | 0, 1 * prior_scale);
  lprior += student_t_lpdf(sigma_mid | 3, 0, 0.15 * prior_scale)
    - 1 * student_t_lccdf(0 | 3, 0, 0.15 * prior_scale);
  lprior += normal_lpdf(b_hind[1] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_hind[2] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_hind[3] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_hind[4] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_hind[5] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_hind[6] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_hind[7] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_hind[8] | 0, 1 * prior_scale);
  lprior += student_t_lpdf(sigma_hind | 3, 0, 0.15 * prior_scale)
    - 1 * student_t_lccdf(0 | 3, 0, 0.15 * prior_scale);
  lprior += normal_lpdf(b_earlin[1] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_earlin[2] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_earlin[3] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_earlin[4] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_earlin[5] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_earlin[6] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_earlin[7] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_earlin[8] | 0, 1 * prior_scale);
  lprior += student_t_lpdf(sigma_earlin | 3, 0, 0.15 * prior_scale)
    - 1 * student_t_lccdf(0 | 3, 0, 0.15 * prior_scale);
  lprior += normal_lpdf(b_eye[1] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_eye[2] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_eye[3] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_eye[4] | 0, 0.5 * prior_scale);
  lprior += normal_lpdf(b_eye[5] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_eye[6] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_eye[7] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_eye[8] | 0, 1 * prior_scale);
  lprior += student_t_lpdf(sigma_eye | 3, 0, 0.15 * prior_scale)
    - 1 * student_t_lccdf(0 | 3, 0, 0.15 * prior_scale);
  lprior += lkj_corr_cholesky_lpdf(Lrescor | 2);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N_fore] mu_fore = rep_vector(0.0, N_fore);
    // initialize linear predictor term
    vector[N_mid] mu_mid = rep_vector(0.0, N_mid);
    // initialize linear predictor term
    vector[N_hind] mu_hind = rep_vector(0.0, N_hind);
    // initialize linear predictor term
    vector[N_earlin] mu_earlin = rep_vector(0.0, N_earlin);
    // initialize linear predictor term
    vector[N_eye] mu_eye = rep_vector(0.0, N_eye);
    // multivariate predictor array
    array[N] vector[nresp] Mu;
    vector[nresp] sigma = transpose([sigma_fore, sigma_mid, sigma_hind, sigma_earlin, sigma_eye]);
    // cholesky factor of residual covariance matrix
    matrix[nresp, nresp] LSigma = diag_pre_multiply(sigma, Lrescor);
    mu_fore += X_fore * b_fore;
    mu_mid += X_mid * b_mid;
    mu_hind += X_hind * b_hind;
    mu_earlin += X_earlin * b_earlin;
    mu_eye += X_eye * b_eye;
    // combine univariate parameters
    for (n in 1:N) {
      Mu[n] = transpose([mu_fore[n], mu_mid[n], mu_hind[n], mu_earlin[n], mu_eye[n]]);
    }
    target += multi_normal_cholesky_lpdf(Y | Mu, LSigma);
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // residual correlations
  corr_matrix[nresp] Rescor = multiply_lower_tri_self_transpose(Lrescor);
  vector<lower=-1,upper=1>[nrescor] rescor;
  // extract upper diagonal of correlation matrix
  for (k in 1:nresp) {
    for (j in 1:(k - 1)) {
      rescor[choose(k - 1, 2) + j] = Rescor[j, k];
    }
  }
}

