// generated with brms 2.21.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_fore;  // number of observations
  vector[N_fore] Y_fore;  // response variable
  int<lower=1> K_fore;  // number of population-level effects
  matrix[N_fore, K_fore] X_fore;  // population-level design matrix
  int<lower=1> Kc_fore;  // number of population-level effects after centering
  int<lower=1> N_mid;  // number of observations
  vector[N_mid] Y_mid;  // response variable
  int<lower=1> K_mid;  // number of population-level effects
  matrix[N_mid, K_mid] X_mid;  // population-level design matrix
  int<lower=1> Kc_mid;  // number of population-level effects after centering
  int<lower=1> N_hind;  // number of observations
  vector[N_hind] Y_hind;  // response variable
  int<lower=1> K_hind;  // number of population-level effects
  matrix[N_hind, K_hind] X_hind;  // population-level design matrix
  int<lower=1> Kc_hind;  // number of population-level effects after centering
  int<lower=1> N_earlin;  // number of observations
  vector[N_earlin] Y_earlin;  // response variable
  int<lower=1> K_earlin;  // number of population-level effects
  matrix[N_earlin, K_earlin] X_earlin;  // population-level design matrix
  int<lower=1> Kc_earlin;  // number of population-level effects after centering
  int<lower=1> N_eye;  // number of observations
  vector[N_eye] Y_eye;  // response variable
  int<lower=1> K_eye;  // number of population-level effects
  matrix[N_eye, K_eye] X_eye;  // population-level design matrix
  int<lower=1> Kc_eye;  // number of population-level effects after centering
  int<lower=1> nresp;  // number of responses
  int nrescor;  // number of residual correlations
  int prior_only;  // should the likelihood be ignored?
  real prior_scale;
}
transformed data {
  matrix[N_fore, Kc_fore] Xc_fore;  // centered version of X_fore without an intercept
  vector[Kc_fore] means_X_fore;  // column means of X_fore before centering
  matrix[N_mid, Kc_mid] Xc_mid;  // centered version of X_mid without an intercept
  vector[Kc_mid] means_X_mid;  // column means of X_mid before centering
  matrix[N_hind, Kc_hind] Xc_hind;  // centered version of X_hind without an intercept
  vector[Kc_hind] means_X_hind;  // column means of X_hind before centering
  matrix[N_earlin, Kc_earlin] Xc_earlin;  // centered version of X_earlin without an intercept
  vector[Kc_earlin] means_X_earlin;  // column means of X_earlin before centering
  matrix[N_eye, Kc_eye] Xc_eye;  // centered version of X_eye without an intercept
  vector[Kc_eye] means_X_eye;  // column means of X_eye before centering
  array[N] vector[nresp] Y;  // response array
  for (i in 2:K_fore) {
    means_X_fore[i - 1] = mean(X_fore[, i]);
    Xc_fore[, i - 1] = X_fore[, i] - means_X_fore[i - 1];
  }
  for (i in 2:K_mid) {
    means_X_mid[i - 1] = mean(X_mid[, i]);
    Xc_mid[, i - 1] = X_mid[, i] - means_X_mid[i - 1];
  }
  for (i in 2:K_hind) {
    means_X_hind[i - 1] = mean(X_hind[, i]);
    Xc_hind[, i - 1] = X_hind[, i] - means_X_hind[i - 1];
  }
  for (i in 2:K_earlin) {
    means_X_earlin[i - 1] = mean(X_earlin[, i]);
    Xc_earlin[, i - 1] = X_earlin[, i] - means_X_earlin[i - 1];
  }
  for (i in 2:K_eye) {
    means_X_eye[i - 1] = mean(X_eye[, i]);
    Xc_eye[, i - 1] = X_eye[, i] - means_X_eye[i - 1];
  }
  for (n in 1:N) {
    Y[n] = transpose([Y_fore[n], Y_mid[n], Y_hind[n], Y_earlin[n], Y_eye[n]]);
  }
}
parameters {
  vector[Kc_fore] b_fore;  // regression coefficients
  real Intercept_fore;  // temporary intercept for centered predictors
  real<lower=0> sigma_fore;  // dispersion parameter
  vector[Kc_mid] b_mid;  // regression coefficients
  real Intercept_mid;  // temporary intercept for centered predictors
  real<lower=0> sigma_mid;  // dispersion parameter
  vector[Kc_hind] b_hind;  // regression coefficients
  real Intercept_hind;  // temporary intercept for centered predictors
  real<lower=0> sigma_hind;  // dispersion parameter
  vector[Kc_earlin] b_earlin;  // regression coefficients
  real Intercept_earlin;  // temporary intercept for centered predictors
  real<lower=0> sigma_earlin;  // dispersion parameter
  vector[Kc_eye] b_eye;  // regression coefficients
  real Intercept_eye;  // temporary intercept for centered predictors
  real<lower=0> sigma_eye;  // dispersion parameter
  cholesky_factor_corr[nresp] Lrescor;  // parameters for multivariate linear models
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(b_fore[1] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_fore[2] | 0, 0.25 * prior_scale);
  lprior += normal_lpdf(b_fore[3] | 0, 0.25 * prior_scale);
  lprior += normal_lpdf(Intercept_fore | 0, 0.5 * prior_scale);
  lprior += student_t_lpdf(sigma_fore | 3, 0, 0.15 * prior_scale)
    - 1 * student_t_lccdf(0 | 3, 0, 0.15 * prior_scale);
  lprior += normal_lpdf(b_mid[1] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_mid[2] | 0, 0.25 * prior_scale);
  lprior += normal_lpdf(b_mid[3] | 0, 0.25 * prior_scale);
  lprior += normal_lpdf(Intercept_mid | 0, 0.5 * prior_scale);
  lprior += student_t_lpdf(sigma_mid | 3, 0, 0.15 * prior_scale)
    - 1 * student_t_lccdf(0 | 3, 0, 0.15 * prior_scale);
  lprior += normal_lpdf(b_hind[1] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_hind[2] | 0, 0.25 * prior_scale);
  lprior += normal_lpdf(b_hind[3] | 0, 0.25 * prior_scale);
  lprior += normal_lpdf(Intercept_hind | 0, 0.5 * prior_scale);
  lprior += student_t_lpdf(sigma_hind | 3, 0, 0.15 * prior_scale)
    - 1 * student_t_lccdf(0 | 3, 0, 0.15 * prior_scale);
  lprior += normal_lpdf(b_earlin[1] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_earlin[2] | 0, 0.25 * prior_scale);
  lprior += normal_lpdf(b_earlin[3] | 0, 0.25 * prior_scale);
  lprior += normal_lpdf(Intercept_earlin | 0, 0.5 * prior_scale);
  lprior += student_t_lpdf(sigma_earlin | 3, 0, 0.15 * prior_scale)
    - 1 * student_t_lccdf(0 | 3, 0, 0.15 * prior_scale);
  lprior += normal_lpdf(b_eye[1] | 0, 1 * prior_scale);
  lprior += normal_lpdf(b_eye[2] | 0, 0.25 * prior_scale);
  lprior += normal_lpdf(b_eye[3] | 0, 0.25 * prior_scale);
  lprior += normal_lpdf(Intercept_eye | 0, 0.5 * prior_scale);
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
    mu_fore += Intercept_fore + Xc_fore * b_fore;
    mu_mid += Intercept_mid + Xc_mid * b_mid;
    mu_hind += Intercept_hind + Xc_hind * b_hind;
    mu_earlin += Intercept_earlin + Xc_earlin * b_earlin;
    mu_eye += Intercept_eye + Xc_eye * b_eye;
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
  // actual population-level intercept
  real b_fore_Intercept = Intercept_fore - dot_product(means_X_fore, b_fore);
  // actual population-level intercept
  real b_mid_Intercept = Intercept_mid - dot_product(means_X_mid, b_mid);
  // actual population-level intercept
  real b_hind_Intercept = Intercept_hind - dot_product(means_X_hind, b_hind);
  // actual population-level intercept
  real b_earlin_Intercept = Intercept_earlin - dot_product(means_X_earlin, b_earlin);
  // actual population-level intercept
  real b_eye_Intercept = Intercept_eye - dot_product(means_X_eye, b_eye);
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

