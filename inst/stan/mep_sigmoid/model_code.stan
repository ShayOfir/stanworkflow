// Hierarchical sigmoid MEP input-output curve (IOcurve9 of ShayOfir/MEP_IOcurve,
// Wonderman Bar Sela et al. 2026), renamed to stanworkflow's data conventions.
// curve: mu = A / (1 + exp(-slope * (intensity - theta))), y ~ lognormal(log(mu), sigma)
// theta, slope, log(A) = population mean + coil effect + side effect + subject effect.
data {
  int<lower=1> N;
  int<lower=1> n_subj;
  array[N] int<lower=1, upper=n_subj> subj_id;
  int<lower=1> n_coil;
  array[N] int<lower=1, upper=n_coil> coil_id;
  int<lower=1> n_side;
  array[N] int<lower=1, upper=n_side> side_id;
  vector[N] intensity;
  vector<lower=0>[N] y;
  vector[n_coil] mso_limit;
  real prior_mu_theta_mean;
  real prior_mu_theta_sd;
  real prior_mu_slope_mean;
  real prior_mu_slope_sd;
  real prior_mu_A_mean;
  real prior_mu_A_sd;
  real prior_sigma_theta_scale;
  real prior_sigma_slope_scale;
  real prior_sigma_A_scale;
  real prior_sigma_scale;
  real prior_delta_theta_coil_scale;
  real prior_delta_theta_side_scale;
  real prior_delta_slope_coil_scale;
  real prior_delta_slope_side_scale;
  real prior_delta_logA_coil_scale;
  real prior_delta_logA_side_scale;
  int<lower=0, upper=1> skip_likelihood;
}
transformed data {
  real sigma_ubound = 2;
  real logA_cap = 2;
}
parameters {
  vector[n_subj] theta_z;
  vector[n_subj] slope_z;
  vector[n_subj] logA_z;
  real mu_theta;
  real mu_slope;
  real mu_logA;
  real<lower=0, upper=sigma_ubound> sigma_theta;
  real<lower=0, upper=sigma_ubound> sigma_slope;
  real<lower=0, upper=sigma_ubound> sigma_logA;
  real<lower=0, upper=sigma_ubound> sigma;
  vector[n_coil] delta_theta_coil_raw;
  vector[n_side] delta_theta_side_raw;
  vector[n_coil] delta_slope_coil_raw;
  vector[n_side] delta_slope_side_raw;
  vector[n_coil] delta_logA_coil_raw;
  vector[n_side] delta_logA_side_raw;
}
transformed parameters {
  // sum-to-zero coil and side effects
  vector[n_coil] delta_theta_coil = delta_theta_coil_raw - mean(delta_theta_coil_raw);
  vector[n_side] delta_theta_side = delta_theta_side_raw - mean(delta_theta_side_raw);
  vector[n_coil] delta_slope_coil = delta_slope_coil_raw - mean(delta_slope_coil_raw);
  vector[n_side] delta_slope_side = delta_slope_side_raw - mean(delta_slope_side_raw);
  vector[n_coil] delta_logA_coil = delta_logA_coil_raw - mean(delta_logA_coil_raw);
  vector[n_side] delta_logA_side = delta_logA_side_raw - mean(delta_logA_side_raw);
  // curve parameters per subject x coil x side
  array[n_subj, n_coil, n_side] real theta;
  array[n_subj, n_coil, n_side] real slope;
  array[n_subj, n_coil, n_side] real logA;
  for (i in 1:n_subj) for (j in 1:n_coil) for (k in 1:n_side) {
    theta[i, j, k] = mu_theta + delta_theta_coil[j] + delta_theta_side[k] + sigma_theta * theta_z[i];
    slope[i, j, k] = mu_slope + delta_slope_coil[j] + delta_slope_side[k] + sigma_slope * slope_z[i];
    logA[i, j, k]  = mu_logA  + delta_logA_coil[j]  + delta_logA_side[k]  + sigma_logA  * logA_z[i];
  }
}
model {
  mu_theta ~ normal(prior_mu_theta_mean, prior_mu_theta_sd);
  mu_slope ~ normal(prior_mu_slope_mean, prior_mu_slope_sd);
  mu_logA  ~ normal(log(prior_mu_A_mean), prior_mu_A_sd);
  sigma_theta ~ student_t(3, 0, prior_sigma_theta_scale);
  sigma_slope ~ student_t(3, 0, prior_sigma_slope_scale);
  sigma_logA  ~ student_t(3, 0, prior_sigma_A_scale);
  sigma       ~ student_t(3, 0, prior_sigma_scale);
  theta_z ~ std_normal();
  slope_z ~ std_normal();
  logA_z  ~ std_normal();
  delta_theta_coil_raw ~ normal(0, prior_delta_theta_coil_scale);
  delta_theta_side_raw ~ normal(0, prior_delta_theta_side_scale);
  delta_slope_coil_raw ~ normal(0, prior_delta_slope_coil_scale);
  delta_slope_side_raw ~ normal(0, prior_delta_slope_side_scale);
  delta_logA_coil_raw  ~ normal(0, prior_delta_logA_coil_scale);
  delta_logA_side_raw  ~ normal(0, prior_delta_logA_side_scale);
  if (!skip_likelihood) {
    for (n in 1:N) {
      // MNAR: intensities above a coil's maximal output cannot be observed
      if (intensity[n] <= mso_limit[coil_id[n]]) {
        int i = subj_id[n]; int j = coil_id[n]; int k = side_id[n];
        real mu = exp(fmin(logA[i, j, k], logA_cap)) /
                  (1 + exp(-slope[i, j, k] * (intensity[n] - theta[i, j, k])));
        y[n] ~ lognormal(log(mu + 1e-9), sigma);
      }
    }
  }
}
generated quantities {
  // mu_rep: the curve (median MEP; the mean is exp(sigma^2 / 2) times larger)
  vector[N] mu_rep;
  vector[N] y_rep;
  for (n in 1:N) {
    int i = subj_id[n]; int j = coil_id[n]; int k = side_id[n];
    mu_rep[n] = exp(fmin(logA[i, j, k], logA_cap)) /
                (1 + exp(-slope[i, j, k] * (intensity[n] - theta[i, j, k])));
    y_rep[n] = lognormal_rng(log(mu_rep[n] + 1e-9), sigma);
  }
}
