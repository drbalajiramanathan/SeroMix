// seromix.stan (Final Version with Re-parameterization)
functions {
  real log_titre(real t, real t_inf, real baseline, real peak, real decay, real t_peak) {
    if (t < t_inf) {
      return baseline;
    }
    real time_since_infection = t - t_inf;
    if (time_since_infection < t_peak) {
      return fmax(baseline, peak * (time_since_infection / t_peak));
    } else {
      return fmax(baseline, (peak - 1.5) * exp(-decay * (time_since_infection - t_peak)) + 1.5);
    }
  }
}
data {
  int<lower=1> N_obs;
  int<lower=1> N_ind;
  array[N_obs] int<lower=1, upper=N_ind> id;
  vector[N_obs] y;
  vector[N_obs] t;
  vector[N_ind] baseline;
}
parameters {
  // Re-parameterized infection time
  array[N_ind] real<lower=-2, upper=2> t_inf_raw; // Sample a standardized value

  array[N_ind] simplex[3] theta;
  real<lower=0, upper=10> peak_primary;
  real<lower=0, upper=10> peak_breakthrough;
  real<lower=0, upper=1> decay_primary;
  real<lower=0, upper=1> decay_breakthrough;
}
transformed parameters {
  // Transform back to the original scale
  array[N_ind] real<lower=0, upper=150> t_inf;
  for (i in 1:N_ind) {
    // Map the standardized t_inf_raw to the study period (approx. day 0-150)
    t_inf[i] = 75 + 37.5 * t_inf_raw[i];
  }
}
model {
  // Prior for the standardized infection time parameter
  t_inf_raw ~ std_normal(); // This is equivalent to normal(0, 1)

  // Tighter priors for stability
  peak_primary ~ normal(2.5, 0.5);
  peak_breakthrough ~ normal(4.0, 0.5);
  decay_primary ~ normal(0.03, 0.005);
  decay_breakthrough ~ normal(0.04, 0.005);

  for (i in 1:N_ind) {
    real log_lik_uninf = 0;
    real log_lik_primary = 0;
    real log_lik_breakthrough = 0;

    for (j in 1:N_obs) {
      if (id[j] == i) {
        log_lik_uninf += normal_lpdf(y[j] | baseline[i], 0.25);
        log_lik_primary += normal_lpdf(y[j] | log_titre(t[j], t_inf[i], baseline[i], peak_primary, decay_primary, 21), 0.25);
        log_lik_breakthrough += normal_lpdf(y[j] | log_titre(t[j], t_inf[i], baseline[i], peak_breakthrough, decay_breakthrough, 14), 0.25);
      }
    }

    row_vector[3] lps = [log(theta[i, 1]) + log_lik_uninf,
                         log(theta[i, 2]) + log_lik_primary,
                         log(theta[i, 3]) + log_lik_breakthrough];
    target += log_sum_exp(lps);
  }
}
