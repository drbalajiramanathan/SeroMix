// seromix_twostate.stan
// A simplified 2-state model (Uninfected vs. Infected) to establish a stable baseline.
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
  array[N_ind] real<lower=-2, upper=2> t_inf_raw;
  array[N_ind] simplex[2] theta; // Simplified to 2 states

  // A single set of kinetic parameters for a generic "infected" state
  real<lower=0, upper=10> peak_infected;
  real<lower=0, upper=1> decay_infected;
}
transformed parameters {
  array[N_ind] real<lower=0, upper=150> t_inf;
  for (i in 1:N_ind) {
    t_inf[i] = 75 + 37.5 * t_inf_raw[i];
  }
}
model {
  // Priors
  t_inf_raw ~ std_normal();
  peak_infected ~ normal(3.0, 1.0); // Prior is an average of the two previous peaks
  decay_infected ~ normal(0.035, 0.01);

  for (i in 1:N_ind) {
    real log_lik_uninf = 0;
    real log_lik_infected = 0;

    for (j in 1:N_obs) {
      if (id[j] == i) {
        log_lik_uninf += normal_lpdf(y[j] | baseline[i], 0.25);
        // Using a single t_peak of 18 days as a compromise
        log_lik_infected += normal_lpdf(y[j] | log_titre(t[j], t_inf[i], baseline[i], peak_infected, decay_infected, 18), 0.25);
      }
    }

    // Simplified 2-state mixture model
    row_vector[2] lps = [log(theta[i, 1]) + log_lik_uninf,
                         log(theta[i, 2]) + log_lik_infected];
    target += log_sum_exp(lps);
  }
}
