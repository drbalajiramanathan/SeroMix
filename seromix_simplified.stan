// seromix_simplified.stan
// A simplified model to test the core mixture framework.
data {
  int<lower=1> N_obs;
  int<lower=1> N_ind;
  array[N_obs] int<lower=1, upper=N_ind> id;
  vector[N_obs] y;
  vector[N_ind] baseline;
}
parameters {
  array[N_ind] simplex[3] theta;
  real<lower=0> peak_primary;
  real<lower=0> peak_breakthrough;
}
model {
  // Priors
  peak_primary ~ normal(2.5, 1);
  peak_breakthrough ~ normal(4.0, 1);

  for (i in 1:N_ind) {
    // Log-likelihood contributions
    real log_lik_uninf = 0;
    real log_lik_primary = 0;
    real log_lik_breakthrough = 0;

    for (j in 1:N_obs) {
      if (id[j] == i) {
        log_lik_uninf += normal_lpdf(y[j] | baseline[i], 0.25);
        log_lik_primary += normal_lpdf(y[j] | peak_primary, 0.25); // Simplified
        log_lik_breakthrough += normal_lpdf(y[j] | peak_breakthrough, 0.25); // Simplified
      }
    }

    // Mixture model
    row_vector[3] lps = [log(theta[i, 1]) + log_lik_uninf,
                         log(theta[i, 2]) + log_lik_primary,
                         log(theta[i, 3]) + log_lik_breakthrough];
    target += log_sum_exp(lps);
  }
}
