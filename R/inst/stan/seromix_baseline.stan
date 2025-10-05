// seromix_baseline.stan
// A baseline Gaussian Mixture Model to test for fundamental separability.
data {
  int<lower=1> N_ind;
  vector[N_ind] avg_titre; // Using average titre per person
}
parameters {
  simplex[2] theta; // Overall proportion in each group
  ordered[2] mu;    // Ordered means for the two groups (low, high)
  vector<lower=0>[2] sigma; // Standard deviations for the two groups
}
model {
  // Priors
  mu[1] ~ normal(0.5, 0.5); // Prior for low titre group
  mu[2] ~ normal(2.5, 1);   // Prior for high titre group
  sigma ~ exponential(1);

  // Mixture Model
  for (i in 1:N_ind) {
    target += log_mix(theta[1],
                      normal_lpdf(avg_titre[i] | mu[1], sigma[1]),
                      normal_lpdf(avg_titre[i] | mu[2], sigma[2]));
  }
}
