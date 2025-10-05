#include <Rcpp.h>
using namespace Rcpp;

// --- Internal C++ helper function for kinetics ---
double biphasic_kinetics_cpp(double t, double t_infection, List params) {
  double time_since_infection = t - t_infection;
  if (R_IsNA(time_since_infection) || time_since_infection < 0) {
    return NA_REAL;
  }
  double t_peak = as<double>(params["t_peak"]);
  double peak = as<double>(params["peak"]);
  double plateau = as<double>(params["plateau"]);
  double decay = as<double>(params["decay"]);
  if (time_since_infection < t_peak) {
    return peak * (time_since_infection / t_peak);
  } else {
    return (peak - plateau) * exp(-decay * (time_since_infection - t_peak)) + plateau;
  }
}

// --- Internal C++ function to calculate log-likelihood for a single person ---
double calculate_person_log_likelihood(
    int person_id,
    DataFrame& observed_data,
    IntegerVector& obs_ids,
    NumericVector& obs_sample_days,
    NumericVector& obs_titre_values,
    CharacterVector& statuses,
    NumericVector& infection_times,
    NumericVector& baselines,
    List& primary_params,
    List& breakthrough_params,
    double noise_sd
) {
  
  double person_log_likelihood = 0.0;
  String status = statuses[person_id - 1];
  
  if (status == "Uninfected") {
    for (int i = 0; i < observed_data.nrows(); ++i) {
      if (obs_ids[i] == person_id) {
        person_log_likelihood += R::dnorm(obs_titre_values[i], baselines[person_id - 1], noise_sd, true);
      }
    }
    return person_log_likelihood;
  }
  
  for (int i = 0; i < observed_data.nrows(); ++i) {
    if (obs_ids[i] == person_id) {
      double true_titre = 0.0;
      List current_params = (status == "Primary") ? primary_params : breakthrough_params;
      
      double titre_from_infection = biphasic_kinetics_cpp(
        obs_sample_days[i], 
                       infection_times[person_id - 1], 
                                      current_params
      );
      
      if (R_IsNA(titre_from_infection)) {
        true_titre = baselines[person_id - 1];
      } else {
        true_titre = std::max(baselines[person_id - 1], titre_from_infection);
      }
      
      person_log_likelihood += R::dnorm(obs_titre_values[i], true_titre, noise_sd, true);
    }
  }
  return person_log_likelihood;
}


// --- The main exported MCMC sampler function ---
//' @export
 // [[Rcpp::export]]
 List run_mcmc_sampler(
     DataFrame observed_data,
     DataFrame initial_state,
     List primary_params,
     List breakthrough_params,
     double noise_sd,
     int iterations,
     double study_duration = 150.0
 ) {
   
   int n_individuals = initial_state.nrows();
   
   NumericMatrix infection_time_chain(iterations, n_individuals);
   CharacterMatrix infection_status_chain(iterations, n_individuals);
   
   NumericVector current_infection_times = initial_state["infection_time"];
   CharacterVector current_statuses = initial_state["infection_status"];
   NumericVector baselines = initial_state["baseline_titre"];
   
   IntegerVector obs_ids = observed_data["id"];
   NumericVector obs_sample_days = observed_data["sample_day"];
   NumericVector obs_titre_values = observed_data["titre_value"];
   
   NumericVector person_log_likelihoods(n_individuals);
   for (int i = 0; i < n_individuals; ++i) {
     person_log_likelihoods[i] = calculate_person_log_likelihood(
       i + 1, observed_data, obs_ids, obs_sample_days, obs_titre_values,
       current_statuses, current_infection_times, baselines, 
       primary_params, breakthrough_params, noise_sd
     );
   }
   
   for (int iter = 0; iter < iterations; ++iter) {
     if ((iter + 1) % 1000 == 0) {
       Rcout << "Iteration " << iter + 1 << "/" << iterations << std::endl;
     }
     
     for (int p = 0; p < n_individuals; ++p) {
       
       double jump_type = R::runif(0, 1);
       
       if (current_statuses[p] == "Uninfected" || jump_type < 0.5) {
         // --- JUMP TYPE 1: Birth/Death Move ---
         CharacterVector proposed_statuses = clone(current_statuses); // CORRECTED: Use clone for a true copy
         NumericVector proposed_times = clone(current_infection_times);   // CORRECTED: Use clone for a true copy
         double log_jacobian = 0.0;
         
         if (current_statuses[p] == "Uninfected") {
           proposed_times[p] = R::runif(1, study_duration);
           proposed_statuses[p] = (baselines[p] < 1.0) ? "Primary" : "Breakthrough";
           log_jacobian = log(study_duration);
         } else {
           proposed_times[p] = NA_REAL;
           proposed_statuses[p] = "Uninfected";
           log_jacobian = -log(study_duration);
         }
         
         double proposed_log_likelihood = calculate_person_log_likelihood(p + 1, observed_data, obs_ids, obs_sample_days, obs_titre_values, proposed_statuses, proposed_times, baselines, primary_params, breakthrough_params, noise_sd);
         double acceptance_ratio = exp(proposed_log_likelihood - person_log_likelihoods[p] + log_jacobian);
         
         if (R::runif(0, 1) < acceptance_ratio) {
           current_statuses = proposed_statuses;
           current_infection_times = proposed_times;
           person_log_likelihoods[p] = proposed_log_likelihood;
         }
         
       } else if (jump_type < 0.85) {
         // --- JUMP TYPE 2: Update infection time (if infected) ---
         double current_time = current_infection_times[p];
         double proposed_time = R::rnorm(current_time, 5.0);
         
         if (proposed_time > 0 && proposed_time < study_duration) {
           NumericVector proposed_times = clone(current_infection_times);
           proposed_times[p] = proposed_time;
           
           double proposed_log_likelihood = calculate_person_log_likelihood(p + 1, observed_data, obs_ids, obs_sample_days, obs_titre_values, current_statuses, proposed_times, baselines, primary_params, breakthrough_params, noise_sd);
           double acceptance_ratio = exp(proposed_log_likelihood - person_log_likelihoods[p]);
           
           if (R::runif(0, 1) < acceptance_ratio) {
             current_infection_times[p] = proposed_time;
             person_log_likelihoods[p] = proposed_log_likelihood;
           }
         }
       } else {
         // --- JUMP TYPE 3: Switch infection type (if infected) ---
         CharacterVector proposed_statuses = clone(current_statuses); // CORRECTED: Use clone for a true copy
         if (current_statuses[p] == "Primary") {
           proposed_statuses[p] = "Breakthrough";
         } else {
           proposed_statuses[p] = "Primary";
         }
         
         double proposed_log_likelihood = calculate_person_log_likelihood(p + 1, observed_data, obs_ids, obs_sample_days, obs_titre_values, proposed_statuses, current_infection_times, baselines, primary_params, breakthrough_params, noise_sd);
         double acceptance_ratio = exp(proposed_log_likelihood - person_log_likelihoods[p]);
         
         if (R::runif(0, 1) < acceptance_ratio) {
           current_statuses = proposed_statuses;
           person_log_likelihoods[p] = proposed_log_likelihood;
         }
       }
     }
     
     infection_time_chain(iter, _) = current_infection_times;
     infection_status_chain(iter, _) = current_statuses;
   }
   
   return List::create(
     _["infection_time_chain"] = infection_time_chain,
     _["infection_status_chain"] = infection_status_chain
   );
 }