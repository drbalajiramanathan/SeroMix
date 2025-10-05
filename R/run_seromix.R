#' Run the Sero-Mix Baseline Model
#'
#' This function takes a dataframe of serological data, calculates the average
#' titre for each individual, and fits a two-component Bayesian mixture model
#' to classify individuals as 'infected' or 'uninfected'.
#'
#' @param data A dataframe with at least two columns: 'id' (patient identifier)
#' and 'titre_value' (the observed log titre).
#' @return A `CmdStanMCMC` object containing the results of the model fit.
#' @export
run_seromix_baseline <- function(data) {
  # 1. Check for required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required. Please install it.")
  }

  # 2. Prepare data: Calculate average titre
  message("Preparing data by calculating average titres...")
  baseline_data <- data %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(avg_titre = mean(titre_value, na.rm = TRUE))

  stan_data_baseline <- list(
    N_ind = nrow(baseline_data),
    avg_titre = baseline_data$avg_titre
  )

  # 3. Find and compile the Stan model included with the package
  message("Compiling the Stan model...")
  stan_file <- system.file("stan", "seromix_baseline.stan", package = "SeroMix")
  model_baseline <- cmdstanr::cmdstan_model(stan_file)

  # 4. Run the sampler
  message("Running the MCMC sampler...")
  fit_baseline <- model_baseline$sample(
    data = stan_data_baseline,
    seed = 123,
    chains = 4,
    parallel_chains = 2,
    refresh = 0
  )
  message("Analysis complete!")

  return(fit_baseline)
}
