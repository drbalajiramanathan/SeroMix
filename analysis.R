# --- PHASE 3: FINAL ANALYSIS WITH STAN ---

# 1. Install and load cmdstanr
if (!require("cmdstanr")) {
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  cmdstanr::install_cmdstan()
}
library(cmdstanr)
library(dplyr)
library(ggplot2)

# 2. Load and prepare data
observed_data <- read.csv("observed_data.csv")
ground_truth <- read.csv("ground_truth.csv")

stan_data <- list(
  N_obs = nrow(observed_data),
  N_ind = nrow(ground_truth),
  id = observed_data$id,
  y = observed_data$titre_value,
  t = observed_data$sample_day,
  baseline = ground_truth$baseline_titre
)

# 3. Compile the Stan model
print("Compiling Stan model (this may take a few minutes)...")
model <- cmdstan_model("seromix.stan")
print("Model compiled.")

# 4. Run the Stan sampler (Corrected for Low Memory)
print("Running Stan sampler (sequentially to save memory)...")
fit <- model$sample(
  data = stan_data,
  seed = 123,
  chains = 2,           # Run 2 chains instead of 4
  parallel_chains = 1,  # IMPORTANT: Run chains one after another
  refresh = 500,
  iter_warmup = 1000,
  iter_sampling = 1000
)
print("Sampler finished.")

# 5. Process and visualize results
print("Processing results...")
posterior_summary <- fit$summary("theta")

# Extract posterior probabilities of being in each group
prob_summary <- posterior_summary %>%
  tibble() %>%
  mutate(
    id = as.integer(gsub("theta\\[([0-9]+),([0-9]+)\\]", "\\1", variable)),
    group = as.integer(gsub("theta\\[([0-9]+),([0-9]+)\\]", "\\2", variable))
  ) %>%
  select(id, group, posterior_prob = mean) %>%
  pivot_wider(names_from = group, values_from = posterior_prob, names_prefix = "prob_group_") %>%
  mutate(prob_infected = prob_group_2 + prob_group_3) # Prob(Primary) + Prob(Breakthrough)

# Join with ground truth for comparison
results_summary <- left_join(ground_truth, prob_summary, by = "id")

print("Posterior Probability of Infection vs. True Status (first 20 individuals):")
print(head(results_summary, 20))

# Create a summary plot
prob_plot <- ggplot(results_summary, aes(x = infection_status, y = prob_infected, fill = infection_status)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Stan Model: Posterior Probability of Infection by True Status",
    x = "True Infection Status",
    y = "Posterior Probability of Infection"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face="bold"))

print(prob_plot)