library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Function to simulate one trainload with fixed parameters
simulate_trainload <- function(nrailcars_per_train, 
                               prob_of_zero,
                               shape_param, 
                               scale_param) {
  
  # Step one: determine if the value is 0 or not
  # (returns list of 0s and 1s)
  value_is_nonzero_indicator <- rbinom(n = nrailcars_per_train, 
                                       size = 1,
                                       p = 1-prob_of_zero)
  
  # Step two: draw from gamma dist
  gamma_component <- rgamma(nrailcars_per_train, 
                            shape = shape_param, 
                            scale = scale_param)  
  
  # Step three: Multiply together to get a single vector
  observed <- value_is_nonzero_indicator*gamma_component
  
  return(observed)
}

# Simple uniform sampling function
simple_sample <- function(ppb, size) {
  sample(length(ppb), size)
}

# Define combinations of parameter values to simulation
param.tb <- expand.grid(prob_of_zero = seq(0, 0.05, by = 0.05), # 2 
                        shape_param = seq(1, 4, by = 1), # 4
                        scale_param = seq(1, 5, by = 1), # 5
                        sample_sizes = seq(25, 200, by = 25))
# Function to simulate 
run_simulation <- function(n_sim_per_param, param.tb, 
                           nrailcars_per_train = 200, 
                           regulatory_threshold = 20) {
  
  # Progress bar setup
  pb <- txtProgressBar(min = 0, max = nrow(param.tb), style = 3)
  
  # Initialize table to save output at each simulation
  sim.tb <- data.frame()
  
  # Loop over parameter combinations
  for (i in 1:nrow(param.tb)) {
    
    parameters <- param.tb[i,]
    
    for (j in 1:n_sim_per_param) {
      
      trainload_sim_ppb <- simulate_trainload(nrailcars_per_train, 
                                              parameters$prob_of_zero,
                                              parameters$shape_param, 
                                              parameters$scale_param)
      
      train_true_safe_prop <- mean(trainload_sim_ppb < regulatory_threshold)
      
      # Draw a uniform sample.
      sample_idx <- simple_sample(trainload_sim_ppb, parameters$sample_sizes)
      sample_ppb <- trainload_sim_ppb[sample_idx]
      
      # Save information on the simulation
      add.tb <- data.frame(parameter_id = i,
                           rep_id = j,
                           
                           # True trainload data
                           train_true_safe_prop = train_true_safe_prop,
                           
                           # Decision rule: the sample is safe if the maximum ppb in the sample is below the regulatory threshold.
                           sample_safe = (max(sample_ppb) < regulatory_threshold))
      sim.tb <- rbind(sim.tb, add.tb)
      
    }
    
    setTxtProgressBar(pb, i)  # Update progress bar
  }
  close(pb)
  return(sim.tb)
}

set.seed(454)
res1 <- run_simulation(n_sim_per_param = 500, param.tb)
save(res1, file = "/Users/changxu/Desktop/NCSU/ST542/simulation1_parameter_grid.RData")

load("/Users/changxu/Desktop/NCSU/ST542/simulation1_parameter_grid.RData")

# A train is considered truly safe if the proportion of railcars below the regulatory threshold is at least tp.
target_safe_proportions <- c(0.9, 0.95, 0.99, 1)
res1 <- res1 %>%
  crossing(target_safe_proportions = target_safe_proportions)
res1 <- res1 %>% 
  mutate(train_is_safe = train_true_safe_prop >= target_safe_proportions)
param.tb$parameter_id <- 1:nrow(param.tb)
res1 <- res1 %>% 
  left_join(param.tb, by = "parameter_id")

# Calculate power, etc
summary1.tb <- res1 %>% 
  group_by(parameter_id, prob_of_zero, shape_param, scale_param, sample_sizes, target_safe_proportions, train_is_safe) %>% 
  reframe(total_trains = n(), 
          num_pass = sum(sample_safe),
          num_fail = sum(!sample_safe),
          frac_pass = sum(sample_safe)/n(),
          frac_pass_moe = 1.96*sqrt(frac_pass*(1-frac_pass)/total_trains),
          frac_pass_ci_lb = frac_pass - frac_pass_moe,
          frac_pass_ci_ub = frac_pass + frac_pass_moe) 
summary1.tb %>% 
  filter(train_is_safe) %>% 
  pull(frac_pass) %>% 
  summary
producer_risk_09 <- summary1.tb %>%
  filter(train_is_safe == TRUE,
         target_safe_proportions == 0.9) %>%
  dplyr::select(prob_of_zero, shape_param, scale_param, sample_sizes, frac_pass)
ggplot(producer_risk_09,
       aes(x = sample_sizes,
           y = 1-frac_pass,
           color = factor(scale_param))) +
  geom_line() +
  geom_point() +
  facet_grid(prob_of_zero ~ shape_param, labeller = label_both) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x     = "Sample Size",
    y     = "Producer’s Risk\n(% ‘unsafe’ samples when train safe)",
    color = "Scale\nParam",
    title = "Producer’s Risk vs. Sample Size\n(target safe = 90%)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  )
producer_risk_095 <- summary1.tb %>%
  filter(train_is_safe == TRUE,
         target_safe_proportions == 0.95) %>%
  dplyr::select(prob_of_zero, shape_param, scale_param, sample_sizes, frac_pass)
ggplot(producer_risk_095,
       aes(x = sample_sizes,
           y = 1-frac_pass,
           color = factor(scale_param))) +
  geom_line() +
  geom_point() +
  facet_grid(prob_of_zero ~ shape_param, labeller = label_both) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x     = "Sample Size",
    y     = "Producer’s Risk\n(% ‘unsafe’ samples when train safe)",
    color = "Scale\nParam",
    title = "Producer’s Risk vs. Sample Size\n(target safe = 95%)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  )
producer_risk_099 <- summary1.tb %>%
  filter(train_is_safe == TRUE,
         target_safe_proportions == 0.99) %>%
  dplyr::select(prob_of_zero, shape_param, scale_param, sample_sizes, frac_pass)
ggplot(producer_risk_099,
       aes(x = sample_sizes,
           y = 1-frac_pass,
           color = factor(scale_param))) +
  geom_line() +
  geom_point() +
  facet_grid(prob_of_zero ~ shape_param, labeller = label_both) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x     = "Sample Size",
    y     = "Producer’s Risk\n(% ‘unsafe’ samples when train safe)",
    color = "Scale\nParam",
    title = "Producer’s Risk vs. Sample Size\n(target safe = 99%)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  )
producer_risk_1 <- summary1.tb %>%
  filter(train_is_safe == TRUE,
         target_safe_proportions == 1) %>%
  dplyr::select(prob_of_zero, shape_param, scale_param, sample_sizes, frac_pass)
ggplot(producer_risk_1,
       aes(x = sample_sizes,
           y = 1-frac_pass,
           color = factor(scale_param))) +
  geom_line() +
  geom_point() +
  facet_grid(prob_of_zero ~ shape_param, labeller = label_both) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x     = "Sample Size",
    y     = "Producer’s Risk\n(% ‘unsafe’ samples when train safe)",
    color = "Scale\nParam",
    title = "Producer’s Risk vs. Sample Size\n(target safe = 100%)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  )
# 4) Extract consumer's risk where train truly unsafe AND target == 0.9
consumer_risk_09 <- summary1.tb %>%
  filter(train_is_safe == FALSE,
         target_safe_proportions == 0.9) %>%
  dplyr::select(prob_of_zero, shape_param, scale_param, sample_sizes, frac_pass)

# Now plot
ggplot(consumer_risk_09,
       aes(x = sample_sizes,
           y = frac_pass,
           color = factor(scale_param))) +
  geom_line() +
  geom_point() +
  facet_grid(prob_of_zero ~ shape_param, labeller = label_both) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x     = "Sample Size",
    y     = "Consumer’s Risk\n(% ‘safe’ samples when train unsafe)",
    color = "Scale\nParam",
    title = "Consumer’s Risk vs. Sample Size\n(target safe = 90%)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  )

consumer_risk_095 <- summary1.tb %>%
  filter(train_is_safe == FALSE,
         target_safe_proportions == 0.95) %>%
  dplyr::select(prob_of_zero, shape_param, scale_param, sample_sizes, frac_pass)

# Now plot
ggplot(consumer_risk_095,
       aes(x = sample_sizes,
           y = frac_pass,
           color = factor(scale_param))) +
  geom_line() +
  geom_point() +
  facet_grid(prob_of_zero ~ shape_param, labeller = label_both) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x     = "Sample Size",
    y     = "Consumer’s Risk\n(% ‘safe’ samples when train unsafe)",
    color = "Scale\nParam",
    title = "Consumer’s Risk vs. Sample Size\n(target safe = 95%)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  )
consumer_risk_099 <- summary1.tb %>%
  filter(train_is_safe == FALSE,
         target_safe_proportions == 0.99) %>%
  dplyr::select(prob_of_zero, shape_param, scale_param, sample_sizes, frac_pass)

# Now plot
ggplot(consumer_risk_099,
       aes(x = sample_sizes,
           y = frac_pass,
           color = factor(scale_param))) +
  geom_line() +
  geom_point() +
  facet_grid(prob_of_zero ~ shape_param, labeller = label_both) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x     = "Sample Size",
    y     = "Consumer’s Risk\n(% ‘safe’ samples when train unsafe)",
    color = "Scale\nParam",
    title = "Consumer’s Risk vs. Sample Size\n(target safe = 99%)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  )
consumer_risk_1 <- summary1.tb %>%
  filter(train_is_safe == FALSE,
         target_safe_proportions == 1) %>%
  dplyr::select(prob_of_zero, shape_param, scale_param, sample_sizes, frac_pass)

# Now plot
ggplot(consumer_risk_1,
       aes(x = sample_sizes,
           y = frac_pass,
           color = factor(scale_param))) +
  geom_line() +
  geom_point() +
  facet_grid(prob_of_zero ~ shape_param, labeller = label_both) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x     = "Sample Size",
    y     = "Consumer’s Risk\n(% ‘safe’ samples when train unsafe)",
    color = "Scale\nParam",
    title = "Consumer’s Risk vs. Sample Size\n(target safe = 100%)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  )
# 1) Compute consumer's risk summary:
consumer_risk.tb <- summary1.tb %>%
  filter(train_is_safe == FALSE) %>%    # only when the train truly isn't safe
  dplyr::select(prob_of_zero, shape_param, scale_param, sample_sizes, frac_pass)

# 2) Plot
ggplot(consumer_risk.tb,
       aes(x = sample_sizes,
           y = frac_pass,
           color = factor(scale_param))) +
  geom_line() +
  geom_point() +
  facet_grid(prob_of_zero ~ shape_param,
             labeller = label_both) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = "Sample Size",
    y     = "Consumer’s Risk\n(% of ‘safe’ samples when train unsafe)",
    color = "Scale\nParam",
    title = "Consumer’s Risk vs. Sample Size\nacross Parameter Combinations"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  )


# 1) Summarize consumer's risk across all combos:
consumer_risk_all <- res1 %>%
  filter(train_is_safe == FALSE) %>%               # only truly-unsafe trains
  group_by(
    prob_of_zero,
    shape_param,
    scale_param,
    sample_sizes,
    target_safe_proportions
  ) %>%
  summarize(
    consumer_risk = mean(sample_safe),            # fraction of “safe” samples
    .groups = "drop"
  )

# 2) Plot it:
ggplot(consumer_risk_all,
       aes(
         x = target_safe_proportions,
         y = consumer_risk,
         color    = factor(sample_sizes),
         linetype = factor(scale_param),
         group    = interaction(sample_sizes, scale_param)
       )) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1, alpha = 0.8) +
  facet_grid(
    prob_of_zero ~ shape_param,
    labeller = label_both
  ) +
  scale_x_continuous(
    breaks = c(0.9, 0.95, 0.99, 1),
    labels = percent_format(accuracy = 1)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x     = "Target Safe Proportion",
    y     = "Consumer’s Risk\n(% “safe” samples when train truly unsafe)",
    color = "Sample\nSize",
    linetype = "Scale\nParam",
    title = "Consumer’s Risk vs. Target Threshold\nacross All Parameter Combinations"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor  = element_blank()
  )

