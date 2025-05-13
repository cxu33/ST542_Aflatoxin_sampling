library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# simulate_trainload 
simulate_trainload <- function(nrailcars_per_train,
                               prob_of_zero,
                               shape_param,
                               scale_param) {
  
  value_is_nonzero_indicator <- rbinom(n = nrailcars_per_train,
                                       size = 1,
                                       p = 1 - prob_of_zero)
  gamma_component <- rgamma(nrailcars_per_train,
                            shape = shape_param,
                            scale = scale_param)
  
  observed <- value_is_nonzero_indicator * gamma_component
  return(observed)
}

# simple_sample function (sample from gamma-distribution)
simple_sample <- function(ppb, size) {
  sample(length(ppb), size)
}

# simulation function: with a simulated "missing the hotspots" procedure AFTER sampling from gamma-distribution
run_simulation <- function(n_sim_per_param, param.tb,
                           nrailcars_per_train = 200,
                           regulatory_threshold = 20) {
  
  pb <- txtProgressBar(min = 0, max = nrow(param.tb), style = 3)
  sim.tb <- data.frame()
  
  for (i in seq_len(nrow(param.tb))) {
    parameters <- param.tb[i, ]
    
    for (j in seq_len(n_sim_per_param)) {
      
      # simulate the full trainload
      trainload_sim_ppb <- simulate_trainload(nrailcars_per_train,
                                              parameters$prob_of_zero,
                                              parameters$shape_param,
                                              parameters$scale_param)
      
      train_true_safe_prop <- mean(trainload_sim_ppb < regulatory_threshold)
      
      # draw your sample
      sample_idx <- simple_sample(trainload_sim_ppb, parameters$sample_sizes)
      sample_ppb <- trainload_sim_ppb[sample_idx]
      
      # randomly drop some sampled values to zero
      # to mimic missing hotspots in the composite
      sample_zero_indicator <- rbinom(length(sample_ppb),
                                      size = 1,
                                      p = 1 - parameters$prob_of_zero)
      sample_ppb <- sample_ppb * sample_zero_indicator
      
      # decision rule
      sample_safe <- (max(sample_ppb) < regulatory_threshold)
      
      # record results
      add.tb <- data.frame(
        parameter_id         = i,
        rep_id               = j,
        train_true_safe_prop = train_true_safe_prop,
        sample_safe          = sample_safe
      )
      sim.tb <- rbind(sim.tb, add.tb)
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(sim.tb)
}

set.seed(454)
param.tb <- expand.grid(
  prob_of_zero  = seq(0, 0.3, by = 0.05),
  shape_param   = seq(1, 4,  by = 1),
  scale_param   = seq(1, 5,  by = 1),
  sample_sizes  = seq(25, 200, by = 25)
)
res3 <- run_simulation(n_sim_per_param = 500, param.tb = param.tb)
save(res3, file = "~/simulation3_parameter_grid.RData")
