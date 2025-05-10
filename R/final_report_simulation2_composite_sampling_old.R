# -------------------------------------------------------------------------# 
# Simulation 2: Composite Sampling Framework
# -------------------------------------------------------------------------# 

library(tidyverse)
library(ggplot2)
sessionInfo()
# R version 4.2.1 (2022-06-23 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5    
# [7] tidyr_1.3.0     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
#   [1] pillar_1.9.0      compiler_4.2.1    tools_4.2.1       digest_0.6.33     timechange_0.3.0 
# [6] evaluate_0.21     lifecycle_1.0.3   gtable_0.3.4      pkgconfig_2.0.3   rlang_1.1.1      
# [11] cli_3.6.1         rstudioapi_0.16.0 yaml_2.3.7        xfun_0.40         fastmap_1.1.1    
# [16] withr_2.5.0       knitr_1.43        generics_0.1.3    vctrs_0.6.5       hms_1.1.3        
# [21] grid_4.2.1        tidyselect_1.2.1  glue_1.6.2        data.table_1.14.2 R6_2.5.1         
# [26] fansi_1.0.3       rmarkdown_2.29    tzdb_0.4.0        magrittr_2.0.3    scales_1.3.0     
# [31] htmltools_0.5.8.1 colorspace_2.0-3  utf8_1.2.2        stringi_1.7.8     munsell_0.5.0   


# Define functions for simulation -------------------------------------

# Function to simulate aflatoxin levels in railcars from one trainload
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


# Function to loop through grid of parameter values and simulate several trainloads
# for each combination of parameters
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
      
      # Ground truth metrics for the trainload
      train_true_safe_prop <- mean(trainload_sim_ppb < regulatory_threshold)
      train_true_avg <- mean(trainload_sim_ppb)
      train_true_avg_of_unsafe <- ifelse(train_true_safe_prop < 1,
                                         mean(trainload_sim_ppb[trainload_sim_ppb >= regulatory_threshold]),
                                         NA_real_)
      
      # Simulate getting composite samples
      num_composites <- nrailcars_per_train / parameters$num_railcars_per_composite
      trainload_sim_tb <- data.frame(values = trainload_sim_ppb,
                                     composite_id = rep(1:num_composites, 
                                                        each = parameters$num_railcars_per_composite ))
      facility_sim_tb <- trainload_sim_tb %>% 
        group_by(composite_id) %>% 
        reframe(afla_in_composite = mean(values))
      
      # Save information on the simulation
      add.tb <- data.frame(parameter_id = i,
                           rep_id = j,
                           
                           # True trainload data
                           train_true_safe_prop = train_true_safe_prop,
                           train_true_avg = train_true_avg,
                           train_true_avg_unsafe = train_true_avg_of_unsafe,
                           
                           # Facility data
                           composites_safe_prop = mean(facility_sim_tb$afla_in_composite < regulatory_threshold),
                           composites_avg = mean(facility_sim_tb$afla_in_composite),
                           composites_sd = sd(facility_sim_tb$afla_in_composite))
      
      # Add to full simulation table
      sim.tb <- rbind(sim.tb, add.tb)
      
    }
    
    setTxtProgressBar(pb, i)  # Update progress bar
  }
  close(pb)
  return(sim.tb)
}


# Run simulation ----------------------------------------

set.seed(123)

# Define combinations of parameter values to simulation
param.tb <- expand.grid(prob_of_zero = seq(0, 0.05, by = 0.05), # 2 values 
                        shape_param = c(0.5, 1, 2, 3, 4), # 5 values
                        scale_param = c(1,2,3,4,5,10,20), # 7 values
                        num_railcars_per_composite = c(2, 5, 10, 20, 40))

res2 <- run_simulation(n_sim_per_param = 500, param.tb)
save(res2, param.tb, file = "sim_data/simulation2_parameter_grid_Apr23_2025.RData")


# Clean results -----------------------------------------------------------

rm(list = setdiff("simulate_trainload",ls()))

load("sim_data/simulation2_parameter_grid_Apr23_2025.RData")

# Add labels
res2 <- res2 %>% 
  mutate(train_is_safe = train_true_safe_prop == 1,
         trainload_pass = composites_safe_prop == 1)

# Add parameter settings
param.tb$parameter_id <- 1:nrow(param.tb)
res2 <- res2 %>% 
  left_join(param.tb, by = "parameter_id")

# Add columns to assist with plotting
res2 <- res2 %>% 
  mutate(shape_plot = paste0("shape: ",shape_param), scale_plot = paste0("scale: ",scale_param),
         parameter_id_sans_nrpc = paste(prob_of_zero, shape_param, scale_param, sep = "_"))
res2 <- res2 %>% 
  mutate(scale_plot = factor(scale_plot, levels = paste0("scale: ",c(1,2,3,4,5,10,20))))

# Add "num unsafe railcars" to res2 for each trainload
res2  <- res2 %>% 
  mutate(num_unsafe_railcars = 200*(1-train_true_safe_prop))

# Calculate power, etc
summary.tb <- res2 %>% 
  group_by(parameter_id, prob_of_zero, shape_param, scale_param, 
           shape_plot, scale_plot, num_railcars_per_composite, parameter_id_sans_nrpc,
           train_is_safe) %>% 
  reframe(total_trains = n(), 
          num_pass = sum(trainload_pass),
          num_fail = sum(!trainload_pass),
          frac_pass = sum(trainload_pass)/n())



# Create figures ----------------------------------------------------------

## Fig S1: Num unsafe trainloads per parameter combination -----------------

# Average number of safe vs unsafe railcars for each parameter combination
# (should be same on average across values of num_railcars_per_composite)
plotdat <- res2 %>% 
  group_by(prob_of_zero, shape_param, shape_plot, scale_plot, train_is_safe) %>% 
  reframe(num_simulated_trains = n(), pct_of_simulated_trains = 100*n()/(500*length(unique(summary.tb$num_railcars_per_composite))))
plotdat <- plotdat %>% 
  mutate(shape = as.factor(shape_param))

# Create plot
ggplot(plotdat %>% filter(prob_of_zero == 0.05), aes(x = shape, y = pct_of_simulated_trains, fill = train_is_safe))+
  scale_fill_viridis_d(option = "magma", name = "train is safe")+
  geom_col()+
  facet_wrap(vars(scale_plot))+
  xlab("shape")+
  ylab("% of simulated trains")+
  theme_bw()+
  theme(text = element_text(size = 12))
ggsave("figs/Figure_S1_pct_simulated_trainloads_safe_vs_unsafe_prob_of_zero5pct.png", width = 6, height = 5)


# compare values when prob of zero is 0%
ggplot(plotdat %>% filter(prob_of_zero == 0), aes(x = shape, y = pct_of_simulated_trains, fill = train_is_safe))+
  scale_fill_viridis_d(option = "magma")+
  geom_col()+
  facet_wrap(vars(scale_plot))+
  xlab("shape")+
  ylab("% of simulated trains")+
  theme_bw()
# Results are very similar, will not show this figure

## Fig S2: Compare consumers risk for prob of zero 0 vs 5% --------------------------------------------

plotdat <- summary.tb %>% 
  filter(!train_is_safe) %>% 
  pivot_wider(id_cols = c(shape_param, scale_param, num_railcars_per_composite), names_from = prob_of_zero, names_prefix  = "p", values_from = frac_pass)

ggplot(plotdat, aes(x  = p0, y = `p0.05`, color = as.factor(num_railcars_per_composite)))+
  geom_point()+
  scale_color_discrete(name = "# railcars per\ncomposite sample")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  theme(text = element_text(size = 12))+
  xlab("Consumer's risk with no zero-inflation")+
  ylab("Consumer's risk with 5% zero-inflation")
ggsave(filename = "figs/Figure_S2_consumers_risk_for_p0_vs_p0.05.png", width = 6, height = 4)

# Okay, so the type 1 error rates do differ in a few cases
plotdat %>% 
  filter(abs(p0 - `p0.05`) > 0.05)
# And for most of the differences, the 5% version has slightly higher type 1 error rates

plotdat %>% 
  filter(is.na(p0) | is.na(`p0.05`))
# 8 occurences, all with scale = 1 or 2
# I'm guessign this is just due to places where one sim had 1 or 2 unsafe trainloads, the other setting did not
# so no real diffeeence here


## Fig 4: Consumer's risk at different combos ----------------------------------------------------------

# Select subset of data plot
plotdat <- summary.tb %>% 
  filter(!train_is_safe, scale_param > 1, prob_of_zero == 0.05)
plotdat <- plotdat %>% 
  mutate(shape = as.factor(shape_param))

# facet by scale, color by shape
ggplot(plotdat, 
       aes(x = num_railcars_per_composite, y = frac_pass, color = shape))+
  geom_point()+
  geom_line(aes(group = parameter_id_sans_nrpc))+
  #scale_color_viridis_b(name = "shape")+
  facet_wrap(vars(scale_plot))+
  ylab("Consumer's risk")+
  xlab("Number of railcars per composite sample")+
  theme_bw()+
  theme(text = element_text(size = 12))

ggsave(filename = "figs/Figure_4_consumers_risk_vs_num_railcars_per_composite_facet_scale.png", width = 8, height = 5)


## Fig 5: Average number of unsafe railcars ----------------------------------------------------------

# Get data to plot
plotdat <- res2 %>% 
  filter(trainload_pass == TRUE, scale_param > 1) %>% 
  group_by(shape_param, scale_param, scale_plot, prob_of_zero, num_railcars_per_composite, parameter_id_sans_nrpc) %>% 
  reframe(avg_num_unsafe_railcars = mean(num_unsafe_railcars))
plotdat$shape <- as.factor(plotdat$shape_param)

# facet by scale, color by shape
ggplot(plotdat %>% filter(prob_of_zero == 0.05), 
       aes(x = num_railcars_per_composite, y = avg_num_unsafe_railcars, color = shape))+
  geom_point()+
  geom_line(aes(group = parameter_id_sans_nrpc))+
  facet_wrap(vars(scale_plot))+
  ylab("Avg # unsafe railcars\nin trainloads labelled as safe")+
  xlab("Number of railcars per composite sample")+
  theme_bw()+
  theme(text = element_text(size = 12))

ggsave(filename = "figs/Figure_5_avg_num_unsafe_railcars_in_accepted_trainloads.png", width = 8, height = 5)

# compare with prob of zero == 0
ggplot(plotdat %>% filter(prob_of_zero == 0), 
       aes(x = num_railcars_per_composite, y = avg_num_unsafe_railcars, color = shape))+
  geom_point()+
  geom_line(aes(group = parameter_id_sans_nrpc))+
  facet_wrap(vars(scale_plot))+
  ylab("Avg # unsafe railcars\nin trainloads labelled as safe")+
  xlab("Number of railcars per composite sample")+
  theme_bw()+
  theme(text = element_text(size = 12))
# Main difference: there are more points where 0 railcars were accepted
# (so there is no points on the graph)
# -> greater conservatism


## Fig 6: Example trainloads -----------------------------------------------

example.dat <- sapply(2:5, simulate_trainload, nrailcars_per_train = 200,
                      prob_of_zero = 0.05, shape_param = 4)
colnames(example.dat) <- paste0("scale: ",2:5)
example.dat <- as.data.frame(example.dat)

example.dat <- example.dat %>% 
  pivot_longer(paste0("scale: ",2:5), values_to = "ppb")

ggplot(example.dat, aes(x = ppb)) +
  geom_histogram()+
  geom_vline(xintercept = 20, color = "blue")+
  facet_wrap(vars(name))+
  xlab("Aflatoxin level (ppb)")+
  #ylab("Number of railcars")+
  theme_bw()+
  theme(text = element_text(size = 12))

ggsave(filename = "figs/Figure_6_example_trainload_dist_prob_0.05_shape_4_scale_facet.png", width = 6, height = 4)



# Generating numbers for results ------------------------------------------

## Number of unsafe vs safe trainloads ---------------------
summary.tb %>% 
  filter(scale_param %in% c(1)) %>% 
  select(train_is_safe, total_trains) %>% 
  arrange(total_trains)
# ONly 4 cases have 1 unsafe trainload

summary.tb %>% 
  filter(scale_param %in% c(10,20)) %>% 
  select(train_is_safe, total_trains) %>% 
  reframe(unique(train_is_safe))
# All FALSE


## safe trainloads identified as safe ----------------------

# See percent of safe trainloads accepted
summary.tb %>% 
  filter(train_is_safe) %>%  # selects safe trainloads
  select(frac_pass) %>% 
  summary
# frac_pass
# Min.   :1  
# 1st Qu.:1  
# Median :1  
# Mean   :1  
# 3rd Qu.:1  
# Max.   :1

# For every parameter combination, 100% of safe trainloads were accepted



## Avg number of unsafe railcars when scale == 2 ----------------


plotdat <- res2 %>% 
  filter(trainload_pass == TRUE) %>% 
  group_by(shape_param, scale_param, scale_plot, prob_of_zero, num_railcars_per_composite) %>% 
  reframe(avg_num_unsafe_railcars = mean(num_unsafe_railcars))

plotdat %>% 
  filter(scale_param == 2, prob_of_zero == 0.05) %>% 
  select(shape_param, scale_param, scale_plot, prob_of_zero, num_railcars_per_composite, avg_num_unsafe_railcars) %>% 
  arrange(-avg_num_unsafe_railcars)
#   shape_param scale_param scale_plot prob_of_zero num_railcars_per_composite avg_num_unsafe_railcars
#         <dbl>       <dbl> <fct>             <dbl>                      <dbl>                   <dbl>
# 1           4           2 scale: 2           0.05                          2                   2.00 
# 2           4           2 scale: 2           0.05                         20                   1.99 
# 3           4           2 scale: 2           0.05                         10                   1.98 
# 4           4           2 scale: 2           0.05                          5                   1.97 
# 5           4           2 scale: 2           0.05                         40                   1.88 
# 6           3           2 scale: 2           0.05                          5                   0.574
