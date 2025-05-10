
# ------------------------------------------------------------------------ # \
# Generating plots to share with Marissa
# May 10, 2025

library(ggplot2)
library(tidyr)
library(dplyr)

# Normal dist -------------------------------------------------------------

xvals <- seq(5, 25, by = 0.1)

norm.tb <- data.frame(x = rep(xvals, times = 3),
                      label = rep(LETTERS[1:3], each = length(xvals)),
                      prob_density = c(dnorm(x = xvals, mean = 15, sd = 1),
                                       dnorm(x = xvals, mean = 10, sd = 1),
                                       dnorm(x = xvals, mean = 15, sd = 2)))

norm.tb <- norm.tb %>% 
  mutate(mean = rep(c(15, 10, 15), each = length(xvals)))
norm.tb<- norm.tb %>% 
  mutate(mean_prob_density = dnorm(x = mean, mean = mean, sd = rep(c(1,1,2), each = length(xvals))))

ggplot(norm.tb, aes(x = x, y= prob_density)) +
  geom_area(fill = "cornflowerblue")+
  geom_line()+
  geom_segment(aes(x = mean, y = 0, yend = mean_prob_density), linetype = "dotted")+
  ylab("probability density")+
  facet_grid(cols = vars(label))

ggsave(filename = "figs/normal_distribution_example.png", width = 6, height = 3)


# Gamma dist -------------------------------------------------------------

xvals <- seq(0, 25, by = 0.1)
shape_vals <- c(1, 2, 3)
scale_vals <- c(1, 2, 3)

gamma.tb <- expand.grid(shape = shape_vals, 
                        scale = scale_vals,
                        x = xvals)
gamma.tb$prob_density <- dgamma(x = gamma.tb$x, 
                               shape = gamma.tb$shape,
                               scale = gamma.tb$scale)

gamma.tb <- gamma.tb %>% 
  mutate(mean = shape*scale,
         variance = shape*scale*scale,
         shape_label = paste0("shape: ",shape),
         scale_label = paste0("scale: ",scale))

# Get the prob density at the max
gamma.tb <- gamma.tb %>% 
  mutate(mean_prob_density = dgamma(x = mean, shape = shape, scale = scale))

ggplot(gamma.tb, aes(x = x, y= prob_density)) +
  geom_area(fill = "cornflowerblue")+
  geom_line()+
  #geom_vline(aes(xintercept = mean), color = "blue")+
  #geom_point(aes(x = mean, y = 0), color = "blue")+
  geom_segment(aes(x = mean, y = 0, yend = mean_prob_density))+
  ylab("probability density")+
  ylim(0, 0.75)+
  facet_grid(cols = vars(shape_label), rows = vars(scale_label))+
  theme(text = element_text(size =15))

ggsave(filename = "figs/gamma_distribution_example.png", width = 6, height = 6)


# GAmma, v2 ---------------------------------------------------------------

xvals <- seq(0.1, 30, by = 0.1)
shape_vals <- c(1, 2, 3)
scale_vals <- c(2, 5)

gamma.tb <- expand.grid(shape = shape_vals, 
                        scale = scale_vals,
                        x = xvals)
gamma.tb$prob_density <- dgamma(x = gamma.tb$x, 
                                shape = gamma.tb$shape,
                                scale = gamma.tb$scale)

gamma.tb <- gamma.tb %>% 
  mutate(mean = shape*scale,
         variance = shape*scale*scale,
         shape_label = paste0("shape: ",shape),
         scale_label = paste0("scale: ",scale))

# Get the prob density at the max
gamma.tb <- gamma.tb %>% 
  mutate(mean_prob_density = dgamma(x = mean, shape = shape, scale = scale))

ggplot(gamma.tb, aes(x = x, y= prob_density)) +
  geom_area(fill = "cornflowerblue")+
  geom_line()+
  #geom_vline(aes(xintercept = mean), color = "blue")+
  #geom_point(aes(x = mean, y = 0), color = "blue")+
  geom_segment(aes(x = mean, y = 0, yend = mean_prob_density), linetype = "dotted")+
  ylab("probability density")+
  xlab("aflatoxin level (ppb)")+
  facet_grid(cols = vars(shape_label), rows = vars(scale_label))
  theme(text = element_text(size =15))

ggsave(filename = "figs/gamma_distribution_example_v2.png", width = 6, height = 4)

