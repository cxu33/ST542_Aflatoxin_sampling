# -------------------------------------------------------------------------# 
# Real data analysis
# -------------------------------------------------------------------------# 

library(tidyverse)
library(ggplot2)
library(MASS)
library(readxl)


# Prepare the data --------------------------------------------------------

dat <- read_xlsx("data/Aflatoxin Testing Results 2024-2025_MHC.xlsx")

# Clean column names
dat <- dat %>% 
  rename(Date = `Sample Date`, ticket_num = `Ticket #`, aflatoxin_ppb = `Aflatoxins, ppb`)


# Figure 1 ----------------------------------------------------------------


ggplot(dat, aes(x = Date, y = aflatoxin_ppb))+
  geom_point(pch = 1)+
  ylab("aflatoxin level (ppb)")+
  theme_bw()+
  theme(text = element_text(size = 12))

ggsave("figs/Figure_1_real_data_analysis_vs_time.png", width = 5, height = 5)


#  Figure 2 ---------------------------------------------------------------

# Fit gamma model to the positive values

dat_pos <- dat %>% 
  filter(aflatoxin_ppb > 0)

fit.gamma <- fitdistr(dat_pos$aflatoxin_ppb, densfun = "gamma")


# Histogram of full data, show gamma fit

# Get values from gamma distribution
xvals <- seq(1e-4, 20, by = 0.1)
gamma_vals <- dgamma(x = xvals, 
        shape = fit.gamma$estimate[1],
       rate = fit.gamma$estimate[2])
gamma.tb <- data.frame(x = xvals, prob_density = gamma_vals)

# Graph
ggplot(dat, aes(x = aflatoxin_ppb)) +
  geom_histogram(mapping = aes(y = after_stat(density)), 
                 breaks = seq(from = -1, to = max(dat$aflatoxin_ppb), by = 0.5), 
                 closed = "right")+
  geom_line(data = gamma.tb, aes(x = x, y = prob_density), col = "cornflowerblue", linewidth = 1.5)+
  theme(text = element_text(size = 12))+
  xlab("aflatoxin level (ppb)")

ggsave(filename = "figs/Figure_2_raw_data_full_histogram_zeros_separate_bin_with_gamma.png", width = 6, height = 6)

# View relevant parameters
data.frame(prob_of_zero = sum(dat$aflatoxin_ppb == 0)/nrow(dat),
           shape = fit.gamma$estimate["shape"],
           scale = 1/fit.gamma$estimate["rate"])

#      prob_of_zero    shape    scale
# shape    0.2197309 1.356312 2.370728
