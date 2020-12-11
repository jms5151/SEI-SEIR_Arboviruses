# run linear regressions for epidemic characteristics

# load data
epid_char_1 <- read.csv("data/Number_outbreaks.csv")
epid_char_2 <- read.csv("data/Duration_intensity_outbreaks.csv")

# Number of outbreaks
summary(lm(epid_char_1$Number_outbreaks_predicted ~ epid_char_1$Number_outbreaks_observed))

# Peak timing of outbreaks
summary(lm(epid_char_2$Peak_timing_predicted ~ epid_char_2$Peak_timing_observed))

# Duration of outbreaks
summary(lm(epid_char_2$Duration_predicted ~ epid_char_2$Duration_observed))

# Maximum infections during outbreaks
summary(lm(epid_char_2$Max_infected_predicted ~ epid_char_2$Max_infected_observed))

# Outbreak size
summary(lm(epid_char_2$Outbreak_size_predicted ~ epid_char_2$Outbreak_size_observed))
