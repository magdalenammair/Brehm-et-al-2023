### Calculate mean proportion and standard deviation of dead Daphnia neonates per treatment

neon_mort <- read.csv("Neon_mort_MM.csv")
str(neon_mort)

neon_mort$Treatments <- as.factor(neon_mort$Treatments)
neon_mort$prop.dead <- neon_mort$dead/5

library(dplyr)

neon <- neon_mort[neon_mort$day == 21,] %>% 
  group_by(Treatments) %>%
  summarise(mean = mean(prop.dead), sd = sd(prop.dead))
neon


### Adult mortality:here, only 1 individual per beaker was tested

adult_mort <- read.csv("adult_mort_MM.csv")
adult_mort <- adult_mort[adult_mort$Time..day.==21, 3:8]
stacked <- stack(adult_mort)

adult <- stacked %>% 
  group_by(ind) %>%
  summarise(mean = mean(values), sd = sd(values), n = n())
adult

