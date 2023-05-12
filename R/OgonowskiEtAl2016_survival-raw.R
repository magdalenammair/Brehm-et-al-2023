library(dplyr)

dat <- read.csv2("OgonowskiEtAl2016_Supplement/Exp_I_reproduction_data.csv")
str(dat)

dat[dat$Type == "Control",]

table(dat$Type, dat$Concentration) # sample size n = 10

# calculate numbers of individuals that were alive and dead after 21 days of exposure:
# for each concentration ("Concentration") and treatment ("Type")

dat %>% group_by(Type, Concentration) %>% summarize(alive = length(surv.days[surv.days == 21]), dead = 10-alive)

