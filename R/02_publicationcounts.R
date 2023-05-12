# load packages
library(dplyr)

# import data
analysis_dat = readRDS("Data/analysis_dat.rds")
str(analysis_dat)

## Total---
dim(analysis_dat)
length(unique(analysis_dat$pub_ID))

###Species--------
spp <- levels(analysis_dat$sp_name)

n.studies= data.frame(
  species = spp,
  n.stud = rep(NA,length(spp))
)

i=1
for(p in spp){
  sub = analysis_dat$pub_ID[analysis_dat$sp_name == p]
  n.studies$n.stud[i] = length(unique(sub))
  i = i+1
} 
n.studies


# number of data points:
analysis_dat %>% group_by(sp_name) %>% summarize(n = n())


####Age-----
# number of data points:
analysis_dat %>% group_by(age_d) %>% summarize(n = n())


###Food-------
#number of studies with food provided
length(unique(analysis_dat$pub_ID[which(analysis_dat$food_present=="yes")]))
# number of samples with food provided
nrow(analysis_dat[which(analysis_dat$food_present=="yes"),])
#number of studies without food provided
length(unique(analysis_dat$pub_ID[which(analysis_dat$food_present=="no")]))
# number of samples without food provided
nrow(analysis_dat[which(analysis_dat$food_present=="no"),])


###Exposure duration--------
analysis_dat$exp_time_d = as.factor(analysis_dat$exp_time_d)
dur <- levels(analysis_dat$exp_time_d)

n.studies= data.frame(
  exp.dur = dur,
  n.stud = rep(NA,length(dur))
)

i=1
for(p in dur){
  sub = analysis_dat$pub_ID[analysis_dat$exp_time_d == p]
  n.studies$n.stud[i] = length(unique(sub))
  i = i+1
} 

n.studies

# number of data points:
analysis_dat %>% group_by(exp_time_d) %>% summarize(n = n())

### Polymer type-------
#number of studies per polymer type:
poly <- levels(analysis_dat$polymer)

n.studies= data.frame(
  polymer = poly,
  n.stud = rep(NA,length(poly))
)

i=1
for(p in poly){
 sub = analysis_dat$pub_ID[analysis_dat$polymer == p]
 n.studies$n.stud[i] = length(unique(sub))
 i = i+1
} 
n.studies

# number of polymer types (data points):
analysis_dat %>% group_by(polymer) %>% summarize(n = n())


### Shape-----
#number of studies for each shape
length(unique(analysis_dat$pub_ID[which(analysis_dat$shape=="Aggregates")]))
length(unique(analysis_dat$pub_ID[which(analysis_dat$shape=="Fibers")]))
length(unique(analysis_dat$pub_ID[which(analysis_dat$shape=="Fragments")]))
length(unique(analysis_dat$pub_ID[which(analysis_dat$shape=="Spheres")]))

# number of samples for each shape
nrow(analysis_dat[which(analysis_dat$shape=="Aggregates"),])
nrow(analysis_dat[which(analysis_dat$shape=="Fibers"),])
nrow(analysis_dat[which(analysis_dat$shape=="Fragments"),])
nrow(analysis_dat[which(analysis_dat$shape=="Spheres"),])

###Size----------
#number of studies with MP vs. NP
length(unique(analysis_dat$pub_ID[which(analysis_dat$MP.NP=="MP")]))
length(unique(analysis_dat$pub_ID[which(analysis_dat$MP.NP=="NP")]))
# number of samples with modified MNP
nrow(analysis_dat[which(analysis_dat$MP.NP=="MP"),])
nrow(analysis_dat[which(analysis_dat$MP.NP=="NP"),])


### Controlled surface and other modifications---------
#number of studies with modifications
length(unique(analysis_dat$pub_ID[which(analysis_dat$modified=="yes")]))
# number of samples with modified NMP
nrow(analysis_dat[which(analysis_dat$modified=="yes"),])

# number of data points for  different modification types:
analysis_dat %>% group_by(mod_type) %>% summarize(n = n())

# number of studies per modification type:
# 1. additives: number of studies
length(unique(analysis_dat$pub_ID[which(analysis_dat$mod_type=="BP-3, melt-blended")]))
length(unique(analysis_dat$pub_ID[which(analysis_dat$mod_type=="13C")]))
#2. fluorescence tags:
length(unique(analysis_dat$pub_ID[which(analysis_dat$fluorescence=="yes")]))# studies
nrow(analysis_dat[which(analysis_dat$fluorescence=="yes"),])# samples
#3. controlled surface modifications: number of studies:
length(unique(analysis_dat$pub_ID[which(analysis_dat$mod_type=="Amidine derivatized")]))
length(unique(analysis_dat$pub_ID[which(analysis_dat$mod_type=="Aminated")]))
length(unique(analysis_dat$pub_ID[which(analysis_dat$mod_type=="Carboxylated")]))
length(unique(analysis_dat$pub_ID[which(analysis_dat$mod_type=="Extracted")]))

### Uncontrolled surface mods------

#number of studies with DOC
length(unique(analysis_dat$pub_ID[which(analysis_dat$DOC_added=="yes")]))
# number of samples with DOC
nrow(analysis_dat[which(analysis_dat$DOC_added=="yes"),])

length(unique(analysis_dat$pub_ID[which(analysis_dat$DOC_type=="fulvic acid")]))
length(unique(analysis_dat$pub_ID[which(analysis_dat$DOC_type=="humic acid")]))
length(unique(analysis_dat$pub_ID[which(analysis_dat$DOC_type=="lake water")]))

#samples per DOC_type:
analysis_dat %>% group_by(DOC_type) %>% summarize(n = n())

# incubation of MP in biotic environment prior to exposure
length(unique(analysis_dat$pub_ID[which(analysis_dat$biotic_env != "no")])) #studies
analysis_dat %>% group_by(biotic_env) %>% summarize(n=n()) #samples

## END --------


