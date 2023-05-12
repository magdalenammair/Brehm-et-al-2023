# load packages
library(dplyr)
library(metafor)

# import data
dat <- read.csv("Data/input_data.csv", na.strings = c("", "NA"), stringsAsFactors = TRUE)

# check variable classes
# str(dat)

# number of studies for which data extraction was tried:
length(unique(dat$pub_ID))

# remove all samples where data extraction was not possible (included = "no"):
dat <- dat[dat$included == "yes",]

# number of studies for which data extraction was successful:
length(unique(dat$pub_ID)) 

## Check variable classes:-------------------------
# change all integer classes to numeric:
dat[sapply(dat, is.integer)] <- lapply(dat[sapply(dat, is.integer)], as.numeric)
# str(dat)


## Prepare response -----
# calculate number dead/alive (immobile/mobile) from mean dead for control and treatment:

sel = is.na(dat$c_dead) & !is.na(dat$c_mean)
dat$c_dead[sel] = round(dat$c_mean[sel] * dat$c_n_individuals[sel])
dat$c_alive[sel] = dat$c_n_individuals[sel] - dat$c_dead[sel]
rm(sel)

sel = is.na(dat$t_dead) & !is.na(dat$t_mean)
dat$t_dead[sel] = round(dat$t_mean[sel] * dat$t_n_individuals[sel])
dat$t_alive[sel] = dat$t_n_individuals[sel] - dat$t_dead[sel]
rm(sel)

# remove lines without survival data
dat <- dat[!is.na(dat$t_dead) & !is.na(dat$c_dead),]


## Prepare predictors: --------------------------
# calculate mg per ml from mg per litre samples:
sel = !is.na(dat$conc_mg_l.1) & is.na(dat$conc_mg_ml.1)
dat$conc_mg_ml.1[sel] <- dat$conc_mg_l.1[sel]/1000
rm(sel)

# calculate mean nm size from mean um size:
sel = is.na(dat$size_nm_mean) & !is.na(dat$size_um_mean)
dat$size_nm_mean[sel] = dat$size_um_mean[sel] * 1000
rm(sel)

# calculate mean um size from mean nm size:
sel = !is.na(dat$size_nm_mean) & is.na(dat$size_um_mean)
dat$size_um_mean[sel] = dat$size_nm_mean[sel] / 1000
rm(sel)

# make new categorical column with MP (more or equal to 100nm) vs. NP (less than 100 nm):

dat$MP.NP = NA
dat$MP.NP[dat$size_nm_mean < 1000] = "NP" 
dat$MP.NP[dat$size_nm_mean >= 1000] = "MP"
dat$MP.NP = as.factor(dat$MP.NP)

# fix shape names:
levels(dat$shape)
levels(dat$shape)[levels(dat$shape) == "fragments"] <-  "Fragments"
levels(dat$shape)[levels(dat$shape) == "fiber"] <- "Fibers"
levels(dat$shape)[levels(dat$shape) == "aggregates"] <- "Aggregates"
levels(dat$shape)[levels(dat$shape) == "spherical"] <- "Spheres"
levels(dat$shape)[levels(dat$shape) == "spheres"] <- "Spheres"
dat$shape[dat$shape == "fragments; size range: 2-60"] <- "Fragments"
dat$shape = factor(dat$shape)

#fix polymer types:
levels(dat$polymer)
dat$polymer[dat$polymer == "PS "] <- "PS"
# remove non-plastic particles
dat <- subset(dat, dat$polymer != "fumed silica")

levels(dat$polymer)[levels(dat$polymer)=="Thermoset amino formaldehyde polymer."] <- "other"
levels(dat$polymer)[levels(dat$polymer)=="unidentified"] <- "other"

# check sample size per polymer type:
poly.samples = dat %>% group_by(polymer) %>% summarize(n = n())
# pool all polymer types with <=2 data points as "other"
other = as.vector(poly.samples$polymer[poly.samples$n <= 2])
levels(dat$polymer)[levels(dat$polymer) %in% other] <- "other"

# change to class factor
dat$polymer <- factor(dat$polymer)
# remove unused factor levels
dat <- droplevels(dat)
# check again
levels(dat$polymer)

#fix species names
dat$sp_name[dat$sp_name == "Daphnia magna "] <- "Daphnia magna"
dat$sp_name <- factor(dat$sp_name)
levels(dat$sp_name)

#fix modification types
levels(dat$mod_type)[levels(dat$mod_type) == "aminified"] <- "Aminated"
levels(dat$mod_type)[levels(dat$mod_type) == "carboxylated"] <- "Carboxylated"
levels(dat$mod_type)[levels(dat$mod_type) == "extracted"] <- "Extracted"
levels(dat$mod_type)[levels(dat$mod_type) == "amidine derivatized"] <- "Amidine derivatized"
levels(dat$mod_type)

# create col_yes_no:
levels(dat$color)
dat$col_yes_no = NA
dat$col_yes_no[dat$color == "no" | dat$color == "clear"] = "no"
dat$col_yes_no[dat$color != "no"] = "yes"
dat$col_yes_no <- as.factor(dat$col_yes_no)

# create sample_ID:
dat$sample_ID = seq(1, nrow(dat),1)

# check again:
# str(dat)
# summary(dat)
# dim(dat)

# Subset dataset to include only columns used in meta-analysis:
datsub <- subset(dat, select = c("sample_ID", "pub_ID", "sp_name", "polymer", "size_um_mean", "size_nm_mean", "shape", "color", "col_yes_no", 
                                 "fluorescence", "modified", "mod_type", "surface_charge", "surfactant", "biotic_env", 
                                 "density_gcm.3", "age_d", "food_present", "DOC_added", "DOC_type", "other_added", 
                                 "weathered","exp_time_d", 
                                 "conc_p_ml.1", "conc_mg_ml.1", "MP.NP", "temperature",
                                 "c_dead", "c_alive", "t_dead", "t_alive"))
str(datsub)
summary(datsub)

## Calculate log(Risk Ratios): -----------
analysis_dat <- escalc(measure="RR", ai = t_dead, bi = t_alive, ci = c_dead, di = c_alive, data = datsub, slab = pub_ID)
str(analysis_dat)

# sampleID should be a factor:
analysis_dat$sample_ID <- as.factor(analysis_dat$sample_ID)

# calculate standard error and inverse standard error for plotting:
analysis_dat$se <- analysis_dat$vi/(analysis_dat$c_dead + analysis_dat$c_alive +analysis_dat$t_dead +analysis_dat$t_alive)
analysis_dat$invse <- 1/analysis_dat$se
analysis_dat$invse.by.10 <- 1/analysis_dat$se/10
analysis_dat$log.invse <- log(analysis_dat$invse)

# Calculate concentrations in mg per ml: estimate mg per ml from particles per ml, size and density
df = analysis_dat
## Samples that do not report concentrations in p per ml or ml per ml
sel = is.na(df$conc_mg_ml.1) & is.na(df$conc_p_ml.1)
sum(sel) # samples report neither mg per ml nor particles per ml; instead they report ppm
# remove these samples from dataset
df = df[!sel,]
rm(sel)

# Calculate particle volume for different shapes according to Thornton-Hampton et al. 2022:
df$p.vol[!is.na(df$shape) & df$shape == "Spheres"] = 4/3 * pi * (df$size_um_mean[!is.na(df$shape) & df$shape == "Spheres"]/2000)^3  #unit:mm^3
df$p.vol[!is.na(df$shape) & df$shape == "Fragments"] = pi/6 * (df$size_um_mean[!is.na(df$shape) & df$shape == "Fragments"]/1000)^3 * 0.4^2 # unit:mm^3
# see: Kooi M, Koelmans AA. Simplifying microplastic via continuous probability distributions for size, shape, and density. Environ Sci Technol Lett. 2019;6(9):551–7. 
# and: Thornton Hampton, Leah M., Heili Lowman, Scott Coffin, Emily Darin, Hannah De Frond, Ludovic Hermabessiere, Ezra Miller, Vera N. de Ruijter, Andrea Faltynkova, Syd Kotar, Laura Monclús, Samreen Siddiqui, Johannes Völker, Susanne Brander, Albert A. Koelmans, Chelsea M. Rochman, Martin Wagner, and Alvine C. Mehinto. 2022. ‘A Living Tool for the Continued Exploration of Microplastic Toxicity’. Microplastics and Nanoplastics 2(1):13. doi: 10.1186/s43591-022-00032-4.
unique(df$pub_ID[df$shape == "Fibers"]) # 
# calculate fiber volumes:
df$p.vol[df$shape=="Fibers" & df$pub_ID == "Kokalj2018"] = pi * (df$size_um_mean[df$shape == "Fibers" & df$pub_ID == "Kokalj2018"]/2000) * 2 * 0.3 #unit:mm^3 length: Kokalj2018: 300µm (mean),
df$p.vol[df$shape=="Fibers" & df$pub_ID == "Zocchi2019"] = pi * (df$size_um_mean[df$shape == "Fibers" & df$pub_ID == "Zocchi2019"]/2000) * 2 * 0.01 #unit:mm^3 length: Zocchi2019: 10µm
df$p.vol[df$shape=="Fibers" & df$pub_ID == "Tourinho2021"] = pi * (df$size_um_mean[df$shape == "Fibers" & df$pub_ID == "Tourinho2021"]/2000) * 2 * 0.366 # unit:mm^3 Tourinho2021: fiber length: 366µm
df$p.vol[df$shape=="Fibers" & df$pub_ID == "Schwarzer2022"] = pi * (df$size_um_mean[df$shape == "Fibers" & df$pub_ID == "Schwarzer2022"]/2000) * 2 * 0.075 # unit:mm^3 Schwarzer2022: fiber length: 75µm
df$p.vol[df$shape=="Fibers" & df$pub_ID == "Kim2021"] = pi * (df$size_um_mean[df$shape == "Fibers" & df$pub_ID == "Kim2021"]/2000) * 2 * 0.120 # unit:mm^3 Kim2021: fiber length: 120µm

# Calculate particle weight in mg
# particle weight(mg) = density(mg/mm^3) * Vol(mm^3)
df$p.weight = df$density_gcm.3 * df$p.vol #unit: mg

# get all samples where concentration in mg ml is missing, but particles per ml and density are available:
sel = is.na(df$conc_mg_ml.1) & !is.na(df$conc_p_ml.1) & !is.na(df$density_gcm.3)
sum(sel)
# Calculate concentration in mg per ml for these samples:
# mg/ml = number of particles per ml * particle weight:
df$conc_mg_ml.1[sel] = df$p.weight[sel] * df$conc_p_ml.1[sel]
rm(sel)

#get all samples where concentration in particles per ml is missing, but mg per ml and particle weight are available
sel = is.na(df$conc_p_ml.1) & !is.na(df$conc_mg_ml.1) & !is.na(df$p.weight)
# Calculate particles per ml from mg per ml and particle weight:
df$conc_p_ml.1[sel] = df$conc_mg_ml.1[sel]/df$p.weight[sel]
rm(sel)
# 
sum(is.na(df$conc_mg_ml.1)) # 4 samples still do not have concentration measures in mg per ml 
df[is.na(df$conc_mg_ml.1),] # 2x Schrank2019, 2xPikuda2019

# remove these data points:
df = df[!is.na(df$conc_mg_ml.1),]
df = droplevels(df)
# summary(df)
# str(df)

#save data frame for use in plotting and modelling:
saveRDS(df, "Data/analysis_dat.rds")


### END-------------------
