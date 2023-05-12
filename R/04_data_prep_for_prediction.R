# load packages
library(missRanger)

# import data:
analysis_dat = readRDS("Data/analysis_dat.rds")
df = as.data.frame(analysis_dat)
str(df)

# 1.Handling NA values-------
# A. Change NA to additional factor level "not reported": food_present, DOC_added, other_added, color, col_yes_no, fluorescence, modified, surface_charge, weathered

NA.to.notreported = function(dframe, colname) {
  dframe[,colname] <- as.character(dframe[,colname])
  dframe[is.na(dframe[,colname]),colname] <- "not reported"
  dframe[,colname] <- as.factor(dframe[,colname])
}

columns = c("food_present", "DOC_added", "other_added", "weathered", "color", "col_yes_no", "fluorescence", "modified", "surface_charge")
for(i in 1:length(columns)){
  df[, columns[i]] = NA.to.notreported(dframe = df, colname = columns[i])
}
summary(df)

# B. Change NA to additional factor "none" for: mod_type, DOC_type
df[,"mod_type"] <- as.character(df[,"mod_type"])
df[is.na(df[,"mod_type"]),"mod_type"] <- "none"
df[,"mod_type"] <- as.factor(df[,"mod_type"])

df[,"DOC_type"] <- as.character(df[,"DOC_type"])
df[is.na(df[,"DOC_type"]),"DOC_type"] <- "none"
df[,"DOC_type"] <- as.factor(df[,"DOC_type"])

summary(df)

# 2. Check concentrations ----
# check again: do all samples have values in mg per ml and particles per ml?

## Samples that do not report concentrations in p per ml or ml per ml have already been removed in "01_data_preparation.R"
sel = is.na(df$conc_mg_ml.1) & is.na(df$conc_p_ml.1)
sum(sel) #0 

# Samples where concentration in particles per ml is missing, but mg per ml and particle weight are available
sel = is.na(df$conc_p_ml.1) & !is.na(df$conc_mg_ml.1) & !is.na(df$p.weight)
sum(sel) #0


# 3. Scale numeric predictors -----
df[,c("age_d", "exp_time_d", "temperature", "conc_p_ml.1", "conc_mg_ml.1", "size_um_mean", "density_gcm.3")] <- scale(df[,c("age_d", "exp_time_d", "temperature", "conc_p_ml.1", "conc_mg_ml.1", "size_um_mean", "density_gcm.3")])
summary(df)

# 4. Select variables for final dataset------
df = subset(df, select = c("pub_ID",
                            "sp_name", "age_d", "food_present", "DOC_added", "DOC_type", "other_added", "weathered", "exp_time_d", "temperature",
                            "conc_p_ml.1", "conc_mg_ml.1",
                            "polymer", "size_um_mean", "MP.NP", "shape", "color", "col_yes_no", "fluorescence", "modified", "mod_type",
                            "surface_charge", "surfactant", "biotic_env", "density_gcm.3",
                            "yi", "vi"))
summary(df)

# Save dataset with NA values
saveRDS(df, "Data/notimp_pred_dat.rds")

# Prep for imputation:-------

# Put response values yi into separate vector:
response = subset(df, select = c("yi"))
summary(response)

# Remove pub_ID and yi for imputation:
df = subset(df, select = -c(pub_ID, yi, vi))
str(df)

## Multiple Imputation:-------
imputed = missRanger::missRanger(df)
summary(imputed)


# add response again:
imputed = cbind(imputed,response)

# save dataset:

saveRDS(imputed, "Data/imputed_dat.rds")

# END----

