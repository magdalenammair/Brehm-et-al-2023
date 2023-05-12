# load packages
library(dplyr)
library(metafor)
### Load packages ----------
library(RColorBrewer)
library(svglite)

### Import helper functions -----
source("R/helper.R")

# Set general color palette:
col_vector<-c('#000075' , '#3cb44b', '#808000', '#4363d8', '#f58231', '#911eb4', '#46f0f0','#ffe119', '#e6194b','#9a6324' ,
               '#008080', '#bcf60c', '#ffd8b1', '#fffac8', '#800000', '#e6beff','#aaffc3' , '#fabebe', '#f032e6', '#808080', '#ffffff', '#000000')
palette(col_vector)

###Import data--------
df = readRDS("Data/analysis_dat.rds")

## Models and plots--------------------

## Code structure for all parameters is similar.##
## Check polymer plot for detailed coding comments.##

### 1. Species---------

# check order of factor levels
levels(df$sp_name)
#change order 
df$sp_name = factor(df$sp_name, levels(df$sp_name)[c(2,3,1)])
# check again
levels(df$sp_name)


svglite("Plots/regplot_species.svg", height = 5, 
        width = 5, pointsize = 20)

par(mar=c(4,4,0.5,2), xpd = TRUE, mfrow = c(1,1), las =1, bty = "l",
    cex = 0.8)

palette(col_vector[c(5,2,4)])

plot(yi ~ log(conc_mg_ml.1), data = df, pch = 21, 
     col = adjustcolor(as.numeric(df$sp_name), alpha.f = 0.1), 
     bg = adjustcolor(as.numeric(df$sp_name), alpha.f = 0.1), cex = log.invse/4,
     ylab = "", xlab = "", ylim = c(-2.2,7))
mtext("log (Risk Ratio)", side = 2, line = 2, cex = 0.8, las = 3)
mtext("Concentration in mg per ml (log-scale)", side = 1, line = 2, cex = 0.8)
lines(c(-17,7),y= c(0,0), col = adjustcolor("grey60", alpha.f= 0.9), lwd = 1.5, lty = 2)
text(x = 6, y = 0.4,"No effect", col = adjustcolor("grey60", alpha.f= 0.9), cex = 0.8)

# Emphasize data points of D. pulex and galeata 
points(yi ~ log(conc_mg_ml.1), data = df[which(df$sp_name == "Daphnia pulex"),], 
       col = adjustcolor(2, alpha.f = 0.8), pch = 21,
       cex = df$log.invse[which(df$sp_name == "Daphnia pulex")]/4,
       bg = adjustcolor(2, alpha.f = 0.2))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$sp_name == "Daphnia galeata"),], 
       col = adjustcolor(3, alpha.f = 0.8), pch = 21,
       cex = df$log.invse[which(df$sp_name == "Daphnia galeata")]/4,
       bg = adjustcolor(3, alpha.f = 0.2))

# fit model and add line for D. magna
m.Dmag <- rma.mv(yi, vi, W = log.invse, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID,
                 data = df[df$sp_name == "Daphnia magna",])
lines(x = c(-15, 6.5), y = c(m.Dmag$b[1]-m.Dmag$b[2]*15 ,m.Dmag$b[1]+m.Dmag$b[2]*6.5), 
      col = as.numeric(df$sp_name[df$sp_name == "Daphnia magna"]), lwd = 1.5)
# fit model and add line for D. pulex
m.Dpul <- rma.mv(yi, vi, W = log.invse, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID,
                 data = df[df$sp_name == "Daphnia pulex",])
lines(x = c(-14, -1), y = c(m.Dpul$b[1]-m.Dpul$b[2]*14 ,m.Dpul$b[1]-m.Dpul$b[2]*1), 
      col = as.numeric(df$sp_name[df$sp_name == "Daphnia pulex"]), lwd = 1.5)

legend(-2, -0.1, levels(df$sp_name), fill = c(1:3), cex = 0.8, text.font = 3, bty = "n")

dev.off()

## 2. Food --------

# remove samples with NA values
df_food <- df[!is.na(df$food_present),]

# check factor levels
levels(df$food_present)

svglite("Plots/regplot_food.svg", height = 5, width = 5, pointsize = 20)

opar <- par(mar=c(4, 3, 0.5, 2), xpd = TRUE, mfrow = c(1,1), bty = "l", las = 1, cex = 0.8, lwd = 1)
doc_pal = col_vector[c(2,4,19, 5, 20)]
palette(doc_pal)

plot(yi ~ log(conc_mg_ml.1), data = df_food, pch = 21, 
     col = adjustcolor(as.numeric(df_food$food_present), alpha.f = 0.3), 
     bg = adjustcolor(as.numeric(df_food$food_present), alpha.f = 0.1), cex = log.invse/4,
     ylab = "", xlab = "", 
     ylim = c(-2.2,7), las = 1)
mtext("log (Risk Ratio)", side = 2, line = 2, cex = 0.8, las = 3)
mtext("Concentration in mg per ml (log-scale)", side = 1, line = 2, cex = 0.8)
lines(c(-17,7),y= c(0,0), col = adjustcolor("grey60", alpha.f= 0.6), lwd = 1.5, lty = 2)
text(x = 6, y = 0.3,"No effect", col = adjustcolor("grey60", alpha.f= 0.8), cex = 0.8)

# fit model and add line for yes, food added 
m.yes <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
                data = df_food[df_food$food_present == "yes",])
lines(x = c(-16, 3.5), y = c(m.yes$b[1]-m.yes$b[2]*16 ,m.yes$b[1]+m.yes$b[2]*3.5), 
      col = 2, lwd = 1.5)
# fit model and add line for no food added
m.no <- rma.mv(yi, vi, W = log.invse, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID,
               data = df[df_food$food_present == "no",])
lines(x = c(-12, 6.5), y = c(m.no$b[1]-m.no$b[2]*12,m.no$b[1]+m.no$b[2]*6.5), 
      col = 1, lwd = 1.5)

legend(-1,0, c("No food present", "Food present"), fill = c(1:2), cex = 0.8, bty = "n")
dev.off()

### 3. Exposure duration  -----

# check factor levels:
df$exp_time_d <- as.factor(as.character(df$exp_time_d))
levels(df$exp_time_d)
# order factor levels:
df$exp_time_d = factor(df$exp_time_d, levels(df$exp_time_d)[c(1,5,8,9,11,2,3,4,6,7,10)])
# check again
levels(df$exp_time_d)

svglite("Plots/regplot_expdur.svg", height= 5, width = 5.5, pointsize = 20)

exp_pal_total = col_vector[c(11,2,3,1,15,4,19,9,5,6,7)] 

opar <- par(mar=c(4,4,0.5,5), xpd = TRUE, mfrow = c(1,1), bty = "l", las = 1, cex = 0.8)

palette(exp_pal_total)
plot(yi ~ log(conc_mg_ml.1), data = df, pch = 21, 
     col = adjustcolor(as.numeric(df$exp_time_d), alpha.f = 0.2), 
     bg = adjustcolor(as.numeric(df$exp_time_d), alpha.f = 0.05), cex = log.invse/4,
     ylab = "", xlab = "", 
     ylim = c(-2.2,7), 
     las = 1)
mtext("log (Risk Ratio)", side = 2, line = 2, cex = 0.8, las = 3)
mtext("Concentration in mg per ml (log-scale)", side = 1, line = 2, cex = 0.8)
lines(c(-17,7),y= c(0,0), col = adjustcolor("grey60", alpha.f= 0.9), lwd = 1.5, lty = 2)
text(x = 3, y = -0.5,"No effect", col = adjustcolor("grey60", alpha.f= 0.9), cex = 0.8)

# Emphasize points for 1,3,8,17,19,63 days 
points(yi ~ log(conc_mg_ml.1), data = df[which(df$exp_time_d == "1"),], 
       col = adjustcolor(1, alpha.f = 0.5), pch = 21,
       cex = df$log.invse[which(df$exp_time_d == "1")]/4,
       bg = adjustcolor(1, alpha.f = 0.2))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$exp_time_d == "3"),], 
       col = adjustcolor(3, alpha.f = 0.8), pch = 21,
       cex = df$log.invse[which(df$exp_time_d == "1")]/4,
       bg = adjustcolor(3, alpha.f = 0.3))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$exp_time_d == "8"),], 
       col = adjustcolor(5, alpha.f = 0.8), pch = 21,
       cex = df$log.invse[which(df$exp_time_d == "1")]/4,
       bg = adjustcolor(5, alpha.f = 0.3))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$exp_time_d == "17"),], 
       col = adjustcolor(7, alpha.f = 0.8), pch = 21,
       cex = df$log.invse[which(df$exp_time_d == "1")]/4,
       bg = adjustcolor(7, alpha.f = 0.3))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$exp_time_d == "19"),], 
       col = adjustcolor(8, alpha.f = 0.8), pch = 21,
       cex = df$log.invse[which(df$exp_time_d == "1")]/4,
       bg = adjustcolor(8, alpha.f = 0.3))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$exp_time_d == "26"),], 
       col = adjustcolor(10, alpha.f = 0.8), pch = 21,
       cex = df$log.invse[which(df$exp_time_d == "1")]/4,
       bg = adjustcolor(10, alpha.f = 0.3))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$exp_time_d == "63"),], 
       col = adjustcolor(11, alpha.f = 0.8), pch = 21,
       cex = df$log.invse[which(df$exp_time_d == "1")]/4,
       bg = adjustcolor(11, alpha.f = 0.3))

# Add regression lines for exposure times with more than 5 data points:
# table(df$exp_time_d)
# >1   2   3   4   8  14  17  19  21  26  63 
# >13 462   1  65   3   4   1   1 136   3   1 

# fit model and add line for 1 day
m.mg.1 <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
                 data = df[which(df$exp_time_d == "1"),])
lines(x = c(-5, -1.5), y = c(m.mg.1$b[1]-m.mg.1$b[2]*5 ,m.mg.1$b[1]-m.mg.1$b[2]*1.5),
      col = 1, lwd = 1.5)
# fit model and add line for 2 days
m.mg.2 <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
                 data = df[which(df$exp_time_d == "2"),])
lines(x = c(-12, 6.5), y = c(m.mg.2$b[1]-m.mg.2$b[2]*12 ,m.mg.2$b[1]+m.mg.2$b[2]*6.5), 
      col = 2, lwd = 1.5)
# fit model and add line for 4 days
m.mg.4 <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
                 data = df[which(df$exp_time_d == "4"),])
lines(x = c(-15.5, 2.3), y = c(m.mg.4$b[1]-m.mg.4$b[2]*15.5 ,m.mg.4$b[1]+m.mg.4$b[2]*2.3),
      col = 4, lwd = 1.5)
# fit model and add line for 21 days
m.mg.21 <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
                  data = df[which(df$exp_time_d == "21"),])
min(log(df$conc_mg_ml.1[df$exp_time_d=="21"]))
max(log(df$conc_mg_ml.1[df$exp_time_d=="21"]))
lines(x = c(-14, 1), y = c(m.mg.21$b[1]-m.mg.21$b[2]*14 ,m.mg.21$b[1]-m.mg.21$b[2]*1), 
      col = 9, lwd = 1.5)

legend(x = 7, y = 4, 
       c("24hours", "48 hours", "72 hours", "96 hours","8 days", "14 days","17 days","19 days","21 days", "26 days", "63 days"), 
       fill = c(1:11), cex = 0.8, bty = "n")
dev.off()

### 4. Polymer types-----------------------

levels(df$polymer)
col_1 = col_vector[c(20, 3, 2, 12, 8, 18, 4, 16, 19, 1, 5, 7, 15, 6, 9)]
leg = data.frame(
  labels = c(levels(df$polymer)),
  labcol = c(col_1),
  colnum = c(20, 3, 2, 12, 8, 18, 4, 16, 19, 1, 5, 7, 15, 6, 9),
  labcol.names= c("grey","olive","forestgreen","lightgreen","yellow","piggish","midblue","lavender",
                  "pink","darkblue","orange","lightblue","darkred","purple","red")
)

# open svg connection to save vector image
svglite("Plots/regplot_polymers.svg", height = 5, 
        width = 5.5, pointsize = 20)

#set graphics paramters
opar <- par(mar=c(4,4,0.5,6), xpd = TRUE, mfrow = c(1,1), bty = "l", las = 1, cex = 0.8)
# set color palette
palette(col_1)

# plot 
plot(yi ~ log(conc_mg_ml.1), data = df, pch = 21, 
     col = adjustcolor(as.numeric(df$polymer), alpha.f = 0.1), 
     bg = adjustcolor(as.numeric(df$polymer), alpha.f = 0.05), cex = log.invse/4, # change 4 to 2, to increase point size relative to inverse standard error
     ylab = "", xlab = "", ylim = c(-2.7,7))
# add y axis label
mtext("log (Risk Ratio)", side = 2, line = 2, cex = 0.8, las = 3)
# add x axis label
mtext("Concentration in mg per ml (log-scale)", side = 1, line = 2, cex = 0.8)
# add horizontal line at log(RR) = 0, i.e. no effect
lines(c(-17,7),y= c(0,0), col = adjustcolor("grey60", alpha.f= 0.9), lwd = 1.5, lty = 2)
# add label to this line
text(x = 3, y = -0.3,"No effect", col = adjustcolor("grey60", alpha.f= 0.9), cex = 0.8)
# use helper function polyplot_mg_ml to fit mixed meta-regression models for each polymer type and add regression lines
lapply(levels(df$polymer), polyplot_mg_ml, d = df)
# emphasize data points for not so frequent polymers
points(log(df$conc_mg_ml.1[df$polymer == "PET/PS/ABS"]), 
       df$yi[df$polymer == "PET/PS/ABS"], 
       col = adjustcolor(as.numeric(df$polymer[df$polymer == "PET/PS/ABS"]), alpha.f = 0.3), 
       bg = adjustcolor(as.numeric(df$polymer[df$polymer == "PET/PS/ABS"]), alpha.f = 0.2),
       cex = df$log.invse[df$polymer == "PET/PS/ABS"]/4,
       pch = 21)
points(log(df$conc_mg_ml.1[df$polymer == "PVC/PE"]), 
       df$yi[df$polymer == "PVC/PE"], 
       col = adjustcolor(as.numeric(df$polymer[df$polymer == "PVC/PE"]), alpha.f = 0.3), 
       bg = adjustcolor(as.numeric(df$polymer[df$polymer == "PVC/PE"]), alpha.f = 0.2),
       cex = df$log.invse[df$polymer == "PVC/PE"]/4,
       pch = 21)
# add legend
legend(x = 8, y = 6, leg$labels,
       #c("other","LDPE","PE","PLA","PS/Latex","PMMA","PMMA-PSMA","PP","PS","PUR","PVC","PVC/PE","PET/PS/ABS"), 
       leg$labcol, cex = 0.8, bty = "n")
# close svg connection
dev.off()


### 5. Size ---------
svglite("Plots/regplot_size.svg", height = 5, 
        width = 5.5, pointsize = 20)

opar <- par(mar=c(4,4,0.5,2), xpd = TRUE, mfrow = c(1,1), bty = "l", las = 1, cex = 0.8)
palette(col_vector[c(11,5)])
m.size <- rma.mv(yi, vi, mods = ~ log10(size_nm_mean), random = ~ 1|pub_ID/sample_ID, data = df)
m.size
# plot regression plot
regplot(m.size,bty = "l", ylab = "", xlab = "  ", psize = df$log.invse/4,
        col = adjustcolor(as.numeric(as.factor(df$MP.NP)), alpha.f = 0.2), 
        bg = adjustcolor(as.numeric(as.factor(df$MP.NP)), alpha.f = 0.1),
        lcol = "grey60", ylim = c(-2.7, 7))
mtext("log (Risk Ratio)", side = 2, line = 2, cex = 0.8, las = 3)
mtext("log (NMP size in nm)", side = 1, line = 2, cex = 0.8)
legend("bottomright", c("MP", "NP"), fill = c(1:2), cex = 0.8, border = c(1:2), bty = "n")
dev.off()


### 6. Shapes -----
# set color palette
shape_pal1 = col_vector[c(2,4,5,6)] #purple, orange, green, blue  for "Aggregates" "Fibers"     "Fragments"  "Spheres"
# check order of factor levels
levels(df$shape)
# change order to spheres, fragments, fibers, aggregates
df$shape = factor(df$shape,levels(df$shape)[c(4:1)])
# check again
levels(df$shape)


svglite("Plots/regplot_shapes.svg", height = 5, 
        width = 5.5, pointsize = 20)

opar <- par(mar=c(4, 4, 0.5, 2), xpd = TRUE, mfrow = c(1,1), bty = "l", las = 1, cex = 0.8)

palette(shape_pal1)

plot(yi ~ log(conc_mg_ml.1), data = df, pch = 21, 
     col = adjustcolor(as.numeric(df$shape), alpha.f = 0.2), 
     bg = adjustcolor(as.numeric(df$shape), alpha.f = 0.1), cex = log.invse/4,
     ylab = "", xlab = "", 
     ylim = c(-2.7,7), 
     las = 1)
mtext("log (Risk Ratio)", side = 2, line = 2, cex = 0.8, las = 3)
mtext("Concentration in mg per ml (log-scale)", side = 1, line = 2, cex = 0.8)
lines(c(-15,7),y= c(0,0), col = adjustcolor("grey60", alpha.f= 0.9), lwd = 1.5, lty = 2)
text(x = 6, y = 0.4,"No effect", col = adjustcolor("grey60", alpha.f= 0.9), cex = 0.8)

# Emphasize points for fibers
points(yi ~ log(conc_mg_ml.1), data = df[which(df$shape == "Fibers"),], 
       col = adjustcolor(3, alpha.f = 0.8), pch = 21,
       cex = df$log.invse[which(df$shape == "Fibers")]/4,
       bg = adjustcolor(3, alpha.f = 0.2))
# Emphasize points for aggregates 
points(yi ~ log(conc_mg_ml.1), data = df[which(df$shape == "Aggregates"),], 
       col = adjustcolor(4, alpha.f = 0.9), pch = 21,
       cex = df$log.invse[which(df$shape == "Aggregates")]/4,
       bg = adjustcolor(4, alpha.f = 0.4))

# fit model and add line for spheres:
m.s1 <- rma.mv(yi, vi, W = log.invse, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID,
                 data = df[which(df$shape == "Spheres"),])
lines(x = c(-12, 6.5), y = c(m.s1$b[1]-m.s1$b[2]*15 ,m.s1$b[1]+m.s1$b[2]*6.5), 
      col = as.numeric(df$shape[df$shape == "Spheres"]), lwd = 1.5)
# fit model and add line for fragments:
m.f1 <- rma.mv(yi, vi, W = log.invse, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID,
               data = df[which(df$shape == "Fragments"),])
lines(x = c(-14, 3), y = c(m.f1$b[1]-m.f1$b[2]*14 ,m.f1$b[1]+m.f1$b[2]*3), 
      col = as.numeric(df$shape[df$shape == "Fragments"]), lwd = 1.5)
# fit model and a line for fibers:
m.fib <- rma.mv(yi, vi, W = log.invse, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID,
               data = df[which(df$shape == "Fibers"),])
lines(x = c(-4, 5), y = c(m.fib$b[1]-m.fib$b[2]*4 ,m.fib$b[1]+m.fib$b[2]*5), 
      col = as.numeric(df$shape[df$shape == "Fibers"]), lwd = 1.5)

legend(x = 0, y= 0, c("Spheres", "Fragments", "Fibers", "Aggregates"), fill = c(1:4), cex = 0.8, bty = "n")
dev.off()


## 7. Additives, fluorescence tags and other polymer modifications --------

# check  factor levels
levels(df$fluorescence)
levels(df$mod_type)
#check whether there are data for MNP with both fluorescence tags and additives:
df[df$fluorescence == "yes" & df$mod_type %in% c("13C", "DiNP(plasticizer)","BP-3, melt-blended"),]
# none

# create grouping variable:
df$fluor.add <- as.character(df$fluorescence)
df$fluor.add[which(df$mod_type == "13C")] = "13C"
df$fluor.add[which(df$mod_type == "BP-3, melt-blended")] = "BP-3, melt-blended"
df$fluor.add[which(df$mod_type == "DiNP(plasticizer)")] = "DiNP(plasticizer)"
df$fluor.add = as.factor(df$fluor.add)
# check again
levels(df$fluor.add)

#define palette:
levels(df$fluor.add) #13C, BP3, no, yes (=fluorescence)

svglite("Plots/regplot_additives.svg", height = 5, 
        width = 5.5, pointsize = 20)

opar <- par(mar=c(4,4,0.5,5), xpd = TRUE, mfrow = c(1,1), bty = "l", las = 1, cex = 0.8)

fluor_tot = col_vector[c(2,19, 20,4)] 
palette(fluor_tot)

plot(yi ~ log(conc_mg_ml.1), data = df, pch = 21, 
     col = adjustcolor(as.numeric(df$fluor.add), alpha.f = 0.1), 
     bg = adjustcolor(as.numeric(df$fluor.add), alpha.f = 0.05), cex = log.invse/4,
     ylab = "", xlab = "", 
     ylim = c(-2.7,7), 
     las = 1)
mtext("log (Risk Ratio)", side = 2, line = 2, cex = 0.8, las = 3)
mtext("Concentration in mg per ml (log-scale)", side = 1, line = 2, cex = 0.8)
lines(c(-17,7),y= c(0,0), col = adjustcolor("grey60", alpha.f= 0.9), lwd = 1.5, lty = 2)
text(x = 6, y = 0.4,"No effect", col = adjustcolor("grey60", alpha.f= 0.9), cex = 0.8)

# Emphasize points for 13C and BP3 
points(yi ~ log(conc_mg_ml.1), data = df[which(df$fluor.add == "13C"),], 
       col = adjustcolor(1, alpha.f = 0.8), pch = 21,
       cex = df$log.invse[which(df$fluor.add == "13C")]/4,
       bg = adjustcolor(1, alpha.f = 0.9))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$fluor.add == "BP-3, melt-blended"),], 
       col = adjustcolor(2, alpha.f = 0.6), pch = 21,
       cex = df$log.invse[which(df$fluor.add == "BP-3, melt-blended")]/4,
       bg = adjustcolor(2, alpha.f = 0.5))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$fluor.add == "yes"),], 
       col = adjustcolor(4, alpha.f = 0.4), pch = 21,
       cex = df$log.invse[which(df$fluor.add == "yes")]/4,
       bg = adjustcolor(4, alpha.f = 0.2))


# fit model and add line for no modifications
m.no <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
               data = df[which(df$fluor.add == "no"),])
lines(x = c(-15, 6), y = c(m.no$b[1]-m.no$b[2]*15 ,m.no$b[1]+m.no$b[2]*6), 
      col = 3, lwd = 1.5)
# fit model and add line for BP-3, melt-blended
m.BP3 <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
                data = df[which(df$fluor.add == "BP-3, melt-blended"),])
lines(x = c(-10, -3.5), y = c(m.BP3$b[1]-m.BP3$b[2]*10 ,m.BP3$b[1]-m.BP3$b[2]*3.5), 
      col = 2, lwd = 1.5)
# fit model and add line for fluorescence tags
m.yes <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
                data = df[which(df$fluor.add == "yes"),])
lines(x = c(-12, 4), y = c(m.yes$b[1]-m.yes$b[2]*12 ,m.yes$b[1]+m.yes$b[2]*4), 
      col = 4, lwd = 1.5)

legend(x = -1, y = 0, 
       c("BP-3 melt-blended", "C13 labelled", "Fluorescence tags", "No tags"), 
       fill = c(2,1,4,3), cex = 0.8, bty = "n")
dev.off()


### 8. Uncontrolled surface modifications / Dissolved Organic Carbon (DOC)  --------

# prepare new variable for plotting categories:
# 1. fulvic acid, 2. humic acid, 3. lake water, 4. no, 5. incubated
# new grouping variable
df$uncontrolled <- as.character(df$DOC_type)
# change all values with DOC_added no to "no"
df$uncontrolled[df$DOC_added == "no"] <- "no"
# change all values for data from experiments with incubated MNP to "incubated"
df$uncontrolled[df$biotic_env != "no"] <- "incubated"
df$uncontrolled <- as.factor(df$uncontrolled)
#check levels
levels(df$uncontrolled)

svglite("Plots/regplot_DOC.svg", height = 5, 
        width = 5.5, pointsize = 20)

opar <- par(mar=c(4, 4, 0.5, 5), xpd = TRUE, mfrow = c(1,1), bty = "l", las = 1, cex = 0.8)

palette(doc_pal)

plot(yi ~ log(conc_mg_ml.1), data = df, pch = 21, 
     col = adjustcolor(as.numeric(df$uncontrolled), alpha.f = 0.1), 
     bg = adjustcolor(as.numeric(df$uncontrolled), alpha.f = 0.1), cex = log.invse/4,
     ylab = "", xlab = "", 
     ylim = c(-2.7,7), 
     las = 1)
mtext("log (Risk Ratio)", side = 2, line = 2, cex = 0.8, las = 3)
mtext("Concentration in mg per ml (log-scale)", side = 1, line = 2, cex = 0.8)
lines(c(-17,7),y= c(0,0), col = adjustcolor("grey60", alpha.f= 0.9), lwd = 1.5, lty = 2)
text(x = 8, y = 0.4,"No effect", col = adjustcolor("grey60", alpha.f= 0.9), cex = 0.8)

# Emphasize points
points(yi ~ log(conc_mg_ml.1), data = df[which(df$uncontrolled == "fulvic acid"),], 
       col = adjustcolor(1, alpha.f = 0.4), pch = 21,
       cex = df$log.invse[which(df$uncontrolled == "fulvic acid")]/4,
       bg = adjustcolor(1, alpha.f = 0.3))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$uncontrolled == "humic acid"),], 
       col = adjustcolor(2, alpha.f = 0.5), pch = 21,
       cex = df$log.invse[which(df$uncontrolled == "humic acid")]/4,
       bg = adjustcolor(2, alpha.f = 0.3))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$uncontrolled == "incubated"),], 
       col = adjustcolor(3, alpha.f = 0.7), pch = 21,
       cex = df$log.invse[which(df$uncontrolled == "incubated")]/4,
       bg = adjustcolor(3, alpha.f = 0.5))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$uncontrolled == "lake water"),], 
       col = adjustcolor(4, alpha.f = 0.4), pch = 21,
       cex = df$log.invse[which(df$uncontrolled == "lake water")]/4,
       bg = adjustcolor(4, alpha.f = 0.2))

# fit model and add line for humic acid
m.HA <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
               data = df[which(df$uncontrolled == "humic acid"),])
lines(x = c(-8, 4), y = c(m.HA$b[1]-m.HA$b[2]*8 ,m.HA$b[1]+m.HA$b[2]*4), 
      col = 2, lwd = 1.5)
# fit model and add line for incubated NMP
m.inc <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
                data = df[which(df$uncontrolled == "incubated"),])
lines(x = c(-3, 3), y = c(m.inc$b[1]-m.inc$b[2]*3 ,m.inc$b[1]+m.inc$b[2]*3), 
      col = 3, lwd = 1.5)
# fit model and add line for exposure in lake water
m.lake <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
                 data = df[which(df$uncontrolled == "lake water"),])
lines(x = c(-8, -3), y = c(m.lake$b[1]-m.lake$b[2]*8 ,m.lake$b[1]-m.lake$b[2]*3), 
      col = 4, lwd = 1.5)
# fit model and add line for no DOC added
m.no <- rma.mv(yi, vi, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID, W = log.invse,
               data = df[which(df$uncontrolled == "no"),])
lines(x = c(-12, 6), y = c(m.no$b[1]-m.no$b[2]*12 ,m.no$b[1]+m.no$b[2]*6), 
      col = 5, lwd = 1.5)

legend(x = 4, y = 0.2, c("+ Fulvic acid","+ Humic acid", "Lake water", "NMP incubated", "No DOC"), fill = c(1,2,4,3,5), cex = 0.8, bty = "n")

dev.off()


### 9. Controlled Surface Modifications ----------
# check factor levels
levels(df$mod_type)

# change to class character to adjust factor levels
df$controlled <- as.character(df$mod_type)
# change NMP samples that have not been modified to modification type "no"
df$controlled[which(df$modified == "no")] <- "no"
# change modification type "13C", "DiNP(plasticizer)", "BP-3, melt-blended" to "no", 
# because these are polymer modifications and not controlled surface modifications
df$controlled[which(df$mod_type == "13C")] <- "no"
df$controlled[which(df$mod_type == "DiNP(plasticizer)")] <- "no"
df$controlled[which(df$mod_type == "BP-3, melt-blended")] <- "no"
# change class back to factor
df$controlled <- as.factor(df$controlled)
# check again
levels(df$controlled)

svglite("Plots/regplot_CSM.svg", height = 5, 
        width = 5.5, pointsize = 20)

opar <- par(mar=c(4,4, 0.5, 5), xpd = TRUE, mfrow = c(1,1), bty = "l", las = 1, cex = 0.8)

cont_pal1 = col_vector[c(2,19, 5,4, 20)] 
palette(cont_pal1)

plot(yi ~ log(conc_mg_ml.1), data = df, pch = 21, 
     col = adjustcolor(as.numeric(df$controlled), alpha.f = 0.2), 
     bg = adjustcolor(as.numeric(df$controlled), alpha.f = 0.1), cex = log.invse/4,
     ylim = c(-2.7,7),
     ylab = "", xlab = "")
mtext("log (Risk Ratio)", side = 2, line = 2, cex = 0.8, las = 3)
mtext("Concentration in mg per ml (log-scale)", side = 1, line = 2, cex = 0.8)
lines(c(-17,6.5),y= c(0,0), col = adjustcolor("grey60", alpha.f= 0.9), lwd = 1.5, lty = 2)
text(x = 5, y = 0.3,"No effect", col = adjustcolor("grey60", alpha.f= 0.9), cex = 0.8)

# Emphasize points for carboxylated and extracted NMP
points(yi ~ log(conc_mg_ml.1), data = df[which(df$controlled == "Extracted"),], 
       col = adjustcolor(4, alpha.f = 0.2), pch = 21,
       cex = df$log.invse[which(df$uncontrolled == "Extracted")]/4,
       bg = adjustcolor(4, alpha.f = 0.1))
points(yi ~ log(conc_mg_ml.1), data = df[which(df$controlled == "Carboxylated"),], 
       col = adjustcolor(3, alpha.f = 0.2), pch = 21,
       cex = df$log.invse[which(df$uncontrolled == "Carboxylated")]/4,
       bg = adjustcolor(3, alpha.f = 0.1))

# fit model and add line for "Carboxylated" :
m.carboxy <- rma.mv(yi, vi,  W = log.invse, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID,
                    data = df[which(df$controlled == "Carboxylated"),])
lines(x = c(-9, 0), y = c(m.carboxy$b[1]-m.carboxy$b[2]*9 ,m.carboxy$b[1]+m.carboxy$b[2]*0), 
      col = as.numeric(df$controlled[which(df$controlled == "Carboxylated")]), lwd = 1.5)

# fit model and add line for "Aminated" :
m.amini <- rma.mv(yi, vi, W = log.invse,  mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID,
                  data = df[which(df$controlled == "Aminated"),])
lines(x = c(-13, 4), y = c(m.amini$b[1]-m.amini$b[2]*13 ,m.amini$b[1]+m.amini$b[2]*4), 
      col = as.numeric(df$controlled[which(df$controlled == "Aminated")]), lwd = 1.5)

# fit model and add line for "Amidine derivatized" :
m.amid <- rma.mv(yi, vi, W = log.invse,  mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID,
                 data = df[which(df$controlled == "Amidine derivatized"),])
lines(x = c(-8, -3.5), y = c(m.amid$b[1]-m.amid$b[2]*8 ,m.amid$b[1]-m.amid$b[2]*3.5), 
      col = as.numeric(df$controlled[which(df$controlled == "Amidine derivatized")]), lwd = 1.5)

# fit model and add line for "Extracted" :
m.ext <- rma.mv(yi, vi,  W = log.invse, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID,
                data = df[which(df$controlled == "Extracted"),])
lines(x = c(-7, -1), y = c(m.ext$b[1]-m.ext$b[2]*7 ,m.ext$b[1]-m.ext$b[2]*1), 
      col = as.numeric(df$controlled[which(df$controlled == "Extracted")]), lwd = 1.5)

# fit model and add line for "no" :
m.notmod <- rma.mv(yi, vi,  W = log.invse, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID/sample_ID,
                   data = df[which(df$controlled == "no"),])
lines(x = c(-14, 7), y = c(m.notmod$b[1]-m.notmod$b[2]*14 ,m.notmod$b[1]+m.notmod$b[2]*7), 
      col = as.numeric(df$controlled[which(df$controlled == "no")]), lwd = 1.5)

legend(0, 0.2, c("Amidine derivatized", "Aminated","Carboxylated","Extracted","Not modified"), fill = c(1:5), cex = 0.8, text.font = 3, bty = "n")
dev.off()


### 10. Age ---------

# fit meta-regression model
m.age <- rma.mv(yi, vi, mods = ~ as.numeric(age_d), random = ~ 1|pub_ID/sample_ID, data = df)

svglite("Plots/regplot_age.svg", height = 5, 
        width = 10, pointsize = 20)

opar <- par(mar=c(4,5,1,1), xpd = TRUE, mfrow = c(1,1), bty = "l", las = 1)
palette(col_vector)
regplot(m.age,bty = "l", ylab = "log(Risk Ratio)", xlab = "Age in days", psize = df$log.invse/4,
        lcol = "grey60", col = adjustcolor(1, alpha.f = 0.3), bg = adjustcolor(1, alpha.f = 0.1))

dev.off()


### 11. Temperature --------

# fit model
m.temperature <- rma.mv(yi, vi, mods = ~ temperature, random = ~ 1|pub_ID, data = df)

svglite("Plots/regplot_temperature.svg", height = 5, 
        width = 10, pointsize = 20)

opar <- par(mar=c(4,5,1,7), xpd = TRUE, mfrow = c(1,1), bty = "l", las = 1)
palette(col_vector)
regplot(m.temperature,bty = "l", ylab = "log(Risk Ratio)", xlab = "Temperature in °C", psize = df$log.invse/4,
         lcol = "grey60", col = adjustcolor(15, alpha.f = 0.3), bg = adjustcolor(15, alpha.f = 0.1))

dev.off()

### END ----------