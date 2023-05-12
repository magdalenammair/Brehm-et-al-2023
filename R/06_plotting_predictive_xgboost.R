## load packages
library(svglite)
library(xgboost)


## import data ----
res = readRDS("Results/importances.RDS")
pred_importance = res$grouped_importances
importance = res$importance

## Gain plot--------
svglite("Plots/group_gain_best_model.svg", height = 6, 
        width = 9, pointsize = 14)


par(las = 2, bty = "l", mar = c(9,5,1,1), xpd = TRUE, lwd = 3)

palette(c( "#FFD966", "#8FAADC", "#A9D18E", "#C9C9C9")) #gelb, blau, grün, grau
barplot(pred_importance$gain , col = pred_importance$boxgroup, border = pred_importance$linegroup, 
        ylim = c(0, 0.25),
        #names.arg = pred_importance$predictor,
        ylab = "Variable importance (gain)"
)
text(pred_importance$label, 
     x = seq(1,25.2,by = 1.2), y = rep(-0.01,21), 
     srt = 45, pos = 2, cex = 0.8)

dev.off()


# SUPPLEMENT - other importance measures and single feature importances - barplots
## Cover plot--------
pred_importance2 = pred_importance[order(pred_importance$cover, decreasing = T),]

svglite("Plots/FigureS8_Cover.svg", height = 6, 
        width = 9, pointsize = 14)


par(las = 2, bty = "l", mar = c(9,5,1,1), xpd = TRUE, lwd = 3)
barplot(pred_importance2$cover , col = pred_importance2$boxgroup, border = pred_importance2$linegroup, 
        ylim = c(0, 0.2),
        #names.arg = pred_importance$predictor,
        ylab = "Variable importance (cover)"
)
text(pred_importance2$label, 
     x = seq(1,25.2,by = 1.2), y = rep(-0.01,21), 
     srt = 45, pos = 2, cex = 0.8)
dev.off()


## Frequency plot--------
pred_importance3 = pred_importance[order(pred_importance$frequency, decreasing = T),]

svglite("Plots/FigureS7_Frequency.svg", height = 6, 
        width = 9, pointsize = 14)

par(las = 2, bty = "l", mar = c(9,5,1,1), xpd = TRUE, lwd = 3)

barplot(pred_importance3$frequency , col = pred_importance3$boxgroup, border = pred_importance3$linegroup, 
        ylim = c(0, 0.25),
        #names.arg = pred_importance$predictor,
        ylab = "Variable importance (frequency)"
)
text(pred_importance3$label, 
     x = seq(1,25.2,by = 1.2), y = rep(-0.01,21), 
     srt = 45, pos = 2, cex = 0.8)

dev.off()

## Single feature importances-----
# gain
svglite("Plots/FigureS6_single_gain.svg", height = 9, 
        width = 9, pointsize = 14)
xgb.plot.importance(importance_matrix = importance)
dev.off()
# cover
xgb.plot.importance(importance_matrix = importance, measure = "Cover")
#frequency
xgb.plot.importance(importance_matrix = importance, measure = "Frequency")

## Plot performances ----

all_performance = readRDS("Results/performance_xgboost.RDS")

# and combine with gain plot:
png("Plots/performance_importances.png", width = 20, height = 12, res = 1000, units = "cm")

lo = layout(mat = matrix(c(1,2,2,2,2),nrow = 1, ncol = 5, byrow = TRUE))
par(bty="l", xaxt = "n", las = 1, mar = c(10,5,1,1), lwd = 1, cex = 0.8, xpd =T)
boxplot(all_performance$rmse ~ all_performance$which.model, 
        xlab = "", ylab = "Predictive error (RMSE)", pch = 20,
        col = c("grey95","grey75"))
mtext(c("Baseline model","True model"), side = 1, line = 1, at = c(1, 2), cex = 0.7, las = 2)
text("A", x = 0.8, y = 2.125)

par(las = 2, lwd = 2.5)
palette(c( "#FFD966", "#8FAADC", "#A9D18E", "#C9C9C9")) #gelb, blau, grün, grau
barplot(pred_importance$gain , col = pred_importance$boxgroup, border = pred_importance$linegroup, 
        ylim = c(0, 0.25),
        #names.arg = pred_importance$predictor,
        ylab = "Variable importance (gain)"
)
text(pred_importance$label, 
     x = seq(1,25.2,by = 1.2), y = rep(-0.01,21), 
     srt = 45, pos = 2, cex = 0.8)
text("B", x = 0.5, y = 0.25)
dev.off()

## END --------
