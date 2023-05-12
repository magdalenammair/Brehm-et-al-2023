polyplot_mg_ml <- function(d, p) {
  d = as.data.frame(d)
  model = rma.mv(yi, vi, W = log.invse, mods = ~ log(conc_mg_ml.1), random = ~ 1|pub_ID, 
                 data = d[d$polymer == p,])
  minx = min(log(d$conc_mg_ml.1[which(d$polymer == p)]), na.rm = T)
  maxx = max(log(d$conc_mg_ml.1[which(d$polymer == p)]), na.rm = T)
  lines(x = c(minx,maxx), y = c(model$b[1]+model$b[2]*minx, model$b[1]+model$b[2]*maxx), 
        col = as.numeric(d$polymer[d$polymer == p]), lwd = 2)
}

