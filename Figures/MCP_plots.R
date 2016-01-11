library(ggplot2)

setwd("~/Dropbox/Economics/Courses/Thesis/Positioning/Figures")



### Mean eccentricity ###
data = read.csv("data/MCP_hunter_mean_eccentricity_20160110_210442_i250_b150.csv")

data[, "polarization"] = NA
data[data$mu <= 0.5, ][, "polarization"] = 1
data[0.5 < data$mu & data$mu < 1, ][, "polarization"] = 2
data[1 <= data$mu, ][, "polarization"] = 3

c1 = ggplot(data, aes(y=MeanEst, x=N, colour=factor(polarization)))
c1 + stat_smooth() + 
  geom_point() + 
  xlim(2,12) + 
  ylim(0, 1.65)

# Subset
datasubset = subset(data, N %in% c(2,3,4,12))

c2 = ggplot(datasubset, aes(y=MeanEst, x=n_ratio, colour=factor(N)))
c2 + stat_smooth(aes(fill = factor(N)))# + 
     #geom_point()



### ENP ###
data = read.csv("data/MCP_hunter_ENP_20160110_210442_i250_b150.csv")

data[, "polarization"] = NA
data[data$mu <= 0.5, ][, "polarization"] = 1
data[0.5 < data$mu & data$mu < 1, ][, "polarization"] = 2
data[1 <= data$mu, ][, "polarization"] = 3

c3 = ggplot(data, aes(y=MeanEst, x=N, colour=factor(polarization)))
c3 + stat_smooth() + 
  geom_abline(intercept = 0, slope = 1) +
  xlim(2,12) + 
  ylim(1, 12)# +
  #geom_point()