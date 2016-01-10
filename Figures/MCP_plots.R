library(ggplot2)

setwd("~/Dropbox/Economics/Courses/Thesis/Positioning/Figures")
#data = read.csv("data/MCP_hunter_mean_eccentricity_20160110_184957_i250_b150.csv")
data = read.csv("data/MCP_hunter_ENP_20160110_184957_i250_b150.csv")

# Groups the polarization variable mu into three groups.
data[, "polarization"] = NA
data[data$mu <= 0.5, ][, "polarization"] = 1
data[0.5 < data$mu & data$mu < 1, ][, "polarization"] = 2
data[1 <= data$mu, ][, "polarization"] = 3


# Mean eccentricity
# Geoms and stats are automatically split by aesthetics that are factors
c = ggplot(data, aes(y=MeanEst, x=N, colour=factor(polarization)))
c + stat_smooth() + 
  geom_point() + 
  xlim(2,12) + 
  ylim(0, 1.65)

# Geoms and stats are automatically split by aesthetics that are factors
c2 = ggplot(data, aes(y=MeanEst, x=n_ratio, colour=factor(N)))
c2 + stat_smooth() + geom_point()


# ENP
# Geoms and stats are automatically split by aesthetics that are factors
c = ggplot(data, aes(y=MeanEst, x=N, colour=factor(polarization)))
c + stat_smooth() + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  xlim(2,12) + 
  ylim(1, 12) 