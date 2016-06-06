# FIGURES
# Sets up figures and formats them for final paper.
# Jonas K. Sekamane


# TO-DO -------------------------------------------------------------------

# stat_smooth ==> Median spline


# 1. SETUP ----------------------------------------------------------------

# 1.1 Libraries
library('ggplot2')
library('plyr') 
library('reshape2') 
library('quantreg')


# 1.2 Working directory
setwd("~/Dropbox/Economics/Courses/Thesis/Positioning/Figures")

# 1.3 Functions

# Return the confidence interval bands for median splines / rqss (Additive Quantile Regression Smoothing).
rqss_ci_group_bands = function(df, y, x, group) {
  pdf(file = NULL) # suppress plot output
  band = plot(rqss(df[, y] ~ qss(df[, x])), bands=TRUE, rug=FALSE, jit=FALSE)
  dev.off() # end suppress plot output
  output = data.frame( band[[1]]$x, band[[1]]$blo, band[[1]]$bhi, unique(df[, group]) )
  names(output) = c(x, paste(y,"BLo",sep=""), paste(y,"BHi",sep=""), group) # set column names
  return(output)
}
rqss_ci_bands = function(dd, y, x, group) {
  group_bands = lapply(split(dd, as.factor(dd[, group])), rqss_ci_group_bands, y=y, x=x, group=group)
  bands = Reduce("rbind", group_bands)
}

# The n_ratio breaking point seperating unimodal and bimodal distributions.
n_ratio_break = function(mu) {
  #if( mu < 0.5 ) {
  #  return(NaN)
  #} else {
    breakpoint = -(sqrt(4*mu^2-1)-2*mu) / (sqrt(4*mu^2-1)+2*mu) * exp(4*sqrt(4*mu^2-1)*mu)
    return(breakpoint)
  #}
}

# Equlibrium condition 1: ENP at the limit of 'No firm has market twice as large as any other firm's market'.
eqcond1_ENP = function(N) {(1+N)^2/(3+N)}

# 2. GRID SWEEP: SAME DECISION RULE MODELS --------------------------------

gs_s1 = read.csv("data/GS_sticker_mean_eccentricity_20160203_211455_r1000.csv")
gs_h1 = read.csv("data/GS_hunter_mean_eccentricity_20160202_194230_i1150_b150.csv")
gs_a1 = read.csv("data/GS_aggregator_mean_eccentricity_20160202_194347_r1000_b50.csv")
gs_m1 = read.csv("data/GS_maxcov_mean_eccentricity_20160203_231850_r1000_b99.csv")

gs_s2 = read.csv("data/GS_sticker_ENP_20160203_211455_r1000.csv")
gs_h2 = read.csv("data/GS_hunter_ENP_20160202_194230_i1150_b150.csv")
gs_a2 = read.csv("data/GS_aggregator_ENP_20160202_194347_r1000_b50.csv")
gs_m2 = read.csv("data/GS_maxcov_ENP_20160203_231850_r1000_b99.csv")

gs_s3 = read.csv("data/GS_sticker_mean_representation_20160203_211455_r1000.csv")
gs_h3 = read.csv("data/GS_hunter_mean_representation_20160202_194230_i1150_b150.csv")
gs_a3 = read.csv("data/GS_aggregator_mean_representation_20160202_194347_r1000_b50.csv")
gs_m3 = read.csv("data/GS_maxcov_mean_representation_20160203_231850_r1000_b99.csv")

# Merge the grid-sweep datasets of every decision rule into one.
gs_1 = as.data.frame(mapply(c, 
                     data.frame(gs_s1[,1:3], 1), 
                     data.frame(gs_h1[,1:3], 2), 
                     data.frame(gs_a1[,1:3], 3),
                     data.frame(gs_m1[,1:3], 4) ))
colnames(gs_1)[4] = "Rule"
gs_1$Rule = factor(gs_1$Rule, labels = c("All-sticker", "All-hunter", "All-aggregator", "All-maxcov"))

gs_2 = as.data.frame(mapply(c, 
                            data.frame(gs_s2[,1:3], 1), 
                            data.frame(gs_h2[,1:3], 2), 
                            data.frame(gs_a2[,1:3], 3),
                            data.frame(gs_m2[,1:3], 4)  ))
colnames(gs_2)[4] = "Rule"
gs_2$Rule = factor(gs_2$Rule, labels = c("All-sticker", "All-hunter", "All-aggregator", "All-maxcov"))

gs_3 = as.data.frame(mapply(c, 
                            data.frame(gs_s3[,1:3], 1), 
                            data.frame(gs_h3[,1:3], 2), 
                            data.frame(gs_a3[,1:3], 3),
                            data.frame(gs_m3[,1:3], 4)  ))
colnames(gs_3)[4] = "Rule"
gs_3$Rule = factor(gs_3$Rule, labels = c("All-sticker", "All-hunter", "All-aggregator", "All-maxcov"))

  # 2.1 Mean eccentricity
  fig21a = ggplot(gs_1, aes(y = MeanEst, x = N, colour = factor(Rule))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig21a + geom_line(size=1) + 
    geom_ribbon(aes(x = N, ymin = MeanEst-StdDev, ymax = pmin(gs_1$MeanEst+gs_1$StdDev, rep(1.7, nrow(gs_1))), group = factor(Rule), fill=Rule, color = NULL), alpha = 0.1) +
    scale_fill_discrete(guide=FALSE) +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
    labs(y = "Mean eccentricity", x = "Number of firms", colour = NULL)
  ggsave("fig21a.pdf", width = 21, height = 16, units = "cm")
  
  # 2.2 ENP
  fig22a = ggplot(gs_2, aes(y=MeanEst, x=N, colour=factor(Rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig22a + geom_line(size=1) + 
    geom_ribbon(aes(x = N, ymin = pmax(gs_2$MeanEst-gs_2$StdDev, rep(1, nrow(gs_2))), ymax = pmin(gs_2$MeanEst+gs_2$StdDev, rep(12, nrow(gs_2))), group = factor(Rule), fill=Rule, color = NULL), alpha = 0.1) +
    scale_fill_discrete(guide=FALSE) +
    geom_abline(intercept = 0, slope = 1, color="gray") +
    stat_function(fun = eqcond1_ENP, colour = "black") + 
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(1, 12), breaks = 2:12, minor_breaks = NULL, expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "Number of firms", colour = NULL)
  ggsave("fig22a.pdf", width = 21, height = 16, units = "cm")
  
  # 2.3 Mean representation
  fig23a_bands = rqss_ci_bands(gs_3, y="MeanEst", x="N", group="Rule")
  fig23a = ggplot(gs_3, aes(y = MeanEst, x = N, colour = factor(Rule))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig23a + geom_line(size=1) + 
    geom_ribbon(aes(x = N, ymin = MeanEst-StdDev, ymax = MeanEst+StdDev, group = factor(Rule), fill = Rule, color = NULL), alpha = 0.1) +
    scale_fill_discrete(guide=FALSE) +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(-4, 0), expand = c(0, 0)) +
    labs(y = "Mean representation", x = "Number of firms", colour = NULL)  
  ggsave("fig23a.pdf", width = 21, height = 16, units = "cm")
  

# 3. MCP: SAME DECISION RULE WITH EXOGENOUS NUMBER OF FIRMS ---------------

# 3.1 ALL-STICKER MODEL
# ...



# 3.2 ALL-HUNTER MODEL
mcp_ex_h1 = read.csv("data/MCP_hunter2_mean_eccentricity_20160420_104505_i151_b150_r100.csv")
mcp_ex_h2 = read.csv("data/MCP_hunter2_ENP_20160420_104505_i151_b150_r100.csv")
mcp_ex_h3 = read.csv("data/MCP_hunter2_mean_representation_20160420_104505_i151_b150_r100.csv")

# Grouping polarization into three categories (less than 0.5. between 0.5-1. above 1).
mcp_ex_h1[, "polarization"] = NA
mcp_ex_h1[mcp_ex_h1$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_h1[0.5 < mcp_ex_h1$mu & mcp_ex_h1$mu < 1, ][, "polarization"] = 2
mcp_ex_h1[1 <= mcp_ex_h1$mu, ][, "polarization"] = 3
mcp_ex_h1$polarization = ordered(mcp_ex_h1$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

mcp_ex_h2[, "polarization"] = NA
mcp_ex_h2[mcp_ex_h2$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_h2[0.5 < mcp_ex_h2$mu & mcp_ex_h2$mu < 1, ][, "polarization"] = 2
mcp_ex_h2[1 <= mcp_ex_h2$mu, ][, "polarization"] = 3
mcp_ex_h2$polarization = ordered(mcp_ex_h2$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

mcp_ex_h3[, "polarization"] = NA
mcp_ex_h3[mcp_ex_h3$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_h3[0.5 < mcp_ex_h3$mu & mcp_ex_h3$mu < 1, ][, "polarization"] = 2
mcp_ex_h3[1 <= mcp_ex_h3$mu, ][, "polarization"] = 3
mcp_ex_h3$polarization = ordered(mcp_ex_h3$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

# Grouping polarization into unimodal and bimodal
mcp_ex_h1$n_ratio_break = n_ratio_break(mcp_ex_h1$mu);
mcp_ex_h1[, "unimodal"] = 0
mcp_ex_h1[is.na(mcp_ex_h1$n_ratio_break) | mcp_ex_h1$n_ratio >= mcp_ex_h1$n_ratio_break, ][, "unimodal"] = 1
mcp_ex_h1$unimodal = ordered(mcp_ex_h1$unimodal, labels = c("Bimodal", "Unimodal"))

mcp_ex_h3$n_ratio_break = n_ratio_break(mcp_ex_h3$mu);
mcp_ex_h3[, "unimodal"] = 0
mcp_ex_h3[is.na(mcp_ex_h3$n_ratio_break) | mcp_ex_h3$n_ratio >= mcp_ex_h3$n_ratio_break, ][, "unimodal"] = 1
mcp_ex_h3$unimodal = ordered(mcp_ex_h3$unimodal, labels = c("Bimodal", "Unimodal"))

# Subset of data only consisting 2, 3, 4 and 12 firms.
mcp_ex_h1_subset = subset(mcp_ex_h1, N %in% c(2,3,4,12))
cfirm4 = c("#253494", "#2c7fb8", "#41b6c4", "#a1dab4") # YlGnBu 5 reversed (but only using first 4). http://colorbrewer2.org


  # 3.2.1 Mean eccentricity
  # All-hunter mean eccentricity as function of number of firms in market.
  fig321a = ggplot(mcp_ex_h1, aes(y = MeanEst, x = N, colour = factor(polarization))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig321a + stat_smooth() + 
    geom_point(size = 1) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
    labs(y = "Mean eccentricity", x = "Number of firms", colour = NULL)
  ggsave("fig321a.pdf", width = 21, height = 16, units = "cm")
  
  fig321b = ggplot(mcp_ex_h1, aes(y = MeanEst, x = N, colour = factor(unimodal))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig321b + stat_smooth() + 
    geom_point(size = 1) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
    labs(y = "Mean eccentricity", x = "Number of firms", colour = NULL)
  
  # All-hunter mean eccentricity as function the relative subpopulation size for a subset number of firms.
  fig321b = ggplot(mcp_ex_h1_subset, aes(y = MeanEst, x = n_ratio, colour = factor(N))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig321b + stat_smooth(aes(fill = factor(N))) + 
    scale_color_manual(values=cfirm4) +
    scale_fill_manual(values = cfirm4, guide="none") +
    scale_x_continuous(limits = c(1, 2), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 1.3), expand = c(0, 0)) +
    labs(y = "Mean eccentricity", x = "Relative size of left subpopulation", colour = "Number of firms")
  ggsave("fig321b.pdf", width = 21, height = 16, units = "cm")
  
  # 3.2.2 ENP
  # All-hunter effective number of firms compared to actual number of firms in market.
  fig322a = ggplot(mcp_ex_h2, aes(y=MeanEst, x=N, colour=factor(polarization))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig322a + stat_smooth() + 
    scale_colour_brewer(palette = "Dark2") +
    geom_abline(intercept = 0, slope = 1, color="gray") +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(1, 12), breaks = 2:12, minor_breaks = NULL, expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "Number of firms", colour = NULL)
  ggsave("fig322a.pdf", width = 21, height = 16, units = "cm")

  # 3.2.3 Mean representation
  # All-hunter mean representation as function of number of firms in market.
  fig323a_bands = rqss_ci_bands(mcp_ex_h3, y="MeanEst", x="N", group="polarization")
  fig323a = ggplot(mcp_ex_h3, aes(y = MeanEst, x = N, colour = factor(polarization))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig323a + stat_quantile(quantiles=0.5, formula = y ~ qss(x), method = "rqss", size = 1) +
    geom_ribbon(data = fig323a_bands, 
                aes(x = N, ymin = MeanEstBLo, ymax = MeanEstBHi, group = factor(polarization), y = NULL, color = NULL), 
                alpha = 0.15) +
    geom_point(size = 1) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(-4, 0), expand = c(0, 0)) +
    labs(y = "Mean representation", x = "Number of firms", colour = NULL)  
  ggsave("fig323a.pdf", width = 21, height = 16, units = "cm")
  
  fig323b_bands = rqss_ci_bands(mcp_ex_h3, y="MeanEst", x="N", group="unimodal")
  fig323b = ggplot(mcp_ex_h3, aes(y = MeanEst, x = N, colour = factor(unimodal))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig323b + stat_quantile(quantiles=0.5, formula = y ~ qss(x), method = "rqss", size = 1) +
    geom_ribbon(data = fig323b_bands, 
                aes(x = N, ymin = MeanEstBLo, ymax = MeanEstBHi, group = factor(unimodal), y = NULL, color = NULL), 
                alpha = 0.15) +
    geom_point(size = 1) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(-4, 0), expand = c(0, 0)) +
    labs(y = "Mean representation", x = "Number of firms", colour = NULL)  

#   fig323a = ggplot(mcp_ex_h3, aes(y = MeanEst, x = N, colour = factor(polarization))) + 
#     theme_minimal() +
#     theme(legend.position = "top", 
#           legend.box = "horizontal", 
#           legend.key = element_rect(fill = NA, colour = NA) )
#   fig323a + stat_smooth( method = "gam", formula = y ~ s(x) ) +
#     geom_point(size = 1) + 
#     scale_colour_brewer(palette = "Dark2") +
#     scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
#     scale_y_continuous(limits = c(-4, 0), expand = c(0, 0)) +
#     labs(y = "Mean representation", x = "Number of firms", colour = NULL)  
  
    
  
# 3.3 ALL-AGGREGATOR MODEL
mcp_ex_a1 = read.csv("data/MCP_aggregator_mean_eccentricity_20160202_195640_i101_b100_r100.csv")
mcp_ex_a2 = read.csv("data/MCP_aggregator_ENP_20160202_195640_i101_b100_r100.csv")
mcp_ex_a3 = read.csv("data/MCP_aggregator_mean_representation_20160202_195640_i101_b100_r100.csv")

# Grouping polarization into three categories (less than 0.5. between 0.5-1. above 1).
mcp_ex_a1[, "polarization"] = NA
mcp_ex_a1[mcp_ex_a1$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_a1[0.5 < mcp_ex_a1$mu & mcp_ex_a1$mu < 1, ][, "polarization"] = 2
mcp_ex_a1[1 <= mcp_ex_a1$mu, ][, "polarization"] = 3
mcp_ex_a1$polarization = ordered(mcp_ex_a1$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

mcp_ex_a2[, "polarization"] = NA
mcp_ex_a2[mcp_ex_a2$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_a2[0.5 < mcp_ex_a2$mu & mcp_ex_a2$mu < 1, ][, "polarization"] = 2
mcp_ex_a2[1 <= mcp_ex_a2$mu, ][, "polarization"] = 3
mcp_ex_a2$polarization = ordered(mcp_ex_a2$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

mcp_ex_a3[, "polarization"] = NA
mcp_ex_a3[mcp_ex_a3$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_a3[0.5 < mcp_ex_a3$mu & mcp_ex_a3$mu < 1, ][, "polarization"] = 2
mcp_ex_a3[1 <= mcp_ex_a3$mu, ][, "polarization"] = 3
mcp_ex_a3$polarization = ordered(mcp_ex_a3$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))


  # 3.3.1 Mean eccentricity
  # All-aggregator mean eccentricity as function of number of firms in market.
  fig331a = ggplot(mcp_ex_a1, aes(y = MeanEst, x = N, colour = factor(polarization))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig331a + stat_smooth() + 
    geom_point(size = 1) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
    labs(y = "Mean eccentricity", x = "Number of firms", colour = NULL)
  ggsave("fig331a.pdf", width = 21, height = 16, units = "cm")
  
  # 3.3.2 ENP
  # All-aggregator effective number of firms compared to actual number of firms in market.
  fig332a = ggplot(mcp_ex_a2, aes(y=MeanEst, x=N, colour=factor(polarization))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig332a + stat_smooth() + 
    scale_colour_brewer(palette = "Dark2") +
    geom_abline(intercept = 0, slope = 1, color="gray") +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(1, 12), breaks = 2:12, minor_breaks = NULL, expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "Number of firms", colour = NULL)
  ggsave("fig332a.pdf", width = 21, height = 16, units = "cm")

  # 3.3.3 Mean representation
  # All-aggregator mean representation as function of number of firms in market.
  fig333a_bands = rqss_ci_bands(mcp_ex_a3, y="MeanEst", x="N", group="polarization")
  fig333a = ggplot(mcp_ex_a3, aes(y = MeanEst, x = N, colour = factor(polarization))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig333a + stat_quantile(quantiles=0.5, formula = y ~ qss(x), method = "rqss", size = 1) +
    geom_ribbon(data = fig333a_bands, 
                aes(x = N, ymin = MeanEstBLo, ymax = MeanEstBHi, group = factor(polarization), y = NULL, color = NULL), 
                alpha = 0.15) +
    geom_point(size = 1) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(-4, 0), expand = c(0, 0)) +
    labs(y = "Mean representation", x = "Number of firms", colour = NULL)
  ggsave("fig333a.pdf", width = 21, height = 16, units = "cm")


# 3.4 ALL-MAXCOV MODEL
mcp_ex_m1 = read.csv("data/MCP_maxcov_mean_eccentricity_20160203_121329_i151_b150_r50.csv")
mcp_ex_m2 = read.csv("data/MCP_maxcov_ENP_20160203_121329_i151_b150_r50.csv")
mcp_ex_m3 = read.csv("data/MCP_maxcov_mean_representation_20160203_121329_i151_b150_r50.csv")

# Grouping polarization into three categories (less than 0.5. between 0.5-1. above 1).
mcp_ex_m1[, "polarization"] = NA
mcp_ex_m1[mcp_ex_m1$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_m1[0.5 < mcp_ex_m1$mu & mcp_ex_m1$mu < 1, ][, "polarization"] = 2
mcp_ex_m1[1 <= mcp_ex_m1$mu, ][, "polarization"] = 3
mcp_ex_m1$polarization = ordered(mcp_ex_m1$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

mcp_ex_m2[, "polarization"] = NA
mcp_ex_m2[mcp_ex_m2$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_m2[0.5 < mcp_ex_m2$mu & mcp_ex_m2$mu < 1, ][, "polarization"] = 2
mcp_ex_m2[1 <= mcp_ex_m2$mu, ][, "polarization"] = 3
mcp_ex_m2$polarization = ordered(mcp_ex_m2$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

mcp_ex_m3[, "polarization"] = NA
mcp_ex_m3[mcp_ex_m3$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_m3[0.5 < mcp_ex_m3$mu & mcp_ex_m3$mu < 1, ][, "polarization"] = 2
mcp_ex_m3[1 <= mcp_ex_m3$mu, ][, "polarization"] = 3
mcp_ex_m3$polarization = ordered(mcp_ex_m3$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

  
  # 3.4.1 Mean eccentricity
  # All-maxcov mean eccentricity as function of number of firms in market.
  fig341a = ggplot(mcp_ex_m1, aes(y = MeanEst, x = N, colour = factor(polarization))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig341a + stat_smooth() + 
    geom_point(size = 1) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
    labs(y = "Mean eccentricity", x = "Number of firms", colour = NULL)
  ggsave("fig341a.pdf", width = 21, height = 16, units = "cm")

  # 3.4.2 ENP
  # All-maxcov effective number of firms compared to actual number of firms in market.
  fig342a = ggplot(mcp_ex_m2, aes(y=MeanEst, x=N, colour=factor(polarization))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig342a + stat_smooth() + 
    scale_colour_brewer(palette = "Dark2") +
    geom_abline(intercept = 0, slope = 1, color="gray") +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(1, 12), breaks = 2:12, minor_breaks = NULL, expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "Number of firms", colour = NULL)
  ggsave("fig342a.pdf", width = 21, height = 16, units = "cm")
  
  # 3.4.3 Mean representation
  # All-maxcov mean representation as function of number of firms in market.
  fig343a_bands = rqss_ci_bands(mcp_ex_m3, y="MeanEst", x="N", group="polarization")
  fig343a = ggplot(mcp_ex_m3, aes(y = MeanEst, x = N, colour = factor(polarization))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig343a + stat_quantile(quantiles=0.5, formula = y ~ qss(x), method = "rqss", size = 1) +
    geom_ribbon(data = fig343a_bands, 
                aes(x = N, ymin = MeanEstBLo, ymax = MeanEstBHi, group = factor(polarization), y = NULL, color = NULL), 
                alpha = 0.15) +
    geom_point(size = 1) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
    scale_y_continuous(limits = c(-4, 0), expand = c(0, 0)) +
    labs(y = "Mean representation", x = "Number of firms", colour = NULL)  
  ggsave("fig343a.pdf", width = 21, height = 16, units = "cm")
    

  
# 4. MCP: SAME DECISION RULE WITH ENDOGENOUS NUMBER OF FIRMS --------------


  
# 5. MCP: MIXED DECISION RULES WITH ENDOGENOUS NUMBER OF FIRMS ------------

# 5.1 MODEL 1
#mcp_en_m1 = read.csv("data/MCP_mixed_20160119_163128_mean_eccentricity_i101_b100_r50.csv")
#mcp_en_m2 = read.csv("data/MCP_mixed_20160119_163128_ENP_i101_b100_r50.csv")
#mcp_en_m3 = read.csv("data/MCP_mixed_20160119_163128_mean_share_i101_b100_r50.csv")
#mcp_en_m4 = read.csv("data/MCP_mixed_20160119_163128_N_i101_b100_r50.csv")
#mcp_en_m5 = read.csv("data/MCP_mixed_20160119_163128_mean_age_death_i101_b100_r50.csv")
mcp_en_m1 = read.csv("data/MCP_mixed2_20160120_154024_mean_eccentricity_i101_b100_r50")
mcp_en_m2 = read.csv("data/MCP_mixed2_20160120_154024_ENP_i101_b100_r50.csv")
mcp_en_m3 = read.csv("data/MCP_mixed2_20160120_154024_mean_share_i101_b100_r50.csv")
mcp_en_m4 = read.csv("data/MCP_mixed2_20160120_154024_N_i101_b100_r50.csv")
mcp_en_m5 = read.csv("data/MCP_mixed2_20160120_154024_mean_age_death_i101_b100_r50.csv")

# Adding decision rule labels.
mcp_en_m1$rule = factor(mcp_en_m1$rule, labels = c("Total", "Sticker", "Hunter", "Aggregator"))
mcp_en_m3$rule = factor(mcp_en_m3$rule, labels = c("Total", "Sticker", "Hunter", "Aggregator"))
mcp_en_m4$rule = factor(mcp_en_m4$rule, labels = c("Total", "Sticker", "Hunter", "Aggregator"))
mcp_en_m5$rule = factor(mcp_en_m5$rule, labels = c("Total", "Sticker", "Hunter", "Aggregator"))
crule3p1 = c("#000000", "#F8766D", "#00BA38", "#619CFF") # Three decisions rules + total (black).

  
  # 5.1.1 Mean eccentricity
  fig511a = ggplot(mcp_en_m1, aes(y=MeanEst, x=tau, colour = factor(rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig511a + stat_smooth() + 
    #geom_point() + 
    scale_color_manual(values = crule3p1) +
    scale_fill_manual(values = crule3p1) +
    scale_x_continuous(limits = c(0.05, 0.3), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 1.65), expand = c(0, 0)) +
    labs(y = "Mean eccentricity", x = "de facto survival threshold (percentage of market)")

  fig511b = ggplot(mcp_en_m1, aes(y=MeanEst, x=psi, colour = factor(rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig511b + stat_smooth() + 
    #geom_point() + 
    scale_color_manual(values = crule3p1) +
    scale_fill_manual(values = crule3p1) +
    scale_x_continuous(limits = c(10, 25), breaks=seq(10, 25, 5), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 1.65), expand = c(0, 0)) +
    labs(y = "Mean eccentricity", x = "Number of ticks per system tick / grace period")

  fig511c = ggplot(mcp_en_m1, aes(y=MeanEst, x=a_f, colour = factor(rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig511c + stat_smooth() + 
    #geom_point() + 
    scale_color_manual(values = crule3p1) +
    scale_fill_manual(values = crule3p1) +
    scale_x_continuous(limits = c(0, 0.9), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 1.65), expand = c(0, 0)) +
    labs(y = "Mean eccentricity", x = "Memory parameter of firms (weight on past fitness)")


  # 5.1.2 ENP
  fig512a = ggplot(mcp_en_m2, aes(y=MeanEst, x=tau)) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig512a + stat_smooth() + 
    geom_point() + 
    scale_x_continuous(limits = c(0.05, 0.3), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(1, 13), minor_breaks = NULL, expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "de facto survival threshold (percentage of market)")
  
  fig512b = ggplot(mcp_en_m2, aes(y=MeanEst, x=psi)) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig512b + stat_smooth() + 
    geom_point() + 
    scale_x_continuous(limits = c(10, 25), breaks=seq(10, 25, 5), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(1, 13), minor_breaks = NULL, expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "Number of ticks per system tick / grace period")
  
  fig512c = ggplot(mcp_en_m2, aes(y=MeanEst, x=a_f)) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig512c + stat_smooth() + 
    geom_point() + 
    scale_x_continuous(limits = c(0, 0.9), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(1, 13), minor_breaks = NULL, expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "Memory parameter of firms (weight on past fitness)")

  
  # 5.1.3 Mean share
  fig513a = ggplot(mcp_en_m3, aes(y=MeanEst, x=tau, colour = factor(rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig513a + stat_smooth() + 
    #geom_point() + 
    scale_color_manual(values = crule3p1) +
    scale_fill_manual(values = crule3p1) +
    scale_x_continuous(limits = c(0.05, 0.3), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
    labs(y = "Mean share of market", x = "de facto survival threshold (percentage of market)", colour="Decision rule")
  
  fig513b = ggplot(mcp_en_m3, aes(y=MeanEst, x=a_f, colour = factor(rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig513b + stat_smooth() + 
    #geom_point() + 
    scale_color_manual(values = crule3p1) +
    scale_fill_manual(values = crule3p1) +
    scale_x_continuous(limits = c(0, 0.9), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
    labs(y = "Mean share of market", x = "Memory parameter of firms (weight on past fitness)", colour="Decision rule")
  
  fig513c = ggplot(mcp_en_m3, aes(y=MeanEst, x=psi, colour = factor(rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig513c + stat_smooth() + 
    #geom_point() + 
    scale_color_manual(values = crule3p1) +
    scale_fill_manual(values = crule3p1) +
    scale_x_continuous(limits = c(10, 25), breaks=seq(10, 25, 5), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
    labs(y = "Mean share of market", x = "Number of ticks per system tick / grace period", colour="Decision rule")
  
  fig513d = ggplot(mcp_en_m3, aes(y=MeanEst, x=mu, colour = factor(rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig513d + stat_smooth() + 
    #geom_point() + 
    scale_color_manual(values = crule3p1) +
    scale_fill_manual(values = crule3p1) +
    scale_x_continuous(limits = c(0, 1.5), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
    labs(y = "Mean share of market", x = "Polarisation", colour="Decision rule")
  
  
  
  # 5.1.4 Number of surviving firms (N)
  fig514a = ggplot(mcp_en_m4, aes(y=MeanEst, x=tau, colour = factor(rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig514a + stat_smooth() + 
    #geom_point() + 
    scale_color_manual(values = crule3p1) +
    scale_fill_manual(values = crule3p1) +
    scale_x_continuous(limits = c(0.05, 0.3), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 15), minor_breaks = NULL, expand = c(0, 0)) +
    labs(y = "Number surviving firms (N)", x = "de facto survival threshold (percentage of market)", colour="Decision rule")

  fig514b = ggplot(mcp_en_m4, aes(y=MeanEst, x=a_f, colour = factor(rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig514b + stat_smooth() + 
    #geom_point() + 
    scale_color_manual(values = crule3p1) +
    scale_fill_manual(values = crule3p1) +
    scale_x_continuous(limits = c(0, 0.9), minor_breaks = NULL) +  
    scale_y_continuous(limits = c(0, 15), minor_breaks = NULL, expand = c(0, 0)) +
    labs(y = "Number surviving firms (N)", x = "Memory parameter of firms (weight on past fitness)", colour="Decision rule")
  
  
  # 5.1.5 Mean age at death
  fig515a = ggplot(mcp_en_m5, aes(y=MeanEst, x=tau, colour = factor(rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig515a + stat_smooth() + 
    #geom_point() + 
    scale_color_manual(values = crule3p1) +
    scale_fill_manual(values = crule3p1) +
    scale_x_continuous(limits = c(0.05, 0.3), minor_breaks = NULL) + 
    #scale_y_continuous(limits = c(0, 15), expand = c(0, 0)) +
    labs(y = "Mean firm age at death", x = "de facto survival threshold (percentage of market)", colour="Decision rule")

  fig515b = ggplot(mcp_en_m5, aes(y=MeanEst, x=a_f, colour = factor(rule))) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig515b + stat_smooth() + 
    #geom_point() + 
    scale_color_manual(values = crule3p1) +
    scale_fill_manual(values = crule3p1) +
    scale_x_continuous(limits = c(0, 0.9), minor_breaks = NULL) +  
    #scale_y_continuous(limits = c(0, 15), expand = c(0, 0)) +
    labs(y = "Mean firm age at death", x = "Memory parameter of firms (weight on past fitness)", colour="Decision rule")
  
  

# 6. MCP: MAXCOV-INDUCTOR & MAXCOV-INDUCTOR-GA MODEL ----------------------

mcp_ex_mi1 = read.csv("data/MCP_maxcov-inductor_mean_eccentricity_20160521_130042_i1001_psi1_b1000_r50.csv")
mcp_ex_mi2 = read.csv("data/MCP_maxcov-inductor_ENP_20160521_130042_i1001_psi1001_b1000_r50.csv")
mcp_ex_mi3 = read.csv("data/MCP_maxcov-inductor_mean_representation_20160521_130042_i1001_psi1001_b1000_r50.csv")
  
# Grouping polarization into three categories (less than 0.5. between 0.5-1. above 1).
mcp_ex_mi1[, "polarization"] = NA
mcp_ex_mi1[mcp_ex_mi1$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_mi1[0.5 < mcp_ex_mi1$mu & mcp_ex_mi1$mu < 1, ][, "polarization"] = 2
mcp_ex_mi1[1 <= mcp_ex_mi1$mu, ][, "polarization"] = 3
mcp_ex_mi1$polarization = ordered(mcp_ex_mi1$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

mcp_ex_mi2[, "polarization"] = NA
mcp_ex_mi2[mcp_ex_mi2$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_mi2[0.5 < mcp_ex_mi2$mu & mcp_ex_mi2$mu < 1, ][, "polarization"] = 2
mcp_ex_mi2[1 <= mcp_ex_mi2$mu, ][, "polarization"] = 3
mcp_ex_mi2$polarization = ordered(mcp_ex_mi2$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

mcp_ex_mi3[, "polarization"] = NA
mcp_ex_mi3[mcp_ex_mi3$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_mi3[0.5 < mcp_ex_mi3$mu & mcp_ex_mi3$mu < 1, ][, "polarization"] = 2
mcp_ex_mi3[1 <= mcp_ex_mi3$mu, ][, "polarization"] = 3
mcp_ex_mi3$polarization = ordered(mcp_ex_mi3$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

# Grouping polarization into unimodal and bimodal
mcp_ex_mi1$n_ratio_break = n_ratio_break(mcp_ex_mi1$mu);
mcp_ex_mi1[, "unimodal"] = 0
mcp_ex_mi1[is.na(mcp_ex_mi1$n_ratio_break) | mcp_ex_mi1$n_ratio >= mcp_ex_mi1$n_ratio_break, ][, "unimodal"] = 1
mcp_ex_mi1$unimodal = ordered(mcp_ex_mi1$unimodal, labels = c("Bimodal", "Unimodal"))

# Grouping polarization into unimodal and bimodal
mcp_ex_mi1$n_ratio_break = n_ratio_break(mcp_ex_mi1$mu);
mcp_ex_mi1[, "unimodal"] = 0
mcp_ex_mi1[is.na(mcp_ex_mi1$n_ratio_break) | mcp_ex_mi1$n_ratio >= mcp_ex_mi1$n_ratio_break, ][, "unimodal"] = 1
mcp_ex_mi1$unimodal = ordered(mcp_ex_mi1$unimodal, labels = c("Bimodal", "Unimodal"))

    # 6.1.1 Mean eccentricity
    # Maxcov-inductor mean eccentricity as function of number of firms in market.
    fig611a = ggplot(mcp_ex_mi1, aes(y = MeanEst, x = N, colour = factor(polarization))) + 
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    fig611a + stat_smooth() + 
      geom_point(size=1) + #geom_point(aes(size = n_ratio)) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
      labs(y = "Mean eccentricity", x = "Number of firms", colour = NULL)
    ggsave("fig611a.pdf", width = 21, height = 16, units = "cm")
    
    fig611b = ggplot(mcp_ex_mi1, aes(y = MeanEst, x = N, colour = factor(unimodal))) + 
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    fig611b + stat_smooth() + 
      geom_point(size = 1) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
      labs(y = "Mean eccentricity", x = "Number of firms", colour = NULL)
    
    fig611c = ggplot(mcp_ex_mi1[which(mcp_ex_mi1$unimodal=="Unimodal"), ], aes(y = MeanEst, x = n_ratio, colour = factor(N))) + 
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    fig611c + stat_smooth(se=FALSE) + 
      geom_point(size = 1) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_x_continuous(limits = c(1, 2), minor_breaks = NULL) + 
      scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
      labs(y = "Mean eccentricity", x = "Relative size of subpopulation", colour = NULL)  
    
    # 6.1.2 ENP
    # Maxcov-inductor effective number of firms compared to actual number of firms in market.
    fig612a = ggplot(mcp_ex_mi2, aes(y=MeanEst, x=N, colour=factor(polarization))) +
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    fig612a + stat_smooth() + 
      scale_colour_brewer(palette = "Dark2") +
      geom_abline(intercept = 0, slope = 1, color="gray") +
      #stat_function(fun = eqcond1_ENP, colour = "black") + 
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(1, 12), breaks = 2:12, minor_breaks = NULL, expand = c(0, 0)) +
      labs(y = "Effective number of firms (ENP)", x = "Number of firms", colour = NULL)
    ggsave("fig612a.pdf", width = 21, height = 16, units = "cm")
    
    # 6.1.3 Mean representation
    # Maxcov-inductor mean representation as function of number of firms in market.
    #fig613a_bands = rqss_ci_bands(mcp_ex_mi3, y="MeanEst", x="N", group="polarization")
    fig613a = ggplot(mcp_ex_mi3, aes(y = MeanEst, x = N, colour = factor(polarization))) + 
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    fig613a + stat_quantile(quantiles=0.5, formula = y ~ qss(x), method = "rqss", size = 1) +
      #geom_ribbon(data = fig613a_bands, 
      #            aes(x = N, ymin = MeanEstBLo, ymax = MeanEstBHi, group = factor(polarization), y = NULL, color = NULL), 
      #            alpha = 0.15) +
      geom_point(size = 1) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(-4, 0), expand = c(0, 0)) +
      labs(y = "Mean representation", x = "Number of firms", colour = NULL)  
    ggsave("fig613a.pdf", width = 21, height = 16, units = "cm")
    
    
mcp_ex_miga1 = read.csv("data/MCP_maxcov-inductor-GA_mean_eccentricity_20160518_095112_i50_psi50_b49_r50.csv")
mcp_ex_miga2 = read.csv("data/MCP_maxcov-inductor-GA_ENP_20160518_095112_i50_psi50_b49_r50.csv")
mcp_ex_miga3 = read.csv("data/MCP_maxcov-inductor-GA_mean_representation_20160518_095112_psi50_i50_b49_r50.csv")

# Grouping polarization into three categories (less than 0.5. between 0.5-1. above 1).
mcp_ex_miga1[, "polarization"] = NA
mcp_ex_miga1[mcp_ex_miga1$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_miga1[0.5 < mcp_ex_miga1$mu & mcp_ex_miga1$mu < 1, ][, "polarization"] = 2
mcp_ex_miga1[1 <= mcp_ex_miga1$mu, ][, "polarization"] = 3
mcp_ex_miga1$polarization = ordered(mcp_ex_miga1$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

mcp_ex_miga2[, "polarization"] = NA
mcp_ex_miga2[mcp_ex_miga2$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_miga2[0.5 < mcp_ex_miga2$mu & mcp_ex_miga2$mu < 1, ][, "polarization"] = 2
mcp_ex_miga2[1 <= mcp_ex_miga2$mu, ][, "polarization"] = 3
mcp_ex_miga2$polarization = ordered(mcp_ex_miga2$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

mcp_ex_miga3[, "polarization"] = NA
mcp_ex_miga3[mcp_ex_miga3$mu <= 0.5, ][, "polarization"] = 1
mcp_ex_miga3[0.5 < mcp_ex_miga3$mu & mcp_ex_miga3$mu < 1, ][, "polarization"] = 2
mcp_ex_miga3[1 <= mcp_ex_miga3$mu, ][, "polarization"] = 3
mcp_ex_miga3$polarization = ordered(mcp_ex_miga3$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))

    # 6.2.1 Mean eccentricity
    # Maxcov-inductor-GA mean eccentricity as function of number of firms in market.
    fig621a = ggplot(mcp_ex_miga1, aes(y = MeanEst, x = N, colour = factor(polarization))) + 
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    fig621a + stat_smooth() + 
      geom_point(size = 1) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
      labs(y = "Mean eccentricity", x = "Number of firms", colour = NULL)
    ggsave("fig621a.pdf", width = 21, height = 16, units = "cm")
    
    fig621c = ggplot(mcp_ex_miga1, aes(y = MeanEst, x = n_ratio, colour = factor(N))) +
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    fig621c + stat_smooth(aes(fill = factor(N)), se=F) + 
      geom_point(size = 1) + 
      #scale_color_manual(values=cfirm4) +
      #scale_fill_manual(values = cfirm4, guide="none") +
      scale_x_continuous(limits = c(1, 2), minor_breaks = NULL) + 
      scale_y_continuous(limits = c(0, 1.3), expand = c(0, 0)) +
      labs(y = "Mean eccentricity", x = "Relative size of left subpopulation", colour = "Number of firms")
    
    # 6.2.2 ENP
    # Maxcov-inductor-GA effective number of firms compared to actual number of firms in market.
    fig622a = ggplot(mcp_ex_miga2, aes(y=MeanEst, x=N, colour=factor(polarization))) +
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    fig622a + stat_smooth() + 
      scale_colour_brewer(palette = "Dark2") +
      geom_abline(intercept = 0, slope = 1, color="gray") +
      stat_function(fun = eqcond1_ENP, colour = "black") + 
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(1, 12), breaks = 2:12, minor_breaks = NULL, expand = c(0, 0)) +
      labs(y = "Effective number of firms (ENP)", x = "Number of firms", colour = NULL)
    ggsave("fig622a.pdf", width = 21, height = 16, units = "cm")
    
    # 6.2.3 Mean representation
    # Maxcov-inductor-GA mean representation as function of number of firms in market.
    #fig623a_bands = rqss_ci_bands(mcp_ex_miga3, y="MeanEst", x="N", group="polarization")
    fig623a = ggplot(mcp_ex_miga3, aes(y = MeanEst, x = N, colour = factor(polarization))) + 
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    fig623a + stat_quantile(quantiles=0.5, formula = y ~ qss(x), method = "rqss", size = 1) +
      #geom_ribbon(data = fig623a_bands, 
      #            aes(x = N, ymin = MeanEstBLo, ymax = MeanEstBHi, group = factor(polarization), y = NULL, color = NULL), 
      #            alpha = 0.15) +
      geom_point(size = 1) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(-4, 0), expand = c(0, 0)) +
      labs(y = "Mean representation", x = "Number of firms", colour = NULL)  
    ggsave("fig623a.pdf", width = 21, height = 16, units = "cm")
    
    
    
# APPENDIX A -- MAXCOVRND COMPARISON
    
    # A.1 Symmetric distribution
    gs_m1 = read.csv("data/GS_maxcov_mean_eccentricity_20160203_231850_r1000_b99.csv")
    gs_mr1 = read.csv("data/GS_maxcovrnd_mean_eccentricity_20160504_221217_r1000_b99.csv")
    gs_m2 = read.csv("data/GS_maxcov_ENP_20160203_231850_r1000_b99.csv")
    gs_mr2 = read.csv("data/GS_maxcovrnd_ENP_20160504_221217_r1000_b99.csv")
    gs_m3 = read.csv("data/GS_maxcov_mean_representation_20160203_231850_r1000_b99.csv")
    gs_mr3 = read.csv("data/GS_maxcovrnd_mean_representation_20160504_221217_r1000_b99.csv")
    
    # Merge the grid-sweep datasets of every decision rule into one.
    gs_1 = as.data.frame(mapply(c, 
                                data.frame(gs_m1[,1:3], 4),
                                data.frame(gs_mr1[,1:3], 5) ))
    colnames(gs_1)[4] = "Rule"
    gs_1$Rule = factor(gs_1$Rule, labels = c("All-maxcov", "All-maxcovrnd"))
    
    gs_2 = as.data.frame(mapply(c, 
                                data.frame(gs_m2[,1:3], 4),
                                data.frame(gs_mr2[,1:3], 5) ))
    colnames(gs_2)[4] = "Rule"
    gs_2$Rule = factor(gs_2$Rule, labels = c("All-maxcov", "All-maxcovrnd"))
    
    gs_3 = as.data.frame(mapply(c, 
                                data.frame(gs_m3[,1:3], 4),
                                data.frame(gs_mr3[,1:3], 5) ))
    colnames(gs_3)[4] = "Rule"
    gs_3$Rule = factor(gs_3$Rule, labels = c("All-maxcov", "All-maxcovrnd"))
    
    
    # A.1.1 Mean eccentricity
    figA11 = ggplot(gs_1, aes(y = MeanEst, x = N, linetype = factor(Rule))) + 
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    figA11 + geom_line(size=1) + 
      #geom_ribbon(aes(x = N, ymin = MeanEst-StdDev, ymax = pmin(gs_1$MeanEst+gs_1$StdDev, rep(1.7, nrow(gs_1))), group = factor(Rule), fill=Rule, color = NULL), alpha = 0.1) +
      scale_fill_discrete(guide=FALSE) +
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
      labs(y = "Mean eccentricity", x = "Number of firms", colour = NULL, linetype = NULL)
    ggsave("figA11.pdf", width = 21, height = 16, units = "cm")
    
    # A.1.2 ENP
    figA12 = ggplot(gs_2, aes(y=MeanEst, x=N, linetype=factor(Rule))) +
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    figA12 + geom_line(size=1) + 
      #geom_ribbon(aes(x = N, ymin = pmax(gs_2$MeanEst-gs_2$StdDev, rep(1, nrow(gs_2))), ymax = pmin(gs_2$MeanEst+gs_2$StdDev, rep(12, nrow(gs_2))), group = factor(Rule), fill=Rule, color = NULL), alpha = 0.1) +
      scale_fill_discrete(guide=FALSE) +
      geom_abline(intercept = 0, slope = 1, color="gray") +
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(1, 12), breaks = 2:12, minor_breaks = NULL, expand = c(0, 0)) +
      labs(y = "Effective number of firms (ENP)", x = "Number of firms", colour = NULL, linetype = NULL)
    ggsave("figA12.pdf", width = 21, height = 16, units = "cm")
    
    # A.1.3 Mean representation
    #fig23a_bands = rqss_ci_bands(gs_3, y="MeanEst", x="N", group="Rule")
    figA13 = ggplot(gs_3, aes(y = MeanEst, x = N, linetype = factor(Rule))) + 
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    figA13 + geom_line(size=1) + 
      #geom_ribbon(aes(x = N, ymin = MeanEst-StdDev, ymax = MeanEst+StdDev, group = factor(Rule), fill = Rule, color = NULL), alpha = 0.1) +
      scale_fill_discrete(guide=FALSE) +
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(-4, 0), expand = c(0, 0)) +
      labs(y = "Mean representation", x = "Number of firms", colour = NULL, linetype = NULL)
    ggsave("figA13.pdf", width = 21, height = 16, units = "cm")
    
    
    # A.2 Asymmetric distribution
    mcp_ex_m1 = read.csv("data/MCP_maxcov_mean_eccentricity_20160203_121329_i151_b150_r50.csv")
    mcp_ex_m2 = read.csv("data/MCP_maxcov_ENP_20160203_121329_i151_b150_r50.csv")
    mcp_ex_m3 = read.csv("data/MCP_maxcov_mean_representation_20160203_121329_i151_b150_r50.csv")
    mcp_ex_mr1 = read.csv("data/MCP_maxcovrnd_mean_eccentricity_20160505_093227_i151_b150_r50.csv")
    mcp_ex_mr2 = read.csv("data/MCP_maxcovrnd_ENP_20160505_093227_i151_b150_r50.csv")
    mcp_ex_mr3 = read.csv("data/MCP_maxcovrnd_mean_representation_20160505_093227_i151_b150_r50.csv")
    
    mcp_ex1 = as.data.frame(mapply(c, 
                                data.frame(mcp_ex_m1[,1:5], 4),
                                data.frame(mcp_ex_mr1[,1:5], 5) ))
    colnames(mcp_ex1)[6] = "Rule"
    mcp_ex1$Rule = factor(mcp_ex1$Rule, labels = c("All-maxcov", "All-maxcovrnd"))
    
    mcp_ex2 = as.data.frame(mapply(c, 
                                   data.frame(mcp_ex_m2[,1:5], 4),
                                   data.frame(mcp_ex_mr2[,1:5], 5) ))
    colnames(mcp_ex2)[6] = "Rule"
    mcp_ex2$Rule = factor(mcp_ex2$Rule, labels = c("All-maxcov", "All-maxcovrnd"))
    
    mcp_ex3 = as.data.frame(mapply(c, 
                                   data.frame(mcp_ex_m3[,1:5], 4),
                                   data.frame(mcp_ex_mr3[,1:5], 5) ))
    colnames(mcp_ex3)[6] = "Rule"
    mcp_ex3$Rule = factor(mcp_ex3$Rule, labels = c("All-maxcov", "All-maxcovrnd"))
    
    # Grouping polarization into three categories (less than 0.5. between 0.5-1. above 1).
    mcp_ex1[, "polarization"] = NA
    mcp_ex1[mcp_ex1$mu <= 0.5, ][, "polarization"] = 1
    mcp_ex1[0.5 < mcp_ex1$mu & mcp_ex1$mu < 1, ][, "polarization"] = 2
    mcp_ex1[1 <= mcp_ex1$mu, ][, "polarization"] = 3
    mcp_ex1$polarization = ordered(mcp_ex1$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))
    
    mcp_ex2[, "polarization"] = NA
    mcp_ex2[mcp_ex2$mu <= 0.5, ][, "polarization"] = 1
    mcp_ex2[0.5 < mcp_ex2$mu & mcp_ex2$mu < 1, ][, "polarization"] = 2
    mcp_ex2[1 <= mcp_ex2$mu, ][, "polarization"] = 3
    mcp_ex2$polarization = ordered(mcp_ex2$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))
    
    mcp_ex3[, "polarization"] = NA
    mcp_ex3[mcp_ex3$mu <= 0.5, ][, "polarization"] = 1
    mcp_ex3[0.5 < mcp_ex3$mu & mcp_ex3$mu < 1, ][, "polarization"] = 2
    mcp_ex3[1 <= mcp_ex3$mu, ][, "polarization"] = 3
    mcp_ex3$polarization = ordered(mcp_ex3$polarization, labels = c("Polarization ≤ 0.5", "0.5 < Polarization < 1", "1 ≤ Polarization"))
    
    # A.2.1 Mean eccentricity
    # All-maxcov mean eccentricity as function of number of firms in market.
    figA21 = ggplot(mcp_ex1, aes(y = MeanEst, x = N, colour = factor(polarization), linetype=factor(Rule))) + 
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    figA21 + stat_smooth(se = FALSE) + 
      #geom_point(size = 1) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
      labs(y = "Mean eccentricity", x = "Number of firms", colour = NULL, linetype = NULL)
    ggsave("figA21.pdf", width = 21, height = 16, units = "cm")
    
    # A.2.2 ENP
    # All-maxcov effective number of firms compared to actual number of firms in market.
    figA22 = ggplot(mcp_ex2, aes(y=MeanEst, x=N, colour=factor(polarization), linetype=factor(Rule))) +
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    figA22 + stat_smooth(se = FALSE) + 
      scale_colour_brewer(palette = "Dark2") +
      geom_abline(intercept = 0, slope = 1, color="gray") +
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(1, 12), breaks = 2:12, minor_breaks = NULL, expand = c(0, 0)) +
      labs(y = "Effective number of firms (ENP)", x = "Number of firms", colour = NULL, linetype = NULL)
    ggsave("figA22.pdf", width = 21, height = 16, units = "cm")
    
    # A.2.3 Mean representation
    # All-maxcov mean representation as function of number of firms in market.
    #fig343a_bands = rqss_ci_bands(mcp_ex3, y="MeanEst", x="N", group="polarization")
    figA23 = ggplot(mcp_ex3, aes(y = MeanEst, x = N, colour = factor(polarization), linetype=factor(Rule))) +
      theme_minimal() +
      theme(legend.position = "top", 
            legend.box = "horizontal", 
            legend.key = element_rect(fill = NA, colour = NA) )
    figA23 + stat_quantile(quantiles=0.5, formula = y ~ qss(x), method = "rqss", size = 1) +
      #geom_ribbon(data = fig343a_bands, 
      #            aes(x = N, ymin = MeanEstBLo, ymax = MeanEstBHi, group = factor(polarization), y = NULL, color = NULL), 
      #            alpha = 0.15) +
      #geom_point(size = 1) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_x_continuous(limits = c(2, 12), breaks = 2:12, minor_breaks = NULL) + 
      scale_y_continuous(limits = c(-4, 0), expand = c(0, 0)) +
      labs(y = "Mean representation", x = "Number of firms", colour = NULL, linetype = NULL)
    ggsave("figA23.pdf", width = 21, height = 16, units = "cm")
    
    