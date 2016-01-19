# FIGURES
# Sets up figures and formats them for final paper.
# Jonas K. Sekamane


# TO-DO -------------------------------------------------------------------

# stat_smooth ==> Median spline


# 1. SETUP ----------------------------------------------------------------

# 1.1 Libraries
library(ggplot2)

# 1.2 Working directory
setwd("~/Dropbox/Economics/Courses/Thesis/Positioning/Figures")



# 2. GRID SWEEP: SAME DECISION RULE MODELS --------------------------------

# ...



# 3. MCP: SAME DECISION RULE WITH EXOGENOUS NUMBER OF FIRMS ---------------

# 3.1 ALL-STICKER MODEL
# ...



# 3.2 ALL-HUNTER MODEL
mcp_ex_h1 = read.csv("data/MCP_hunter_mean_eccentricity_20160110_210442_i250_b150.csv")
mcp_ex_h2 = read.csv("data/MCP_hunter_ENP_20160110_210442_i250_b150.csv")

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

# Subset of data only consisting 2, 3, 4 and 12 firms.
mcp_ex_h1_subset = subset(mcp_ex_h1, N %in% c(2,3,4,12))


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
  
  # All-hunter mean eccentricity as function the relative subpopulation size for a subset number of firms.
  cfirm4 = c("#253494", "#2c7fb8", "#41b6c4", "#a1dab4") # YlGnBu 5 reversed (but only using first 4). http://colorbrewer2.org
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
    scale_y_continuous(limits = c(1, 12), breaks = 2:12, expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "Number of firms", colour = NULL)

  
# 3.3 ALL-AGGREGATOR MODEL
mcp_ex_a1 = read.csv("data/MCP_aggregator_mean_eccentricity_20160111_185750_i101_b100_r100.csv")
mcp_ex_a2 = read.csv("data/MCP_aggregator_ENP_20160111_185750_i101_b100_r100.csv")

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
    scale_y_continuous(limits = c(1, 12), breaks = 2:12, expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "Number of firms", colour = NULL)


  
# 4. MCP: SAME DECISION RULE WITH ENDOGENOUS NUMBER OF FIRMS --------------


  
# 5. MCP: MIXED DECISION RULES WITH ENDOGENOUS NUMBER OF FIRMS ------------

# 5.1 MODEL 1
mcp_en_m1 = read.csv("data/MCP_mixed_mean_eccentricity_20160115_015809_i101_b100_r50.csv")
mcp_en_m2 = read.csv("data/MCP_mixed_ENP_20160115_015809_i101_b100_r50.csv")
  
  # 5.1.1 Mean eccentricity
  fig511a = ggplot(mcp_en_m1, aes(y=MeanEst, x=tau)) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig511a + stat_smooth() + 
    geom_point() + 
    scale_x_continuous(limits = c(0.05, 0.3), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 1.65), expand = c(0, 0)) +
    labs(y = "Mean eccentricity", x = "de facto survival threshold (percentage of market)")

  fig511b = ggplot(mcp_en_m1, aes(y=MeanEst, x=psi)) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig511b + stat_smooth() + 
    geom_point() + 
    scale_x_continuous(limits = c(10, 25), breaks=seq(10, 25, 5), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(0, 1.65), expand = c(0, 0)) +
    labs(y = "Mean eccentricity", x = "Number of ticks per system tick / grace period")

  fig511c = ggplot(mcp_en_m1, aes(y=MeanEst, x=a_f)) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig511c + stat_smooth() + 
    geom_point() + 
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
    scale_y_continuous(limits = c(1, 13), expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "de facto survival threshold (percentage of market)")
  
  fig512b = ggplot(mcp_en_m2, aes(y=MeanEst, x=psi)) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig512b + stat_smooth() + 
    geom_point() + 
    scale_x_continuous(limits = c(10, 25), breaks=seq(10, 25, 5), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(1, 13), expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "Number of ticks per system tick / grace period")
  
  fig512c = ggplot(mcp_en_m2, aes(y=MeanEst, x=a_f)) +
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA) )
  fig512c + stat_smooth() + 
    geom_point() + 
    scale_x_continuous(limits = c(0, 0.9), minor_breaks = NULL) + 
    scale_y_continuous(limits = c(1, 13), expand = c(0, 0)) +
    labs(y = "Effective number of firms (ENP)", x = "Memory parameter of firms (weight on past fitness)")




