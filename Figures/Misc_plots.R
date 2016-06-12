# FIGURES
# Sets up figures and formats them for final paper.
# Jonas K. Sekamane


# TO-DO -------------------------------------------------------------------


# 1. SETUP ----------------------------------------------------------------

# 1.1 Libraries
library('ggplot2')
#library('plyr') 
#library('reshape2') 
#library('quantreg')
library(grid)
library(gridExtra)
library('quantreg')
require(scales)
library(deldir)



# 1.2 Working directory
setwd("~/Dropbox/Economics/Courses/Thesis/Positioning/Figures")

# 1.3 Functions

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# 2. DETERMINING BURN IN --------------------------------------------------

# 2.1 Determenistic time-homogenous Markov chain (single state) - All-hunter
burn_det = read.csv("data/All-aggregator_N12_mu15_nratio2_i100.csv")

figb21a = ggplot(burn_det, aes(y = mean_eccentricity, x = iteration)) + 
  theme_minimal() +
  theme(legend.position = "top", 
        legend.box = "horizontal", 
        legend.key = element_rect(fill = NA, colour = NA) )
figb21a + geom_line(size=1, colour="#f35e5a") + 
  scale_x_continuous(limits = c(1, 100), minor_breaks = NULL) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0), minor_breaks = NULL) +
  labs(y = "Mean eccentricity", x = "Iteration")
ggsave("figb21a.pdf", width = 12, height = 9.6, units = "cm")


# 2.2 Stocastic time-homogenous Markov chain - All-hunter
burn_stc = read.csv("data/All-hunter_N12_mu0_nratio1_i500.csv")
burn_stc_2nd_mean = mean(burn_stc[c(251:500),3])
burn_stc_2nd_std = sd(burn_stc[c(251:500),3])

#burn_stc$mean_eccentricity > burn_stc_2nd_mean+burn_stc_2nd_std

figb22a = ggplot(burn_stc, aes(y = mean_eccentricity, x = iteration)) + 
  theme_minimal() +
  theme(legend.position = "top", 
        legend.box = "horizontal", 
        legend.key = element_rect(fill = NA, colour = NA) )
figb22a + geom_line(size=1, colour="#f35e5a") + 
  scale_x_continuous(limits = c(1, 500), minor_breaks = NULL) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0), minor_breaks = NULL) +
  labs(y = "Mean eccentricity", x = "Iteration") +
  geom_hline(yintercept=burn_stc_2nd_mean, colour="black") +
  geom_hline(yintercept=burn_stc_2nd_mean-burn_stc_2nd_std, colour="black", linetype="dashed") +
  geom_hline(yintercept=burn_stc_2nd_mean+burn_stc_2nd_std, colour="black", linetype="dashed")
ggsave("figb22a.pdf", width = 12, height = 9.6, units = "cm")


# 2.3 Stocastic time-homogenous Markov chain (time average not representative) - All-hunter
# data_mean_eccentricity_20160420_180651.csv
# filter-time tool
# figb23a.pdf, figb23b.pdf



# 3. FIRM MOVEMENT --------------------------------------------------------

# 3.1 -- ALL-HUNTER
moves = read.csv("data/xy_all-hunter_N5_mu15_nratio2_i100.csv")

move_panel = function(df, start, end, tick) {
  # Number of firms
  N = 5
  # Hide or show the tick labes on the x and y-axis.
  if (tick) ticklabelcolor = "black" else ticklabelcolor = "white"
  
  figure = ggplot(df[((start-1)*N+1):(end*N),], aes(y = y, x = x, colour = factor(firm))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA),
          axis.text.x = element_text(colour = ticklabelcolor), 
          axis.text.y = element_text(colour = ticklabelcolor),
          plot.title = element_text(size=10, hjust = 0),
          plot.margin = unit(c(0.2, 0, -0.4, 0), "lines") ) + 
    geom_point(data = subset(df[((start-1)*N+1):(end*N),], iteration == end), size = 2) +
    geom_line() +
    scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1.5), minor_breaks = NULL, name="") + 
    scale_y_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2, 2, by = 1), minor_breaks = NULL, name="") +
    scale_colour_discrete(guide = FALSE) +
    geom_point(aes(x = -0.5, y = 0), shape = 3, colour = "black")+
    labs(title = paste("Iteration", sprintf("%1.0f", end), sep =" "))
  
  return(figure)
}

#figmfull = move_panel(moves, 1, 100, TRUE)
#figmfull

p1 = move_panel(moves, 1, 10, TRUE)
p2 = move_panel(moves, 11, 20, FALSE)
p3 = move_panel(moves, 21, 30, FALSE)
p4 = move_panel(moves, 31, 40, FALSE)
p5 = move_panel(moves, 41, 50, FALSE)
p6 = move_panel(moves, 51, 60, FALSE)
p7 = move_panel(moves, 61, 70, FALSE)
p8 = move_panel(moves, 71, 80, FALSE)
p9 = move_panel(moves, 81, 90, FALSE)
p10 = move_panel(moves, 91, 100, TRUE)

lay = t(matrix(seq(1,10), 2, 5, byrow = T))
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol=2, layout_matrix = lay)

figm = arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol=2, layout_matrix = lay)
ggsave(file="figm.pdf", figm, width = 22, height = 25.5, units = "cm") #saves g


# 3.2 -- MAXCOV VS MAXCOV-INDUCTOR

moves_m = read.csv("data/xy_maxcov_N4_mu15_nratio15_pbi49.csv")

move_panel2 = function(df, start, end, tick, burnin) {
  # Number of firms
  N = 4
  # Hide or show the tick labes on the x and y-axis.
  if (tick) ticklabelcolor = "black" else ticklabelcolor = "white"
  
  endpoints = subset(df[((start-1)*N+1):(end*N),], iteration == burnin+end)
  #voronoi = deldir(endpoints$x, endpoints$y, rw = c(-3, 3, -2, 2))
  
  figure = ggplot(df[((start-1)*N+1):(end*N),], aes(y = y, x = x, colour = factor(firm))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA),
          axis.text.x = element_text(colour = ticklabelcolor), 
          axis.text.y = element_text(colour = ticklabelcolor),
          plot.title = element_text(size=10, hjust = 0),
          plot.margin = unit(c(0.2, 0, -0.4, 0), "lines") ) + 
    geom_point(data = endpoints, size = 2) +
    geom_line() +
    scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1.5), minor_breaks = NULL, name="") + 
    scale_y_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2, 2, by = 1), minor_breaks = NULL, name="") +
    scale_colour_discrete(guide = FALSE) +
    geom_point(aes(x = -0.3, y = 0), shape = 3, colour = "black") +
    #geom_segment(data = voronoi$dirsgs, aes(x = x1, y = y1, xend = x2, yend = y2), size = 0.3, linetype = 1, color= "#000000") + 
    labs(title = paste("Iteration", sprintf("%1.0f", burnin+end), sep =" "))
  
  return(figure)
}

moves_m$firm = rep(c(1, 4, 3, 2), 50) # matching colors
figm1 = move_panel2(moves_m, 1, 5, TRUE, 150) # line is location at the last 4 iterations.
figm1
ggsave(file="figm_maxcov.pdf", figm1, width = 10, height = 6.6, units = "cm") #saves g


moves_mi = read.csv("data/xy_maxcov-inductor_N4_mu15_nratio15_pbi49.csv")
moves_mi$firm = rep(c(3, 4, 1, 2), 50) # match colors with MAXCOV
figmi1 = move_panel2(moves_mi, 1, 7, TRUE, 1000) # line is location at the last 6 iterations.
figmi1
ggsave(file="figm_maxcov-inductor.pdf", figmi1, width = 10, height = 6.6, units = "cm") #saves g


# 3.3 -- CLUSTERING

moves_cmu = read.csv("data/xy_clustering_maxcov_N12_mu0_nratio1_pbi50.csv")
moves_cmb = read.csv("data/xy_clustering_maxcov_N12_mu15_nratio2_pbi50.csv")
moves_cmiu = read.csv("data/xy_clustering_maxcov-inductor_N12_mu0_nratio1_pbi50.csv")
moves_cmib = read.csv("data/xy_clustering_maxcov-inductor_N12_mu15_nratio2_pbi50.csv")

move_panel3 = function(df, start, end, tick, burnin) {
  # Number of firms
  N = 12
  # Hide or show the tick labes on the x and y-axis.
  if (tick) ticklabelcolor = "black" else ticklabelcolor = "white"
  
  endpoints = subset(df[((start-1)*N+1):(end*N),], iteration == burnin+end)
  voronoi = deldir(endpoints$x, endpoints$y, rw = c(-3, 3, -2, 2))
  
  figure = ggplot(df[((start-1)*N+1):(end*N),], aes(y = y, x = x, colour = factor(firm))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA),
          axis.text.x = element_text(colour = ticklabelcolor), 
          axis.text.y = element_text(colour = ticklabelcolor),
          plot.title = element_text(size=10, hjust = 0),
          plot.margin = unit(c(0.2, 0, -0.4, 0), "lines") ) + 
    geom_segment(data = voronoi$dirsgs, aes(x = x1, y = y1, xend = x2, yend = y2), size = 0.3, linetype = 1, color= "#000000") + 
    geom_line(alpha = 0.5) +
    geom_point(data = endpoints, size = 2) +
    scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1.5), minor_breaks = NULL, name="") + 
    scale_y_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2, 2, by = 1), minor_breaks = NULL, name="") +
    scale_colour_discrete(guide = FALSE) +
    #geom_point(aes(x = -0.3, y = 0), shape = 3, colour = "black") +
    labs(title = paste("Iteration", sprintf("%1.0f", burnin+end), sep =" "))
  
  return(figure)
}

figm_cmu = move_panel3(moves_cmu, 1, 10, TRUE, 150) # line is location at the last 6 iterations.
figm_cmu
ggsave(file="figm_cmu.pdf", figm_cmu, width = 10, height = 6.6, units = "cm")

figm_cmb = move_panel3(moves_cmb, 1, 10, TRUE, 150) # line is location at the last 6 iterations.
figm_cmb
ggsave(file="figm_cmb.pdf", figm_cmb, width = 10, height = 6.6, units = "cm")

figm_cmiu = move_panel3(moves_cmiu, 1, 10, TRUE, 1000) # line is location at the last 6 iterations.
figm_cmiu
ggsave(file="figm_cmiu.pdf", figm_cmiu, width = 10, height = 6.6, units = "cm")

figm_cmib = move_panel3(moves_cmib, 1, 10, TRUE, 1000) # line is location at the last 6 iterations.
figm_cmib
ggsave(file="figm_cmib.pdf", figm_cmib, width = 10, height = 6.6, units = "cm")

#width = 10, height = 6.6


# 4. Median spline - Eccentricity vs. share -------------------------------

allhunter = read.csv("data/All-hunter_20160502_231508_N5_mu0_n1_b150_i6150.csv")
cfirm5 = c("#000000", "#000000", "#000000", "#000000", "#000000")#, "#000000")

fig3ms = ggplot(allhunter, aes(y = share, x = eccentricity, colour = factor(firm))) + 
  theme_minimal() +
  theme(legend.position = "none", 
        legend.box = "horizontal", 
        legend.key = element_rect(fill = NA, colour = NA) )
fig3ms + geom_vline(xintercept=0.4416208, colour="red") +
  stat_smooth(se = FALSE) + 
  scale_color_manual(values=cfirm5) +
  scale_x_continuous(limits = c(0, 1.55), minor_breaks = NULL) + 
  scale_y_continuous(expand = c(0, 0), labels=percent) +
  labs(y = "Share of market", x = "Eccentricity", colour = NULL)
ggsave("fig3ms.pdf", width = 14, height = 12, units = "cm")

#geom_segment(aes(x = 0.4416208, y = 0, xend = 0.4416208, yend = 0.216), color = "red") +
#geom_quantile(quantiles=0.5, formula = y ~ qss(x, lambda=5), method = "rqss", size = 1) +
#geom_point(size = 0.1) + 



# 5. Long-run market share (three firms) -------------------------------

s_m = read.csv("data/three_maxcov_shares_N3_mu15_nratio2_i250.csv", header = FALSE)
s_mi = read.csv("data/three_maxcov-inductor_shares_N3_mu15_nratio2_i2500.csv", header = FALSE)
s_miga = read.csv("data/three_maxcov-inductor-ga_shares_N3_mu15_nratio2_i2500.csv", header = FALSE)
cf_mi_1 = read.csv("data/three_maxcov-inductor_cf_N3_mu15_nratio2_i2500_1.csv", header = FALSE)
cf_mi_2 = read.csv("data/three_maxcov-inductor_cf_N3_mu15_nratio2_i2500_2.csv", header = FALSE)
cf_mi_3 = read.csv("data/three_maxcov-inductor_cf_N3_mu15_nratio2_i2500_3.csv", header = FALSE)

colnames(s_m) = c("iteration","firm","share")
colnames(s_mi) = c("iteration","firm","share")
colnames(s_miga) = c("iteration","firm","share")
colnames(cf_mi_1) = c("iteration","cf", "used")
colnames(cf_mi_2) = c("iteration","cf", "used")
colnames(cf_mi_3) = c("iteration","cf", "used")

color3 = c("#f35e5a", "#17b12b", "#5086ff")

# MAXCOV
s_m$firm = rep(c(2, 3, 1), each=250) # match colors with MAXCOV-INDUCTOR
fig5s_m = ggplot(s_m, aes(y = share, x =iteration, colour = factor(firm)))  + 
  theme_minimal() +
  theme(legend.position = "none", 
        legend.box = "horizontal", 
        legend.key = element_rect(fill = NA, colour = NA) )
fig5s_m + geom_line() +
  scale_x_continuous(limits = c(0, 250), expand = c(0, 5), minor_breaks = NULL) + 
  scale_y_continuous(limits = c(0, 0.7), expand = c(0, 0), labels=percent) +
  labs(y = "Share of market", x = "Iteration", colour = NULL)
ggsave("fig5s_m.pdf", width = 24, height = 5.3, units = "cm")

# MAXCOV-INDUCTOR
fig5s_mi = ggplot(s_mi, aes(y = share, x = iteration, colour = factor(firm)))  + 
  theme_minimal() +
  theme(legend.position = "none", 
        legend.box = "horizontal", 
        legend.key = element_rect(fill = NA, colour = NA) )
fig5s_mi + geom_line() +
  scale_x_continuous(limits = c(0, 1000), expand = c(0, 20), minor_breaks = NULL) + 
  scale_y_continuous(limits = c(0, 0.7), expand = c(0, 0), labels=percent) +
  labs(y = "Share of market", x = "Iteration", colour = NULL)
ggsave("fig5s_mi.pdf", width = 24, height = 5.3, units = "cm")

fig5cf_mi_1 = ggplot(cf_mi_1, aes(x = iteration, y = cf, fill = used)) + 
  theme_minimal() +
  theme(legend.position = "none", 
        legend.box = "horizontal", 
        legend.key = element_rect(fill = NA, colour = NA),
        panel.grid.major.y = element_blank(),
        plot.margin=unit(c(0.2,0.4,-0.4,0), "cm"),
        axis.line = element_line(colour = "grey90") ) +
  geom_raster() +
  scale_fill_gradientn(colours=rep(color3[1],3)) +
  scale_x_continuous(limits = c(0, 1000), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(1, 100), expand = c(0, 0), minor_breaks = NULL) +
  labs(y = " ", x = " ", colour = NULL)
fig5cf_mi_1

fig5cf_mi_2 = ggplot(cf_mi_2, aes(x = iteration, y = cf, fill = used)) + 
  theme_minimal() +
  theme(legend.position = "none", 
        legend.box = "horizontal", 
        legend.key = element_rect(fill = NA, colour = NA),
        panel.grid.major.y = element_blank(),
        plot.margin=unit(c(0.2,0.4,-0.4,0), "cm"),
        axis.line = element_line(colour = "grey90") ) +
  geom_raster() +
  scale_fill_gradientn(colours=rep(color3[2],3)) +
  scale_x_continuous(limits = c(0, 1000), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(1, 100), expand = c(0, 0), minor_breaks = NULL) +
  labs(y = "Condition/forecast rule", x = " ", colour = NULL)
fig5cf_mi_2

fig5cf_mi_3 = ggplot(cf_mi_3, aes(x = iteration, y = cf, fill = used)) + 
  theme_minimal() +
  theme(legend.position = "none", 
        legend.box = "horizontal", 
        legend.key = element_rect(fill = NA, colour = NA),
        panel.grid.major.y = element_blank(),
        plot.margin=unit(c(0.2,0.4,0,0), "cm"),
        axis.line = element_line(colour = "grey90") ) +
  geom_raster() +
  scale_fill_gradientn(colours=rep(color3[3],3)) +
  scale_x_continuous(limits = c(0, 1000), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(1, 100), expand = c(0, 0), minor_breaks = NULL) +
  labs(y = " ", x = "Iteration", colour = NULL)
fig5cf_mi_3

cf_lay = t(matrix(seq(1,3), 1, 3, byrow = T))
grid.arrange(fig5cf_mi_1, fig5cf_mi_2, fig5cf_mi_3, ncol=3, layout_matrix = cf_lay)

fig5cf_mi = arrangeGrob(fig5cf_mi_1, fig5cf_mi_2, fig5cf_mi_3, ncol=3, layout_matrix = cf_lay)
ggsave(file="fig5cf_mi.pdf", fig5cf_mi, width = 21, height = 12, units = "cm") #saves g

# MAXCOV-INDUCTOR-GA
#s_miga$firm = rep(c(1, 3, 2), each=2500) # match colors with MAXCOV-INDUCTOR
fig5s_miga = ggplot(s_miga, aes(y = share, x = iteration, colour = factor(firm)))  + 
  theme_minimal() +
  theme(legend.position = "none", 
        legend.box = "horizontal", 
        legend.key = element_rect(fill = NA, colour = NA) )
fig5s_miga + geom_line() +
  scale_x_continuous(limits = c(0, 2500), expand = c(0, 20), minor_breaks = NULL) + 
  scale_y_continuous(limits = c(0, 0.7), expand = c(0, 0), labels=percent) +
  labs(y = "Share of market", x = "Iteration", colour = NULL)
ggsave("fig5s_miga.pdf", width = 24, height = 5.3, units = "cm")

#library(R.matlab)
#path = system.file("data", package="R.matlab")
#pathname = file.path(path, "three_maxcov-inductor-ga_xy_N3_mu15_nratio2_i2500.mat")
#xy_miga_mat = readMat("data/three_maxcov-inductor-ga_xy_N3_mu15_nratio2_i2500.mat")
#library (plyr)
#xy_miga = t(ldply(xy_miga_mat, data.frame))
#xy_miga = xy_miga[-1,]

xy_miga = read.csv("data/three_maxcov-inductor-ga_xy_N3_mu15_nratio2_i2500.csv", header = FALSE)
colnames(xy_miga) = c("iteration","firm","x", "y")



# 6. Accuracy (five firms) -------------------------------

a_moves_miga = read.csv("data/xy_accuracy_maxcov-inductor-ga_N5_mu13_nratio15_psi50_i50.csv")
a_miga = read.csv("data/accuracy_maxcov-inductor-ga_N5_mu13_nratio15_psi50_i50.csv")
a_mi_same = read.csv("data/accuracy_maxcov-inductor_N5_mu13_nratio15_psi50_i50_same.csv")
a_miga$forecasterror = sqrt(a_miga$accuracy)
a_mi_same$forecasterror = sqrt(a_miga$accuracy)
a_miga_500 = a_miga[ which(a_miga$iteration <= 500), ]
a_mi_same_500 = a_mi_same[ which(a_mi_same$iteration <= 500), ]


fig6a_miga = ggplot(a_miga_500, aes(y = accuracy, x = iteration, colour = factor(firm)))  + 
  theme_minimal() +
  theme(legend.position = "none", 
        legend.box = "horizontal", 
        legend.key = element_rect(fill = NA, colour = NA) )
fig6a_miga + geom_line() +
  annotate("rect", xmin=271, xmax=320, ymin=0, ymax=0.004, alpha=0.1, fill="black")  + 
  scale_x_continuous(limits = c(1, 500), expand = c(0, 5), breaks = seq(0, 500, 50), minor_breaks = NULL) + 
  scale_y_continuous(limits = c(0, 0.004), expand = c(0, 0), breaks = seq(0, 0.004, 0.002)) +
  labs(y = "Average forecast error variance", x = "Iteration", colour = NULL)
ggsave(file="fig6a_miga.pdf", width = 12, height = 14, units = "cm")

fig6a_mi_same = ggplot(a_mi_same_500, aes(y = accuracy, x = iteration, colour = factor(firm)))  + 
  theme_minimal() +
  theme(legend.position = "none", 
        legend.box = "horizontal", 
        legend.key = element_rect(fill = NA, colour = NA) )
fig6a_mi_same + geom_line() +
  scale_x_continuous(limits = c(1, 500), expand = c(0, 5), breaks = seq(0, 500, 50), minor_breaks = NULL) + 
  scale_y_continuous(limits = c(0, 0.008), expand = c(0, 0)) +
  labs(y = "Average forecast error variance", x = "Iteration", colour = NULL)
ggsave(file="fig6a_mi_same.pdf", width = 12, height = 14, units = "cm")


move6 = function(df, start, end, tick) {
  # Number of firms
  N = 5
  # Hide or show the tick labes on the x and y-axis.
  if (tick) ticklabelcolor = "black" else ticklabelcolor = "white"
  
  #endpoints = subset(df[((start-1)*N+1):(end*N),], iteration == end)
  #voronoi = deldir(endpoints$x, endpoints$y, rw = c(-3, 3, -2, 2))
  
  figure = ggplot(df[((start-1)*N+1):(end*N),], aes(y = y, x = x, colour = factor(firm))) + 
    theme_minimal() +
    theme(legend.position = "top", 
          legend.box = "horizontal", 
          legend.key = element_rect(fill = NA, colour = NA),
          axis.text.x = element_text(colour = ticklabelcolor), 
          axis.text.y = element_text(colour = ticklabelcolor),
          plot.title = element_text(size=10, hjust = 0),
          plot.margin = unit(c(0.2, 0, -0.4, 0), "lines") ) + 
    geom_line(alpha = 0.5) +
    geom_point(data = subset(df[((start-1)*N+1):(end*N),], iteration == end), size = 2) +
    #geom_segment(data = voronoi$dirsgs, aes(x = x1, y = y1, xend = x2, yend = y2), size = 0.3, linetype = 1, color= "#000000") + 
    scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1.5), minor_breaks = NULL, name="") + 
    scale_y_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2, 2, by = 1), minor_breaks = NULL, name="") +
    scale_colour_discrete(guide = FALSE) +
    geom_point(aes(x = -0.26, y = 0), shape = 3, colour = "black")+
    labs(title = paste("Iteration", sprintf("%1.0f", start), "to", sprintf("%1.0f", end), sep =" "))
  
  return(figure)
}

fig6m_miga1 = move6(a_moves_miga, 50, 270, TRUE)
fig6m_miga1
ggsave(file="fig6m_miga1.pdf", fig6m_miga1, width = 7.5, height = 6.6, units = "cm")
fig6m_miga2 = move6(a_moves_miga, 271, 320, TRUE)
fig6m_miga2
ggsave(file="fig6m_miga2.pdf", fig6m_miga2, width = 7.5, height = 6.6, units = "cm")
fig6m_miga3 = move6(a_moves_miga, 321, 2500, TRUE)
fig6m_miga3
ggsave(file="fig6m_miga3.pdf", fig6m_miga3, width = 7.5, height = 6.6, units = "cm")

#width = 10, height = 6.6