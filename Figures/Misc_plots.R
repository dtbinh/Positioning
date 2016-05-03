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
ggsave("figb21a.pdf", width = 21, height = 16, units = "cm")


# 2.2 Stocastic time-homogenous Markov chain - All-hunter
burn_stc = read.csv("data/All-hunter_N12_mu0_nratio1_i500.csv")
burn_stc_2nd_mean = mean(burn_stc[c(251:500),3])
burn_stc_2nd_std = sd(burn_stc[c(251:500),3])

burn_stc$mean_eccentricity > burn_stc_2nd_mean+burn_stc_2nd_std

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
ggsave("figb22a.pdf", width = 21, height = 16, units = "cm")


# 2.3 Stocastic time-homogenous Markov chain (time average not representative) - All-hunter
# data_mean_eccentricity_20160420_180651.csv
# filter-time tool
# figb23a.pdf, figb23b.pdf



# 3. FIRM MOVEMENT --------------------------------------------------------

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
    geom_point(data = subset(moves[((start-1)*N+1):(end*N),], iteration == end), size = 2) +
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
ggsave(file="figm.pdf", figm, width = 16, height = 21, units = "cm") #saves g



# 3. Median spline - Eccentricity vs. share -------------------------------

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
ggsave("fig3ms.pdf", width = 21, height = 12, units = "cm")

#geom_segment(aes(x = 0.4416208, y = 0, xend = 0.4416208, yend = 0.216), color = "red") +
#geom_quantile(quantiles=0.5, formula = y ~ qss(x, lambda=5), method = "rqss", size = 1) +
#geom_point(size = 0.1) + 

