# Script for determining the differences between glucose consumption rates
# over various glucose and oxygen levels over a culture period of 10 days.

# Author <- Lukas Jaworski, l.jaworski@umiami.edu
# Institution <- University of Miami
# PI <- Alicia Jackson, a.jackson2@miami.edu

# Copyright University of Miami 2017
###############################################################################
## Library imports
###############################################################################

library(ggplot2)
library(ggsignif)
library(cowplot)

###############################################################################
## Functions
###############################################################################
## Summarizes data.This function was copied from cookbook-r.com
## Gives count, mean, standard deviation, standard error of the mean, and 
## confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be 
##               summariezed
##   groupvars: a vector containing names of columns that contain grouping 
##              variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval(default is 95%)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

facet_signif <- function(plot, mapping){
  #add in significance bars on a faceted plot.
  #plot: a ggplot2 object
  #mapping: which of the facets each of the significance bars goes to
  #         it is a list containing a vector for each facet with the
  #         position of the annotations in ggsignif
  
  # example: 
  # my.plot <- ggplot(test.data, aes(time, dependant.variable))+
  #   geom_boxplot()+
  #   facet_grid(. ~ grp) +
  #   ## it takes these significance annotations ##
  #   geom_signif(annotations = c('*','*', '**', '*'),
  #               y_position = c(-4.5, -1.5, -2, -1), 
  #               ## ie y_position[1] = -4.5, y_position[2] = -1.5, etc.
  #               xmin = c(1, 1, 1, 1),
  #               xmax = c(3, 3, 2, 3), 
  #               vjust = 0.5) +
  #
  ## the significance mapping shows where each annotation will go ##
  #
  # my.plot <- facet_signif(my.plot, list(c(1),c(2),c(3,4)))
  #
  ## The plot has three facets, the first on recieves annotation 1,
  ## the second facet receives the annotation 2 and the thrid facet
  ## receives annotations 3 and 4. 
  
  
  #pull out the plot data
  plot.info <- ggplot_build(plot)
  signif.drawing <- plot.info$data[[2]]
  new.drawing <- data.frame()
  new.drawing <- merge(signif.drawing, new.drawing)
  # took too long to comment this so there's some black magic below
  # it adds to each facet the number position of the annotations added 
  # for ggsignif
  for (i in 1:length(mapping)){
    pv <- signif.drawing[((signif.drawing$group %in% mapping[[i]]) &
                         (signif.drawing$PANEL == i)), ]
    
    new.drawing <- merge(new.drawing, 
                         pv, 
                         all.x = T, 
                         all.y = T)
  }
  plot.info$data[[2]] <- new.drawing
  return(ggplot_gtable(plot.info))
}


###############################################################################
## Script
###############################################################################
# Retrieving and arranging the data

# Set the working directory to where the data is kept, 
# make it workstation agnostic
if (Sys.info()[['sysname']]== 'Linux'){
  working.directory <- paste0('/home/',Sys.info()[['user']],
                              '/Dropbox/Lab/Experiments/GCR/',
                              'a2-gcr(gluc.ox~days)-NP/Script')
} else {
  working.directory <- paste0('C:/Users/',Sys.info()[['user']],
                              '/Dropbox/Lab/Experiments/GCR/',
                              'a2-gcr(gluc.ox~days)-NP/Script')
}
working.directory <- paste0('P:/Dropbox/Lab/Experiments/GCR/',
                              'a2-gcr(gluc.ox~days)-NP/Script')
setwd(working.directory)

# for the glucose concentration groups
gluc.gcr.vals <- read.csv('a2-GCR-gluc.csv')

# for the oxygen tension groups
ox.gcr.vals <- read.csv('a2-GCR-ox.csv')

# assigning names to the data columns: 
# Sample-Sample Name; GCR-glucose consumption rate

# for the glucose concentration groups
names(gluc.gcr.vals)[1] <- 'Sample'
names(gluc.gcr.vals)[2] <- 'GCR'

# for the oxygen tension groups
names(ox.gcr.vals)[1] <- 'Sample'
names(ox.gcr.vals)[2] <- 'GCR'

# changing the sample names to characters instead of factors
gluc.gcr.vals$Sample <- as.character(gluc.gcr.vals$Sample)
ox.gcr.vals$Sample   <- as.character(ox.gcr.vals$Sample)

# adding a column to identify from which pig each sample came from
# this will be the 'blocking' variable in the ANOVA

# for the glucose concentration groups
gluc.gcr.vals$pig <- NA
gluc.gcr.vals$pig[grep('P7', gluc.gcr.vals$Sample)]  <- 'P7'
gluc.gcr.vals$pig[grep('P8', gluc.gcr.vals$Sample)]  <- 'P8'
gluc.gcr.vals$pig[grep('P9', gluc.gcr.vals$Sample)]  <- 'P9'
gluc.gcr.vals$pig[grep('P10', gluc.gcr.vals$Sample)] <- 'P10'
gluc.gcr.vals$pig <- as.factor(gluc.gcr.vals$pig)

# for the oxygen tension groups
ox.gcr.vals$pig <- NA
ox.gcr.vals$pig[grep('P11', ox.gcr.vals$Sample)] <- 'P11'
ox.gcr.vals$pig[grep('P12', ox.gcr.vals$Sample)] <- 'P12'
ox.gcr.vals$pig[grep('P13', ox.gcr.vals$Sample)] <- 'P13'
ox.gcr.vals$pig[grep('P14', ox.gcr.vals$Sample)] <- 'P14'
ox.gcr.vals$pig[grep('P15', ox.gcr.vals$Sample)] <- 'P15'
ox.gcr.vals$pig <- as.factor(ox.gcr.vals$pig)

# adding a column identifying the first independant variable:
# glucose concentration or oxygen tension

# for the glucose concentration groups
gluc.gcr.vals$glucGroup <- NA
gluc.gcr.vals$glucGroup[grep(' 1s', gluc.gcr.vals$Sample)] <- '1mM'
gluc.gcr.vals$glucGroup[grep(' 2s', gluc.gcr.vals$Sample)] <- '2.5mM'
gluc.gcr.vals$glucGroup[grep(' 5s', gluc.gcr.vals$Sample)] <- '5mM'
gluc.gcr.vals$glucGroup[grep(' 1-', gluc.gcr.vals$Sample)] <- '1mM'
gluc.gcr.vals$glucGroup[grep(' 2-', gluc.gcr.vals$Sample)] <- '2.5mM'
gluc.gcr.vals$glucGroup[grep(' 5-', gluc.gcr.vals$Sample)] <- '5mM'
gluc.gcr.vals$glucGroup <- factor(gluc.gcr.vals$glucGroup, 
                                  levels=c("1mM", "2.5mM", "5mM"))

# for the oxygen tension groups
ox.gcr.vals$O2group <- NA
ox.gcr.vals$O2group[grep(' 1%-', ox.gcr.vals$Sample)]  <- '1%'
ox.gcr.vals$O2group[grep(' 5%-', ox.gcr.vals$Sample)]  <- '5%'
ox.gcr.vals$O2group[grep(' 21%-', ox.gcr.vals$Sample)] <- '21%'
ox.gcr.vals$O2group <- factor(ox.gcr.vals$O2group,
                              levels = c('1%', '5%', '21%'))

# adding a column identifying the second independant variable: time

# for the glucose concentration groups
gluc.gcr.vals$Day <- NA
gluc.gcr.vals$Day[grep('D1 ', gluc.gcr.vals$Sample)]  <- 'Day 1'
gluc.gcr.vals$Day[grep('D5 ', gluc.gcr.vals$Sample)]  <- 'Day 5'
gluc.gcr.vals$Day[grep('D10 ', gluc.gcr.vals$Sample)] <- 'Day 10'
gluc.gcr.vals$Day <- factor(gluc.gcr.vals$Day,
                            levels = c('Day 1', 'Day 5', 'Day 10'))

# for the oxygen tension groups
ox.gcr.vals$Day <- NA
ox.gcr.vals$Day[grep('D1 ', ox.gcr.vals$Sample)]  <- 'Day 1'
ox.gcr.vals$Day[grep('D5 ', ox.gcr.vals$Sample)]  <- 'Day 5'
ox.gcr.vals$Day[grep('D10 ', ox.gcr.vals$Sample)] <- 'Day 10'
ox.gcr.vals$Day <- factor(ox.gcr.vals$Day, 
                          levels = c('Day 1', 'Day 5', 'Day 10'))


###############################################################################
# running the ANOVA and the post hoc tests 

# for the glucose concentration groups
gluc.GCR.anova <- aov(GCR ~ glucGroup * Day + pig,
                      data = gluc.gcr.vals)
summary.aov(gluc.GCR.anova)
TukeyHSD(gluc.GCR.anova)
comparisons.gluc <- TukeyHSD(gluc.GCR.anova)$'glucGroup:Day'
comparisons.gluc[comparisons.gluc[, 'p adj'] <= 0.05 &
                 comparisons.gluc[, 'p adj'] >  0.01, ]
comparisons.gluc[comparisons.gluc[, 'p adj'] <= 0.01 &
                 comparisons.gluc[, 'p adj'] >  0.001, ]
comparisons.gluc[comparisons.gluc[, 'p adj'] <= 0.001, ]

# for the oxygen tension groups
ox.GCR.anova <- aov(GCR ~ O2group * Day + pig,
                    data = ox.gcr.vals)
summary.aov(ox.GCR.anova)
TukeyHSD(ox.GCR.anova)
comparisons.ox <- TukeyHSD(ox.GCR.anova)$'O2group:Day'
comparisons.ox[comparisons.ox[, 'p adj'] <= 0.05 &
               comparisons.ox[, 'p adj'] >  0.01, ]
comparisons.ox[comparisons.ox[, 'p adj'] <= 0.01 &
               comparisons.ox[, 'p adj'] >  0.001, ]
comparisons.ox[comparisons.ox[, 'p adj'] <= 0.001, ]
drop1(ox.GCR.anova, ~., test = "F")

# Getting summray statisics to be used in plotting

# for the glucose concentration groups
gluc.summary <- summarySE(gluc.gcr.vals, 
                          measurevar = "GCR", 
                          groupvars = c("Day", "glucGroup"))

# for the oxygen tension groups
ox.summary <- summarySE(ox.gcr.vals,
                        measurevar = "GCR", 
                        groupvars = c("Day", "O2group"))

###############################################################################
# Graphing the results
theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(), 
             panel.background = element_blank(),
             axis.line = element_line(colour = "black"),
             text = element_text(size = 22))

# for the glucose concentration groups
# plot comparing differences between glucose concentrations within each day
gluc.plot.color <- ggplot(gluc.summary, aes(Day, GCR, fill = glucGroup)) + 
  geom_col(position = position_dodge(),
           color = 'black') + 
  geom_signif(annotations = c('*', '**', '***', '***', '***', '*', '**'),
              y_position = c(473, 527, 488, 560, 508, 590, 425),
              xmin = c(0.7, 0.7, 1.7, 1, 2, 1.3, 2.3),
              xmax = c(1.7, 2.7, 2.7, 3, 3, 3.3, 3.3),
              vjust = .5) +
  geom_errorbar(aes(ymin = GCR, ymax = GCR + se),
                position = position_dodge()) +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(bquote('GCR (nmol/10'^6*' cells/hr)')) +
  scale_fill_manual(name = 'Glucose\nConcentration', values = c('#005030', '#F47321', '#FFFFFF'))

plot(gluc.plot.color)

gluc.plot.days <- ggplot(gluc.summary, aes(glucGroup, GCR)) +
  geom_col() + 
  facet_grid(. ~ Day) +
  geom_errorbar(aes(ymin = GCR, ymax = GCR + se)) +
  xlab('Glucose Concentration (mM)') +
  ylab(bquote('GCR (nmol/10'^6*' cells/hr)')) 

# plot comparing differences between days within each glucose concentration
gluc.plot.groups <- ggplot(gluc.summary, aes(Day, GCR)) +
  geom_col() + 
  facet_grid(. ~ glucGroup) +
  geom_signif(annotations = c('*', '**', '***', '***', '***', '*', '**'),
              y_position = c(490, 540, 590, 500, 450, 500, 450),
              xmin = c(1, 1, 2, 1, 2, 1, 2),
              xmax = c(2, 3, 3, 3, 3, 3, 3),
              vjust = .5) +
  geom_errorbar(aes(ymin = GCR, ymax = GCR + se)) +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(bquote('GCR (nmol/10'^6*' cells/hr)'))
gluc.plot.groups <- facet_signif(gluc.plot.groups, 
                                 list(c(1, 2, 3), c(4, 5), c(6, 7)))


# plot both sets of comparisons to highlight either glucose or time differences
plot_grid(gluc.plot.days, 
          gluc.plot.groups, 
          labels = c('A', 'B'), 
          ncol = 1, 
          nrow = 2)

# for the oxygen tension groups
# plot comparing differences between O2 tensions within each day
ox.plot.days <- ggplot(ox.summary, aes(O2group, GCR)) +
  geom_col() + 
  facet_grid(. ~ Day) +
  geom_signif(annotations = c('**', '*', '**'),
              y_position = c(310, 320, 270),
              xmin = c(1, 1, 2),
              xmax = c(3, 3, 3), 
              vjust = 0.5) +
  geom_errorbar(aes(ymin = GCR, ymax = GCR + se)) +
  xlab('Oxygen Tension') +
  ylab(bquote('GCR (nmol/10'^6*' cells/hr)')) +
  scale_y_continuous(limits = c(0, 350)) 

ox.plot.days <- facet_signif(ox.plot.days, list(c(), c(1), c(2, 3)))

# plot comparing differences between days within each O2 tension
ox.plot.groups <- ggplot(ox.summary, aes(Day, GCR)) +
  geom_col() + 
  facet_grid(. ~ O2group) +
  geom_signif(annotations = c('*', '***'),
              y_position = c(330, 310),
              xmin = c(1, 1),
              xmax = c(3, 3), 
              vjust = 0.5) +
  geom_errorbar(aes(ymin = GCR, ymax = GCR + se)) +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(bquote('GCR (nmol/10'^6*' cells/hr)')) +
  scale_y_continuous(limits = c(0, 350)) 

ox.plot.groups <- facet_signif(ox.plot.groups, list(c(), c(1), c(2)))

# plot both sets of comparisons to highlight ethier oxygen or time differences
plot_grid(ox.plot.days, 
          ox.plot.groups, 
          labels = c('a', 'b'), 
          ncol = 1, 
          nrow = 2)

