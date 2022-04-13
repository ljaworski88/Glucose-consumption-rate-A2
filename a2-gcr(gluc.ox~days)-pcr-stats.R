# Script for determining the differences between genes over various glucose and
# oxygen levels over a culture period of 10 days.

# Author <- Lukas Jaworski, l.jaworski@umiami.edu
# Institution <- University of Miami
# PI <- Alicia Jackson, a.jackson2@miami.edu

# Copyright University of Miami 2017
###############################################################################
## Library imports
###############################################################################

library(plyr)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(ReadqPCR)
library(NormqPCR)

###############################################################################
## Functions
###############################################################################

aov.results <- function(genes, data, formula){
  # A function to get the statistically significant results without any results,
  # which are not significant
  
  
  anova.results <- list()
  for (x in genes){
    anova.results[[x]] <- aov(as.formula(gsub('{}', x, formula, fixed = T)),
                              data = data)
  }
  # ANOVA results and perform post hoc tests
  aov.summary<-list()
  for (x in genes){
    temp.summary <- summary.aov(anova.results[[x]])[[1]]
    aov.summary[[x]][[1]] <- temp.summary[temp.summary['Pr(>F)'] <= 0.05 &
                                            !is.na(temp.summary['Pr(>F)']),
                                          'Pr(>F)',
                                          drop = F]
    if (nrow(aov.summary[[x]][[1]])>0){
      temp.tukey <- TukeyHSD(anova.results[[x]])
      aov.summary[[x]][[2]] <- lapply(temp.tukey,
                                      function(x) x[x[, 'p adj'] <= 0.05, 'p adj'])
    }
  }
  return(aov.summary)
}

facet_signif <- function(plot, mapping){
  #add in significance bars on a faceted plot.
  #plot: a ggplot2 object
  #mapping: which of the facets each of the significance bars goes to
  #         it is a list containing a vector for each facet with the
  #         position of the annotations in ggsignif
  
  # example: 
  # my.plot <- ggplot(test.data, aes(time, dependant.variable)) +
  #   geom_boxplot() +
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
  ## The plot has three facets, the first on recieves the first annotation,
  ## the second facet receives the second annotation and the thrid facet
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
    new.drawing <- merge(new.drawing, pv, all.x = T, all.y = T)
  }
  plot.info$data[[2]] <- new.drawing
  return(ggplot_gtable(plot.info))
}

###############################################################################
## Script
###############################################################################
# Retrieving and arranging the data

# Set the working directory [workstation agnostic]
if (Sys.info()[['sysname']] == 'Linux'){
  working.directory <- paste0('/home/', Sys.info()[['user']], 
                              '/Dropbox/Lab/Experiments/PCR/',
                              'a2-gcr(gluc.ox~days)-NP/cleaned_data')
  
} else {
  working.directory <- paste0('C:/Users/', Sys.info()[['user']],
                              '/Dropbox/Lab/Experiments/PCR/',
                              'a2-gcr(gluc.ox~days)-NP/cleaned_data')
} 
  working.directory <- paste0('P:/Dropbox/Lab/Experiments/PCR/',
                              'a2-gcr(gluc.ox~days)-NP/cleaned_data')
setwd(working.directory)

# Load the Phenotype data from the text file
pheno.data <- read.csv("PCR_phenodata.txt", sep = "\t")

# Order the phenotype data according to sample name so it aligns with the 
# assayData
# (This should do a hash lookup so the sample names don't have to be ordered)
pheno.data <- pheno.data[with(pheno.data, order(Sample)), ]

#change day number into a factor
pheno.data$Day <- factor(pheno.data$Day)

pigs <- c(rep('8', 9), rep('9', 23), rep('10', 34), rep('11', 15), rep('12', 17), 
          rep('14', 19), rep('15', 51), rep('7', 11), rep('8', 9))
pheno.data$Pig <-factor(pigs)
# Copy the sample names into the rownames 
# (needs to be done otherwise read.qPCR won't be able to match the data)
rownames(pheno.data) <- pheno.data$Sample

# Metadata for our annotated dataframe
meta <- data.frame(c('Days cultured', 'exp group', 'sample number', 'pig'))

#Gotta match that eSet format so boom
names(meta) <- 'labelDescription'

#Finally making the annotated dataframe
annotated.pheno <- AnnotatedDataFrame(pheno.data,
                                      meta, 
                                      c('sampleNames','sampleColumns'))

# Getting the Cqs from a text file and adding the phenotype data 
Cq.vals <- read.qPCR('PCR_cts.txt', annotated.pheno)

# Subsetting the data

# Glucose concentration experiment
# The control genes were verified using another script
hk.gene.gluc <- c('YWHAZ', 
                  'RPL4',
                  'HPRT_I', 
                  'GAPDH', 
                  'X18s',
                  'HMBS', 
                  'ACTB', 
                  'TBP')

dCq.vals <- deltaCq(qPCRBatch = Cq.vals, 
                    hkgs = hk.gene.gluc, 
                    calc="geom")

dCq.vals.gluc <- dCq.vals[, dCq.vals@phenoData@data$Exp_Grp == '1mM'   | 
                            dCq.vals@phenoData@data$Exp_Grp == '2.5mM' |
                            dCq.vals@phenoData@data$Exp_Grp == '5mM']

# Oxygen tension experiment
# The control genes were verified using another script
hk.gene.o2 <- c('YWHAZ', 
                'RPL4',
                'HPRT_I', 
                'GAPDH', 
                'X18s', 
                'HMBS')

dCq.vals <- deltaCq(qPCRBatch = Cq.vals, 
                    hkgs = hk.gene.o2, 
                    calc = "geom")

dCq.vals.O2 <- dCq.vals[, dCq.vals@phenoData@data$Exp_Grp == '21%' |
                          dCq.vals@phenoData@data$Exp_Grp == '5%'  |
                          dCq.vals@phenoData@data$Exp_Grp == '1%']

###############################################################################
## Run ANOVAs on all of the genes and display the results

# In the O2 group 

# unpack qPCRbatch values into a dataframe for statistical analysis
dCq.O2 <- data.frame(t(exprs(dCq.vals.O2)))
dCq.O2$day <- revalue(factor(dCq.vals.O2@phenoData@data$Day), 
                      c("1"  = "Day 1",
                        "5"  = "Day 5", 
                        "10" = "Day 10"))

dCq.O2$grp <- factor(dCq.vals.O2@phenoData@data$Exp_Grp, 
                     levels = c('1%', '5%', '21%'))

dCq.O2$pig <- factor(dCq.vals.O2@phenoData@data$Pig)

# The list of genes to test
genes.to.test <- c('Agg',
                   'Col_IA1',
                   'Col_IIA1',
                   'MMP_1',
                   'MMP_3', 
                   'MMP_13',
                   'MT_1',
                   'T', 
                   'TIMP_1',
                   'TIMP_2',
                   'TIMP_3')

# run the anovas
anova.summary.o2 <- aov.results(genes.to.test,
                                dCq.O2,
                                '{} ~ day * grp + pig ')

## In the glucose group

# unpack qPCRbatch values into a dataframe for statistical analysi
dCq.gluc <- data.frame(t(exprs(dCq.vals.gluc)))
dCq.gluc$day <- revalue(factor(dCq.vals.gluc@phenoData@data$Day), 
                      c("1"  = "Day 1",
                        "5"  = "Day 5", 
                        "10" = "Day 10"))

dCq.gluc$grp <- dCq.vals.gluc@phenoData@data$Exp_Grp

dCq.gluc$pig <- dCq.vals.gluc@phenoData@data$Pig

# run the ANOVAs
aov.summary.gluc <- aov.results(genes.to.test, 
                                dCq.gluc,
                                '{} ~ day * grp + pig')

## wilcoxon rank sum and t.test comparisons between all pairs of variables both
## were done to do a sensitivity analysis and there was a question to the 
## normality of gene expression data. After furthur research Cq/Ct vbalues are
## generally normally distrubuted and the ANOVA with tukey post hoc is ok,
## as far as I am aware. This code will remain here incase I am wrong. 
## the results of each analysis will need to have their p values adjusted for
## multiple comparisons though.


# days <- levels(dCq.O2$day)
# grps <- c('1%','5%','21%')
# pvals.wilcox <- data.frame(row.names = genes.to.test)
# for (z in 1:(length(days))){
#   for (x in 1:(length(grps)-1)){
#     for (y in (x+1):length(grps)){
#       for(i in row.names(pvals.wilcox)) {
#         case <- dCq.O2[dCq.O2$day==days[z] & dCq.O2$grp==grps[x], i]
#         control <- dCq.O2[dCq.O2$day==days[z] & dCq.O2$grp==grps[y], i]
#         if(sum(is.na(case)) != 0 | sum(is.na(control)) != 0) {
#           # gets rid of NaN values (mostly the housekeeping gene)
#           pvals.wilcox[i, paste(days[z],
#                          grps[x],
#                          grps[y],
#                          sep = '_')] <- NA
#         }else{
#           pvals.wilcox[i, paste(days[z],
#                          grps[x],
#                          grps[y],
#                          sep = '_')] <- wilcox.test(case,control)$p.value
#         }
#       }
#     }
#   }
# }
# for (z in 1:(length(grps))){
#   for (x in 1:(length(days)-1)){
#     for (y in (x+1):length(days)){
#       for(i in row.names(pvals.wilcox)) {
#         case <- dCq.O2[dCq.O2$grp==grps[z] & dCq.O2$day==days[x], i]
#         control <- dCq.O2[dCq.O2$grp==grps[z] & dCq.O2$day==days[y], i]
#         if(sum(is.na(case)) != 0 | sum(is.na(control)) != 0) {
#           # gets rid of NaN values (mostly the housekeeping gene)
#           pvals.wilcox[i, paste(grps[z],
#                          days[x],
#                          days[y],
#                          sep = '_')] <- NA
#         } else {
#           pvals.wilcox[i, paste(grps[z],
#                          days[x],
#                          days[y],
#                          sep = '_')] <- wilcox.test(case,control)$p.value
#         }
#       }
#     }
#   }
# }
# 
# pvals.t.test <- data.frame(row.names = genes.to.test)
# for (z in 1:(length(days))){
#   for (x in 1:(length(grps)-1)){
#     for (y in (x+1):length(grps)){
#       for(i in row.names(pvals.t.test)) {
#         case <- dCq.O2[dCq.O2$day==days[z] & dCq.O2$grp==grps[x], i]
#         control <- dCq.O2[dCq.O2$day==days[z] & dCq.O2$grp==grps[y], i]
#         if(sum(is.na(case)) != 0 | sum(is.na(control)) != 0) {
#           # gets rid of NaN values (mostly the housekeeping gene)
#           pvals.t.test[i, paste(days[z],
#                                 grps[x],
#                                 grps[y],
#                                 sep = '_')] <- NA
#         } else {
#           pvals.t.test[i, paste(days[z],
#                                 grps[x],
#                                 grps[y],
#                                 sep = '_')] <- t.test(case,control)$p.value
#         }
#       }
#     }
#   }
# }
# 
# for (z in 1:(length(grps))){
#   for (x in 1:(length(days)-1)){
#     for (y in (x+1):length(days)){
#       for(i in row.names(pvals.t.test)) {
#         case <- dCq.O2[dCq.O2$grp==grps[z] & dCq.O2$day==days[x], i]
#         control <- dCq.O2[dCq.O2$grp==grps[z] & dCq.O2$day==days[y], i]
#         if(sum(is.na(case)) != 0 | sum(is.na(control)) != 0) {
#           # gets rid of NaN values (mostly the housekeeping gene)
#           pvals.t.test[i, paste(grps[z],
#                                 days[x],
#                                 days[y],
#                                 sep = '_')] <- NA
#         }else{
#           pvals.t.test[i, paste(grps[z],
#                                 days[x],
#                                 days[y],
#                                 sep = '_')] <- t.test(case,control)$p.value
#         }
#       }
#     }
#   }
# }

# finding correlations between Brachyury transcription factor and 
# other genes studied

dCq <- exprs(dCq.vals)

dCq.t <- t(dCq)
p.correlations <- NULL
genes.remaining <- setdiff(genes.to.test, 'Agg')
for (y in genes.to.test){
  genes.remaining <- setdiff(genes.remaining, y)
  for (x in genes.remaining){
    p.correlations <- rbind(p.correlations,
                            c(y, 
                              x, 
                              cor.test(dCq.t[, y],
                                       dCq.t[, x], 
                                       method = 'pearson')$estimate,
                              cor.test(dCq.t[, y],
                                       dCq.t[, x], 
                                       method = 'pearson')$p.value))
  }
}
adjPVals.correl <- p.adjust(as.numeric(p.correlations[, 4]),
                            method = "bonferroni")
p.correlations <- cbind(p.correlations, adjPVals.correl)
p.correlations.possible <- p.correlations[as.numeric(p.correlations[, 4])
                                          < (10^-6), ]

###############################################################################
# Graphing the results

## for the O2 groups

# plot brachyury expression
theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             axis.line = element_line(colour = "black"),
             text = element_text(size = 22))

plot.T.o2.grps  <-  ggplot(dCq.O2, aes(grp, -T)) +
  geom_boxplot() +
  scale_x_discrete(limits = c('1%', '5%', '21%')) +
  facet_grid(. ~ day) +
  geom_signif(annotations = c('*', '***'),
              y_position = c(-3.8, -3.3),
              xmin = c(1, 1),
              xmax = c(2, 3),
              vjust = 0.5) +
  ggtitle('Brachyury') +
  xlab('Oxygen Tension') +
  ylab(expression(paste('-', Delta, 'Cq'))) 
plot.T.o2.grps <- facet_signif(plot.T.o2.grps, list(c(), c(), c(1, 2)))

plot.T.o2.days  <-  ggplot(dCq.O2, aes(day, -T)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('**', '*'),
              y_position = c(-4.8, -4.3),
              xmin = c(2, 1),
              xmax = c(3, 3),
              vjust = 0.5) +
  ggtitle('Brachyury') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  xlab('Day') +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.T.o2.days <- facet_signif(plot.T.o2.days, list(c(), c(), c(1, 2)))

plot_grid(plot.T.o2.grps, 
          plot.T.o2.days,
          ncol = 1,
          nrow = 2,
          labels = 'AUTO')

# Plot the anabolic genes
plot.agg.o2.grps  <-  ggplot(dCq.O2, aes(grp, -Agg)) +
  geom_boxplot() +
  scale_x_discrete(limits = c('1%', '5%', '21%')) +
  facet_grid(. ~ day) +
  ggtitle('Aggrecan') +
  xlab('Oxygen Tension') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.agg.o2.days  <-  ggplot(dCq.O2, aes(day, -Agg)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  ggtitle('Aggrecan') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.col1.o2.grps  <-  ggplot(dCq.O2, aes(grp, -Col_IA1)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1%', '5%', '21%')) +
  ggtitle('Collagen I') +
  xlab('Oxygen Tension') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.col1.o2.days  <-  ggplot(dCq.O2, aes(day, -Col_IA1)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('*'),
              y_position = c(-7.5),
              xmin = c(1),
              xmax = c(2), 
              vjust = 0.5) +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ggtitle('Collagen I') +
  xlab('Day') +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.col1.o2.days <- facet_signif(plot.col1.o2.days, list(c(), c(1), c()))

plot.col2.o2.grps  <-  ggplot(dCq.O2, aes(grp, -Col_IIA1)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  geom_signif(annotations = c('**'),
              y_position = c(9.5),
              xmin = c(1),
              xmax = c(3), 
              vjust = 0.5) +
  scale_x_discrete(limits = c('1%', '5%', '21%')) +
  ggtitle('Collagen II') +
  xlab('Oxygen Tension') +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.col2.o2.grps <- facet_signif(plot.col2.o2.grps, list(c(), c(), c(1)))

plot.col2.o2.days  <-  ggplot(dCq.O2, aes(day, -Col_IIA1)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('*'),
              y_position = c(9.5),
              xmin = c(1),
              xmax = c(3), 
              vjust = 0.5) +
  ggtitle('Collagen II') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.col2.o2.days <- facet_signif(plot.col2.o2.days, list(c(1), c(), c()))

plot_grid(plot.agg.o2.grps, 
          plot.col1.o2.grps, 
          plot.col2.o2.grps, 
          plot.agg.o2.days, 
          plot.col1.o2.days, 
          plot.col2.o2.days,
          labels = 'AUTO', 
          ncol = 3, 
          nrow = 2)

# Plot the catabolic genes
plot.mmp1.o2.grps <-  ggplot(dCq.O2, aes(grp, -MMP_1)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1%', '5%', '21%')) +
  ggtitle('MMP 1') +
  xlab('Oxygen Tension') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.mmp1.o2.days  <-  ggplot(dCq.O2, aes(day, -MMP_1)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('**','***', '**', '*'),
              y_position = c(-1.5, 2.5, 3, 4),
              xmin = c(1, 1, 2, 1),
              xmax = c(3, 3, 3, 3), 
              vjust = 0.5) +
  ggtitle('MMP 1') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.mmp1.o2.days <- facet_signif(plot.mmp1.o2.days, list(c(1),c(2),c(3,4)))

plot.mmp3.o2.grps <-  ggplot(dCq.O2, aes(grp, -MMP_3)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1%', '5%', '21%')) +
  geom_signif(annotations = c('*'),
              y_position = c(-3),
              xmin = c(2),
              xmax = c(3), 
              vjust = 0.5) +
  ggtitle('MMP 3') +
  xlab('Oxygen Tension') +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.mmp3.o2.grps <- facet_signif(plot.mmp3.o2.grps, list(c(1), c(), c()))

plot.mmp3.o2.days  <-  ggplot(dCq.O2, aes(day, -MMP_3)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('***', '***'),
              y_position = c(.25, -.25),
              xmin = c(1, 2),
              xmax = c(3, 3), 
              vjust = 0.5) +
  ggtitle('MMP 3') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.mmp3.o2.days <- facet_signif(plot.mmp3.o2.days, list(c(), c(1, 2), c()))

plot.mmp13.o2.grps <-  ggplot(dCq.O2, aes(grp, -MMP_13)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1%', '5%', '21%')) +
  ggtitle('MMP 13') +
  xlab('Oxygen Tension') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.mmp13.o2.days  <-  ggplot(dCq.O2, aes(day, -MMP_13)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('*','*', '**', '*'),
              y_position = c(-4.5, -1.5, -2, -1),
              xmin = c(1, 1, 1, 1),
              xmax = c(3, 3, 2, 3), 
              vjust = 0.5) +
  ggtitle('MMP 13') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.mmp13.o2.days <- facet_signif(plot.mmp13.o2.days, list(c(1), c(2), c(3, 4)))

plot.mt1.o2.grps <-  ggplot(dCq.O2, aes(grp, -MT_1)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1%', '5%', '21%')) +
  ggtitle('MT 1') +
  xlab('Oxygen Tension') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.mt1.o2.days  <-  ggplot(dCq.O2, aes(day, -MT_1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  geom_boxplot() +
  facet_grid(. ~ grp) +
  ggtitle('MT 1') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))

plot_grid(plot.mmp1.o2.grps,
          plot.mmp3.o2.grps, 
          plot.mmp13.o2.grps, 
          plot.mt1.o2.grps, 
          plot.mmp1.o2.days,
          plot.mmp3.o2.days, 
          plot.mmp13.o2.days,
          plot.mt1.o2.days,
          labels = 'AUTO', 
          ncol = 4, 
          nrow = 2)

# Plot the catabolic inhibitor genes
plot.timp1.o2.grps <-  ggplot(dCq.O2, aes(grp, -TIMP_1)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1%', '5%', '21%')) +
  ggtitle('TIMP 1') +
  xlab('Oxygen Tension') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.timp1.o2.days  <-  ggplot(dCq.O2, aes(day, -TIMP_1)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  ggtitle('TIMP 1') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.timp2.o2.grps <-  ggplot(dCq.O2, aes(grp, -TIMP_2)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1%', '5%', '21%')) +
  ggtitle('TIMP 2') +
  xlab('Oxygen Tension') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.timp2.o2.days  <-  ggplot(dCq.O2, aes(day, -TIMP_2)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('**'),
              y_position = c(10),
              xmin = c(1),
              xmax = c(3), 
              vjust = 0.5) +
  ggtitle('TIMP 2') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.timp2.o2.days <- facet_signif(plot.timp2.o2.days, list(c(1), c(), c()))

plot.timp3.o2.grps <-  ggplot(dCq.O2, aes(grp, -TIMP_3)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1%', '5%', '21%')) +
  ggtitle('TIMP 3') +
  xlab('Oxygen Tension') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.timp3.o2.days  <-  ggplot(dCq.O2, aes(day, -TIMP_3)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('*'),
              y_position = c(0),
              xmin = c(1),
              xmax = c(2), 
              vjust = 0.5) +
  ggtitle('TIMP 3') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.timp3.o2.days <- facet_signif(plot.timp3.o2.days, list(c(), c(), c(1)))

plot_grid(plot.timp1.o2.grps,
          plot.timp2.o2.grps, 
          plot.timp3.o2.grps, 
          plot.timp1.o2.days, 
          plot.timp2.o2.days, 
          plot.timp3.o2.days,
          labels = 'AUTO', 
          ncol = 3, 
          nrow = 2)

## For the glucose groups
# plot brachyury expression
plot.T.gluc.grps  <-  ggplot(dCq.gluc, aes(grp, -T)) +
  geom_boxplot() +
  scale_x_discrete(limits = c('1mM', '2.5mM', '5mM')) +
  facet_grid(. ~ day) +
  ggtitle('Brachyury') +
  xlab('Glucose Concentration') +
  ylab(expression(paste('-', Delta, 'Cq'))) 

plot.T.gluc.days  <-  ggplot(dCq.gluc, aes(day, -T)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  ggtitle('Brachyury') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  xlab('Day') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot_grid(plot.T.gluc.grps, 
          plot.T.gluc.days,
          ncol = 1,
          nrow = 2,
          labels = 'AUTO')

# Plot the anabolic genes
plot.agg.gluc.grps  <-  ggplot(dCq.gluc, aes(grp, -Agg)) +
  geom_boxplot() +
  scale_x_discrete(limits = c('1mM', '2.5mM', '5mM')) +
  facet_grid(. ~ day) +
  ggtitle('Aggrecan') +
  xlab('Glucose Concentration') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.agg.gluc.days  <-  ggplot(dCq.gluc, aes(day, -Agg)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  ggtitle('Aggrecan') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.col1.gluc.grps  <-  ggplot(dCq.gluc, aes(grp, -Col_IA1)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1mM', '2.5mM', '5mM')) +
  ggtitle('Collagen I') +
  xlab('Glucose Concentration') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.col1.gluc.days  <-  ggplot(dCq.gluc, aes(day, -Col_IA1)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('**', '***', '**', '***'),
              y_position = c(-6, -4, -4, -4),
              xmin = c(1, 1, 1, 1),
              xmax = c(2, 3, 3, 3), 
              vjust = 0.5) +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ggtitle('Collagen I') +
  xlab('Day') +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.col1.gluc.days <- facet_signif(plot.col1.gluc.days, list(c(1, 2), c(3), c(4)))

plot.col2.gluc.grps  <-  ggplot(dCq.gluc, aes(grp, -Col_IIA1)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1mM', '2.5mM', '5mM')) +
  ggtitle('Collagen II') +
  xlab('Glucose Concentration') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.col2.gluc.days  <-  ggplot(dCq.gluc, aes(day, -Col_IIA1)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('***','***', '*'),
              y_position = c(4, 4, 4),
              xmin = c(1, 1, 2),
              xmax = c(2, 2, 3), 
              vjust = 0.5) +
  ggtitle('Collagen II') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.col2.gluc.days <- facet_signif(plot.col2.gluc.days, list(c(1), c(2), c(3)))

plot_grid(plot.agg.gluc.grps, 
          plot.col1.gluc.grps, 
          plot.col2.gluc.grps, 
          plot.agg.gluc.days, 
          plot.col1.gluc.days, 
          plot.col2.gluc.days,
          labels = 'AUTO',
          ncol = 3, 
          nrow = 2)

# Plot the catabolic genes
plot.mmp1.gluc.grps <-  ggplot(dCq.gluc, aes(grp, -MMP_1)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1mM', '2.5mM', '5mM')) +
  ggtitle('MMP 1') +
  xlab('Glucose Concentration') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.mmp1.gluc.days  <-  ggplot(dCq.gluc, aes(day, -MMP_1)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('**','**', '*'),
              y_position = c(0, -1, 2),
              xmin = c(1, 1, 1),
              xmax = c(3, 3, 3), 
              vjust = 0.5) +
  ggtitle('MMP 1') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.mmp1.gluc.days <- facet_signif(plot.mmp1.gluc.days, list(c(1),c(2),c(3)))

plot.mmp3.gluc.grps <-  ggplot(dCq.gluc, aes(grp, -MMP_3)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1mM', '2.5mM', '5mM')) +
  ggtitle('MMP 3') +
  xlab('Glucose Concentration') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.mmp3.gluc.days  <-  ggplot(dCq.gluc, aes(day, -MMP_3)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('***','***','***','***','***'),
              y_position = c(1, 0, 1, 0, -1),
              xmin = c(1, 2, 1, 2, 1),
              xmax = c(3, 3, 3, 3, 3), 
              vjust = 0.5) +
  ggtitle('MMP 3') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.mmp3.gluc.days <- facet_signif(plot.mmp3.gluc.days, list(c(1, 2), c(3, 4), c(5)))

plot.mmp13.gluc.grps <-  ggplot(dCq.gluc, aes(grp, -MMP_13)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1mM', '2.5mM', '5mM')) +
  ggtitle('MMP 13') +
  xlab('Glucose Concentration') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.mmp13.gluc.days  <-  ggplot(dCq.gluc, aes(day, -MMP_13)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('***','***','***','***','***', '**'),
              y_position = c(1, 0, 1, 0, 1, 0),
              xmin = c(1, 2, 1, 2, 1, 2),
              xmax = c(3, 3, 3, 3, 3, 3), 
              vjust = 0.5) +
  ggtitle('MMP 13') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.mmp13.gluc.days <- facet_signif(plot.mmp13.gluc.days, list(c(1, 2), c(3, 4), c(5, 6)))

plot.mt1.gluc.grps <-  ggplot(dCq.gluc, aes(grp, -MT_1)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1mM', '2.5mM', '5mM')) +
  ggtitle('MT 1') +
  xlab('Glucose Concentration') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.mt1.gluc.days  <-  ggplot(dCq.gluc, aes(day, -MT_1)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('*'),
              y_position = c(5),
              xmin = c(2),
              xmax = c(3), 
              vjust = 0.5) +
  ggtitle('MT 1') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.mt1.gluc.days <- facet_signif(plot.mt1.gluc.days, list(c(1), c(), c()))

plot_grid(plot.mmp1.gluc.grps,
          plot.mmp3.gluc.grps, 
          plot.mmp13.gluc.grps, 
          plot.mt1.gluc.grps, 
          plot.mmp1.gluc.days,
          plot.mmp3.gluc.days, 
          plot.mmp13.gluc.days,
          plot.mt1.gluc.days,
          labels = 'AUTO',
          ncol = 4, 
          nrow = 2)

# Plot the catabolic inhibitor genes
plot.timp1.gluc.grps <-  ggplot(dCq.gluc, aes(grp, -TIMP_1)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1mM', '2.5mM', '5mM')) +
  ggtitle('TIMP 1') +
  xlab('Glucose Concentration') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.timp1.gluc.days  <-  ggplot(dCq.gluc, aes(day, -TIMP_1)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('***','***','***','***','***'),
              y_position = c(5, 4, 5, 4, 3.5),
              xmin = c(1, 2, 1, 2, 1),
              xmax = c(3, 3, 3, 3, 3), 
              vjust = 0.5) +
  ggtitle('TIMP 1') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.timp1.gluc.days <- facet_signif(plot.timp1.gluc.days, list(c(1, 2), c(3, 4), c(5)))

plot.timp2.gluc.grps <-  ggplot(dCq.gluc, aes(grp, -TIMP_2)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  scale_x_discrete(limits = c('1mM', '2.5mM', '5mM')) +
  ggtitle('TIMP 2') +
  xlab('Glucose Concentration') +
  ylab(expression(paste('-', Delta, 'Cq')))

plot.timp2.gluc.days  <-  ggplot(dCq.gluc, aes(day, -TIMP_2)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('***','*','***','**','***'),
              y_position = c(4, 3, 4, 3, 2),
              xmin = c(1, 2, 1, 2, 1),
              xmax = c(3, 3, 3, 3, 3), 
              vjust = 0.5) +
  ggtitle('TIMP 2') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.timp2.gluc.days <- facet_signif(plot.timp2.gluc.days, list(c(1, 2), c(3, 4), c(5)))

plot.timp3.gluc.grps <-  ggplot(dCq.gluc, aes(grp, -TIMP_3)) +
  geom_boxplot() +
  facet_grid(. ~ day) +
  geom_signif(annotations = c('*'),
              y_position = c(4.5),
              xmin = c(2),
              xmax = c(3), 
              vjust = 0.5) +
  scale_x_discrete(limits = c('1mM', '2.5mM', '5mM')) +
  ggtitle('TIMP 3') +
  xlab('Glucose Concentration') +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.timp3.gluc.grps <- facet_signif(plot.timp3.gluc.grps, list(c(), c(1), c()))

plot.timp3.gluc.days  <-  ggplot(dCq.gluc, aes(day, -TIMP_3)) +
  geom_boxplot() +
  facet_grid(. ~ grp) +
  geom_signif(annotations = c('***','***','***','***','***', '**'),
              y_position = c(6, 5, 5, 4, 6, 5),
              xmin = c(1, 2, 1, 2, 1, 2),
              xmax = c(3, 3, 3, 3, 3, 3), 
              vjust = 0.5) +
  ggtitle('TIMP 3') +
  xlab('Day') +
  scale_x_discrete(breaks = c("Day 1", "Day 5", "Day 10"),
                   labels = c("1", "5", "10")) +
  ylab(expression(paste('-', Delta, 'Cq')))
plot.timp3.gluc.days <- facet_signif(plot.timp3.gluc.days, list(c(1, 2), c(3, 4), c(5, 6)))

plot_grid(plot.timp1.gluc.grps,
          plot.timp2.gluc.grps, 
          plot.timp3.gluc.grps, 
          plot.timp1.gluc.days, 
          plot.timp2.gluc.days, 
          plot.timp3.gluc.days,
          labels = 'AUTO',
          ncol = 3, 
          nrow = 2)
