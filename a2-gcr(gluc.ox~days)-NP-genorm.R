# This Script runs two gene normalization strategies on A select group of 
# genes to determine the most stable ones to use later in the determination of
# deltaCqs.

# Author <- Lukas Jaworski, l.jaworski@umiami.edu
# Institution <- University of Miami
# PI <- Alicia Jackson, a.jackson2@miami.edu

# Copyright University of Miami 2017
###############################################################################
## Library imports
###############################################################################
library(NormqPCR)
library(ggplot2)
library(cowplot)
###############################################################################
## Functions
###############################################################################

###############################################################################
## Script
###############################################################################
# Retrieving and arranging the data

# Set the working directory to where the data is kept, make it workstation agnostic
if (Sys.info()[['sysname']]=='Linux'){
  working.directory <- paste0('/home/', Sys.info()[['user']], 
                              '/Dropbox/Lab/Experiments/PCR/',
                              'a2-gcr(gluc.ox~days)-NP-ref-genes/Cleaned Data')
                              
} else {
  working.directory <- paste0('C:/Users/', Sys.info()[['user']],
                              '/Dropbox/Lab/Experiments/PCR/',
                              'a2-gcr(gluc.ox~days)-NP-ref-genes/Cleaned Data')
}
working.directory <- paste0('P:/Dropbox/Lab/Experiments/PCR/',
                              'a2-gcr(gluc.ox~days)-NP-ref-genes/Cleaned Data')
setwd(working.directory)

# Load the Phenotype data from the text file
pheno.data <- read.csv("GeneNorm_PhenoData.txt", sep="\t")

# Order the phenotype data according to sample name so it aligns with the assayData
# (This should do a hash lookup so the sample names don't have to be ordered)
pheno.data <- pheno.data[with(pheno.data, order(Sample)), ]

# change day number into a factor
pheno.data$Day <- factor(pheno.data$Day)

# Copy the sample names into the rownames 
# (needs to be done otherwise read.qPCR won't be able to match the data)
rownames(pheno.data) <- pheno.data$Sample

# Metadata for our annotated dataframe
meta <- data.frame(c('Days cultured', 'exp group', 'sample number'))
 
# Need to at a label to the metadata column names, to match the eSet format
names(meta) <- 'labelDescription'

# Finally making the annotated dataframe
annotated_pheno <- AnnotatedDataFrame(pheno.data, meta, c('sampleNames','sampleColumns'))

# Getting the Cqs from a text file and adding the phenotype data 
qPCRBatch.Cq <- read.qPCR('GeneNorm_Cqs.txt',annotated_pheno)

###############################################################################
# Picking the housekeeping Genes using GeNorm and Normfinder

# Since all the data is housekeeping genes we can just take the rownames
hkgn <- rownames(exprs(qPCRBatch.Cq))

# select the experimental groups that were present in the study:
# this is used in the Normfinder analysis
groups <- pData(qPCRBatch.Cq)[,'Exp_Grp']

# Run the GeNorm analysis on the whole dataset
genorm.result <- selectHKs(qPCRBatch.Cq, 
                           method = "geNorm", 
                           Symbols = hkgn, 
                           minNrHK = 2,
                           log = TRUE)

# Rerun the analysis using the NormFinder algorithm on the whole data set
normfinder.result <- selectHKs(qPCRBatch.Cq, 
                               method = "NormFinder",
                               Symbols = hkgn, 
                               minNrHK = (length(hkgn)-1), 
                               log = TRUE, 
                               group = groups)

# Compare ranks of the two algorithms
genorm.result$ranking
normfinder.result$ranking

###############################################################################
# Rerunning the analysis on the two pooled experiments seperately

# Subset the data
# Selecting the O2 tension experiment data
qPCRBatch.Cq.O2 <- qPCRBatch.Cq[,qPCRBatch.Cq@phenoData@data$Exp_Grp == '21%' |
                                  qPCRBatch.Cq@phenoData@data$Exp_Grp == '5%' |
                                  qPCRBatch.Cq@phenoData@data$Exp_Grp == '1%']

# Selecting the glucose concentration experiment data
qPCRBatch.Cq.gluc <- qPCRBatch.Cq[,qPCRBatch.Cq@phenoData@data$Exp_Grp == '1 mM' | 
                                    qPCRBatch.Cq@phenoData@data$Exp_Grp == '2.5 mM' |
                                    qPCRBatch.Cq@phenoData@data$Exp_Grp == '5 mM']

# Run the GeNorm analysis on the data subsets
# Run the GeNorm algorithm on the O2 tension experiment data
genorm.result.O2 <- selectHKs(qPCRBatch.Cq.O2,
                              method = "geNorm", 
                              Symbols = hkgn, 
                              minNrHK = 2, 
                              log = TRUE, 
                              na.rm=TRUE)

# Run the GeNorm algorithm on the glucose concentration experiment data
genorm.result.gluc <- selectHKs(qPCRBatch.Cq.gluc,
                                method = "geNorm",
                                Symbols = hkgn, 
                                minNrHK = 2,
                                log = TRUE, 
                                na.rm=TRUE)

# Rerun the analysis using the NormFinder algorithm
# This analysis doesn't work as some of the genes return NaN values
# 
# # Selecting the subgroups for the Normfinder analysis
# groups.O2 <- as.factor(as.character(pData(qPCRBatch.Cq.O2)[,'Exp_Grp']))
# groups.gluc <- as.factor(as.character(pData(qPCRBatch.Cq.gluc)[,'Exp_Grp']))
# # Run the Normfinder algorithm on the O2 tension experiment data
# normfinder.result.O2 <- selectHKs(qPCRBatch.Cq.O2, 
#                                   method = "NormFinder", 
#                                   Symbols = colnames(t(exprs(qPCRBatch.Cq.O2))),
#                                   minNrHK = (length(hkgn)-1), 
#                                   log = TRUE, 
#                                   group = groups.O2,
#                                   na.rm=TRUE)
# # Run the Normfinder algorithm on the O2 tension experiment data
# normfinder.result.gluc <- selectHKs(qPCRBatch.Cq.gluc, 
#                                     method = "NormFinder", 
#                                     Symbols = hkgn,
#                                     minNrHK = (length(hkgn)-1), 
#                                     log = TRUE, 
#                                     group = groups.gluc,
#                                     na.rm=TRUE)

# The normfinder algorithm does not seem to work with only a subset of the data,
# will look into it later. For now only GeNorm results will be graphed and used

# Compare the rankings of the different experiments
genorm.result.gluc$ranking
# normfinder.result.gluc$ranking
genorm.result.O2$ranking
# normfinder.result.O2$ranking


###############################################################################
# Graphing the results

# Graphing the overall stability and variability of the housekeeping genes 
# across the entire dataset

# Making a dataframe with only the mean stability ranking M of the whole dataset
genorm.for.graph <- cbind(genorm.result$ranking[8:2],genorm.result$meanM)
genorm.for.graph[7,1] <- paste0(genorm.result$ranking[2],'\n', genorm.result$ranking[1])
colnames(genorm.for.graph) <- c('ranking', 'meanM')
genorm.df <- data.frame(genorm.for.graph)
genorm.df$meanM <- as.numeric(as.character(genorm.df$meanM))

# Making a dataframe with only the Pairwaise Variation of the whole dataset
genorm.for.graph.var <- cbind(names(genorm.result$variation), as.numeric(as.character(genorm.result$variation)))
colnames(genorm.for.graph.var) <- c('gene.number', 'variation')
genorm.df.var <- data.frame(genorm.for.graph.var)
genorm.df.var$variation <- as.numeric(as.character(genorm.df.var$variation))

# Making the Mean stability graph of the whole dataset
meanM.graph <- ggplot(genorm.df, aes(ranking, meanM)) +
  theme_classic(base_size = 22) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=genorm.df$ranking) +
  scale_y_continuous(breaks=seq(0, round(max(genorm.df$meanM), digits=2), by=0.5)) +
  xlab('Genes') +
  ylab('Mean Stability M')
# plot(meanM.graph)

# Making the Pairwise Variation graph of the whole dataset
var.graph <- ggplot(genorm.df.var, aes(gene.number, variation)) +
  theme_classic(base_size = 22) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=genorm.df.var$gene.number) +
  scale_y_continuous(breaks=seq(0, round(max(genorm.df.var$variation), digits=2), by=0.1)) +
  xlab(bquote('Normalization Factors ('*NF[n]/NF[n+1]*')')) +
  ylab('Pairwise Variation')
# plot(var.graph)

# Plotting the Mean Stability and Pairwaise Variation in one figure
plot_grid(meanM.graph, var.graph,nrow = 2, ncol = 1, labels = c('A', 'B'))

# Making a dataframe with only the mean stability ranking M of the O2 tension data
genorm.O2.for.graph <- cbind(genorm.result.O2$ranking[8:2],genorm.result.O2$meanM)
genorm.O2.for.graph[7,1] <- paste0(genorm.result.O2$ranking[2],'\n', genorm.result.O2$ranking[1])
colnames(genorm.O2.for.graph) <- c('ranking', 'meanM')
genorm.df.O2 <- data.frame(genorm.O2.for.graph)
genorm.df.O2$meanM <- as.numeric(as.character(genorm.df.O2$meanM))

# Making a dataframe with only the Pairwaise Variation of the O2 tension data
genorm.O2.for.graph.var <- cbind(names(genorm.result.O2$variation), as.numeric(as.character(genorm.result.O2$variation)))
colnames(genorm.O2.for.graph.var) <- c('gene.number', 'variation')
genorm.df.O2.var <- data.frame(genorm.O2.for.graph.var)
genorm.df.O2.var$variation <- as.numeric(as.character(genorm.df.O2.var$variation))

# Making the Mean stability graph of the O2 tension data
O2.meanM.graph <- ggplot(genorm.df.O2, aes(ranking, meanM)) +
  theme_classic(base_size = 22) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=genorm.df.O2$ranking) +
  scale_y_continuous(breaks=seq(0, round(max(genorm.df.O2$meanM), digits=2), by=0.5)) +
  xlab('Genes') +
  ylab('Average stability M')
# plot(O2.meanM.graph)

# Making the Pairwise Variation graph of the O2 tension data
O2.var.graph <- ggplot(genorm.df.O2.var, aes(gene.number, variation)) +
  theme_classic(base_size = 22) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=genorm.df.O2.var$gene.number) +
  scale_y_continuous(breaks=seq(0, round(max(genorm.df.O2.var$variation), digits=2), by=0.1)) +
  xlab(bquote('Normalization Factors ('*NF[n]/NF[n+1]*')')) +
  ylab('Pairwise Variation')
# plot(O2.var.graph)\

# Plotting the Mean Stability and Pairwaise Variation in one figure
plot_grid(O2.meanM.graph, O2.var.graph,nrow = 2, ncol = 1, labels = c('A', 'B'))

# Making a dataframe with only the mean stability ranking M of the glucose concentration data
genorm.gluc.for.graph <- cbind(genorm.result.gluc$ranking[8:2],genorm.result.gluc$meanM)
genorm.gluc.for.graph[7,1] <- paste0(genorm.result.gluc$ranking[2],'\n', genorm.result.gluc$ranking[1])
colnames(genorm.gluc.for.graph) <- c('ranking', 'meanM')
genorm.df.gluc <- data.frame(genorm.gluc.for.graph)
genorm.df.gluc$meanM <- as.numeric(as.character(genorm.df.gluc$meanM))

# Making a dataframe with only the Pairwaise Variation of the glucose concentration data
genorm.gluc.for.graph.var <- cbind(names(genorm.result.gluc$variation), as.numeric(as.character(genorm.result.gluc$variation)))
colnames(genorm.gluc.for.graph.var) <- c('gene.number', 'variation')
genorm.df.gluc.var <- data.frame(genorm.gluc.for.graph.var)
genorm.df.gluc.var$variation <- as.numeric(as.character(genorm.df.gluc.var$variation))

# Making the Mean stability graph of the glucose concentration data
gluc.meanM.graph <- ggplot(genorm.df.gluc, aes(ranking, meanM)) +
  theme_classic(base_size = 22) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=genorm.df.gluc$ranking) +
  scale_y_continuous(breaks=seq(0, round(max(genorm.df.gluc$meanM), digits=2), by=0.5)) +
  xlab('Genes') +
  ylab('Mean Stability M')
# plot(gluc.meanM.graph)

# Making the Pairwise Variation graph of the glucose concentration data
gluc.var.graph <- ggplot(genorm.df.gluc.var, aes(gene.number, variation)) +
  theme_classic(base_size = 22) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=genorm.df.gluc.var$gene.number) +
  scale_y_continuous(breaks=seq(0, round(max(genorm.df.gluc.var$variation), digits=2), by=0.1)) +
  xlab(bquote('Normalization Factors ('*NF[n]/NF[n+1]*')')) +
  ylab('Pairwise Variation')
# plot(gluc.var.graph)

# Plotting the Mean Stability and Pairwaise Variation in one figure
plot_grid(gluc.meanM.graph, gluc.var.graph,nrow = 2, ncol = 1, labels = c('A', 'B'))

