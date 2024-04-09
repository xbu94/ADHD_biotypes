### Script to run resampling stability analysis

library(corrplot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)

setwd('D:/Helab/ADHDsubtyping/Stats/subtyping_nm_gender/spectral/CT')
source("D:/Helab/ADHDsubtyping/Stats/subtyping_nm_gender/spectral/CT/CT_bootstrapping.R")

subjects <- read.csv("./ids.csv", header=TRUE) #list of participants
CT <- read.csv("./ct_z.csv", header=TRUE) # cortical thickness (n=219)
directory <- ("D:/Helab/ADHDsubtyping/Stats/subtyping_nm_gender/spectral/CT/clustering_stability/")

numboot=1000
nsub=351
K=30
alpha=0.8
bootsize=0.8
clusters=2

# Setting up which participants will be included in each permutation
permutation_matrix <- bootstrapping_SNF(numboot=numboot, nsub=nsub, bootsize=bootsize)
# Getting clustering solutions for all the permutations of sampled participants using SNF 
clus_sil <- clustering(perms=permutation_matrix, bootsize=bootsize, K=K, alpha=alpha, clusters=clusters, CT=CT)
# Dividing output matrix into clusters and silhouette widths
clus_out <- clus_sil[1:(numboot), ]
silhouette_width <- clus_sil[(numboot+1):(numboot*2), ]
## getting NMI scores for each permutation
All_NMI_scores <- NMI_scores(perms=permutation_matrix, bootsize=bootsize, K=K, alpha=alpha, clusters=clusters, CT=CT)
# getting the adjusted rand index between all clustering solutions
list_randindex <- stability(clus_out=clus_out, perms=permutation_matrix) # returns b rand index:adjusted rand index
list_adjustedrandindex <- list_randindex[ ,1001:2000]
# Calculate how often each participant is clustered together and the probability that they will be clustered together
percent_agree <- percent_agree(clus_out=clus_out)

write.csv(permutation_matrix, file=file.path(directory, paste("permutation_matrix_1000perms.csv", sep="")))
write.csv(percent_agree, file=file.path(directory, paste("Percent_agree_1000perms.csv", sep="")))
write.csv(list_adjustedrandindex, file=file.path(directory, paste("adjrandindex_1000perms.csv", sep="")))
write.csv(list_randindex, file=file.path(directory, paste("Rand_indices_1000perms.csv", sep="")))
write.csv(clus_out, file=file.path(directory, paste("clu_solution_1000perms.csv", sep="")))
write.csv(silhouette_width, file=file.path(directory, paste("silhouette_1000perms.csv", sep="")))
write.csv(All_NMI_scores, file=file.path(directory, paste("NMI_scores_0.8_1000perms.csv", sep="")))

##### Sorting out measures and getting the top 20
measures <- read.csv("D:/Helab/ADHDsubtyping/Stats/subtyping_nm_gender/roi_anat.csv", header = F)
measures = data.frame(roi = measures[c(1:219),])
library(stringr)
measures$roi = gsub("CT_", "", measures$roi) # remove "CT_" from roi names
#All_NMI_scores <- read.csv(file=file.path(directory, paste("NMI_scores_0.8_1000perms.csv", sep="")))
#All_NMI_scores$X <- NULL

top20_NMI_measures <- data.frame(matrix(0, nrow = 20, ncol = numboot))
for (i in seq(1:length(All_NMI_scores))) {
  intermediate <- All_NMI_scores[ , i]
  intermediate <- cbind(measures, intermediate)
  intermediate$intermediate <- as.numeric(as.character(intermediate$intermediate))
  intermediate <- intermediate[order(-intermediate$intermediate), ]
  top20_NMI_measures[ ,i] <- intermediate[1:20, 1]
}
rm(intermediate)
write.csv(top20_NMI_measures, file=file.path(directory, paste("Top20_NMI_scores_0.8_1000perms.csv", sep="")))

#### counting how often measures are in the top 20
top20_measures <- top20_NMI_measures
top20_measures$Column <- seq(1:length(top20_measures[,1]))
top20_measures <- melt(top20_measures, id.vars=c('Column'),var='Index')
top20_measures$Column <- NULL
top20_measures$Index <- NULL

measure_count <- as.data.frame(table(top20_measures$value))
names(measure_count)[names(measure_count) == "Var1"] <- "Measures"
measure_count <- measure_count[order(-measure_count$Freq),]
write.csv(measure_count, file=file.path(directory, paste("Count_top20_NMI_0.8_1000perms.csv", sep="")))

