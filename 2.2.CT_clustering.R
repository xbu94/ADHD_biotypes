### Running spectral clustering on the CT data
## This script is run after Running_parameter_iterations.r in order to determine the optimal 
# cluster number, hyperparameter (alpha), and nearest neighbors parameter (K)

library(SNFtool)
library(cluster)
library(MASS)
library(ggplot2)
library(corrplot)
library(dplyr)
library(broom)
library(tidyr)
library(psych)
library(dunn.test)
library(fossil)
setwd('D:/Helab/ADHDsubtyping/Stats/subtyping_nm_gender/spectral/CT')

# source code for function to determine data integration and clustering across resampling using SNF and spectral clustering
source("D:/Helab/ADHDsubtyping/Stats/subtyping_nm_gender/spectral/CT/robust_core-clustering_function_CT.R")

# set up for clustering analysis
# importing the different data types as individual participant matrices 
subjects <- read.csv("./sex_validation/ids.csv", header=TRUE) #list of participants
CT <- read.csv("./sex_validation/Z_ct_allsex.csv", header=TRUE) # cortical thickness (n=219)

# normalizing measures within each data type using a function from the SNF package
CT = standardNormalization(CT)

# setting the parameters (finalized after comparisons using Running_parameter_iterations.r )
K =30;		# number of neighbors, usually (10~30), usually sample size/10
alpha = 0.8;  	# hyperparameter, usually (0.3~0.8)

# creating participant distance matrices using euclidean distances
Dist_CT = dist2(as.matrix(CT),as.matrix(CT))

# creating participant affinity matrices within each data type
AM_CT = affinityMatrix(Dist_CT,K,alpha) 

### Calculating similarity matrix and spectral clustering groups
# setting output directory
directory <- ("D:/Helab/ADHDsubtyping/Stats/subtyping_nm_gender/spectral/CT/sex_validation/")

# Clustering
C = 2 #setting cluster number
subgroup = spectralClustering(AM_CT,C)
displayClustersWithHeatmap(AM_CT, subgroup)
# calling function to cluster participants using spectral clustering across resampling 80% of participants 1000 times
robust.W = RobustCoreClusteringMatrix(feature.affinity.mat.list = list(AM_CT),exp.num.samples = 1000, num.clusts = C)
#Two matrices - Dense Core Cluster Matrix and Sparse Core Cluster Matrix
dense <- robust.W[1]
dense <- matrix(unlist(dense), ncol = 351, byrow = TRUE)
sparse <- robust.W[2]
sparse <- matrix(unlist(sparse), ncol = 351, byrow = TRUE)

# displaying clusters
displayClustersWithHeatmap(dense, spectralClustering(dense, C))
displayClustersWithHeatmap(sparse, spectralClustering(sparse, C))

# saving an image of the cluster heatmap
png('./SNF_analysis_newsa/SN_dense_0.8_30_1000perms.png',width = 1500, height = 1500, units = "px", bg = "white", res = 300)
displayClustersWithHeatmap(dense, spectralClustering(dense, C))
dev.off()

png('./SNF_analysis_newsa/SN_sparse_0.8_30_1000perms.png',width = 1500, height = 1500, units = "px", bg = "white", res = 300)
displayClustersWithHeatmap(sparse, spectralClustering(sparse, C))
dev.off()

# calculating normalized mutual information (NMI) based off of the original data types and the clustering similarity matrix
CT_NMIScores <-rankFeaturesByNMI(list(CT), dense) 

# separating and organizing scores
CT_scores <- as.data.frame(CT_NMIScores[[1]][1])
names(CT_scores) <- c("NMI")

CT_rank <- as.data.frame(CT_NMIScores[[2]][1])
names(CT_rank) <- c("rank")

#saving csv of NMI scores and clustering similarity matrix
write.csv(CT_scores, file=file.path(directory, paste("CT_scores_k30_0.8_1000perms.csv", sep="")))
write.csv(CT_rank, file=file.path(directory, paste("CT_rank_k30_0.8_1000perms.csv", sep="")))
write.matrix(dense, file=file.path(directory, paste("dense_k30_0.8_matrix.csv", sep="")),sep = ",")
write.matrix(sparse, file=file.path(directory, paste("sparse_k30_0.8_matrix.csv", sep="")),sep = ",")

# Find cluster labels of individuals using the robust clustering similarity matrix
robust.groups.df = RobustCoreClusteringClusters(core.clustering.list = robust.W,num.clusts = C,verbose = TRUE)
clusters <- cbind(subjects, robust.groups.df$groups)
colnames(clusters)[2] = 'groups'
table(clusters$groups)
# rename cluster to align with main results
clusters$subtype = ifelse(clusters$groups == "1", "2", NA)
clusters$subtype = ifelse(clusters$groups == "2", "1", clusters$subtype)

## calculating silouette width for each participant and the silhouette plot
dissim <- 1 - dense
dissim <- as.matrix(dissim)
clusters$subtype <- as.integer(clusters$subtype)
sil <- silhouette(clusters$subtype, dmatrix = dissim)
summary(sil)

# saving the silhouette plot
plot(sil, col = c("palevioletred", "sandybrown"))
tiff("Silhouette_CT.tiff", width = 18, height = 18, units = "cm", bg="white",res = 200)
plot(sil, col = c("palevioletred", "sandybrown"))
dev.off()

# adding silhouette widths for each individual
clusters$silhouette_width <- 0
for (i in 1:351){
  clusters[i, 4] <- sil[i ,3]
}

write.csv(clusters, file=file.path(directory, paste("2clust_groups_CT.csv", sep="")))


## Check if anyone doesn't reliably cluster 
unreliable.patients = UnrelClustPatientsFullCohort(core.clustering.list = robust.W,verbose = TRUE)

## Check if anyone doesn't reliably cluster within their given cluster group
unreliable.patients.by.grp = UnrelClustPatientsByGrp(core.clustering.list = robust.W,id.group.df = robust.groups.df,verbose = TRUE)

# subtypes visualization
library(qgraph)
subtypes = list(clusters$groups)
#dense_k30_0_8_matrix <- read_csv("~/BNU_Helab/data/nm_subtyping_adhd/snf_gender_newsa/dense_k30_0.8_matrix.csv", 
#                                 col_names = FALSE)
node_color = lapply(subtypes, function(x) replace(x, x == 2, 'pink'))
node_color = lapply(node_color, function(x) replace(x, x == 1, 'sandybrown'))
node_color = node_color[[1]]
qgraph(dense, layout = 'spring', groups = subtypes,
       minimum = 0.4, color = node_color, vsize=3,labels = FALSE,
       posCol = c('burlywood4'), filetype = 'pdf', filename = 'subtypes_qgraph2',
       width = 10, height = 6)




