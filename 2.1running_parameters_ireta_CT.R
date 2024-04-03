library(dplyr)
library(SNFtool)

setwd('D:/Helab/ADHDsubtyping/Stats/subtyping_nm_gender/spectral/CT')
dir.create('./parameter_iterations/alpha_0.3',recursive = TRUE)
dir.create("./parameter_iterations/alpha_0.4",recursive = TRUE)
dir.create("./parameter_iterations/alpha_0.5",recursive = TRUE)
dir.create("./parameter_iterations/alpha_0.6",recursive = TRUE)
dir.create("./parameter_iterations/alpha_0.7",recursive = TRUE)
dir.create("./parameter_iterations/alpha_0.8",recursive = TRUE)

# setting up empty data-frames to compare nearest neighbors (K) and alpha parameters
ktests=seq(10,30)
resultdf_1 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))
resultdf_2 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))
resultdf_3 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))
resultdf_4 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))

# setting up empty dataframes to save NMI scores across parameters
CT_NMI = data.frame(matrix(0, ncol = 21, nrow = 219))

# range of parameters tested
alphas_list <- c("0.3", "0.4", "0.5", "0.6", "0.7", "0.8")

# calculating number of clusters per alpha and K value
# calculating each K for each alpha and saving them

# importing the different data types as individual participant matrices 
subjects <- read.csv("ids.csv", header=TRUE)
CT <- read.csv("ct_z.csv", header=TRUE)

# creating participant distance matrices using euclidean distances
CT = standardNormalization(CT)
Dist_CT = dist2(as.matrix(CT),as.matrix(CT));

# for each parameter combination - SNF is used to create a similarity matrix and (1) the optimal number of clusters as well as
# (2) the normalized mutual information (NMI) between each measure and the matrix are calculated and recorded
for(param in 1:length(alphas_list)){
  alpha <- as.numeric(as.character(alphas_list[param]))
  print(alpha)
  for(test in 1:nrow(resultdf_1)){
    K <- resultdf_1$K[test]
    print(test)
    
    # creating participant affinity matrices for each data type
    AM_CT = affinityMatrix(Dist_CT,K,alpha)
    
    # calculating optimal number of clusters given the created similarity matrix
    EstClust_1 <-estimateNumberOfClustersGivenGraph(AM_CT, NUMC=2:10)[[1]]
    EstClust_2 <-estimateNumberOfClustersGivenGraph(AM_CT, NUMC=2:10)[[2]]
    EstClust_3 <-estimateNumberOfClustersGivenGraph(AM_CT, NUMC=2:10)[[3]] 
    EstClust_4 <-estimateNumberOfClustersGivenGraph(AM_CT, NUMC=2:10)[[4]]
    
    # adding cluster number to tracking file for the given parameter
    resultdf_1[test,(param +1)]<-EstClust_1
    resultdf_2[test,(param + 1)]<-EstClust_2
    resultdf_3[test,(param +1)]<-EstClust_3
    resultdf_4[test,(param +1)]<-EstClust_4
    
    #calculating NMI scores based on similarity matrix
    CT_NMIScores <-rankFeaturesByNMI(list(CT), AM_CT)
    # organizing NMI scores
    CT_NMI[, test] <- CT_NMIScores[[1]][1]
  }
  # setting directory and saving NMI score files
  directory <- (paste("./parameter_iterations/alpha_", alpha, sep=""))
  write.csv(CT_NMI, file=file.path(directory, paste("/NMIscores_CT.csv", sep="")))
}

# saving estimated number of clusters files
directory <- ("./parameter_iterations/")

write.csv(resultdf_1, file=file.path(directory, paste("Bstnumclus_eigengaps.csv", sep="")))
write.csv(resultdf_2, file=file.path(directory, paste("2nd_Bstnumclus_eigengaps.csv", sep="")))
write.csv(resultdf_3, file=file.path(directory, paste("Bstnumclus_rotationcost.csv", sep="")))
write.csv(resultdf_4, file=file.path(directory, paste("2nd_Bstnumclus_rotationcost.csv", sep="")))
