# **Normative modeling of brain morphology reveals ADHD biological subtypes**

This repository provides main code, relevant toolbox, and data for manuscript "Bu, X., Zhao, Y., Zheng, X., Fu, Z., Zhang, K., Sun, X., Cui, Z., Xia, M., Ma, L., Liu, N., Lu, J., Zhao, G., Ding, Y., Deng, Y., Wang, J., Chen, R., Zhang, H., Men, W., Wang, Y., Gao, J., Tan, S., Sun, L., Qin, S., Tao, S., Wang, Y., Dong, Q., Cao, Q., Yang, L., & He, Y. (2024). Normative growth modeling of brain morphology reveals neuroanatomical heterogeneity and biological subtypes in children with ADHD. bioRxiv https://www.biorxiv.org/content/10.1101/2024.03.16.582202v1."

## **Overview**
Following contents include software, source code, and demo data. All data necessary to replicate our results have been made publicly available, including regional cortical thickness, cortical thickness deviations from normative modeling, and demographic information.

## **Code**
Our analyses were carried out using following open source packages:
1. **Multiscale Desikan-Kiliany parcellation** files [1] were downloaded using the netneurotools toolbox (https://github.com/netneurolab/netneurotools). We use 219 parcellations for our analysis.

2. **Normative modeling** was implemented with the PCNtoolkit (https://pcntoolkit.readthedocs.io/en/latest/) [2].

3. **Spectral clustering** was applied to the similarity matrix of cortical thickness deviation to cluster the ADHD into subgroups. This analysis was performed using an open R package SNFtool (http://compbio.cs.toronto.edu/SNF/SNF/Software.html) [3].

4. **ComBat** was used to harmonize raw cortical thickness across sites for group-level analyses (https://github.com/Jfortin1/ComBatHarmonization) [4].

5. The **Allen Human Brain Atlas** (AHBA) datasets were preprocessed using the Python toolbox abagen (https://github.com/rmarkello/abagen) [5].

6. The **Gene Ontology enrichment analysis** was conducted using online tool Metascape (https://metascape.org) [6].

7. **Visualization**: Python nilearn and R ggplot.

The main analyses were carried out using following codes:
### **1. Normative modeling analysis:**
   running normative modeling: 1.structure_hbr_sex.py

### **2. Clustering analysis:**
   determining best parameters and cluster number: 2.1 running_parameters_ireta_CT.R
   
   running spectral clustering: 2.2 CT_clustering.R
   
   evaluating clutering stability by bootstrapping: 2.3 CT_clustering_stability.R

### **3. Subgroup analysis:**
   univariate analysis: 3.1 group_analysis.R
   
   multivariate analysis (brain-behavior relationship): 3.2 pls_symptom_predict.m
   
   medication response: 3.3 repeated_measure_drugeffect.R
   
   gene expression analysis: 3.4 PLS_gene_brain.m

   spin test: https://github.com/frantisekvasa/rotate_parcellation.

Data analysis was performed using MATLAB R2021b, Python v3.8, and R v4.2.2.
   
## **Data**
1. Basic demographic information for all participants and ADHD subtypes: demographic_ADHD.csv, demographic_TDC.csv, and ADHDsubtypes_info.csv

2. Raw regional cortical thickness for each participant: cortical_thickness_ADHD.csv and cortical_thickness_TDC.csv

3. Cortical thickness devations estimated form normative modeling: cortical thickness deviation.csv

## **References**
1. Cammoun L, Gigandet X, Meskaldji D, et al., Mapping the human connectome at multiple scales with diffusion spectrum MRI. J Neurosci Meth 2012; 203: 386-97.

2. Rutherford S, Kia SM, Wolfers T, et al., The normative modeling framework for computational psychiatry. Nat Protoc 2022; 17: 1711-34.

3. Wang B, Mezlini AM, Demir F, et al., Similarity network fusion for aggregating data types on a genomic scale. Nat Methods 2014; 11: 333-U19.

4. Fortin JP, Cullen N, Sheline YI, et al., Harmonization of cortical thickness measurements across scanners and sites. Neuroimage 2018; 167: 104-20.

5. Markello RD, Arnatkeviciute A, Poline JB, et al., Standardizing workflows in imaging transcriptomics with the abagen toolbox. Elife 2021; 10.

6. Zhou YY, Zhou B, Pache L, et al., Metascape provides a biologist-oriented resource for the analysis of systems-level datasets. Nat Commun 2019; 10.
