import abagen
import pandas as pd
import numpy as np

# fetch gene data from AHBA
files = abagen.fetch_microarray(donors='all', verbose=0, data_dir='D:/Helab/ADHDsubtyping/AHBA/data_gene/')

# get gene expression matrix (region * gene)
expression_out_left = abagen.get_expression_data(atlas='D:/Helab/ADHDsubtyping/AHBA/atlas/ADHD2vsHC_ct_roi.nii.gz',
                                               atlas_info='D:/Helab/ADHDsubtyping/AHBA/atlas/ct_atlas_adhd2vscon_L.csv',
                                               return_donors=True,
                                               ibf_threshold=0.5,
                                               data_dir='D:/Helab/ADHDsubtyping/AHBA/data_gene')
# Removes genes in expression with differential stability < threshold (0.1)
expression_left_filter,stability_left = abagen.keep_stable_genes(list(expression_out_left.values()),
                                                             threshold=0.1,percentile=False,return_stability=True)
expression_left_filter = pd.concat(expression_left_filter).groupby('label').mean()
np.savetxt('D:/Helab/ADHDsubtyping/AHBA/gene_expression_on_brain/left_stability.csv', stability_left, delimiter=',')
expression_left_filter.to_csv('D:/Helab/ADHDsubtyping/AHBA/gene_expression_on_brain/left_DS01_expression.csv')
