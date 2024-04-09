import os
import pandas as pd
import numpy as np
import pickle
import pcntoolkit as pcn

# Step 1 grab data files
os.chdir('.\\anat_HBR_update')
processing_dir = os.getcwd()
norm_ct = pd.read_csv('ct_norm14y.csv', index_col=0)
adhd_ct = pd.read_csv('ct_adhd14y.csv', index_col=0)
roi_cor = norm_ct.columns[3:223]

# Step 2: Prepare training norms and testing ADHD sets
# preparing training dataset
cov_train = norm_ct[['age', 'sex', 'site']]
x_train = cov_train[['age']]
# batch effect for train set
batch_effects_train = cov_train[['sex', 'site']]
batch_effects_train = batch_effects_train.join(pd.get_dummies(batch_effects_train.site))
batch_effects_train.drop('site', axis=1, inplace=True)
# select training imaging features
y_train_ct = norm_ct.loc[:, roi_cor]

# preparing testing dataset
cov_test = adhd_ct[['age', 'sex', 'site']]
x_test = cov_test[['age']]
# batch effect for test set
batch_effects_test = cov_test[['sex', 'site']]
batch_effects_test = batch_effects_test.join(pd.get_dummies(batch_effects_test.site))
batch_effects_test.drop('site', axis=1, inplace=True)
# since there are no data of site 1/2 in train set, add two columns for site1&2 with 0 value\n
# to make sure the dimension of batch-effect files are same between train and test set
batch_effects_test.insert(loc=1, column='1', value=np.full((batch_effects_test.shape[0], 1), 0))
batch_effects_test.insert(loc=2, column='2', value=np.full((batch_effects_test.shape[0], 1), 0))
# imaging features for test set
y_test_ct = adhd_ct.loc[:, roi_cor]

# Step 3 Create binary files for model estimation
with open('x_train.pkl', 'wb') as file:
    pickle.dump(x_train, file)
with open('y_train_ct.pkl', 'wb') as file:
    pickle.dump(y_train_ct, file)
with open('trbefile.pkl', 'wb') as file:
    pickle.dump(batch_effects_train, file)

with open('x_test.pkl', 'wb') as file:
    pickle.dump(x_test, file)
with open('y_test_ct.pkl', 'wb') as file:
    pickle.dump(y_test_ct, file)
with open('tsbefile.pkl', 'wb') as file:
    pickle.dump(batch_effects_test, file)

# Step 4: Files and Folders grooming
output_path = os.path.join(processing_dir, 'Models/')
log_dir = os.path.join(processing_dir, 'log/')
if not os.path.isdir(output_path):
    os.mkdir(output_path)
if not os.path.isdir(log_dir):
    os.mkdir(log_dir)

# Step 5: Estimate and apply the normative model to patients
pcn.normative.estimate(covfile=os.path.join(processing_dir, 'x_train.pkl'),
                       respfile=os.path.join(processing_dir, 'y_train_ct.pkl'),
                       trbefile=os.path.join(processing_dir, 'trbefile.pkl'),
                       testcov=os.path.join(processing_dir, 'x_test.pkl'),
                       testresp=os.path.join(processing_dir, 'y_test_ct.pkl'),
                       tsbefile=os.path.join(processing_dir, 'tsbefile.pkl'),
                       alg='hbr',
                       log_path=os.path.join(processing_dir, 'log/'),
                       output_path=os.path.join(processing_dir, 'Models/'),
                       outputsuffix='_hbrct',
                       savemodel=True,
                       saveoutput=True)

# Step 6 Check results
EV_ct = pd.read_pickle('EXPV_hbrct.pkl')
MSLL_ct = pd.read_pickle('MSLL_hbrct.pkl')
NLL_ct = pd.read_pickle('NLL_hbrct.pkl')
pRho_ct = pd.read_pickle('pRho_hbrct.pkl')
Rho_ct = pd.read_pickle('Rho_hbrct.pkl')
RMSE_ct = pd.read_pickle('RMSE_hbrct.pkl')
SMSE_ct = pd.read_pickle('SMSE_hbrct.pkl')
Z_ct = pd.read_pickle('Z_hbrct.pkl')

# Step 7 save z-score
x_test = pd.read_pickle('x_test.pkl')
Z_ct.columns = roi_cor
Z_ct.index = x_test.index
Z_ct.to_csv('Z_ct.csv', header=True, index=True)

# Step 8 10-fold within training set to evaluate the robustness of the models
os.chdir('.\\10foldCV')
processing_dir = os.getcwd()
x_train = pd.read_pickle('x_train.pkl')
batch_effects_train = pd.read_pickle('trbefile.pkl')
y_train_ct = pd.read_pickle('y_train_ct.pkl')

pcn.normative.estimate(covfile=os.path.join(processing_dir, 'x_train.pkl'),
                       respfile=os.path.join(processing_dir, 'y_train_ct.pkl'),
                       trbefile=os.path.join(processing_dir, 'trbefile.pkl'),
                       cvfolds=10,
                       alg='hbr',
                       log_path=os.path.join(processing_dir, 'log/'),
                       output_path=os.path.join(processing_dir, 'Models/'),
                       outputsuffix='_10foldct',
                       savemodel=True,
                       saveoutput=True)
