# Install necessary libraries
import os
import pandas as pd
import pickle
import pcntoolkit as pcn
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1 grab data files
os.chdir('D:\\Helab\\normative_model\\anat_HBR_sex')
processing_dir = os.getcwd()

x_train = pd.read_pickle('x_train.pkl')
batch_effects_train = pd.read_pickle('trbefile.pkl')
# add sex to training set
x_train = pd.concat([x_train, batch_effects_train[['sex']]], axis=1)
batch_effects_train.drop('sex', axis=1, inplace=True)

x_test = pd.read_pickle('x_test.pkl')
batch_effects_test = pd.read_pickle('tsbefile.pkl')
# add sex to testing set
x_test = pd.concat([x_test, batch_effects_test[['sex']]], axis=1)
batch_effects_test.drop('sex', axis=1, inplace=True)

# save files
with open('x_train_sex.pkl', 'wb') as file:
    pickle.dump(x_train, file)
with open('x_test_sex.pkl', 'wb') as file:
    pickle.dump(x_test, file)
with open('trbefile_site.pkl', 'wb') as file:
    pickle.dump(batch_effects_train, file)
with open('tsbefile_site.pkl', 'wb') as file:
    pickle.dump(batch_effects_test, file)

# Files and Folders grooming
output_path = os.path.join(processing_dir, 'Models/')
log_dir = os.path.join(processing_dir, 'log/')
if not os.path.isdir(output_path):
    os.mkdir(output_path)
if not os.path.isdir(log_dir):
    os.mkdir(log_dir)

y_train_sa = pd.read_pickle('y_train_sa.pkl')
y_test_sa = pd.read_pickle('y_test_sa.pkl')
y_train_sa_scale100 = y_train_sa/100
y_test_sa_scale100 = y_test_sa/100
with open('y_train_sa_scale100.pkl', 'wb') as file:
    pickle.dump(y_train_sa_scale100, file)
with open('y_test_sa_scale100.pkl', 'wb') as file:
    pickle.dump(y_test_sa_scale100, file)

total_sa_tr = pd.DataFrame()
total_sa_tr['LH'] = pd.DataFrame(y_train_sa.iloc[:, 0:112]).apply(lambda x: x.sum(), axis=1)
total_sa_tr['RH'] = pd.DataFrame(y_train_sa.iloc[:, 112:220]).apply(lambda x: x.sum(), axis=1)
total_sa_ts = pd.DataFrame()
total_sa_ts['LH'] = pd.DataFrame(y_test_sa.iloc[:, 0:112]).apply(lambda x: x.sum(), axis=1)
total_sa_ts['RH'] = pd.DataFrame(y_test_sa.iloc[:, 112:220]).apply(lambda x: x.sum(), axis=1)


# Estimate and apply the normative model to patients
pcn.normative.estimate(covfile=os.path.join(processing_dir, 'x_train_sex.pkl'),
                       respfile=os.path.join(processing_dir, 'y_train_sa_scale100.pkl'),
                       trbefile=os.path.join(processing_dir, 'trbefile_site.pkl'),
                       testcov=os.path.join(processing_dir, 'x_test_sex.pkl'),
                       testresp=os.path.join(processing_dir, 'y_test_sa_scale100.pkl'),
                       tsbefile=os.path.join(processing_dir, 'tsbefile_site.pkl'),
                       alg='hbr',
                       log_path=os.path.join(processing_dir, 'log/'),
                       output_path=os.path.join(processing_dir, 'Models/'),
                       outputsuffix='_hbrsa100',
                       savemodel=True,
                       saveoutput=True)

# Check results
EV_ct = pd.read_pickle('EXPV_hbrctgender.pkl')
MSLL_ct = pd.read_pickle('MSLL_hbrct.pkl')
NLL_ct = pd.read_pickle('NLL_hbrct.pkl')
pRho_ct = pd.read_pickle('pRho_hbrct.pkl')
Rho_ct = pd.read_pickle('Rho_hbrct.pkl')
RMSE_ct = pd.read_pickle('RMSE_hbrct.pkl')
SMSE_ct = pd.read_pickle('SMSE_hbrct.pkl')

EV_sa = pd.read_pickle('EXPV_hbrsa.pkl')
MSLL_sa = pd.read_pickle('MSLL_hbrsa.pkl')
NLL_sa = pd.read_pickle('NLL_hbrsa.pkl')
pRho_sa = pd.read_pickle('pRho_hbrsa.pkl')
Rho_sa = pd.read_pickle('Rho_hbrsa.pkl')
RMSE_sa = pd.read_pickle('RMSE_hbrsa.pkl')
SMSE_sa = pd.read_pickle('SMSE_hbrsa.pkl')

EV_sb = pd.read_pickle('EXPV_hbrsb.pkl')
MSLL_sb = pd.read_pickle('MSLL_hbrsb.pkl')
NLL_sb = pd.read_pickle('NLL_hbrsb.pkl')
pRho_sb = pd.read_pickle('pRho_hbrsb.pkl')
Rho_sb = pd.read_pickle('Rho_hbrsb.pkl')
RMSE_sb = pd.read_pickle('RMSE_hbrsb.pkl')
SMSE_sb = pd.read_pickle('SMSE_hbrsb.pkl')

mod_sa1 = pd.read_pickle('./Models/NM_0_0_hbrsa.pkl')
Z_ct = pd.read_pickle('./results_ct/Z_hbrctgender.pkl')
Z_sa = pd.read_pickle('./results_sa/Z_hbrsagender.pkl')
Z_sa100 = pd.read_pickle('./results_sa100/Z_hbrsa100.pkl')
Z_sb = pd.read_pickle('./results_sb/Z_hbrsbgender.pkl')

# Step 7 save z-score
x_test = pd.read_pickle('x_test.pkl')
y_test_ct = pd.read_pickle('y_test_ct.pkl')
y_test_sb = pd.read_pickle('y_test_sb.pkl')

Z_ct.columns = y_test_ct.columns
Z_ct.index = x_test.index
Z_sa100.columns = y_test_ct.columns
Z_sa100.index = x_test.index
Z_sb.columns = y_test_sb.columns
Z_sb.index = x_test.index

Z_ct.to_csv('Z_ct.csv', header=True, index=True)
Z_sa100.to_csv('Z_sa100.csv', header=True, index=True)
Z_sb.to_csv('Z_sb.csv', header=True, index=True)

# 10-fold within training set to evaluate the robustness of the models
os.chdir('D:\\Helab\\normative_model\\anat_HBR_sex\\10foldCV')
processing_dir = os.getcwd()

x_train = pd.read_pickle('x_train_sex.pkl')
batch_effects_train = pd.read_pickle('trbefile_site.pkl')
y_train_ct = pd.read_pickle('y_train_ct.pkl')
y_train_sa_scale100 = pd.read_pickle('y_train_sa_scale100.pkl')
y_train_sb = pd.read_pickle('y_train_sb.pkl')

pcn.normative.estimate(covfile=os.path.join(processing_dir, 'x_train_sex.pkl'),
                       respfile=os.path.join(processing_dir, 'y_train_ct.pkl'),
                       trbefile=os.path.join(processing_dir, 'trbefile_site.pkl'),
                       cvfolds=10,
                       alg='hbr',
                       log_path=os.path.join(processing_dir, 'log/'),
                       output_path=os.path.join(processing_dir, 'Models/'),
                       outputsuffix='_10CV',
                       savemodel=True,
                       saveoutput=True)

# Check results
EV = pd.read_pickle('./EXPV_10CV.pkl')
Rho = pd.read_pickle('./Rho_10CV.pkl')
SMSE = pd.read_pickle('./SMSE_10CV.pkl')
MSLL = pd.read_pickle('./MSLL_10CV.pkl')
EV.mean(axis=0)
Rho.mean(axis=0)
SMSE.mean(axis=0)
MSLL.mean(axis=0)
z = pd.read_pickle('./Z_10CV.pkl')
EV.to_csv('./EV_10CV.csv')
Rho.to_csv('./Rho_10CV.csv')
SMSE.to_csv('./SMSE_10CV.csv')
MSLL.to_csv('./MSLL_10CV.csv')
# Plotting evaluation measures
SMSE.max(axis=0)
MSLL.max(axis=0)

sns.set(style='white', context='paper', font_scale=1)
sns.set_style({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans']})

f, ax = plt.subplots(1, 2)
f.set_figwidth(3.5)
f.set_figheight(1.75)
sns.histplot(SMSE, ax=ax[0], kde=True)
ax[0].set_xlabel('Standardized \nmean squared error')
ax[0].tick_params(pad=-2)
ax[0].set_ylabel('Count', labelpad=-1)
ax[0].set_ylim([0, 60])
ax[0].get_legend().remove()
sns.histplot(EV, ax=ax[1], kde=True)
ax[1].set_xlabel('Explained variance')
ax[1].tick_params(pad=-2)
ax[1].set_ylabel('')
ax[1].set_ylim([0, 60])
ax[1].set_yticklabels('')
ax[1].get_legend().remove()
f.savefig('./model_evaluation.png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.close()




