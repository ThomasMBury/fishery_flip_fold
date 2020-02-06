#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:58:49 2020

Compute roc curves for EWS in the Ricker model (without bootstrapping)
Export data in a df for plotting in MMA

@author: tbury
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.metrics as metrics


# Function to compute ROC data from truth and indicator vals
# and return a df.
def roc_compute(truth_vals, indicator_vals):
    
    # Compute ROC curve and threhsolds using sklearn
    fpr, tpr, thresholds = metrics.roc_curve(truth_vals,indicator_vals)
    
    # Compute AUC (area under curve)
    auc = metrics.auc(fpr, tpr)
    
    # Put into a DF
    dic_roc = {'fpr':fpr, 'tpr':tpr, 'thresholds':thresholds, 'auc':auc}
    df_roc = pd.DataFrame(dic_roc)

    return df_roc


# Import Kendall tau values for Fold, Flip and Null trajectories
columns=['Variance','Lag-1 AC',
         'Lag-2 AC','Lag-3 AC',
         'Kurtosis','Skewness',
         'Coefficient of variation',
         'Standard deviation',
         'Smax']
df_ktau_fold = pd.read_csv('../data_export/ricker_fold_sigma0p04/ktau.csv',usecols=columns)[columns]
df_ktau_flip = pd.read_csv('../data_export/ricker_flip_sigma0p04/ktau.csv',usecols=columns)[columns]
df_ktau_null = pd.read_csv('../data_export/ricker_null_sigma0p04_noAIC/ktau.csv',usecols=columns)[columns]

# Add column for truth value
df_ktau_fold['Truth']=1
df_ktau_flip['Truth']=1
df_ktau_null['Truth']=0

# Combine dataframes to provide roc data
df_roc_fold = pd.concat([df_ktau_fold,df_ktau_null],axis=0,ignore_index=True)
df_roc_flip = pd.concat([df_ktau_flip,df_ktau_null],axis=0,ignore_index=True)




## ROC curve for Variance
df_roc_var = roc_compute(df_roc_fold['Truth'],df_roc_fold['Variance'])




## Plot ROC curve for ML
tpr = df_roc_var['tpr']
fpr = df_roc_var['fpr']
auc = df_roc_var['auc'][0]
plt.title('Receiver Operating Characteristic')
plt.plot(fpr, tpr, 'b', label = 'DL prediction, AUC = %0.2f' % auc)
plt.legend(loc = 'lower right')
plt.plot([0, 1], [0, 1],'r--')
plt.xlim([-0.02, 1])
plt.ylim([0, 1.02])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()





### Concatenate ROC dataframes for export
#df_roc = pd.concat([df_roc_ml, df_roc_var, df_roc_ac], axis=0)
#df_roc.to_csv('data/sims1/roc_data/roc_bif_gen_t300_400.csv')

    






