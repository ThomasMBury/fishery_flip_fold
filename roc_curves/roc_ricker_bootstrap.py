#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:58:49 2020

Compute roc curves for EWS in the Ricker model with bootstrapping
Is there improved statistical performance with bootstrapping?
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




## ROC curve data for EWS prior to Fold
df_roc_var = roc_compute(df_roc_fold['Truth'],df_roc_fold['Variance'])
df_roc_ac = roc_compute(df_roc_fold['Truth'],df_roc_fold['Lag-1 AC'])
df_roc_smax = roc_compute(df_roc_fold['Truth'],df_roc_fold['Smax'])

## Plot ROC curve for EWS prior to Fold bif

plt.title('Prior to Fold bifurcation')
plt.plot(df_roc_var['fpr'], df_roc_var['tpr'], 'b', label = 'Variance, AUC = %0.2f' % df_roc_var['auc'][0])
plt.plot(df_roc_ac['fpr'], df_roc_ac['tpr'], 'g', label = 'Lag-1 AC, AUC = %0.2f' % df_roc_ac['auc'][0])
plt.plot(df_roc_smax['fpr'], df_roc_smax['tpr'], 'k', label = 'Smax, AUC = %0.2f' % df_roc_smax['auc'][0])

plt.legend(loc = 'lower right')
plt.plot([0, 1], [0, 1],'r--')
plt.xlim([-0.02, 1])
plt.ylim([0, 1.02])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()



## ROC curve data for EWS prior to Flip
df_roc_var = roc_compute(df_roc_flip['Truth'],df_roc_flip['Variance'])
df_roc_ac = roc_compute(df_roc_flip['Truth'],df_roc_flip['Lag-2 AC'])
df_roc_smax = roc_compute(df_roc_flip['Truth'],df_roc_flip['Smax'])

## Plot ROC curve for EWS prior to Flip bif

plt.title('Prior to Flip bifurcation')
plt.plot(df_roc_var['fpr'], df_roc_var['tpr'], 'b', label = 'Variance, AUC = %0.2f' % df_roc_var['auc'][0])
plt.plot(df_roc_ac['fpr'], df_roc_ac['tpr'], 'g', label = 'Lag-2 AC, AUC = %0.2f' % df_roc_ac['auc'][0])
plt.plot(df_roc_smax['fpr'], df_roc_smax['tpr'], 'k', label = 'Smax, AUC = %0.2f' % df_roc_smax['auc'][0])

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

    






