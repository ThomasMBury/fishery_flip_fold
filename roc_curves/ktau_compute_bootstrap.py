#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 11:32:28 2020


Compute the kendall tau values of the bootstrapped EWS

@author: tbury
"""

import pandas as pd
import numpy as np
import scipy.stats as stats



# Import bootstrapped EWS data for Fold bifurcation
job_nums = np.concatenate((np.arange(52602,52639),np.arange(52640,52703)))
cols = ['Time','Variance','Lag-1 AC','Lag-2 AC','Smax']

df_list=[]
for job_num in job_nums:
    filepath = '../hagrid/ricker_bootstrap/Jobs/job-'+str(job_num)+'/flip/ews_intervals.csv'
    df_temp = pd.read_csv(filepath)
    df_temp2 = df_temp[df_temp['Unnamed: 1']=='Mean'][cols]
    df_temp2['job_num']=job_num
    df_list.append(df_temp2)
    
# Concatenate list of dataframes
df_ews_boot = pd.concat(df_list,axis=0).set_index(['job_num','Time'])


# Function to compute kendall tau values from dataframe of EWS
def ktau_compute(df):
    ews_names = df.columns
    ktau_dic = {}
    for ews in ews_names:
        ktau, pval = stats.kendalltau(df.index,df[ews])
        ktau_dic[ews] = ktau
    
    return ktau_dic


# Compute kendall tau values for all job numbers
list_dic_ktau = []
for job in job_nums:
    df_temp = df_ews_boot.loc[job]
    dic_ktau = ktau_compute(df_temp)
    dic_ktau['job_num']=job
    list_dic_ktau.append(dic_ktau)
    
# Combine kendall tau dictionaries into single dataframe
dic_full_ktau = {}
for k in list_dic_ktau[0].keys():
    dic_full_ktau[k] = np.array([d[k] for d in list_dic_ktau])

# Create dataframe for ktua vlaues
    
df_ktau = pd.DataFrame(dic_full_ktau)

# Export df
#df_ktau.to_csv('ktau_vals/ktau_bootstrap_flip.csv')




