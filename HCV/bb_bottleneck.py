# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 13:41:48 2020

@author: Noam
"""
import pandas as pd
import os
import re
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns


infection_order = {'HCV-P8':1, 'HCV-P1':2, 'HCV-P11':3, 'HCV-P7':5, 'HCV-P4':6, 'HCV-P3':7, 'HCV-P5':8, 'HCV-P2':9, 'HCV-P6':10, 'HCV-P9':11, 'HCV-P10':12, 'HCV-PS2':0}
infection_order = pd.DataFrame.from_dict(infection_order, orient='index').reset_index().rename(columns={'index':'Sample', 0:'infection_order'})


def prepare_for_bb(donor_freq, recipient_freq, out_path):
    donor = pd.read_csv(donor_freq, '\t')
    recipient = pd.read_csv(recipient_freq, '\t')
    merged = pd.merge(donor, recipient, on=['Ref', 'Pos', 'Base'], suffixes=['_donor', '_recipient'])
    merged['Var_read_count_recipient'] = merged.Freq_recipient * merged.Read_count_recipient
    merged['Var_read_count_recipient'] = merged['Var_read_count_recipient'].round()
    merged = merged[merged.Pos.isin(range(738,2758))]
    # filtering - change this
    # remove hvrs
    merged = merged[~(merged.Pos.isin(range(1479,1561))) & ~(merged.Pos.isin(range(1707,1783))) & ~(merged.Pos.isin(range(2037,2068)))]
    #merged = merged[(merged.Ref != '-') & (merged.Var_read_count_recipient > 0) & (merged.Rank_donor != 0)]
    merged = merged[(merged.Ref != '-') & (merged.Var_read_count_recipient > 0) & (merged.Rank_donor == 0)]
    merged['Freq_donor'] = 1 - merged['Freq_donor']
    merged['Freq_recipient'] = 1 - merged['Freq_recipient']
    merged.loc[merged.Freq_recipient < 0.001, 'Freq_recipient'] = 0
    #merged.loc[merged.Freq_donor < 0.01, 'Freq_donor'] = 0
    merged = merged[merged.Freq_donor >= 0.001]
    #
    merged.to_csv(out_path + '.full.csv', index=False)
    merged[['Freq_donor', 'Freq_recipient', 'Read_count_donor', 'Var_read_count_recipient']].to_csv(out_path, header=False, index=False, sep='\t')
    return merged

for f in ['Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/freqs/' + f for f in os.listdir('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/freqs/') if f.count('-') < 2]:
    prepare_for_bb('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/freqs/HCV-PS2.freqs', 
               f,
               'X:/volume2/noam/hcv/BB_bottleneck_output/7/' + f.split('/')[-1].split('.')[0] + '.csv')

#for f in ['/sternadi/home/volume2/noam/hcv/BB_bottleneck_output/' + f for f in os.listdir('/sternadi/home/volume2/noam/hcv/BB_bottleneck_output/') if f.endswith('')]:
#    os.system('module load R/3.6.1 & Rscript /sternadi/home/volume2/noam/hcv/BB_bottleneck/Bottleneck_size_estimation_approx.r --file ' + f + ' > ' + f + '.bb.txt')
    
 

def parse_bb_output(out_file):
    with open(out_file) as f:
        f = f.read()
    pattern = re.compile('"Bottleneck size"\n\[1\] (\d+)\n\[1\] "confidence interval left bound"\n\[1\] (\d+)\n\[1\] "confidence interval right bound"\n\[1\] (\d+)')
    results = pattern.findall(f)[0]
    return [int(d) for d in results]

all_results = []
for f in ['X:/volume2/noam/hcv/BB_bottleneck_output/7/' + f for f in os.listdir('X:/volume2/noam/hcv/BB_bottleneck_output/7/') if f.endswith('.bb.txt')]:

    try:
        results = parse_bb_output(f)
        results = [f.split('/')[-1].split('.')[0]] + results
        all_results.append(results)
    except:
        print(f)
all_results = pd.DataFrame(all_results, columns=['Sample', 'bottleneck', 'ci_left', 'ci_right'])
all_results = pd.merge(all_results, infection_order, on=['Sample'], how='left')
all_results = all_results[~(all_results.infection_order.isna()) & (all_results.Sample != 'HCV-PS2')].sort_values('infection_order')
fig, ax = plt.subplots(nrows=1, ncols=1)
slope, intercept, r_value, p_value, std_err = stats.linregress(all_results.infection_order, all_results.bottleneck)
sns.regplot(all_results.infection_order, all_results.bottleneck, ax = ax)
ax.set_title("y=%fx+%f, r=%f, p=%f" % (slope,intercept, r_value, p_value))

