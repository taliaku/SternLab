import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import sys
sys.path.append('X:/volume2/noam/SternLab')
from freqs_utilities import compare_positions_between_freqs, estimate_insertion_freq
from blast_utilities import blast_to_df
import re
from blast_utilities import get_unaligned_reads

#sys.path.append('X:/volume2/noam/AccuNGS/')
#import co-occurs_to_stretches

infection_order = {'HCV-P8':1, 'HCV-P1':2, 'HCV-P11':3, 'HCV-P7':5, 'HCV-P4':6, 'HCV-P3':7, 'HCV-P5':8, 'HCV-P2':9, 'HCV-P6':10, 'HCV-P9':11, 'HCV-P10':12, 'HCV-PS2':0}
infection_order = pd.DataFrame.from_dict(infection_order, orient='index').reset_index().rename(columns={'index':'Sample', 0:'infection_order'})
infection_order['paper_name'] = 'Patient ' + infection_order.infection_order.astype(str)
infection_order.loc[infection_order['Sample'] == 'HCV-PS2', 'paper_name'] = 'Source'




def unite_all_freq_files(freqs_dir, out=None):
    """
    unites all frequency file into one file. add the following:
    Replica, Degree, Sample, Time
    :param freqs_dir: a directory with frequency files in a format of pTIME-DEGREEREPLICA.freqs
    :param out: output directory to save results
    :return: a merged data frame
    """
    freqs_df = []
    # read all freq files and add the sample name to the data frame
    files = [os.path.join(freqs_dir, f) for f in os.listdir(freqs_dir)]
    for f in files:
        print(f)
        curr_df = pd.read_csv(f, sep='\t')
        sample = os.path.basename(f).split('.')[0]
        #sample = os.path.basename(f)
        curr_df['Sample'] = sample
        freqs_df.append(curr_df)
    df = pd.concat(freqs_df)

    if out != None:
        df.to_csv(os.path.join(out, 'all_freqs.csv'), index=False)
    return df

unite_all_freq_files('Z:/volume1/noam/hcv_data/180503_OST_FINAL_03052018_pipeline/freqs', 'Z:/volume1/noam/hcv_data/180503_OST_FINAL_03052018_pipeline/')
unite_all_freq_files('Z:/volume1/noam/hcv_data/180503_OST_FINAL_03052018_pipeline_optimized1/freqs', 'Z:/volume1/noam/hcv_data/180503_OST_FINAL_03052018_pipeline_optimized1/')
unite_all_freq_files('Z:/volume1/noam/hcv_data/180503_OST_FINAL_03052018_pipeline_optimized2/freqs', 'Z:/volume1/noam/hcv_data/180503_OST_FINAL_03052018_pipeline_optimized2/')
df = unite_all_freq_files('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/freqs/', 'Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/')


df = pd.read_csv('Z:/volume1/noam/hcv_data/180503_OST_FINAL_03052018_pipeline_optimized2/all_freqs.csv')

colors={'HCV-WG-P1':'blue', 'HCV-WG-P6':'red', 'HCV-WG-P8':'green', 'HCV-WG-S2':'orange'}
fig, ax = plt.subplots(nrows=1, ncols=1)
for key, group in df[(df.Base != df.Ref) & (df.Ref == '-')].groupby('Sample'):
    group.plot(kind='scatter', x='Pos', y='estimated_freq', ax=ax, color=colors[key])

fig, ax = plt.subplots(nrows=1, ncols=1)
for key, group in df[(df.Base != df.Ref) & (df.Ref != '-')].groupby('Sample'):
    group.plot(kind='scatter', x='Pos', y='Freq', ax=ax, label=key)

def compare_two_samples(freq_path1, freq_path2):
    freq_path1 = pd.read_csv(freq_path1, '\t')
    freq_path2 = pd.read_csv(freq_path2, '\t')
    df = pd.merge(freq_path1, freq_path2, on=['Pos', 'Base', 'Ref'])
    fig, ax = plt.subplots(nrows=1, ncols=1)
    df[(df.Base != df.Ref) & (df.Ref != '-')].plot(x='Freq_x', y='Freq_y', kind='scatter', ax=ax)
    return fig
    



#### all mutations 
df = pd.read_csv('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/all_freqs.csv')
mutation_type = pd.read_csv('X:/volume2/noam/hcv/references/D90208.1_optimized4.mutation_type.peptides.csv')
df = df[(df.Sample != 'HCV-P12')]
df = df[df.Pos.isin(range(738,2758))]
df = pd.merge(df, mutation_type, on=['Pos', 'Base', 'Ref'], how='left')
df.loc[df['mutation_type'] == '[]', 'mtype'] = 'synonymous'
df.loc[((df['mutation_type'] != '[]') & ~(df['mutation_type'].isna())), 'mtype'] = 'non-synonymous'
df.loc[((df['mutation_type'].isna()) & (df.Base != df.Ref)), 'mtype'] = 'non-coding'
df.loc[df['Base'] == '-', 'mtype'] = 'frameshift'
df['rb'] = df.Ref + df.Base
df.loc[df['rb'].isin(['AG', 'GA', 'CT', 'TC']), 'transition_transversion'] = 'transition'
df.loc[df['rb'].isin(['AC', 'AT', 'CA', 'CG', 'TA', 'TG', 'GT', 'GC']), 'transition_transversion'] = 'transversion'


## new polyprotein

df = pd.read_csv('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/all_freqs.csv')
mutation_type = pd.read_csv('X:/volume2/noam/hcv/references/D90208.1_optimized4.mutation_type.csv')
mutation_type = mutation_type[['Pos', 'Ref', 'Base', 'protein', 'mutation_type']]
mutation_type['aa'] = (((mutation_type.Pos - 330) / 3) + 1).astype(int).astype(float)
df = df[(df.Sample != 'HCV-P12')]
df = df[df.Pos.isin(range(738,2758))]
df = pd.merge(df, mutation_type, on=['Pos', 'Base', 'Ref'], how='left')
df.loc[df['mutation_type'] == '[]', 'mtype'] = 'synonymous'
df.loc[((df['mutation_type'] != '[]') & ~(df['mutation_type'].isna())), 'mtype'] = 'non-synonymous'
df.loc[((df['mutation_type'].isna()) & (df.Base != df.Ref)), 'mtype'] = 'non-coding'
df.loc[df['Base'] == '-', 'mtype'] = 'frameshift'
df['rb'] = df.Ref + df.Base
df.loc[df['rb'].isin(['AG', 'GA', 'CT', 'TC']), 'transition_transversion'] = 'transition'
df.loc[df['rb'].isin(['AC', 'AT', 'CA', 'CG', 'TA', 'TG', 'GT', 'GC']), 'transition_transversion'] = 'transversion'







### synonymous non-synonymous plot
fig, axes = plt.subplots(nrows=4, ncols=6)
fig.set_size_inches(25,18)
axes = axes.flatten()
fig.subplots_adjust(hspace=0.3)
i = 0
for key, group in df[(df.Base != df.Ref) & (df.Ref != '-') & (df.Read_count > 5)].groupby('Sample'):
    group[group.mtype == 'synonymous'].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='blue', label='syn')
    group[group.mtype == 'non-synonymous'].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='red', label='non-syn')
    group[group.mtype == 'non-coding'].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='black', label='non-coding')
    group[group.mtype == 'frameshift'].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='grey', label='frameshift')
    axes[i].set_title(key)
    axes[i].legend().remove()
    i += 1
axes[-4].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axes[-3].set_visible(False)
axes[-2].set_visible(False)
axes[-1].set_visible(False)
fig.savefig('X:/volume2/noam/hcv/variant_frequencies.png' , bbox_inches='tight', dpi=800)


#### transition trnasversions plot
fig, axes = plt.subplots(nrows=4, ncols=6)
fig.set_size_inches(25,18)
axes = axes.flatten()
fig.subplots_adjust(hspace=0.3)
i = 0
for key, group in df[(df.Base != df.Ref) & (df.Ref != '-') & (df.Read_count > 5)].groupby('Sample'):
    group[group.transition_transversion == 'transition'].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='blue', label='transition')
    group[group.transition_transversion == 'transversion'].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='red', label='transversion')
    group[group.mtype == 'frameshift'].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='grey', label='frameshift')
    axes[i].set_title(key)
    axes[i].legend().remove()
    i += 1
axes[-4].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axes[-3].set_visible(False)
axes[-2].set_visible(False)
axes[-1].set_visible(False)
fig.savefig('X:/volume2/noam/hcv/variant_frequencies2.png' , bbox_inches='tight', dpi=800)


###### syn not-syn and transition transversion plots




fig, axes = plt.subplots(nrows=3, ncols=4)
fig.set_size_inches(14,10)
axes = axes.flatten()
fig.subplots_adjust(hspace=0.3)
i = 0
sample_order = ['HCV-PS2',
 'HCV-P8',
 'HCV-P1',
 'HCV-P11',
 'HCV-P7',
 'HCV-P4',
 'HCV-P3',
 'HCV-P5',
 'HCV-P2',
 'HCV-P6',
 'HCV-P9',
 'HCV-P10',]
for key in sample_order:
    group = df[(df.Base != df.Ref) & (df.Ref != '-') & (df.Read_count > 5) & (df.Sample.str.count('-') == 1) & (df.Freq >= 0.01) & (df.Sample == key)]
    group[(group.mtype == 'synonymous') & (group.transition_transversion == 'transition')].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='#4285BB', label='synonymous transition')
    #group[(group.mtype == 'non-coding') & (group.transition_transversion == 'transition')].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='black', label='non-coding transition')
    try:
        group[(group.mtype == 'synonymous') & (group.transition_transversion == 'transversion')].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='#558B2F', label='synonymous transversion')
    except:
        pass
    group[(group.mtype == 'non-synonymous') & (group.transition_transversion == 'transition')].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='#FDA318', label='non-synonymous transition')
    try:
        group[(group.mtype == 'non-synonymous') & (group.transition_transversion == 'transversion')].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='#D51313', label='non-synonymous transversion')
    except:
        pass
    #group[(group.mtype == 'non-coding') & (group.transition_transversion == 'transversion')].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='black', marker='x', label='non-coding transversion')
    #group[group.mtype == 'frameshift'].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='grey', label='frameshift')
    axes[i].set_title( infection_order.loc[infection_order.Sample == key, 'paper_name'].values[0], fontsize=14)
    #axes[i].set_yscale('log')
    axes[i].legend().remove()
    axes[i].set_ylim(-0.05, 1.05)
    axes[i].set_xlabel('')
    axes[i].set_ylabel('')
    i += 1
axes[9].legend(loc='center', bbox_to_anchor=(1, -0.7), fontsize=16, ncol=2)
fig.text(0.5, 0.04, 'Position (Base)', ha='center', va='center', fontsize=16)
fig.text(0.06, 0.5, 'Mutation Frequency', ha='center', va='center', rotation='vertical', fontsize=16)
fig.savefig('X:/volume2/noam/hcv/variant_frequencies3.png' , bbox_inches='tight', dpi=800)





############ quantify diversity
sys.path.append('X:/volume2/noam/SternLab/RG_HIVC_analysis/')
from diversity_calculation import pis_calc
diversities = pis_calc(df, ['Sample'])
diversities = pd.merge(diversities, infection_order, on='Sample')[['paper_name', 'infection_order', 'Pi']]
diversities.to_csv('X:/volume2/noam/hcv/diversities.csv', index=False)

fig, ax = plt.subplots(nrows=1, ncols=1)
diversities = diversities[diversities.paper_name != 'Source']
sns.regplot(diversities.infection_order, diversities.Pi, ax=ax)
c, pvalue = stats.pearsonr(diversities.infection_order, diversities.Pi)
ax.set_title('Order of Infection vs. Pi Diversity\n' + 'r = ' + str(round(c, 3)) + ', p-value = ' + str(round(pvalue, 3)))
fig.set_size_inches(5,3)
ax.set_xlabel('Order of infection')
ax.set_ylabel('Pi Diversity')
fig.savefig('X:/volume2/noam/hcv/diversities_regression.png' , bbox_inches='tight', dpi=800)






## frequencies compared to source - correlations
compare_df = pd.merge(df, df[df.Sample=='HCV-PS2'], on=['Pos', 'Base', 'Ref'], suffixes=['', '_source'])
fig, axes = plt.subplots(nrows=4, ncols=6)
fig.set_size_inches(25,18)
axes = axes.flatten()
fig.subplots_adjust(hspace=0.4)
i = 0
correlations = {}
for key, group in compare_df[(compare_df.Ref != '-') & (compare_df.Read_count > 5) & (compare_df.Read_count_source > 5)].groupby('Sample'):
    sns.scatterplot(data=group, x='Freq', y='Freq_source', ax=axes[i], color='purple')
    #slope, intercept, r_value, p_value, std_err = stats.linregress(group['Freq'], group['Freq_source'])
    #axes[i].set_title(key + '\n' + "y=%fx+%f, r2=%f" % (slope,intercept, r_value**2))
    axes[i].set_title(key)
    correlations[key] = r_value
    axes[i].legend().remove()
    axes[i].set_xlabel('')
    axes[i].set_ylabel('')
    i += 1  
axes[-4].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axes[-3].set_visible(False)
axes[-2].set_visible(False)
axes[-1].set_visible(False)
correlations = pd.DataFrame.from_dict(correlations, orient='index').reset_index().rename(columns={'index':'Sample', 0:'correlation'})
correlations = pd.merge(correlations, infection_order, on='Sample', how='left')
correlations = correlations.sort_values('correlation', ascending=False)
correlations['correlation'] = correlations.correlation.round(3)
fig.text(0.5,0.07, "Sample frequencies", ha="center", va="center", fontsize=20)
fig.text(0.07,0.5, "Source frequencies", ha="center", va="center", rotation=90, fontsize=20)
fig.savefig('X:/volume2/noam/hcv/correlations.png' , bbox_inches='tight', dpi=800)

#fig, ax = plt.subplots(nrows=1, ncols=1)
#correlations2 = correlations[~(correlations.infection_order.isna())]
#slope, intercept, r_value, p_value, std_err = stats.linregress(correlations2.infection_order, correlations2.correlation)
#sns.regplot(x=correlations2.infection_order, y=correlations2.correlation, ax = ax)
#ax.set_title("y=%fx+%f, r2=%f" % (slope,intercept, r_value**2))
#

## color coded
fig, axes = plt.subplots(nrows=3, ncols=4)
fig.set_size_inches(14,10)
axes = axes.flatten()
fig.subplots_adjust(hspace=0.3)
i = 0
sample_order = ['HCV-PS2',
 'HCV-P8',
 'HCV-P1',
 'HCV-P11',
 'HCV-P7',
 'HCV-P4',
 'HCV-P3',
 'HCV-P5',
 'HCV-P2',
 'HCV-P6',
 'HCV-P9',
 'HCV-P10',]
compare_df = pd.merge(df, df[(df.Sample == 'HCV-PS2')], how='right', on=['Pos', 'Base','Ref','protein', 'mutation_type', 'aa', 'mtype', 'rb','transition_transversion'], suffixes=['', '_source'])
variants_to_look_at = df[(df.Base != df.Ref) & (df.Ref != '-') & (df.Base != '-') & (df.Read_count > 5) & (df.Freq > 0.01)][['Pos', 'Base']].drop_duplicates()
compare_df = pd.merge(compare_df, variants_to_look_at, how='right', on=['Pos', 'Base'])
correlations = []
for key in sample_order:
    group = compare_df[(compare_df.Sample == key)]
    group[(group.mtype == 'synonymous') & (group.transition_transversion == 'transition')].plot(kind='scatter', x='Freq_source', y='Freq', ax=axes[i], color='#4285BB', label='synonymous transition')
    try:
        group[(group.mtype == 'synonymous') & (group.transition_transversion == 'transversion')].plot(kind='scatter', x='Freq_source', y='Freq', ax=axes[i], color='#558B2F', label='synonymous transversion')
    except:
        pass
    group[(group.mtype == 'non-synonymous') & (group.transition_transversion == 'transition')].plot(kind='scatter', x='Freq_source', y='Freq', ax=axes[i], color='#FDA318', label='non-synonymous transition')
    try:
        group[(group.mtype == 'non-synonymous') & (group.transition_transversion == 'transversion')].plot(kind='scatter', x='Freq_source', y='Freq', ax=axes[i], color='#D51313', label='non-synonymous transversion')
    except:
        pass
    #group[(group.mtype == 'non-coding') & (group.transition_transversion == 'transversion')].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='black', marker='x', label='non-coding transversion')
    #group[group.mtype == 'frameshift'].plot(kind='scatter', x='Pos', y='Freq', ax=axes[i], color='grey', label='frameshift')
    axes[i].legend().remove()
    axes[i].set_ylim(-0.05, 1.05)
    #axes[i].set_xlim(0.005, 1.005)
    sns.regplot(x=group.Freq_source, y=group.Freq, ax=axes[i], scatter=False, color='grey')
    axes[i].set_title( infection_order.loc[infection_order.Sample == key, 'paper_name'].values[0], fontsize=14)
    axes[i].set_xlabel('')
    axes[i].set_ylabel('')
    i += 1
    c, pvalue = stats.pearsonr(group.Freq_source, group.Freq)
    correlations.append((infection_order.loc[infection_order.Sample == key, 'paper_name'].values[0], c, pvalue))
correlations = pd.DataFrame(correlations, columns=['Sample', 'correlation_coefficient', 'pvalue'])
axes[9].legend(loc='center', bbox_to_anchor=(1, -0.7), fontsize=16, ncol=2)
fig.text(0.5, 0.04, 'Mutation Frequency in Source', ha='center', va='center', fontsize=16)
fig.text(0.06, 0.5, 'Mutation Frequency in Sample', ha='center', va='center', rotation='vertical', fontsize=16)
fig.savefig('X:/volume2/noam/hcv/correlations.png' , bbox_inches='tight', dpi=800)
correlations.to_csv('X:/volume2/noam/hcv/correlations.csv', index=False)



fig, ax = plt.subplots(nrows=1, ncols=1)
correlations = correlations[correlations.Sample != 'Source']
correlations = pd.merge(correlations, infection_order[['infection_order', 'paper_name']], left_on='Sample', right_on='paper_name')
sns.regplot(correlations[correlations.Sample != 'Source'].infection_order, correlations[correlations.Sample != 'Source'].correlation_coefficient, ax=ax)
c, pvalue = stats.pearsonr(correlations[correlations.Sample != 'Source'].infection_order, correlations[correlations.Sample != 'Source'].correlation_coefficient)
ax.set_title('Order of Infection vs. Correlation Coefficient\n' + 'r = ' + str(round(c, 3)) + ', p-value = ' + str(round(pvalue, 3)))
fig.set_size_inches(5,3)
ax.set_xlabel('Order of infection')
ax.set_ylabel('Correlation Coefficient to Source')
fig.savefig('X:/volume2/noam/hcv/correlations_regression.png' , bbox_inches='tight', dpi=800)




compare_df = pd.merge(df, df[df.Sample=='HCV-PS2'], on=['Pos', 'Base', 'Ref'], suffixes=['', '_source'])
positions_to_look_at = compare_df[(compare_df.Ref != '-') & (compare_df.Read_count > 5) & (compare_df.Read_count_source > 5) & (compare_df.Ref != compare_df.Base) & (compare_df.Freq > 0.01)].Pos.drop_duplicates().tolist()
compare_df = compare_df[(compare_df.Pos.isin(positions_to_look_at)) & (compare_df.Ref != '-') & (compare_df.Read_count > 5) & (compare_df.Read_count_source > 5)]
correlations = []
for key, group in compare_df.groupby('Sample'):
    c, pvalue = stats.pearsonr(group.Freq, group.Freq_source)
    correlations.append((key, c, pvalue))
correlations = pd.DataFrame(correlations, columns=['Sample', 'correlation_coefficient', 'pvalue'])
correlations2 = pd.merge(correlations, infection_order, on='Sample', how='left')
correlations2 = correlations2.sort_values('correlation_coefficient', ascending=False)
correlations3 = correlations2[correlations2.infection_order.isna() == False]

fig, ax = plt.subplots(nrows=1, ncols=1)
slope, intercept, r_value, p_value, std_err = stats.linregress(correlations3.infection_order, correlations3.correlation_coefficient)
sns.regplot(x=correlations3.infection_order, y=correlations3.correlation_coefficient, ax = ax)
ax.set_title("y=%fx+%f, r2=%f, pvalue=%f" % (slope,intercept, r_value**2, p_value))





############################ corelation between time points
## color coded
fig, axes = plt.subplots(nrows=3, ncols=3)
fig.set_size_inches(10,10)
axes = axes.flatten()
fig.subplots_adjust(hspace=0.3)
i = 0
pair_order = [('HCV-P8','HCV-P8-1'),
              ('HCV-P8','HCV-P8-2'),
              ('HCV-P11','HCV-P11-1'),
              ('HCV-P7','HCV-P7-2'),
              ('HCV-P3','HCV-P3-1'),
              ('HCV-P5','HCV-P5-1'),
              ('HCV-P2','HCV-P2-1'),
              ('HCV-P6','HCV-P6-1'),
              ('HCV-P9','HCV-P9-2'),]

#compare_df = pd.merge(df, df[(df.Sample == 'HCV-PS2')], how='right', on=['Pos', 'Base','Ref','protein', 'mutation_type', 'aa', 'mtype', 'rb','transition_transversion'], suffixes=['', '_source'])
variants_to_look_at = df[(df.Base != df.Ref) & (df.Ref != '-') & (df.Base != '-') & (df.Read_count > 5) & (df.Freq > 0.01)][['Pos', 'Base']].drop_duplicates()
#compare_df = pd.merge(compare_df, variants_to_look_at, how='right', on=['Pos', 'Base'])
df = pd.merge(df, variants_to_look_at, how='right', on=['Pos', 'Base'])
df['full_mutation'] = df.Ref + df.Pos.astype(int).astype(str) + df.Base
correlations = []
for pair in pair_order:
    group = pd.merge(df[df.Sample == pair[0]], df[df.Sample == pair[1]], on=['Pos', 'Base','Ref','protein', 'mutation_type', 'aa', 'mtype', 'rb','transition_transversion', 'full_mutation'], suffixes=['_0', '_1'])
    #sns.regplot(x=group.Freq_0, y=group.Freq_1, ax=axes[i], scatter=False, color='red')
    group[~group.full_mutation.isin(strain_grey + strain_blue + ['A1479G'])].plot(kind='scatter', x='Freq_0', y='Freq_1', ax=axes[i], color='black', label='other')
    group[group.full_mutation.isin(strain_blue)].plot(kind='scatter', x='Freq_0', y='Freq_1', ax=axes[i], color='#468CC1', label='blue strain')
    group[group.full_mutation.isin(strain_grey)].plot(kind='scatter', x='Freq_0', y='Freq_1', ax=axes[i], color='grey', label='grey strain')
    group[group.full_mutation.isin(['G1479A'])].plot(kind='scatter', x='Freq_0', y='Freq_1', ax=axes[i], color='#B268DC', label='A1479G')
    axes[i].legend().remove()
    axes[i].set_ylim(-0.05, 1.05)
    axes[i].set_title( infection_order.loc[infection_order.Sample == pair[0], 'paper_name'].values[0] + ', time point ' + pair[1].split('-')[-1], fontsize=14)
    axes[i].set_xlabel('')
    axes[i].set_ylabel('')
    i += 1
    c, pvalue = stats.pearsonr(group.Freq_0, group.Freq_1)
    correlations.append((infection_order.loc[infection_order.Sample == pair[0], 'paper_name'].values[0] + ', time point ' + pair[1].split('-')[-1], c, pvalue))
correlations = pd.DataFrame(correlations, columns=['Sample', 'correlation_coefficient', 'pvalue'])
axes[7].legend(loc='center', bbox_to_anchor=(0.4, -0.55), fontsize=14, ncol=2)
fig.text(0.5, 0.07, 'Mutation Frequency in First Sample', ha='center', va='center', fontsize=16)
fig.text(0.06, 0.5, 'Mutation Frequency in Current Sample', ha='center', va='center', rotation='vertical', fontsize=16)
fig.savefig('X:/volume2/noam/hcv/time_points_correlations.png' , bbox_inches='tight', dpi=800)
#correlations.to_csv('X:/volume2/noam/hcv/time_points_correlations.csv', index=False)




def check_correlation(x_data, y_data):
   c, pvalue = stats.pearsonr(x_data, y_data)
   sns.regplot(x=x_data, y=y_data)
   print('correaltion: ', c, '\npvalue: ', pvalue)

#### save mutations over time with mutation types

df = pd.read_csv('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/all_freqs.csv')
#mutation_type = pd.read_csv('X:/volume2/noam/hcv/references/D90208.1_optimized4.mutation_type.peptides.csv')
mutation_type = pd.read_csv('X:/volume2/noam/hcv/references/D90208.1_optimized4.mutation_type.polypeptides.csv')
df = df[(df.Sample != 'HCV-P12')]
df = df[df.Pos.isin(range(738,2758))]
mutations_to_keep = df[(df.Ref != df.Base) & (df.Ref != '-') & (df.Freq > 0.01)][['Pos', 'Base']].drop_duplicates()
mutations_to_keep = pd.merge(mutations_to_keep, df, on=['Pos', 'Base'])
to_pivot = mutations_to_keep.pivot_table(values='Freq', index=['Pos', 'Ref', 'Base'], columns='Sample').reset_index()
to_pivot = pd.merge(to_pivot, mutation_type, how='left', on=['Pos', 'Base', 'Ref'])

to_pivot.loc[to_pivot['mutation_type'] == '[]', 'mtype'] = 'synonymous'
to_pivot.loc[((to_pivot['mutation_type'] != '[]') & ~(to_pivot['mutation_type'].isna())), 'mtype'] = 'non-synonymous'
#to_pivot.loc[((to_pivot['mutation_type'].isna()) & (to_pivot.Base != to_pivot.Ref)), 'mtype'] = 'non-coding'
to_pivot.loc[to_pivot['Base'] == '-', 'mtype'] = 'frameshift'

to_pivot = to_pivot[['Pos', 'Ref', 'Base', 'mutation_type', 'aa', 'protein', 'mtype', 'HCV-P1', 'HCV-P10', 'HCV-P11', 'HCV-P11-1',
       'HCV-P2', 'HCV-P2-1', 'HCV-P3', 'HCV-P3-1', 'HCV-P4', 'HCV-P5',
       'HCV-P5-1', 'HCV-P6', 'HCV-P6-1', 'HCV-P7', 'HCV-P7-2', 'HCV-P8',
       'HCV-P8-1', 'HCV-P8-2', 'HCV-P9', 'HCV-P9-2', 'HCV-PS2',]]
#to_pivot.to_excel('X:/volume2/noam/hcv/mutations_over_0.01.xlsx', index=False)
#to_pivot.to_excel('X:/volume2/noam/hcv/mutations_over_0.01.polyprotein.xlsx', index=False)


to_pivot['full_mutation'] = to_pivot.Ref + to_pivot.Pos.astype(int).astype(str) + to_pivot.Base
columns_order = ['full_mutation', 'Patient 1', 'Patient 2', 'Patient 3', 
                 'Patient 5', 'Patient 6', 'Patient 7', 'Patient 8',
                 'Patient 9', 'Patient 10', 'Patient 11', 'Patient 12', 'Source']
to_pivot = to_pivot.rename(columns=pd.Series(infection_order.paper_name.values,index=infection_order.Sample).to_dict())
to_pivot[to_pivot['Source'] >= 0.09][columns_order].round(4).sort_values('Source', ascending=False).to_excel('X:/volume2/noam/hcv/mutations_over_0.09.order_names.compact.xlsx', index=False)
to_pivot[to_pivot['Source'] >= 0.01][columns_order].round(4).sort_values('Source', ascending=False).to_excel('X:/volume2/noam/hcv/mutations_over_0.01.order_names.compact.xlsx', index=False)




# cluster
df = pd.read_csv('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/all_freqs.csv')
#mutation_type = pd.read_csv('X:/volume2/noam/hcv/references/D90208.1_optimized4.mutation_type.peptides.csv')
df = df[(df.Sample != 'HCV-P12')]
df = df[df.Pos.isin(range(738,2758))]
df['full_mutation'] = df.Ref + df.Pos.astype(int).astype(str) + df.Base
mutations_to_keep = df[(df.Ref != df.Base) & (df.Ref != '-') & (df.Freq > 0.09) & (df['Sample'] == 'HCV-PS2')][['full_mutation']].drop_duplicates()
mutations_to_keep = pd.merge(mutations_to_keep, df, on=['full_mutation'])
mutations_to_keep = mutations_to_keep[~mutations_to_keep.Sample.isin(['HCV-P2-1', 'HCV-P3-1', 'HCV-P5-1', 'HCV-P6-1', 'HCV-P7-2', 'HCV-P8-1', 'HCV-P8-2', 'HCV-P9-2', 'HCV-P11-1'])]
to_pivot = mutations_to_keep.pivot_table(values='Freq', index=['full_mutation'], columns='Sample')
#to_pivot = to_pivot.rename(columns=pd.Series(infection_order.paper_name.values,index=infection_order.Sample).to_dict())
sns.clustermap(to_pivot)









########## intrastrain correlation ############
strain_grey = ['T905C',
 'C2720T',
 'G1666A',
 'A1559G',
 'G1830A',
 'G1211A',
 'C2390T',
 'C2021T',
 'T1058C']
strain_blue = ['G1493A',
 'G1527A',
 'C2594T',
 'T1580C',
 'G1630A',
 'T1244C',
 'C1181T',
 'G2582A',
 'T1748C',
 'C1190T',
 'C2408T',
 'G2286A',
 'T1946C',
 'G1659A',
 'C1938T',
 'T2715C',
 'T2050C',
 'C2381T']
to_pivot = to_pivot.reset_index()
to_pivot = to_pivot.rename(columns={'HCV-PS2':'Source'})
strain_blue = to_pivot[to_pivot.full_mutation.isin(strain_blue)]
strain_grey = to_pivot[to_pivot.full_mutation.isin(strain_grey)]
#columns_order = ['Patient 1', 'Patient 2', 'Patient 3', 
#                 'Patient 5', 'Patient 6', 'Patient 7', 'Patient 8',
#                 'Patient 9', 'Patient 10', 'Patient 11', 'Patient 12', 'Source']
columns_order = ['Source',
 'HCV-P8',
 'HCV-P1',
 'HCV-P11',
 'HCV-P7',
 'HCV-P4',
 'HCV-P3',
 'HCV-P5',
 'HCV-P2',
 'HCV-P6',
 'HCV-P9',
 'HCV-P10',]
#strain_blue.loc[strain_blue['Source'] < 0.04, c] = 0
#strain_grey.loc[strain_grey['Source'] < 0.04, c] = 0


fig, axes = plt.subplots(nrows=3, ncols=4)
fig.set_size_inches(14,10)
axes = axes.flatten()
fig.subplots_adjust(hspace=0.4)
i = 0
correlations = []
for c in columns_order:
    #strain_grey.loc[strain_grey[c] < 0.04, c] = 0
    #strain_blue.loc[strain_blue[c] < 0.04, c] = 0
    slope, intercept, r_value, p_value, std_err = stats.linregress(strain_grey['Source'], strain_grey[c])
    correlations.append((c, 'grey', r_value, p_value, slope))
    sns.regplot(x=strain_grey['Source'], y=strain_grey[c], ax = axes[i], color='grey')
    title = c + '\ngrey r: ' + str(round(r_value, 3)) + ', p: ' + "{:.1e}".format(p_value) + ', slope: ' + str(round(slope, 3))
    slope, intercept, r_value, p_value, std_err = stats.linregress(strain_blue['Source'], strain_blue[c])
    correlations.append((c, 'blue', r_value, p_value, slope))
    sns.regplot(x=strain_blue['Source'], y=strain_blue[c], ax = axes[i], color='#468CC1')
    title +='\nblue r: ' + str(round(r_value, 3)) + ', p: ' + "{:.1e}".format(p_value) + ', slope: ' + str(round(slope, 3))
    axes[i].legend().remove()
    axes[i].set_title(title, fontsize=10)
    axes[i].set_xlabel('')
    axes[i].set_ylabel('')
    axes[i].set_ylim(0,1)
    i += 1
correlations = pd.DataFrame(correlations, columns=['Sample', 'strain', 'r_value', 'pvalue', 'slope'])
axes[9].legend(loc='center', bbox_to_anchor=(1, -0.7), fontsize=16, ncol=2)
fig.text(0.5, 0.04, 'Mutation Frequency in Source', ha='center', va='center', fontsize=16)
fig.text(0.06, 0.5, 'Mutation Frequency in Sample', ha='center', va='center', rotation='vertical', fontsize=16)
fig.savefig('X:/volume2/noam/hcv/intra_strain_correlations.png' , bbox_inches='tight', dpi=800)
correlations.to_csv('X:/volume2/noam/hcv/intra_strain_correlations.csv', index=False)




##### trying to rebuild source haplotypes using patients
df = unite_all_freq_files('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/freqs/', 'Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/')
df['mutation'] = df['Ref'] + df['Pos'].astype(str) + df['Base']

haplotypes = {}
haplotypes['HCV-P10'] = df[(df.Sample == 'HCV-P10') & (df.Ref != df.Base) & (df.Ref != '-') & (df.Pos.isin(range(738,2757))) & (df.Freq > 0.8)][['Pos', 'Ref', 'Base', 'mutation']]
haplotypes['HCV-P6'] = df[(df.Sample == 'HCV-P6') & (df.Ref != df.Base) & (df.Ref != '-') & (df.Pos.isin(range(738,2757))) & (df.Freq > 0.8)][['Pos', 'Ref', 'Base', 'mutation']]

pd.merge(pd.merge(haplotypes['HCV-P10'], haplotypes['HCV-P6'], on=['Pos', 'Ref', 'Base', 'mutation']), df[df.Sample == 'HCV-PS2']).Freq.max()
h = pd.merge(haplotypes['HCV-P10'], haplotypes['HCV-P6'], on=['Pos', 'Ref', 'Base', 'mutation']).mutation.tolist()
df_source = df[(df.Sample == 'HCV-PS2') & (df.mutation.isin())]









# create haplotypes from patient samples, than try to explain the most of the 
# source sample using these haplotypes (try to use them as building blocks )

def break_up_source_into_haplotypes(source_freq, haplotype_list):
    source_freq = pd.read_csv(source_freq, '\t')
    source_freq = source_freq[(source_freq.Pos.isin(range(738,2757))) & (source_freq.Read_count > 100)]
    source_freq['mutation'] = source_freq.Ref + source_freq.Pos.astype(str) + source_freq.Base
    for haplotype in haplotype_list: # does order of haplotypes matter? 
        # find maximum value for strain, both variants and reference
        haplotype_positions = [float(i[1:-1]) for i in haplotype]
        min1 = source_freq[source_freq.mutation.isin(haplotype)].Freq.min()
        min2 = source_freq[~(source_freq.Pos.isin(haplotype_positions)) & (source_freq.Ref == source_freq.Base) & (source_freq.Ref != '-')].Freq.min()
        haplotype_max_freq = min([min1, min2])
        # decrease haplotype_max_freq from frequencies, treating frequencies as counts here, not touching read_count
        source_freq['Freq'] = np.where(source_freq.mutation.isin(haplotype), source_freq.Freq - haplotype_max_freq, source_freq.Freq)
        source_freq['Freq'] = np.where((~(source_freq.Pos.isin(haplotype_positions)) & (source_freq.Ref == source_freq.Base) & (source_freq.Ref != '-')), source_freq.Freq - haplotype_max_freq, source_freq.Freq)
    return source_freq # problem - decreases frequency twice from positions with the haplotype

def try_haplotype_combos():
    pass

haplotypes = {}
haplotypes['HCV-P10'] = df[(df.Sample == 'HCV-P10') & (df.Ref != df.Base) & (df.Ref != '-') & (df.Pos.isin(range(738,2757))) & (df.Freq > 0.8)].mutation.tolist()
haplotypes['HCV-P6'] = df[(df.Sample == 'HCV-P6') & (df.Ref != df.Base) & (df.Ref != '-') & (df.Pos.isin(range(738,2757))) & (df.Freq > 0.8)].mutation.tolist()
haplotypes['HCV-P6-edited'] = df[(df.Sample == 'HCV-P6') & (df.Ref != df.Base) & (df.Ref != '-') & (df.Pos.isin(range(738,2757))) & (df.Freq > 0.8) & (df.Pos != 1799)].mutation.tolist()


build_source_from_haplotypes('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/HCV-PS2/HCV-PS2.freqs', list(haplotypes.values()))





##### accungs haplotyping
com = obtain_comutations(load_file('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_haplotyping/HCV-PS2/all_results.csv', 'HCV-PS2'), distance=1.5)
com['Stretch_str'] = 's' + com.Stretch.astype(str)
f = pd.read_csv('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/HCV-PS2/HCV-PS2.freqs', '\t')
f = f[(f.Ref != f.Base) & (f.Ref != '-')]
f = f[f.groupby(['Pos'])['Freq'].transform(max) == f['Freq']]
f = f[['Pos', 'Ref', 'Base', 'Freq']]
f.Pos = f.Pos.astype(int)
comm = pd.merge(com, f, right_on='Pos', left_on='Pos1', how='left')
comm = comm.sort_values('Freq_x')

fig, ax = plt.subplots(nrows=1, ncols=1)
sns.scatterplot(data=comm, x='Pos1', y='Freq_y',hue='Stretch_str', legend=False, palette = sns.color_palette("bright", 47), ax=ax)

    
    
##### diffferent blast ids
def plot_variants(freq_file):
    df = pd.read_csv(freq_file, '\t')
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.set_yscale('log')
    #ax.set_ylim(10**-3, 10**0)
    df[(df.Base != df.Ref) & (df.Ref != '-') & (df.Pos.isin(range(738,2757)))].plot(x='Pos', y='Read_count', kind='scatter', ax=ax)
    return

def bar_plot_sample

plot_variants('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_blast95/HCV-PS2/HCV-PS2.freqs')
plot_variants('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_blast90/HCV-PS2/HCV-PS2.freqs')
plot_variants('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_blast80/HCV-PS2/HCV-PS2.freqs')
plot_variants('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_blast75/HCV-PS2/HCV-PS2.freqs')

dfs = []
for i in [95,90, 80,75]:
    d = pd.read_csv('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_blast%d/HCV-PS2/HCV-PS2.freqs' % (i), '\t')
    d['blast_id'] = i
    dfs.append(d)
dfs = pd.concat(dfs)

compare_blast = compare_positions_between_freqs({'b95':'Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_blast95/HCV-PS2/HCV-PS2.freqs', 'b90':'Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_blast90/HCV-PS2/HCV-PS2.freqs', 'b85':'Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/HCV-PS2/HCV-PS2.freqs', 'b80':'Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_blast80/HCV-PS2/HCV-PS2.freqs', 'b75':'Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_blast75/HCV-PS2/HCV-PS2.freqs'})
# compare read counts
compare_blast[(compare_blast.Ref == compare_blast.Base) & (compare_blast.Pos.isin(range(738,2757)))].plot(x='Read_count_b95', y='Read_count_b85', kind='scatter')
# compare variant freqs
fig, ax = plt.subplots(nrows=1, ncols=1)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(10**-3, 10**0)
ax.set_xlim(10**-3, 10**0)
compare_blast[(compare_blast.Ref != compare_blast.Base) & (compare_blast.Ref != '-') &  (compare_blast.Pos.isin(range(738,2757)))].plot(x='Freq_b95', y='Freq_b85', kind='scatter', ax=ax)





#
df = pd.read_csv('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/all_freqs.csv')
df['mutation'] = df['Ref'] + df['Pos'].astype(str) + df['Base']

df = estimate_insertion_freq(df, extra_columns=['Sample'])


##### long indels




def get_aligned_more_than_once(tmp_dir):
    files = [os.path.join(tmp_dir, f) for f in os.listdir(tmp_dir) if f.endswith('.blast')]
    dfs = []
    for f in files:
        print(f)
        b = blast_to_df(f)
        dfs.append(b)
    blast_df = pd.concat(dfs)
    alignment_count = blast_df.groupby('read').strand.count().reset_index()
    aligned_more_than_once = blast_df[blast_df.read.isin(alignment_count[alignment_count.strand > 2].read.tolist())]
    print(len(blast_df.read.unique()))
    print(len(aligned_more_than_once.read.unique()))
    print(100 * len(aligned_more_than_once.read.unique()) / len(blast_df.read.unique()))
    return aligned_more_than_once


blast = get_aligned_more_than_once('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/HCV-P1-2/tmp')
blast = blast.sort_values('start_ref')
blast = blast[blast.strand == 'plus']
blast['seq'] = blast.groupby('read')['start_ref'].rank()

fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches(10,16)
blast = blast.sort_values(['seq', 'start_ref'])
sns.swarmplot(x='start_ref', y='read', data=blast, hue='seq', ax=ax, alpha=0.5)
sns.swarmplot(x='end_ref', y='read', data=blast, hue='seq', ax=ax, alpha=0.5)
ax.set_yticks([])
ax.set_title('HCV-P11')




#### get mutation types
mutation_type = pd.read_csv('X:/volume2/noam/hcv/references/D90208.1_optimized4.mutation_type.csv')
mutation_type = mutation_type[['Pos', 'Ref', 'Base', 'protein', 'mutation_type']]
mutation_type['aa'] = (((mutation_type.Pos - 330) / 3) + 1).astype(int).astype(float)
proteins = create_protein_dataframe('X:/volume2/noam/hcv/references/BAA14233.1_features_table.txt')
proteins2 = []
for row in proteins.itertuples():
    for i in range(row.start, row.end+1):
        proteins2.append([i, row.protein])
proteins2 = pd.DataFrame(proteins2, columns=['aa', 'protein'])   
mutation_type = pd.merge(mutation_type.drop(columns=['protein']), proteins2, on='aa')
mutation_type.to_csv('X:/volume2/noam/hcv/references/D90208.1_optimized4.mutation_type.peptides.csv', index=False)

def create_protein_dataframe(ncbi_feature_table):
    genes = []
    with open(ncbi_feature_table) as f:
        fi = f.read()
    # get cds features
    features = [s for s in re.compile('^(?![\t]).+?^\t\t\t.+?^(?![\t])', re.MULTILINE | re.DOTALL).findall(fi) if 'Region' in s]
    fi = ''.join(features)
    for row in fi.splitlines():
        if not row.startswith('\t') and not row.startswith('>') and row != '':
            try:
                protein = fi.split(row)[1].split('region\t')[1].split('\n')[0]
                #protein = fi.split(row)[1].split('gene\t')[1].split('\n')[0]
                start = int(row.split('\t')[0].replace('<', ''))
                end = int(row.split('\t')[1].replace('<', ''))
                complement = start > end
                genes.append((protein, start, end, complement))
            except:
                print('problem with row: ', row)
        fi = fi.replace(row, '', 1)
    gene_df = pd.DataFrame(genes, columns=['protein', 'start', 'end', 'complement'])
    return gene_df

### get mutation types with non-protein areas of polyprotein
mutation_type = pd.read_csv('X:/volume2/noam/hcv/references/D90208.1_optimized4.mutation_type.csv')
mutation_type = mutation_type[['Pos', 'Ref', 'Base', 'protein', 'mutation_type']]
mutation_type['aa'] = (((mutation_type.Pos - 330) / 3) + 1).astype(int).astype(float)
proteins = create_protein_dataframe('X:/volume2/noam/hcv/references/BAA14233.1_features_table.txt')
proteins2 = []
for row in proteins.itertuples():
    for i in range(row.start, row.end+1):
        proteins2.append([i, row.protein])
proteins2 = pd.DataFrame(proteins2, columns=['aa', 'protein'])   
mutation_type = pd.merge(mutation_type.drop(columns=['protein']), proteins2, on='aa', how='left')
mutation_type.to_csv('X:/volume2/noam/hcv/references/D90208.1_optimized4.mutation_type.polypeptides.csv', index=False)












############### read mapping percent
with open('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/mapping_perecent.txt') as f:
    f = f.read()
reads_list = []
mapped_list = []
for row in f.split('\n'):
    sample = row.split('/')[0]
    if 'reads:' in row:
        reads = int(row.split(' ')[-1])
        reads_list.append((sample, reads))
    if 'reference:' in row:
        mapped = int(row.split(' ')[-1])
        mapped_list.append((sample, mapped))
mapped_percent = pd.merge(pd.DataFrame(mapped_list, columns=['sample', 'mapped']), pd.DataFrame(reads_list, columns=['sample', 'reads']), on='sample')
mapped_percent['percent'] = mapped_percent.mapped * 100 / mapped_percent.reads      
   mapped_percent[mapped_percent['sample'] != 'HCV-P12'].sort_values('percent').to_clipboard(index=False)     

# get unaligned reads
   
def get_aligned_read_ids(blast_dir, out_file):
    '''
    Gets a fastq file and a directory with blast file(s) (the tmp dir
    created by pipeline can be used here) and writes to out_file a fastq 
    file only with the unaligned reads.
    '''
    blasts = []
    for b in [blast_dir + '/' + f for f in os.listdir(blast_dir) if f.endswith('.blast')]:
        blasts.append(blast_to_df(b))
    blast = pd.concat(blasts)
    blast.read.drop_duplicates().to_csv(out_file, index=False, header=False)
    return

get_aligned_read_ids('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/HCV-P9/tmp', 'Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/unaligned/HCV-P9/HCV-P9_S9_L001_001.fastq.aligned.txt')




######  phylogeny
for f in ['Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_individual_concensus/' + f for f in os.listdir('Z:/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4_individual_concensus')]:
    with open(f) as text:
        text = text.read()
    text = text.replace('D90208.1 Hepatitis C virus ORF gene, complete cds', f.split('/')[-1])
    with open(f, 'w') as new_text:
        new_text.write(text)