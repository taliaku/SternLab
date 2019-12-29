# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 21:17:36 2018

@author: Noam
"""
import re
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np
import os
import sys
sys.path.append('/sternadi/home/volume2/noam/SternLab')
sys.path.append(r'X:\volume2\noam\Sternlab')
from blast_utilities import blast_to_mutations_list, blast_to_df
from freqs_utilities import estimate_insertion_freq
import matplotlib.pyplot as plt

def big_freqs_file():
    dfs = []
    freqs_list = ['Z:/volume1/noam/cmv/cmv_medicine/1_pipeline_output/1.freqs', 
                  'Z:/volume1/noam/cmv/cmv_medicine/2_pipeline_output/2.freqs', 
                  'Z:/volume1/noam/cmv/cmv_medicine/3_pipeline_output/3.freqs',
                  'Z:/volume1/noam/cmv/cmv_medicine/4_pipeline_output/4.freqs',
                  'Z:/volume1/noam/cmv/cmv_medicine/5_pipeline_output/5.freqs',
                  'Z:/volume1/noam/cmv/cmv_medicine/6_pipeline_output/6.freqs',
                  'Z:/volume1/noam/cmv/cmv_medicine/7_pipeline_output/7.freqs',
                  'Z:/volume1/noam/cmv/cmv_medicine/8_pipeline_output/8.freqs',
                  'Z:/volume1/noam/cmv/cmv_no_medicine/P16LT/P16LT.freqs',
                  'Z:/volume1/noam/cmv/cmv_no_medicine/P16ND/P16ND.freqs',
                  'Z:/volume1/noam/cmv/cmv_no_medicine/P28ND/P28ND.freqs',
                  ]
    for f in freqs_list:
        df = pd.read_csv(f, sep='\t')
        df['file'] = f.split('/')[-1].split('.')[0]
        dfs.append(df)
    df = pd.concat(dfs)
    return df

df = big_freqs_file()
df.to_csv('Z:/volume1/noam/cmv/all_freqs.with_without_medicine.csv', index=False)





def compare_freqs(df, mutation_type_df, file1, file2, output_dir):
    # snps and deletions
    dfx = df[df.file == file1]
    dfy = df[df.file == file2]
    new_df = pd.merge(dfx, dfy, on=['Pos', 'Base', 'Ref'], suffixes=['_' + str(file1), '_' + str(file2)])
    deletions = new_df[(((new_df['Freq_' + str(file1)] - new_df['Freq_' + str(file2)]).abs() > 0.25) & (new_df.Base != new_df.Ref) & (new_df.Ref != '-') & (new_df.Base == '-'))]
    mutations = new_df[(((new_df['Freq_' + str(file1)] - new_df['Freq_' + str(file2)]).abs() > 0.25) & (new_df.Base != new_df.Ref) & (new_df.Ref != '-') & (new_df.Base != '-'))]
    mutations = pd.merge(mutations, f_df, on=['Pos', 'Ref', 'Base'], how='left')
    #mutations = mutations[['Ref', 'Pos', 'Base', 'Freq_' + str(file1), 'Freq_' + str(file2), 'Read_count_' + str(file1), 'Read_count_' + str(file2), 'file_' + str(file1), 'file_' + str(file2), 'protein', 'mutation_type']]
    mutations = mutations[['Ref', 'Pos', 'Base', 'Freq_' + str(file1), 'Freq_' + str(file2), 'Read_count_' + str(file1), 'Read_count_' + str(file2), 'protein', 'mutation_type']]
    mutations.rename(columns={'mutation_type':'protein_mutation'}, inplace=True)
    conditions = [(mutations.protein_mutation == '[]'), (mutations.protein_mutation.isna()), ((mutations.protein_mutation != '[]') & ~(mutations.protein_mutation.isna()))]
    choices = ['synonymous', 'non_coding', 'non_synonymous']
    mutations['mutation_type'] = np.select(conditions, choices)
    
    deletions = pd.merge(deletions, mutation_type_df[['Pos', 'protein']].drop_duplicates(), on='Pos', how='left')
    #deletions = deletions[['Ref', 'Pos', 'Base', 'Freq_' + str(file1), 'Freq_' + str(file2), 'Read_count_' + str(file1), 'Read_count_' + str(file2), 'file_' + str(file1), 'file_' + str(file2), 'protein']]
    deletions = deletions[['Ref', 'Pos', 'Base', 'Freq_' + str(file1), 'Freq_' + str(file2), 'Read_count_' + str(file1), 'Read_count_' + str(file2), 'protein']]
    # insertions
    #insertions = new_df[((new_df.Ref == '-') & ((new_df.Read_count_x - new_df.Read_count_y).abs() > 20))]    
    # positions appearing in one and not the other
    #weird_pos = pd.concat([dfx[(~(dfx.Pos.isin(dfy.Pos.tolist())))], dfy[(~(dfy.Pos.isin(dfx.Pos.tolist())))]])
    mutations.to_csv(output_dir + 'mutations_' + str(file1) + '_' + str(file2) + '.csv', index=False)
    deletions.to_csv(output_dir + 'deletions_' + str(file1) + '_' + str(file2) + '.csv', index=False)
    #insertions.to_csv(output_dir + 'insertions_' + file1 + '_' + file2 + '.csv', index=False)
    #weird_pos.to_csv(output_dir + 'weird_positions_' + file1 + '_' + file2 + '.csv', index=False)
    return
 
for i in range(1,9):
    for j in range(i+1,9):
        if i != j:
            compare_freqs(freqs, mutation_type, i, j, 'Z:/volume1/noam/cmv/comparisons2/')
    







#####
            
to_pivot = pd.merge(freqs, non_synonymous, how='right', on=['Ref', 'Pos', 'Base'])
to_pivot = pd.merge(to_pivot, mutation_type[['Base', 'Ref', 'Pos', 'protein', 'mutation_type']], how='left', on=['Ref', 'Pos', 'Base'])
to_pivot['mutation'] = to_pivot.Ref + to_pivot.Pos.astype(str) + to_pivot.Base
to_pivot = to_pivot.pivot_table(values='Freq', index=['protein', 'mutation'], columns='file')
to_pivot = to_pivot[[3,4,7]]
to_pivot = to_pivot[~(to_pivot.isnull().any(axis=1))]


## mutations
t = to_pivot[(to_pivot[7] > 0.8)].reset_index()
conditions = [(t[3] < 0.5), (t[3] >= 0.5) & (t[3] <= 0.7), (t[3] > 0.7) & (t[3] <= 0.8), (t[3] > 0.8)]
choices = ['less than 0.5', 'between 0.5 and 0.7', 'between 0.7 and 0.8', 'more than 0.8']
t['3_grouped'] = np.select(conditions, choices)

# deletions
freqs[(freqs.Pos.isin(freqs[(freqs.Base == '-') & (freqs.Freq > 0.5) & (freqs.Read_count > 5)].Pos.drop_duplicates().tolist())) & (freqs.Base == '-')].pivot_table(values='Freq', index='Pos', columns='file')


# mapped and unmapped
b = freqs[(freqs.Base == freqs.Ref)].pivot_table(values='Read_count', index='Pos', columns='file')
a = freqs[(freqs.Base == freqs.Ref) & (freqs.Ref != '-')].pivot_table(values='Read_count', index='mutation', columns='file')
a = a[(a).isnull().any(axis=1)]
a = a[a.max(axis=1) > 1]

c = pd.merge(a, mutation_type[['Pos', 'protein']].drop_duplicates(), left_on='mutation', right_on='Pos')


d = pd.merge(c.groupby('protein').Pos.count().reset_index().rename(columns={'Pos':'disappearance_count'}), mutation_type[['protein', 'Pos']].drop_duplicates().groupby('protein').count().reset_index().rename(columns={'Pos':'protein_length'}), on='protein')

# get all mutations for protein RL1
pd.merge(mutation_type[mutation_type.protein =='RL1'], non_synonymous, on=['Pos', 'Ref', 'Base']).sort_values('Pos')
# with frequencies from p16
pd.merge(pd.merge(mutation_type[mutation_type.protein =='RL1'], non_synonymous, on=['Pos', 'Ref', 'Base']).sort_values('Pos'), freqs[freqs.file==3], on=['Pos', 'Ref', 'Base'])







######## controls
nd_p16 = pd.read_csv('Z:/volume1/noam/cmv/cmv_no_medicine/P16ND/P16ND.freqs', '\t')
nd_p28 = pd.read_csv('Z:/volume1/noam/cmv/cmv_no_medicine/P28ND/P28ND.freqs', '\t')
lt_p16 = pd.read_csv('Z:/volume1/noam/cmv/cmv_no_medicine/P16Lt/P16LT.freqs', '\t')

def estimate_insertion_freq(df):
    #read_counts = df[(df.Ref != '-')][['Pos', 'Read_count']].drop_duplicates()
    read_counts = df[(df.Ref != '-')][['Pos', 'Read_count', 'file']].drop_duplicates()
    read_counts.rename(columns={'Read_count':'estimated_read_count', 'Pos':'rounded_pos'}, inplace=True)
    insertions = df[(df.Ref == '-')]
    not_insertions = df[(df.Ref != '-')]
    insertions['rounded_pos'] = insertions.Pos.astype(int).astype(float)
    insertions = pd.merge(insertions, read_counts, how='left', on=['rounded_pos', 'file'])
    insertions['estimated_freq'] = insertions.Freq * insertions.Read_count / insertions.estimated_read_count
    df = pd.concat([insertions, not_insertions])
    return df.sort_values(['Pos'])

a = pd.read_csv('Z:/volume1/noam/cmv/cmv_medicine/all_freqs.csv')
a['file'] = a.file.astype(str)

nd_p16['file'] = 'nd_p16'
nd_p28['file'] = 'nd_p28'
lt_p16['file'] = 'lt_p16'
a = pd.concat([a, nd_p16, nd_p28, lt_p16])
a = estimate_insertion_freq(a)
a['full_mutation'] = a.Ref + a.Pos.astype(str) + a.Base

a.to_csv('Z:/volume1/noam/cmv/all_freqs.csv', index=False)



######

def stats(df):
    df = df[df.Read_count > 5]
    substitutions = df[(df.Base != df.Ref) & (df.Ref != '-') & (df.Base != '-') & (df.Freq > 0.2)]
    deletions = df[(df.Base != df.Ref) & (df.Ref != '-') & (df.Base == '-') & (df.Freq > 0.2)]
    insertions = df[(df.Base != df.Ref) & (df.estimated_freq > 0.2)]
    return len(substitutions), len(deletions), len(insertions)

# substitutions in common, inner merge
# deletions in common, inner merge
# insertions in common, inner merge
pd.merge(a[(a.Read_count > 5) & (a.Base != a.Ref) & (a.estimated_freq > 0.2) & (a.file.str.contains('nd_p16'))][['Ref', 'Pos', 'Base']].drop_duplicates(), a[(a.Read_count > 5) & (a.Base != a.Ref) & (a.estimated_freq > 0.2) & (a.file.str.contains('lt_p16'))][['Ref', 'Pos', 'Base']].drop_duplicates(), how='inner', on=['Ref', 'Pos', 'Base'])    

######

a = pd.read_csv('Z:/volume1/noam/cmv/all_freqs.csv')
mutation_type = pd.read_csv('Z:/volume1/noam/cmv/mutation_type_position_df.csv')

# mutations exceeding 0.2 in one of the three samples and merging them with mutation type
pd.merge(a[(a.file.str.contains('p')) & (a.Base != a.Ref) & (a.Ref != '-') & (a.Base != '-') & (a.Freq > 0.2)][['Pos', 'Ref', 'Base', 'file']].groupby(['Pos', 'Base', 'Ref']).file.apply(','.join).reset_index(), mutation_type[mutation_type.mutation_type != '[]'][['Pos', 'Ref', 'Base', 'protein', 'mutation_type']], how='left', on=['Pos', 'Ref', 'Base'])

b = a[(a.Base != a.Ref) & (a.Ref != '-') & (a.Read_count > 5)].pivot_table(values='Freq', index='full_mutation', columns='file')

b[(b['nd_p16'] < 0.2) & (b.nd_p28 < 0.2) & (b.lt_p16 > 0.2)]


to_pivot = pd.merge(a, non_synonymous, how='right', on=['Ref', 'Pos', 'Base'])
to_pivot = pd.merge(to_pivot, mutation_type[['Base', 'Ref', 'Pos', 'protein', 'mutation_type']], how='left', on=['Ref', 'Pos', 'Base'])
to_pivot['mutation'] = to_pivot.Ref + to_pivot.Pos.astype(str) + to_pivot.Base
to_pivot = to_pivot.pivot_table(values='Freq', index=['protein', 'mutation'], columns='file')
to_pivot = to_pivot[[3,4,7]]
to_pivot = to_pivot[~(to_pivot.isnull().any(axis=1))]


## mutations
t = to_pivot[(to_pivot[7] > 0.8)].reset_index()
conditions = [(t[3] < 0.5), (t[3] >= 0.5) & (t[3] <= 0.7), (t[3] > 0.7) & (t[3] <= 0.8), (t[3] > 0.8)]
choices = ['less than 0.5', 'between 0.5 and 0.7', 'between 0.7 and 0.8', 'more than 0.8']
t['3_grouped'] = np.select(conditions, choices)






###### context - didn't find context dependency

refs = a[a.Ref != '-'][['Pos', 'Ref']].drop_duplicates().sort_values('Pos')
refs['ref_previous_1'] = refs.Ref.shift(1)
refs['ref_previous_2'] = refs.Ref.shift(2)
refs['ref_previous_3'] = refs.Ref.shift(3)
refs['ref_previous_4'] = refs.Ref.shift(4)
refs['ref_next_1'] = refs.Ref.shift(-1)
refs['ref_next_2'] = refs.Ref.shift(-2)
refs['ref_next_3'] = refs.Ref.shift(-3)










################## ad169 samples
p28 = pd.read_csv('Z:/volume1/noam/cmv/cmv_AD169/AD169P28_pipeline/AD169P28.freqs', '\t')
wt = pd.read_csv('Z:/volume1/noam/cmv/cmv_AD169/AD169wt_pipeline/AD169wt.freqs', '\t')
p28['Sample'] = 'p28'
wt['Sample'] = 'wt'
a = pd.concat([p28,wt])
a = estimate_insertion_freq(a, extra_columns=['Sample'])
a['full_mutation'] = a.Ref + a.Pos.astype(str) + a.Base

merged = pd.merge(wt, p28, suffixes=['_wt', '_p28'], on=['Pos', 'Ref', 'Base'], how='outer')
merged['full_mutation'] = merged.Ref + merged.Pos.astype(str) + merged.Base

fig, ax = plt.subplots(nrows=1, ncols=1)
merged[(merged.Base != merged.Ref) & (merged.Ref != '-') & (merged.Read_count_wt > 5) & (merged.Read_count_p28 > 5)].plot(x='Freq_wt', y='Freq_p28', kind='scatter', ax=ax)
fig.set_size_inches(5,3)
ax.set_xlabel('Frequency WT')
ax.set_ylabel('Frequency p28')
ax.set_title('Variant Frequency - WT vs p28')


mutation_type = pd.read_csv('X:/volume2/noam/cmv/references/FJ527563.1.mutation_type.csv')
mutations_in_ad169 = a[(a.Read_count > 5) & (a.Ref != a.Base) & (a.Ref != '-') & (a.Freq > 0.2)].full_mutation.drop_duplicates().tolist()
merged[merged.full_mutation.isin(mutations_in_ad169)]




########### tb40 samples
df = pd.read_csv('Z:/volume1/noam/cmv/all_freqs.with_without_medicine.csv')
df['full_mutation'] = df.Ref + df.Pos.astype(str) + df.Base
mutation_type = pd.read_csv('Z:/volume1/noam/cmv/mutation_type_position_df.csv')
mutation_type = mutation_type.drop(columns=['Unnamed: 0'])


## look at UL54 mutations
ul_54 = pd.merge(df, mutation_type[(mutation_type.protein == 'UL54')][['Pos', 'Ref', 'Base', 'mutation_type']], on=['Pos', 'Ref', 'Base'])
conditions = [ul_54.mutation_type == '[]', ul_54.mutation_type != '[]']
choices = ['syn', 'nonsyn']
ul_54['mtype'] = np.select(conditions, choices)

# plot mutation scatter all samples
fig, axes = plt.subplots(nrows=4, ncols=3,figsize=(8, 10))
axes = axes.flatten()
fig.subplots_adjust(hspace=0.5)
for (i,(key, group)) in enumerate(ul_54.groupby('file')):
    group = group[(group.Ref != group.Base) & (group.Ref != '-') & (group.Read_count >5)]
    group[group.mtype == 'syn'].plot(x='Pos', y='Freq', ax=axes[i], color='blue', kind='scatter')
    group[group.mtype == 'nonsyn'].plot(x='Pos', y='Freq', ax=axes[i], color='red', kind='scatter')
    axes[i].set_title(key)
    axes[i].set_ylim(-0.1,1.1)
    axes[i].set_xlabel('')
    axes[i].set_ylabel('')
fig.text(0.5, 0.04, 'Position (bp)', ha='center')
fig.text(0.04, 0.5, 'Frequency', va='center', rotation='vertical')
fig.savefig('Z:/volume1/noam/cmv/UL_54/syn_not_syn_mutations.csv',bbox_layout='tight', dpi=800)

# get excel of mutation over time
to_pivot = ul_54[ul_54.full_mutation.isin(ul_54[(ul_54.Ref != ul_54.Base) & (ul_54.Ref != '-') & (ul_54.Freq > 0.2) & (ul_54.Read_count > 5)].full_mutation.tolist())]
to_pivot.pivot_table(values='Freq', index=['full_mutation', 'mutation_type'], columns='file').to_excel('Z:/volume1/noam/cmv/UL_54/mutation_over_0.2.xlsx')

# check that all positions were seqeunced
ul_54[(ul_54.Ref != '-') & (ul_54.Read_count > 5)][['Pos', 'file']].drop_duplicates().groupby('file').count()





###### compare 7+8
df1 = df[df.file.isin(['7','8'])].copy()
merged = pd.merge(df1[df1.file == '7'], df1[df1.file == '8'], suffixes=['_7', '_8'], on=['Pos', 'Ref', 'Base'], how='outer')
fig, ax = plt.subplots(nrows=1, ncols=1)
merged[(merged.Base != merged.Ref) & (merged.Ref != '-') & (merged.Read_count_7 > 5) & (merged.Read_count_8 > 5)].plot(x='Freq_7', y='Freq_8', kind='scatter', ax=ax)
fig.set_size_inches(5,3)
ax.set_xlabel('Frequency 7')
ax.set_ylabel('Frequency 8')
ax.set_title('Variant Frequency - 7 vs 8')

# only mutations that show change in frequency
changed = merged[(merged.Base != merged.Ref) & (merged.Ref != '-') & (merged.Freq_7 > 0.5) & (merged.Freq_8 < 0.2) & (merged.Read_count_7 > 5) & (merged.Read_count_8 > 5)]
changed = pd.merge(changed, mutation_type[['Pos', 'Ref', 'Base', 'mutation_type', 'protein']], on=['Pos', 'Ref', 'Base'], how='left')


to_pivot = df[df.full_mutation.isin(changed.full_mutation_7)]
to_pivot = pd.merge(to_pivot, mutation_type[['Base', 'Ref', 'Pos', 'protein', 'mutation_type']], how='left', on=['Ref', 'Pos', 'Base'])
to_pivot['mutation'] = to_pivot.Ref + to_pivot.Pos.astype(str) + to_pivot.Base
to_pivot = to_pivot.pivot_table(values='Freq', index=['protein', 'mutation', 'mutation_type'], columns='file')