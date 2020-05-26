
import pandas as pd
import json
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
from tqdm import tqdm
import numpy as np
import math

#
##########
df1 = pd.read_csv('Z:/volume1/noam/covid_data/coronaTech1_20200415/python_pipeline_x1_c0_v10-9_after_ptrimmer/freqs/technion1_1e-9_mutations_all.csv')
df2 = pd.read_csv('Z:/volume1/noam/covid_data/coronaTech2_20200427/python_pipeline_x1_c0_v10-9_after_ptrimmer/freqs/technion2_1e-9_mutations_all.csv')
df3 = pd.read_csv('Z:/volume1/noam/covid_data/coronaTech3_20200504/python_pipeline_x1_c0_v10-9_after_ptrimmer/freqs/technion3_1e-9_mutations_all.csv')
sh16 = pd.read_csv('Z:/volume1/noam/covid_data/old/SH16.fasta.mutations_list.csv')
sh16 = sh16.rename(columns={'position':'ref_position', 'ref':'ref_base', 'file':'sample'})
sh16['sample'] = 'SH16'

df = pd.concat([df1, df2, df3, sh16])
df['sample'] = df['sample'].astype(str)
df = df.drop(columns=['Unnamed: 0'])

include_list = pd.read_csv('Z:/volume1/noam/covid_data/include.txt', header=None).rename(columns={0:'sample'})['sample'].tolist()
include_list = [i.split('/')[1] for i in include_list]
df = df[df['sample'].isin(include_list)]

df.to_csv('Z:/volume1/noam/covid_data/all_mutations_in_included_cons.csv', index=False)

# read data
df = pd.read_csv('Z:/volume1/noam/covid_data/all_mutations_in_included_cons.csv')
df = pd.read_csv('/Volumes/STERNADILABTEMP$/volume1/noam/covid_data/all_mutations_in_included_cons.csv')


## snps pre sample graph
snps_per_sample = df[(df.base != '-') & (df.ref_base != '-')].groupby('sample').base.count().reset_index().rename(columns={'base':'snps'})
samples_with_snps = snps_per_sample.groupby('snps').count().reset_index()
samples_with_snps = samples_with_snps.append({'snps':14, 'sample':0}, ignore_index=True).sort_values('snps')
fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches(5,3)
sns.barplot(data=samples_with_snps, x='snps', y='sample', ax=ax, color='#307B91')
ax.set_xlabel('Number of SNPs in sample', fontsize=14)
ax.set_ylabel('Count of samples', fontsize=14)
plt.savefig('Z:/volume1/noam/covid_data/snps_per_sample.png' , bbox_inches='tight', dpi=800)



## snps over dates
metadata = pd.read_excel('X:/volume2/noam/covid/mutation_type/all_sequening_run_metadata.xlsx')
metadata['sample'] = metadata['sample'].astype(str)
snps_per_date = pd.merge(snps_per_sample, metadata[['date', 'sample']], on='sample', how='left').sort_values('date')
snps_per_date['date'] = pd.to_datetime(snps_per_date.date)
snps_per_date.to_csv('Z:/volume1/noam/covid_data/snps_per_sample_date.csv', index=False)

fig, ax = plt.subplots(nrows=1, ncols=1)
ax.scatter(x=snps_per_date.date, y=snps_per_date.snps)
ax.set_xlim('2020-02-28', '2020-04-25')
fig.set_size_inches(10,3)

counts = snps_per_date.groupby(['date']).snps.count().reset_index().rename(columns={'snps':'count'})
averages = snps_per_date.groupby(['date']).snps.mean().reset_index().rename(columns={'snps':'average_snps'})
a = pd.merge(counts, averages, on='date')
fig, ax = plt.subplots(nrows=1, ncols=1)
sns.scatterplot(data=a, x='date', y='average_snps', size='count', ax=ax)
ax.set_xlim('2020-02-28', '2020-04-25')
ax.legend(bbox_to_anchor=(1.1, 1.05))
fig.set_size_inches(10,3)

counts = snps_per_date.groupby(['date', 'snps']).count().reset_index().rename(columns={'sample':'count'})
fig, ax = plt.subplots(nrows=1, ncols=1)
sns.scatterplot(data=counts, x='date', y='snps', size='count', ax=ax)
ax.set_xlim('2020-02-28', '2020-04-25')
ax.legend(bbox_to_anchor=(1.1, 1.05))
fig.set_size_inches(10,3)




#### count ns in each fasta seq
metadata = pd.read_excel('Z:/volume1/noam/covid_data/all sequening run metadata.xlsx')
metadata['sample'] = metadata['sample'].astype(str)

with open('Z:/volume1/noam/covid_data/all_israel_concensuses.fasta') as f:
    f = f.read()
n_count = []
for seq in f.split('>'):
    if seq != '':
        n_count.append([seq.split('\n')[0], seq.split('\n')[1].count('N')])

n_count = pd.DataFrame(n_count, columns=['sample', 'n_count'])
n_count['precent_of_bases_sequenced'] = ((29892 - n_count.n_count) / 29892)
n_count['sample'] = n_count['sample'].str.split('/').str[1]
n_count = n_count.drop(columns=['n_count'])
pd.merge(metadata, n_count, on='sample', how='outer').to_excel('Z:/volume1/noam/covid_data/all sequening run metadata.percent_mapped.xlsx', index=False)



##### syn, non-synonymous
df = pd.read_csv('/Volumes/STERNADILABTEMP$//volume1/noam/covid_data/all_mutations_in_included_cons.csv')

mutation_type = pd.read_csv('/Volumes/STERNADILABHOME$/volume2/noam/covid/mutation_type/mutation_type.short.csv')
mutation_type = mutation_type.rename(columns={'Base':'base', 'Ref':'ref_base', 'Pos':'ref_position'})
df = pd.merge(df, mutation_type, on=['ref_position', 'ref_base', 'base'], how='left')

# non sysn
print(len(df[(df.base != '-') & (df.ref_base != '-') & ~(df.mutation_type.astype(str) == '[]') & (df.mutation_type.isna() == False)].groupby('ref_position').count()))
# synonymous
print(len(df[(df.base != '-') & (df.ref_base != '-') & (df.mutation_type.astype(str) == '[]')].groupby('ref_position').count()))
   # non oding
print(len(df[(df.base != '-') & (df.ref_base != '-') & (df.mutation_type.isna() == True)].groupby('ref_position').count()))

fig, ax = plt.subplots(nrows=1, ncols=1)
df[(df.base != '-') & (df.ref_base != '-') & ~(df.mutation_type.astype(str) == '[]') & (df.mutation_type.isna() == False)].groupby('ref_position').ref_base.count().reset_index().plot(x='ref_position', y='ref_base', kind='scatter', color='red', ax=ax, label='non-synonymous')
df[(df.base != '-') & (df.ref_base != '-') & (df.mutation_type.astype(str) == '[]')].groupby('ref_position').ref_base.count().reset_index().plot(x='ref_position', y='ref_base', kind='scatter', color='#22DA05', ax=ax, label='synonymous')
df[(df.base != '-') & (df.ref_base != '-') & (df.mutation_type.isna() == True)].groupby('ref_position').ref_base.count().reset_index().plot(x='ref_position', y='ref_base', kind='scatter', color='blue', ax=ax, label='non-coding')
fig.set_size_inches(20,4)
ax.set_xlabel('Position (nucleotide base)', fontsize=20)
ax.set_ylabel('Count of samples with mutation',fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.set_xlim(0,30000)
ax.set_yscale('log')
ax.legend(bbox_to_anchor=(0.8, -0.25), fontsize=20, ncol=3)

colors={'orf1ab polyprotein':'blue',
 'surface glycoprotein':'red',
 'ORF3a protein':'purple',
 'envelope protein':'yellow',
 'membrane glycoprotein':'orange',
 'ORF6 protein':'green',
 'ORF7a protein':'blue',
 'ORF8 protein':'purple',
 'nucleocapsid phosphoprotein':'yellow',
 'ORF10 protein':'orange'}
for i in mutation_type[['start', 'end', 'protein']].drop_duplicates().itertuples():
    ax.axvspan(i[1], i[2], facecolor=colors[i[3]], alpha=0.2)
    if i[3] != 'orf1ab polyprotein':
        if i[3] == 'surface glycoprotein':
            ax.text(((((i[2] - i[1])/2)+i[1])/29892), 1.06, i[3] + '\n(spike)', transform=ax.transAxes, fontsize=14, rotation=90, horizontalalignment='center')
        else:
            ax.text(((((i[2] - i[1])/2)+i[1])/29892), 1.06, i[3], transform=ax.transAxes, fontsize=14, rotation=90, horizontalalignment='center')
    else:
        ax.text(10000/29892, 1.06, i[3], transform=ax.transAxes, fontsize=14, rotation=90, horizontalalignment='center')
plt.ylim(bottom=0.8)
#plt.savefig('Z:/volume1/noam/covid_data/cumulative_syn_nonsyn.png' , bbox_inches='tight', dpi=1600)
plt.savefig('/Volumes/STERNADILABTEMP$/volume1/noam/covid_data/cumulative_syn_nonsyn.pdf' , bbox_inches='tight', dpi=1600)



## df top mutations chart
df = pd.read_csv('Z:/volume1/noam/covid_data/all_mutations_in_included_cons.csv')
df['full_mutation'] = df.ref_base + df.ref_position.astype(int).astype(str) + df.base
count_mutations = df.groupby(['full_mutation', 'ref_base', 'ref_position', 'base']).sample.count().sort_values().reset_index().rename(columns={'sample':'count_of_samples'})
mutation_type = pd.read_csv('X:/volume2/noam/covid/mutation_type/mutation_type.short.csv')
mutation_type = mutation_type.rename(columns={'Base':'base', 'Ref':'ref_base', 'Pos':'ref_position'})
count_mutations = pd.merge(count_mutations, mutation_type, on=['ref_position', 'ref_base', 'base'], how='left')
count_mutations.to_csv('Z:/volume1/noam/covid_data/mutations_counted_details.csv', index=False)

# how many unique per type
len(count_mutations[(count_mutations.base != '-') & (count_mutations.ref_base != '-') & ~(count_mutations.mutation_type.astype(str) == '[]') & (count_mutations.mutation_type.isna() == False)])
len(count_mutations[(count_mutations.base != '-') & (count_mutations.ref_base != '-') & (count_mutations.mutation_type.astype(str) == '[]')])
len(count_mutations[(count_mutations.base != '-') & (count_mutations.ref_base != '-') & (count_mutations.mutation_type.isna() == True)])


#### heatmap
df = pd.read_csv('Z:/volume1/noam/covid_data/all_mutations_in_included_cons.csv')
df['temp'] = 1

#count_mutations = df.groupby('ref_position').sample.count().sort_values().reset_index()

#to_pivot = df[df.ref_position.isin(count_mutations[count_mutations['sample'] > 10].ref_position.tolist())].pivot_table(values='temp', index=['ref_position'], columns='sample')
to_pivot = df.pivot_table(values='temp', index=['ref_position'], columns='sample')

to_pivot = to_pivot.fillna(0)
clustergrid = sns.clustermap(to_pivot, mask=True)
to_pivot = to_pivot[[to_pivot.columns[s] for s in clustergrid.dendrogram_col.reordered_ind]]
to_pivot['mutation_order'] = pd.Categorical(to_pivot.index, [to_pivot.index[s] for s in clustergrid.dendrogram_row.reordered_ind])
to_pivot.sort_values('mutation_order').to_excel('Z:/volume1/noam/covid_data/pivot_mutations_cons.xlsx')



############

df = pd.read_csv('/Volumes/STERNADILABHOME$/volume2/noam/covid/database_research/msa_0520/msa_0520.compare.remove_edges.csv')

# check if deletions are associated with mutation
chi2_data = []
#samples_with_del = df[(df.base == '-')]['sample'].drop_duplicates().tolist()
for snp in tqdm(df[(df.base != 'N')][['position']].drop_duplicates().position.tolist()[:13]):
    temp_matrix = pd.DataFrame(np.zeros(shape=(2,2)), columns=['snp',0], index=['del',0])
    temp_matrix.at['del','snp'] = len(df[(df['sample'].isin(samples_with_del)) & (df.position == snp) & (df.base != 'N')][['sample']].drop_duplicates())
    temp_matrix.at['del',0] = len(df[(df['sample'].isin(samples_with_del))][['sample']].drop_duplicates()) - len(df[(df['sample'].isin(samples_with_del)) & (df.position == snp) & (df.base != 'N')][['sample']].drop_duplicates())
    temp_matrix.at[0,'snp'] = len(df[~(df['sample'].isin(samples_with_del)) & (df.position == snp) & (df.base != 'N')][['sample']].drop_duplicates())
    temp_matrix.at[0,0] = len(df[~(df['sample'].isin(samples_with_del))][['sample']].drop_duplicates()) - len(df[~(df['sample'].isin(samples_with_del)) & (df.position == snp) & (df.base != 'N')][['sample']].drop_duplicates())
    chi2, pvalue, dof, expected = scipy.stats.chi2_contingency(temp_matrix)
    chi2_data.append(('del',snp,pvalue,chi2))
chi2_data = pd.DataFrame(chi2_data, columns=['pos1', 'pos2', 'pvalue', 'chi2'])