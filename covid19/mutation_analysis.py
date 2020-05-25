
import pandas as pd
import json
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import math


#### organize mutations types - syn and non-syn


def flatten_json(nested_json):
    """
        Flatten json object with nested keys into a single level.
        Args:
            nested_json: A nested json object.
        Returns:
            The flattened json object if successful, None otherwise.
    """
    out = {}
    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + '_')
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i) + '_')
                i += 1
        else:
            out[name[:-1]] = x
    flatten(nested_json)
    return out


with open('/Volumes/STERNADILABTEMP$/volume1/covid/ncov_concensus_ISR/auspice/ncov.json', 'r') as o:
    data = json.load(o)
    
data = flatten_json(data)


with open('/Volumes/STERNADILABHOME$/volume2/noam/covid/mutation_type/tmp_json.txt', 'w') as outfile:
    json.dump(data, outfile)
    
    

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
df = pd.read_csv('Z:/volume1/noam/covid_data/all_mutations_in_included_cons.csv')

mutation_type = pd.read_csv('X:/volume2/noam/covid/mutation_type/mutation_type.short.csv')
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
plt.savefig('Z:/volume1/noam/covid_data/cumulative_syn_nonsyn.png' , bbox_inches='tight', dpi=1600)


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




####### fasta differences
def compare_fastas_to_ref_unaligned(fastas, ref_fasta, output_excel):
    with open(fastas) as f:
        f = f.read()
    f = f.split('>')
    #f = [i for i in f[1:] if i.split('/')[1] in ['2089839', '2089852', '13077726', '2086033', ]]
    f = {i.split('\n')[0]:''.join(i.split('\n')[1:]) for i in f}

    with open(ref_fasta) as g:
        g = g.read()
    g = ['reference', ''.join(g.replace('>', '').split('\n')[1:])]
    f[g[0]] = g[1]

    diffs = []
    for i in range(len(f['reference'])):
        for sample in f:
            if sample != 'reference':
                try:
                    if f[sample][i+1] != f['reference'][i]:
                        diffs.append((i+1, sample, f['reference'][i], f[sample][i+1]))
                except:
                    pass
    return pd.DataFrame(diffs, columns=['position', 'sample', 'ref_base', 'base']).to_excel(output_excel, index=False)

compare_fastas_to_ref_unaligned('/Volumes/STERNADILABTEMP$/volume1/noam/covid_data/all_israel_concensuses.fasta', '/Volumes/STERNADILABHOME$/volume2/noam/covid/MN908947.fasta', '/Volumes/STERNADILABTEMP$/volume1/noam/covid_data/all_israel_concensuses.all_diffs.xlsx')

df = pd.read_excel('/Volumes/STERNADILABTEMP$/volume1/noam/covid_data/all_israel_concensuses.all_diffs.xlsx')
all = df.groupby(['position', 'ref_base', 'base']).sample.count().sort_values()
clade = df[df['sample'].str.split('/').str[1].isin(['2089839', '2089852', '13077726', '2086033', ])].groupby(['position', 'ref_base', 'base']).sample.count().sort_values()
counted = pd.merge(all, clade, on=['position', 'ref_base', 'base'], how='outer', suffixes=['_all', '_clade']).fillna(0)
counted['all_freq'] = counted.sample_all / 213
counted['clade_freq'] = counted.sample_clade / 4





######## find deletions in all sequences

with open('/Volumes/STERNADILABHOME$/volume3/COVID19/data/gisaid_hcov-19_2020_05_05_11.fasta') as f:
    f = f.read()

f = f.split('>')
f = {i.split('\n')[0]:''.join(i.split('\n')[1:]) for i in f}

f_with_20755 = {}
for i in f:
    if 'AAAATTATGGTGATCGTGCAACATTA' in f[i]:
        f_with_20755[i] = f[i]

with open('/Volumes/STERNADILABHOME$/volume2/noam/covid/deletions_outside_israel/fastas_with_20755.fasta', 'w') as w:
    w.write(''.join(['>' + i + '\n' + f_with_20755[i] + '\n' for i in f_with_20755]))


###############################

def compare_fastas_to_ref(fastas, ref_seq_name, output_excel):
    # fasta needs to be aligned
    with open(fastas) as f:
        f = f.read()
    f = f.split('>')
    f = {i.split('\n')[0]:''.join(i.split('\n')[1:]) for i in f if i != ''}
    diffs = []
    ref_pos = 0
    for i in tqdm(range(len(f[ref_seq_name]))):
        if f[ref_seq_name][i] != '-':
            ref_pos = math.floor(ref_pos) + 1
            for sample in f:
                if f[sample][i] != f[ref_seq_name][i]:
                    diffs.append((ref_pos, sample, f[ref_seq_name][i], f[sample][i]))
        else:
            ref_pos += 0.001
            for sample in f:
                if f[sample][i] != f[ref_seq_name][i]:
                    diffs.append((ref_pos, sample, f[ref_seq_name][i], f[sample][i]))
    return pd.DataFrame(diffs, columns=['position', 'sample', 'ref_base', 'base']).to_csv(output_excel, index=False)

compare_fastas_to_ref('/Volumes/STERNADILABHOME$/volume2/noam/covid/database_research/msa_0520/msa_0520.fasta', 'hCoV-19/Wuhan-Hu-1/2019|EPI_ISL_402125|2019-12-31|Asia', '/Volumes/STERNADILABHOME$/volume2/noam/covid/database_research/msa_0520/msa_0520.compare.csv')