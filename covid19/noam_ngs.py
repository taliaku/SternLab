
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append('/sternadi/home/volume1/shared/SternLab')
sys.path.append('X:/volume2/noam/SternLab')
from blast_utilities import blast_to_df
from freqs_utilities import estimate_insertion_freq

def unite_all_freq_files(freqs_dir, out_path=None):
    """
    unites all frequency file into one file. add the following:
    Replica, Degree, Sample, Time
    :param freqs_dir: a directory with frequency files in a format of pTIME-DEGREEREPLICA.freqs
    :param out: output directory to save results
    :return: a merged data frame
    """
    freqs_df = []
    # read all freq files and add the sample name to the data frame
    files = [os.path.join(freqs_dir, f) for f in os.listdir(freqs_dir) if not f.startswith('all')]
    for f in files:
        print(f)
        #curr_df = pd.read_csv(f, sep='\t')
        curr_df = pd.read_csv(f)
        sample = os.path.basename(f).split('.')[0].split('_')[0]
        #sample = os.path.basename(f)
        curr_df['Sample'] = sample
        freqs_df.append(curr_df)
    df = pd.concat(freqs_df)

    if out_path != None:
        df.to_csv(out_path, index=False)
    return df



####### TMN samples
df = pd.read_csv('Z:/volume1/noam/covid_data/TMNcorona_20200410/TMNcorona_20200410_old/pipelines/all_freqs.csv')

for sample in df.Sample.drop_duplicates().tolist():
    fig, ax = plt.subplots(nrows=1, ncols=1)
    df[(df.Sample == sample) & (df.Base == df.Ref) & (df.Ref != '-')].plot(x='Pos', y='Read_count', ax=ax, kind='scatter')
    ax.set_yscale('log')
    ax.set_title(sample)
    fig

primer_positions = [[192,213],
                    [7629,7653],
                    [7290,7312],
                    [14776,14794],
                    [14718,14743],
                    [22419,22439],
                    [22182,22203],
                    [29683,29706],
                    ]
primer_positions = sum([list(range(pair[0], pair[1])) for pair in primer_positions], [])
unamplified_positions = list(range(192)) + list(range(29706,30000))

# make concensuses
df[(df.Sample == 'SH15') & (df.Ref != df.Base) & (df.Ref != '-') & (df.Freq > 0.3)].sort_values('Pos')
df[(df.Sample == 'SH15') & (df.Ref != df.Base) & (df.Ref != '-') & (df.Freq > 0.3) & (df.Pos.isin(primer_positions + unamplified_positions))].sort_values('Pos')
# we have some positions in the 5' that should not be mapped, we can drop them with using
# a cutoff of 50 read counts for cconcesus without loosing other SNPs
# python SternLab/NGS_analysis/make_reference_from_consensus.py -f covid/MN908947.fasta -p /sternadi/nobackup/volume1/noam/covid_data/TMNcorona_20200410_old/pipelines/SH15/SH15.freqs -o /sternadi/nobackup/volume1/noam/covid_data/TMNcorona_20200410_concensus/SH15.fasta -c 50

df[(df.Sample == 'SH16') & (df.Ref != df.Base) & (df.Ref != '-') & (df.Freq > 0.3)].sort_values('Pos')
df[(df.Sample == 'SH16') & (df.Ref != df.Base) & (df.Ref != '-') & (df.Freq > 0.3) & (df.Pos.isin(primer_positions + unamplified_positions))].sort_values('Pos')
# cutoff of 50 as well works here



# look at blast around SH15 deletion
mutations = pd.read_csv('Z:/volume1/noam/covid_data/TMNcorona_20200410_concensus/SH15_blasts/mutations.csv')
blasts = pd.read_csv('Z:/volume1/noam/covid_data/TMNcorona_20200410_concensus/SH15_blasts/blasts.csv')

a = blasts[blasts.read.isin(mutations[(mutations.position == 28254) & (mutations.base == '-')].read.tolist())]


b = blasts.groupby('read').start_ref.count().reset_index()






######## technion1
# before cleaning
df = pd.read_csv('Z:/volume1/noam/covid_data/coronaTech1_20200415/python_pipeline_x1_c0/freqs/all.fixed.csv')

# after cleaning
df = pd.read_csv('Z:/volume1/noam/covid_data/coronaTech1_20200415/python_pipeline_x1_c0_after_ptrimmer/freqs/all.freqs.csv')


##########################################

fig, ax = plt.subplots(nrows=6, ncols=6)
ax = ax.flatten()
i = 0
for sample in df.Sample.drop_duplicates().tolist():
    df[(df.Sample == sample) & (df.base == df.ref_base) & (df.ref_base != '-')].plot(x='ref_position', y='coverage', ax=ax[i], kind='scatter')
    ax[i].set_yscale('log')
    ax[i].set_title(sample.split('_')[0])
    i += 1
fig
fig.subplots_adjust(hspace=0.2)
fig.set_size_inches(20,20)


df[(df.ref_position.isin(range(31,29866))) &(df.base != df.ref_base) & (df.ref_base != '-') & (df['rank'] == 0) & (df.coverage > 10)].groupby('ref_position').base.count().reset_index().to_clipboard()


fig, ax = plt.subplots(nrows=6, ncols=6)
ax = ax.flatten()
i = 0
for sample in df.Sample.drop_duplicates().tolist():
    if sample not in ['2047927_S16_merge', '990333263_S4_merge']:
        df[(df.Sample == sample) & (df.base != df.ref_base) & (df.ref_base != '-') & (df.frequency > 0.01)].plot(x='ref_position', y='frequency', ax=ax[i], kind='scatter')
    ax[i].set_title(sample.split('_')[0])
    ax[i].set_yscale('log')
    ax[i].set_ylim(10**-3,1)
    i += 1
fig
fig.subplots_adjust(hspace=0.2)
fig.set_size_inches(20,20)


##### find problematic amplicons

fig, ax = plt.subplots(nrows=2, ncols=2)
ax = ax.flatten()
i = 0
for sample in df.Sample.drop_duplicates().tolist():
    if i < 4:
        df[(df.ref_position.isin(range(31,29866))) & (df.Sample == sample) & (df.base == df.ref_base) & (df.ref_base != '-')].plot(x='ref_position', y='coverage', ax=ax[i], kind='scatter')
        ax[i].set_yscale('log')
        ax[i].set_title(sample.split('_')[0])
    i += 1
fig
fig.subplots_adjust(hspace=0.2)
fig.set_size_inches(20,20)


###### coverage difference for qpcr genes
# N, E, RDR, S
genes = {'E':(26245,26472),
         'S':(21563,25384),
         'N':(28274,29533),
         'RDR':(15102,15283)}
genes_df = []
for g in genes:
    for i in range(genes[g][0], genes[g][1] + 1):
        genes_df.append([g, float(i)])
genes_df = pd.DataFrame(genes_df, columns=['gene', 'ref_position'])     

genes_df = pd.merge(df, genes_df, on='ref_position')

sns.boxplot(x='Sample', hue='gene', y='coverage',data=genes_df[genes_df.Sample.isin(['990430265_S1_merge', '1639953_S11_merge'])])


####### cts and coverag statistics

df = pd.read_csv('Z:/volume1/noam/covid_data/coronaTech1_20200415/python_pipeline_x1_c0_after_ptrimmer/freqs/all.freqs.csv')

coverage_over_10 = df[(df.ref_position.isin(range(31,29866))) & (df.base == df.ref_base) & (df.ref_base != '-') & (df.coverage >= 10)].groupby('Sample').base.count().reset_index()
average_coverage = df[(df.ref_position.isin(range(31,29866))) & (df.base == df.ref_base) & (df.ref_base != '-') & (df.coverage >= 10)].groupby('Sample').coverage.mean().reset_index()
summary = pd.merge(coverage_over_10, average_coverage, on='Sample').rename(columns={'base':'bases_with_coverage_over_10', 'coverage':'average_coverage'})
summary['Sample'] = summary.Sample.str.split('_').str[0]
summary['bases_with_coverage_over_10_percent'] = 100 * (summary.bases_with_coverage_over_10 / (29866 - 31))

SH15 = pd.read_csv('Z:/volume1/noam/covid_data/TMNcorona_20200410/TMNcorona_20200410_python_pipeline_c0/SH15/SH15_S14_merge.freqs.csv')
SH15['Sample'] = '13075832'
SH16 = pd.read_csv('Z:/volume1/noam/covid_data/TMNcorona_20200410/TMNcorona_20200410_python_pipeline_c0/SH16/SH16_S15_merge.freqs.csv')
SH16['Sample'] = '13075879'
df = pd.concat([SH15, SH16])
coverage_over_10 = df[(df.ref_position.isin(range(192,29706))) & (df.base == df.ref_base) & (df.ref_base != '-') & (df.coverage >= 10)].groupby('Sample').base.count().reset_index()
average_coverage = df[(df.ref_position.isin(range(192,29706))) & (df.base == df.ref_base) & (df.ref_base != '-') & (df.coverage >= 10)].groupby('Sample').coverage.mean().reset_index()
summary_tmn = pd.merge(coverage_over_10, average_coverage, on='Sample').rename(columns={'base':'bases_with_coverage_over_10', 'coverage':'average_coverage'})
summary_tmn['bases_with_coverage_over_10_percent'] = 100 * (summary.bases_with_coverage_over_10 / (29706 - 192))

summary = pd.concat([summary_tmn, summary])
cts = pd.read_excel('X:/volume2/noam/covid/samples_and_cts.xlsx')
cts['Sample'] = cts['Sample'].astype(str)
summary = pd.merge(summary, cts, on='Sample')

#df[df.Sample.str.split('_').str[0].isin(['2047927','990333263'])]




######## cluster
df = pd.read_csv('Z:/volume1/noam/covid_data/coronaTech1_20200415/python_pipeline_x1_c0_after_ptrimmer/freqs/all.freqs.csv')
df = df[(df.ref_position.isin(range(31,29866)))]
df = df[~df['Sample'].isin(['2047927_S16_merge', '990333263_S4_merge'])]

SH15 = pd.read_csv('Z:/volume1/noam/covid_data/TMNcorona_20200410/TMNcorona_20200410_python_pipeline_c0/SH15/SH15_S14_merge.freqs.csv')
SH15['Sample'] = '13075832'
SH16 = pd.read_csv('Z:/volume1/noam/covid_data/TMNcorona_20200410/TMNcorona_20200410_python_pipeline_c0/SH16/SH16_S15_merge.freqs.csv')
SH16['Sample'] = '13075879'
df = pd.concat([df, SH15, SH16])

df['Sample'] = df['Sample'].str.split('_').str[0]
df['full_mutation'] = df.ref_base + df.ref_position.astype(int).astype(str) + df.base
#mutations_to_keep = df[(df.ref_base != df.base) & (df.ref_base != '-') & (df['rank'] == 0) & (df.coverage >= 5)][['full_mutation']].drop_duplicates()
mutations_to_keep = df[(df.ref_base != df.base) & (df.ref_base != '-') & (df.frequency >= 0.2) & (df.coverage >= 5)][['full_mutation']].drop_duplicates()
mutations_to_keep = pd.merge(mutations_to_keep, df, on=['full_mutation'], how='left')
to_pivot = mutations_to_keep.pivot_table(values='frequency', index=['full_mutation'], columns='Sample')
to_pivot = to_pivot.dropna()

clustergrid = sns.clustermap(to_pivot, mask=True)


to_pivot = to_pivot[[to_pivot.columns[s] for s in clustergrid.dendrogram_col.reordered_ind]]
to_pivot['mutation_order'] = pd.Categorical(to_pivot.index, [to_pivot.index[s] for s in clustergrid.dendrogram_row.reordered_ind])
to_pivot.sort_values('mutation_order').to_excel('Z:/volume1/noam/covid_data//all_mutations_over0.2_pivot.xlsx')





##### snps in mutations
israel_seqs = pd.read_csv('X:/volume2/noam/covid/technion1/technion1.fasta.mutations_list.csv')
gisaid_seqs = pd.read_csv('X:/volume2/noam/covid/gisaid_cov2020_sequences.fasta.mutations_list.csv')

# filter out pangolin sequences and very far sequences, filters about 1 percent
gisaid_seq_mutation_count = gisaid_seqs[~(gisaid_seqs.file.str.contains('pangolin')) & (gisaid_seqs.base != 'N')].groupby('file').position.count().reset_index()
gisaid_seqs = gisaid_seqs[gisaid_seqs.file.isin(gisaid_seq_mutation_count[gisaid_seq_mutation_count.position <= 30].file.tolist())]


gisaid_seqs[gisaid_seqs.base != 'N'][['position', 'ref', 'base']].drop_duplicates().groupby(['ref', 'base']).count().to_clipboard()
israel_seqs[israel_seqs.base != 'N'][['position', 'ref', 'base']].drop_duplicates().groupby(['ref', 'base']).count().to_clipboard()











######### technion2
# after cleaning
df = pd.read_csv('Z:/volume1/noam/covid_data/coronaTech2_20200427/python_pipeline_x1_c0_after_ptrimmer/freqs/all.freqs.csv')
df['Sample'] = df.Sample.astype(str)

coverage_over_10 = df[(df.ref_position.isin(range(31,29866))) & (df.base == df.ref_base) & (df.ref_base != '-') & (df.coverage >= 5)].groupby('Sample').base.count().reset_index()
average_coverage = df[(df.ref_position.isin(range(31,29866))) & (df.base == df.ref_base) & (df.ref_base != '-') & (df.coverage >= 5)].groupby('Sample').coverage.mean().reset_index()
summary = pd.merge(coverage_over_10, average_coverage, on='Sample').rename(columns={'base':'bases_with_coverage_over_10', 'coverage':'average_coverage'})
summary['Sample'] = summary.Sample.str.split('_').str[0]
summary['bases_with_coverage_over_10_percent'] = 100 * (summary.bases_with_coverage_over_10 / (29866 - 31))

#cts = pd.read_excel('X:/volume2/noam/covid/samples_and_cts.xlsx')
#cts['Sample'] = cts['Sample'].astype(str)
#summary = pd.merge(summary, cts, on='Sample')

#df[df.Sample.str.split('_').str[0].isin(['2047927','990333263'])]



### look at deletions

df[(df.ref_base != df.base) & (df.base == '-') & (df['rank'] == 0) & (df.ref_position.isin(range(31,29866))) & (df.coverage >= 5)].to_excel('Z:/volume1/noam/covid_data/investigate_deleltions.xlsx', index=False)
df = estimate_insertion_freq(df, ['Sample',])


def estimate_insertion_freq(df, extra_columns=[]):
    '''
    This function gets a freqs file(s) dataframe, calculates the frequency of insertions by using the read count of the
    previous base and returns a dataframe including this.
    :param df: a dataframe of freqs file(s).
    :param extra_columns: if df contains more than the basic freqs columns, for example a column of Sample_id etc., 
    provide a list of the extra columns to be included.
    :return: df with extra columns describing insertion frequency.
    '''
    read_counts = df[(df.ref_base != '-')][ extra_columns + ['ref_position', 'coverage']].drop_duplicates()
    read_counts.rename(columns={'coverage':'estimated_read_count', 'ref_position':'rounded_pos'}, inplace=True)
    insertions = df[(df.ref_base == '-')]
    not_insertions = df[(df.ref_base != '-')]
    insertions['rounded_pos'] = insertions.ref_position.astype(int).astype(float)
    insertions = pd.merge(insertions, read_counts, how='left', on= extra_columns + ['rounded_pos'])
    insertions['estimated_freq'] = insertions.frequency * insertions.coverage / insertions.estimated_read_count
    df = pd.concat([insertions, not_insertions])
    return df.sort_values(extra_columns + ['ref_position'])



df[(df.Sample == '2089839') & (df.ref_position.isin(range(3878, 3905))) & (df['rank'] == 0)]
df[(df.Sample == '2089852') & (df.ref_position.isin(range(3878, 3905))) & (df['rank'] == 0)]


df[(df.Sample == '13077726') & (df.ref_position.isin(range(27385, 27405))) & (df['rank'] == 0)]



######## big freqs csv 
df = pd.read_csv('Z:/volume1/noam/covid_data/coronaTech1_20200415/python_pipeline_x1_c0_after_ptrimmer/freqs/all.freqs.csv')

df2 = pd.read_csv('Z:/volume1/noam/covid_data/coronaTech2_20200427/python_pipeline_x1_c0_after_ptrimmer/freqs/all.freqs.csv')

SH15 = pd.read_csv('Z:/volume1/noam/covid_data/TMNcorona_20200410/TMNcorona_20200410_python_pipeline_c0/SH15/SH15_S14_merge.freqs.csv')
SH15['Sample'] = '13075832'
SH16 = pd.read_csv('Z:/volume1/noam/covid_data/TMNcorona_20200410/TMNcorona_20200410_python_pipeline_c0/SH16/SH16_S15_merge.freqs.csv')
SH16['Sample'] = '13075879'
df = pd.concat([df, SH15, SH16, df2])
df['Sample'] = df['Sample'].astype(str)
df['Sample'] = df['Sample'].str.split('_').str[0]
df.to_csv('Z:/volume1/noam/covid_data/all_freqs.upto_tech2.csv', index=False)

##### read all freqs
df =  pd.read_csv('Z:/volume1/noam/covid_data/all_freqs.upto_tech2.csv')

# cluster

df = df[(df.ref_position.isin(range(31,29866)))]

# all but none clustered
df['full_mutation'] = df.ref_base + df.ref_position.astype(int).astype(str) + df.base
#mutations_to_keep = df[(df.ref_base != df.base) & (df.ref_base != '-') & (df['rank'] == 0) & (df.coverage >= 5)][['full_mutation']].drop_duplicates()
mutations_to_keep = df[(df.ref_base != df.base) & (df.ref_base != '-') & (df.frequency >= 0.2) & (df.coverage >= 5)][['full_mutation']].drop_duplicates()
mutations_to_keep = pd.merge(mutations_to_keep, df, on=['full_mutation'], how='left')
to_pivot = mutations_to_keep.pivot_table(values='frequency', index=['full_mutation'], columns='Sample')
to_pivot.to_excel('Z:/volume1/noam/covid_data//all_mutations_over0.2_pivot.including_tech2.xlsx')

# drop 3 samples and mutations that do not appear in all samples
df = df[~df['Sample'].isin(['2047927', '990333263', '13077494'])]
mutations_to_keep = df[(df.ref_base != df.base) & (df.ref_base != '-') & (df.frequency >= 0.2) & (df.coverage >= 5)][['full_mutation']].drop_duplicates()
mutations_to_keep = pd.merge(mutations_to_keep, df, on=['full_mutation'], how='left')
to_pivot = mutations_to_keep.pivot_table(values='frequency', index=['full_mutation'], columns='Sample')
to_pivot = to_pivot.dropna()
clustergrid = sns.clustermap(to_pivot, mask=True)
to_pivot = to_pivot[[to_pivot.columns[s] for s in clustergrid.dendrogram_col.reordered_ind]]
to_pivot['mutation_order'] = pd.Categorical(to_pivot.index, [to_pivot.index[s] for s in clustergrid.dendrogram_row.reordered_ind])
to_pivot.sort_values('mutation_order').to_excel('Z:/volume1/noam/covid_data//all_mutations_over0.2_pivot.including_tech2.clustered.xlsx')





### saved all mapped more than once
def blast_mapped_too_much(in_blast_df, out_blast_df):
    a = pd.read_csv(in_blast_df)
    reads_to_keep = a.groupby('read').btop.count().reset_index()
    a[a.read.isin(reads_to_keep[reads_to_keep['btop'] > 2].read.tolist())].to_csv(out_blast_df, index=False)
    return 



a['seq'] = a.sort_values('start_read').groupby('read').cumcount()
b = pd.merge(a[a.seq == 0], a[a.seq == 1], on='read')
b = b.sort_values(['start_read_x', 'end_read_x', 'start_read_y', 'end_read_y'])
fig, ax = plt.subplots(nrows=1, ncols=1)
sns.swarmplot(x='start_read_x', y='read', data=b, color='red', ax=ax)
sns.swarmplot(x='end_read_x', y='read', data=b, color='red', ax=ax)
sns.swarmplot(x='start_read_y', y='read', data=b, color='blue', ax=ax)
sns.swarmplot(x='end_read_y', y='read', data=b, color='blue', ax=ax)
fig.set_size_inches(12, 48)
