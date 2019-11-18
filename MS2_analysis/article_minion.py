# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 09:19:39 2019

@author: Noam
"""
import sys
sys.path.append('/sternadi/home/volume2/noam/SternLab')
sys.path.append(r'X:\volume2\noam\Sternlab')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
from freqs_utilities import compare_positions_between_freqs
from mpl_toolkits.mplot3d import Axes3D



COLORS = {'T1764.0-':'#F50202', 'A1664.0G':'#F49D09', 'A535.0G':'#5EC0D2', 'T1440.0C':'#F1F87A', 'T1440.0G':'#C4A0E4', 'A1443.0G':'#FE1D1D', 'A1611.0G':'#327CFD', 'C1724.0T':'#8FD95A', 'A1744.0G':'#FBB3DD', 'G1560.0A':'#A7F0A2', 'G1906.0A':'#A3A3A3', 'C3358.0T':'#26451C', 'G3114.0A':'#B37A42', 'A1770.0G':'#7603BA', 'G2310.0A':'#033E86', 'A2626.0G':'#8FD95A', 'C3299.0T':'#211785', 'C1718.0T':'#DFC236', 'T862.0C':'#880E05', 'A2790.0T':'#DF36C6', 'G1736.0A':'#CFFD2F', 'C1549.0T':'#2CA403', 'G531.0A':'#972FFE', 'C1050.0T':'#13B908'}
#n = [17, 18, 19, 20, 21, 22, 3462, 3463, 3524, 3542, 3543, 3544, 3545, 3546, 3547, 3548, 3564, 3566]
n = list(range(30)) + list(range(3540, 3570))


#################### perrcent of variants per read graph ###########
os.environ['QT_QPA_PLATFORM']='offscreen'
#from blast_utilities import blast_to_mutations_list, blast_to_df
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def percent_variants_per_read_plot(input_blast_df, input_mutations_df, output_path, title):
    blast_df = pd.read_csv(input_blast_df)
    mutations_df = pd.read_csv(input_mutations_df)
    # only keep reads that were mapped only once
    blast_df['read_count'] = blast_df.groupby('read')['start_ref'].transform('count')
    blast_df = blast_df[(blast_df.read_count == 1)]
    blast_df['ref_length'] = blast_df.end_ref - blast_df.start_ref
    mutations_df = mutations_df.groupby('read').position.count().reset_index().rename(columns={'position':'variant_count'})
    df = pd.merge(blast_df, mutations_df, on='read', how='left')

    df['variant_percent'] = (df.variant_count / df.ref_length) * 100
    plot = df.variant_percent.hist(grid=False, bins=20)
    plot.set_title(title, fontsize=16, x=0.5, y=1.05)
    plot.set_xlabel('Percent of variants on the same read', fontsize=14)
    plot.set_ylabel('Number of reads', fontsize=14)
    fig = plot.get_figure()
    fig.savefig(output_path, dpi=800, bbox_inches='tight')

#percent_variants_per_read_plot('X:/volume2/noam/minion/p15-37A/association_tests/blasts.csv', 'X:/volume2/noam/minion/p15-37A/association_tests/mutations.csv', 'X:/volume2/noam/minion/p15-37A/percent_variants_per_read.png', '37A')
#percent_variants_per_read_plot('X:/volume2/noam/minion/p15-37B/association_tests/blasts.csv', 'X:/volume2/noam/minion/p15-37B/association_tests/mutations.csv', 'X:/volume2/noam/minion/p15-37B/percent_variants_per_read.png', '37B')
percent_variants_per_read_plot('/sternadi/home/volume2/noam/minion/p15-37A/association_tests/blasts.csv', '/sternadi/home/volume2/noam/minion/p15-37A/association_tests/mutations.csv', '/sternadi/home/volume2/noam/minion/p15-37A/percent_variants_per_read.png', 'p15A')
percent_variants_per_read_plot('/sternadi/home/volume2/noam/minion/p15-37B/association_tests/blasts.csv', '/sternadi/home/volume2/noam/minion/p15-37B/association_tests/mutations.csv', '/sternadi/home/volume2/noam/minion/p15-37B/percent_variants_per_read.png', 'p15B')


########## coverage ############3


def coverage_plot(freqs_path, output_path, title='Read Count Per Base'):
    df = pd.read_csv(freqs_path, sep='\t')
    plot = df[(df.Base == df.Ref) & (df.Ref != '-')].plot(x='Pos', y='Read_count', legend=False, color='#079418')
    plot.set_title(title, fontsize=16, x=0.5, y=1.05)
    plot.set_xlabel('Position in the genome (basepairs)', fontsize=14)
    plot.set_ylabel('Number of reads', fontsize=14)
    fig = plot.get_figure()
    fig.savefig(output_path, dpi=800, bbox_inches='tight')
    
coverage_plot('X:/volume2/noam/minion/p15-37A/pipeline_q0_b60/fastq.freqs', 'X:/volume2/noam/minion/p15-37A/coverage_plot.png', title='p15-37A')
coverage_plot('X:/volume2/noam/minion/p15-37B/pipeline_q0_b60/fastq.freqs', 'X:/volume2/noam/minion/p15-37B/coverage_plot.png', title='p15-37B')
coverage_plot('X:/volume2/noam/minion/p1-37A/pipeline_q0_b60/FAK41568.freqs', 'X:/volume2/noam/minion/p1-37A/coverage_plot.png', title='p1-37A')


def coverage_plot2(freqs_dict, output_path, title='Read Count Per Base'):
    fig, ax = plt.subplots()
    for f in freqs_dict:
        df = pd.read_csv(freqs_dict[f], sep='\t')
        df = df[(df.Base == df.Ref) & (df.Ref != '-')]
        df.plot(x='Pos', y='Read_count', ax=ax, label=f)
    ax.set_title(title, fontsize=16, x=0.5, y=1.05)
    ax.set_xlabel('Position in the genome (basepairs)', fontsize=14)
    ax.set_ylabel('Number of reads', fontsize=14)
    ax.set_yscale('log')
    fig.savefig(output_path, dpi=800, bbox_inches='tight')

freqs_dict = {'p15A':'X:/volume2/noam/minion/p15-37A/pipeline_q0_b60/fastq.freqs','p15B':'X:/volume2/noam/minion/p15-37B/pipeline_q0_b60/fastq.freqs', 'control':'X:/volume2/noam/minion/p1-37A/pipeline_q0_b60/FAK41568.freqs'}
coverage_plot2(freqs_dict, 'X:/volume2/noam/minion/coverage_plot_log.png')
#######  error boxplot ######


def estimate_insertion_freq(df, extra_columns=[]):
    read_counts = df[(df.Ref != '-')][ extra_columns + ['Pos', 'Read_count']].drop_duplicates()
    read_counts.rename(columns={'Read_count':'estimated_read_count', 'Pos':'rounded_pos'}, inplace=True)
    insertions = df[(df.Ref == '-')]
    not_insertions = df[(df.Ref != '-')]
    insertions['rounded_pos'] = insertions.Pos.astype(int).astype(float)
    insertions = pd.merge(insertions, read_counts, how='left', on= extra_columns + ['rounded_pos'])
    insertions['estimated_freq'] = insertions.Freq * insertions.Read_count / insertions.estimated_read_count
    insertions['Freq'] = insertions.estimated_freq
    df = pd.concat([insertions, not_insertions])
    return df.sort_values(extra_columns + ['Pos'])

def prepare_freqs_for_error_analysis(df):
    df2 = estimate_insertion_freq(df)
    df2['rounded_pos_max_read_count'] = df2.groupby('rounded_pos').Read_count.transform('max')
    insertions_1 = df2[(df2.rounded_pos_max_read_count == df2.Read_count) & (df2.Pos.astype(str).str.endswith('.1'))]
    natural_pos = insertions_1.rounded_pos.drop_duplicates().tolist()
    fake_insertions = insertions_1[['Base', 'Ref',]].drop_duplicates()
    fake_positions = []
    for i in range(1,3570):
        if float(i) not in natural_pos:
            fake_positions.append(float(i))
    fake_positions = pd.DataFrame(fake_positions).rename(columns={0:'Pos'})
    fake_positions['Freq'] = 0.0
    fake_insertions['Freq'] = 0.0
    fake_insertions = pd.merge(fake_insertions, fake_positions, how='outer', on='Freq').sort_values('Pos')
    substitutions_and_deletions = df2[(df2.Base != df2.Ref) & (df2.Ref != '-')].sort_values('Ref')
    return pd.concat([substitutions_and_deletions, insertions_1, fake_insertions, ]).sort_values(['Pos'])

    

def boxplot_mutations(df, out_path, title):
    import matplotlib.pyplot as plt
    import seaborn as sns
    df = df[~(df.Pos.isin(n))]
    df['mutation'] = df.Ref + df.Base
    insertions = df[df.Ref == '-'].sort_values('Base')
    deletions = df[df.Base == '-'].sort_values('Base')
    others = df[(df.Ref != '-') & (df.Base != '-')].sort_values(['Ref', 'Base'])
    df = pd.concat([deletions, insertions, others])
    df = df[(df.Ref != df.Base)]
    sns.boxplot(data=df, x='mutation', y='Freq')
    plt.title(title, fontsize=16, x=0.5, y=1.05)
    plt.xlabel('Variant', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.ylim(0,0.8)
    plt.savefig(out_path, dpi=800)
    return


df = pd.read_csv('X:/volume2/noam/minion/p15-37B/pipeline_q0_b60/fastq.freqs', sep='\t')
df2 = prepare_freqs_for_error_analysis(df)
df2['mutation'] = df2.Ref + df2.Base
boxplot_mutations(df2, 'X:/volume2/noam/minion/p15-37B/variant_boxplot.png', 'MinION-p15B')

df = pd.read_csv('X:/volume2/noam/minion/p15-37A/pipeline_q0_b60/fastq.freqs', sep='\t')
df2 = prepare_freqs_for_error_analysis(df)
df2['mutation'] = df2.Ref + df2.Base
boxplot_mutations(df2, 'X:/volume2/noam/minion/p15-37A/variant_boxplot.png', 'MinION-p15A')

df = pd.read_csv('X:/volume2/noam/passages/freqs/p15-37B.freqs', sep='\t')
df2 = prepare_freqs_for_error_analysis(df)
df2['mutation'] = df2.Ref + df2.Base
boxplot_mutations(df2, 'X:/volume2/noam/minion/p15-37B/variant_boxplot_miseq.png', 'Illumina-p15B')

df = pd.read_csv('X:/volume2/noam/passages/freqs/p15-37A.freqs', sep='\t')
df2 = prepare_freqs_for_error_analysis(df)
df2['mutation'] = df2.Ref + df2.Base
boxplot_mutations(df2, 'X:/volume2/noam/minion/p15-37A/variant_boxplot_miseq.png', 'Illumina-p15B')

df = pd.read_csv('X:/volume2/noam/minion/p1-37A/pipeline_q0_b60/FAK41568.freqs', sep='\t')
df2 = prepare_freqs_for_error_analysis(df)
df2['mutation'] = df2.Ref + df2.Base
boxplot_mutations(df2, 'X:/volume2/noam/minion/p1-37A/variant_boxplot_minion.png', 'MinION-control')

df = pd.read_csv('X:/volume2/noam/passages/freqs/p1-37A.freqs', sep='\t')
df2 = prepare_freqs_for_error_analysis(df)
df2['mutation'] = df2.Ref + df2.Base
boxplot_mutations(df2, 'X:/volume2/noam/minion/p1-37A/variant_boxplot_miseq.png', 'Illumina-control')

##### error rates calculation ###########


def error_rates_from_freq(df, out_path):
    '''
    out_string = 'Accuracy: %s\nSubstitution rate: %s\nInsertion rate: %s\nDeletion rate: %s\nSimple substitution frequency average:%s\nSimple deletion frequency average:%s\n\nSimple substitution frequency: %s\n\nSimple deletion frequency:%s'
    
    # Substitutions: sum of incorrectly called aligned bases divided by 
    # all alligned bases (without deletions and insertions)
    a = df[(df.Ref != df.Base) & (df.Ref != '-') & (df.Base != '-')]
    a1 = (a.Read_count * a.Freq).sum()
    b = df[(df.Ref != '-') & (df.Base != '-')]
    b1 = (b.Read_count * b.Freq).sum()
    substitution_rate = a1 / b1

    # Accuracy: sum of correctly called and aligned bases divided by 
    # all aligned bases + insertions + deletions
    h = df[(df.Ref == df.Base) & (df.Ref != '-') & (df.Base != '-')]
    h1 = (h.Read_count * h.Freq).sum()
    c1 = (df.Read_count * df.Freq).sum()
    accuracy = h1 / c1

    # insertions out of all called bases (insertions + correct aligned)
    d = df[(df.Ref ==  '-')]
    d1 = (d.Read_count * d.Freq).sum()
    e = df[(df.Base != '-')]
    e1 = (e.Read_count * e.Freq).sum()
    insertion_rate = d1 / e1

    # deletions out of all alligned bases + deletions (without insertions)
    f = df[(df.Base == '-')]
    f1 = (f.Read_count * f.Freq).sum()
    g = df[(df.Ref != '-')] 
    g1 = (g.Read_count * g.Freq).sum()
    deletion_rate = f1 / g1
    '''
    df = df[~(df.Pos.isin(n))]
    # simple freqeuncy average for all substitution mutations
    simple_substitution_average = df[(df.Ref != df.Base) & (df.Ref != '-') & (df.Base != '-')].groupby('Pos').Freq.sum().reset_index().Freq.describe(percentiles=[.25, .5, .75, .9, .95])
    #simple_substitution_average2 = df[(df.Ref != df.Base) & (df.Ref != '-') & (df.Base != '-')].groupby('Pos').Freq.sum().reset_index().Freq.mean()
    
    # simple frequency average for all deletion mutations
    simple_deletion_average = df[(df.Ref != df.Base) & (df.Ref != '-') & (df.Base == '-')].Freq.describe(percentiles=[.25, .5, .75, .9, .95])
    #simple_deletion_average2 = df[(df.Ref != df.Base) & (df.Ref != '-') & (df.Base == '-')].Freq.mean()
    
    # simple_frequency average for insertions .1 mutations
    simple_insertion_average = df[(df.Base != df.Ref) & (df.Ref == '-')].groupby('Pos').Freq.sum().describe(percentiles=[.25, .5, .75, .9, .95])
    out_string = 'Substitutions:\n' + str(simple_substitution_average) + '\nDeletions:\n' + str(simple_deletion_average) + '\nInsertions:\n' + str(simple_insertion_average)
    with open(out_path, 'w') as f:
        f.write(out_string)
        #f.write(out_string%(accuracy, substitution_rate, insertion_rate, deletion_rate, simple_substitution_average2, simple_deletion_average2, simple_substitution_average, simple_deletion_average))
 

df = pd.read_csv('X:/volume2/noam/minion/p15-37B/pipeline_q0_b60/fastq.freqs', sep='\t')
df = prepare_freqs_for_error_analysis(df)
df['mutation'] = df.Ref + df.Base
error_rates_from_freq(df, 'X:/volume2/noam/minion/p15-37B/error_rate.txt')

df = pd.read_csv('X:/volume2/noam/minion/p15-37A/pipeline_q0_b60/fastq.freqs', sep='\t')
df = prepare_freqs_for_error_analysis(df)
df['mutation'] = df.Ref + df.Base
error_rates_from_freq(df, 'X:/volume2/noam/minion/p15-37A/error_rate.txt')


df = pd.read_csv('X:/volume2/noam/passages/freqs/p15-37B.freqs', sep='\t')
df = prepare_freqs_for_error_analysis(df)
df['mutation'] = df.Ref + df.Base
error_rates_from_freq(df, 'X:/volume2/noam/minion/p15-37B/error_rate_miseq.txt')


df = pd.read_csv('X:/volume2/noam/passages/freqs/p15-37A.freqs', sep='\t')
df = prepare_freqs_for_error_analysis(df)
df['mutation'] = df.Ref + df.Base
error_rates_from_freq(df, 'X:/volume2/noam/minion/p15-37A/error_rate_miseq.txt')

df = pd.read_csv('X:/volume2/noam/minion/p1-37A/pipeline_q0_b60/FAK41568.freqs', sep='\t')
df = prepare_freqs_for_error_analysis(df)
df['mutation'] = df.Ref + df.Base
error_rates_from_freq(df, 'X:/volume2/noam/minion/p1-37A/error_rate_minion.txt')

df = pd.read_csv('X:/volume2/noam/passages/freqs/p1-37A.freqs', sep='\t')
df = prepare_freqs_for_error_analysis(df)
df['mutation'] = df.Ref + df.Base
error_rates_from_freq(df, 'X:/volume2/noam/minion/p1-37A/error_rate_miseq.txt')


######## compare freuqencies of familiar mutations ###

FREQS_DICT_37B = {'MinION':'X:/volume2/noam/minion/p15-37B/pipeline_q0_b60/fastq.freqs',
              'MiSeq':'X:/volume2/noam/passages/freqs/p15-37B.freqs'}
c = compare_positions_between_freqs(FREQS_DICT_37B, positions_to_compare=[float(i[1:-1]) for i in mutations_37B]) 
c['Full_mutation'] = c.Ref + c.Pos.astype(str) + c.Base
c = c[['Full_mutation', 'Freq_MiSeq', 'Freq_MinION']]
c[c.Full_mutation.isin(mutations_37B)].to_csv('X:/volume2/noam/minion/p15-37B/compare_freqs.csv', index=False)


FREQS_DICT_37A = {'MinION':'X:/volume2/noam/minion/p15-37A/pipeline_q0_b60/fastq.freqs',
              'MiSeq':'X:/volume2/noam/passages/freqs/p15-37A.freqs'}
c = compare_positions_between_freqs(FREQS_DICT_37A, positions_to_compare=[float(i[1:-1]) for i in mutations_37A])
c['Full_mutation'] = c.Ref + c.Pos.astype(str) + c.Base
c = c[['Full_mutation', 'Freq_MiSeq', 'Freq_MinION']]
c[c.Full_mutation.isin(mutations_37A)].to_csv('X:/volume2/noam/minion/p15-37A/compare_freqs.csv', index=False)




############# association graphs #############33


def association_scatter(df, out_png, title, proximity_limit=15, start_pos=False, end_pos=False):
    '''
    For every position, plot all the chi square values, except for values of tests 
    between that position and positions [proximity_limit] bases away or closer (default 15).
    Can use start_pos and end_pos to zoom in on part of the genome.
    '''
    df = df[(df.pos1 - df.pos2).abs() > proximity_limit]
    if start_pos and end_pos:
        df = df[df.pos1.isin(range(start_pos, end_pos))]
    plot = df.plot(x='pos1', y='chi2', kind='scatter')
    plot.set_ylabel('Chi square statistic', fontsize=14)
    #plot.set_xlabel('Position in genome (bp)', fontsize=14)
    plot.set_xlabel('', fontsize=14)
    plot.set_title(title, fontsize=16, x=0.5, y=1.05)
    fig = plot.get_figure()
    fig.savefig(out_png, dpi=800, bbox_inches='tight')
    return

df = pd.read_csv('X:/volume2/noam/minion/p15-37B/association_tests/association_results.csv')
association_scatter(df, 'X:/volume2/noam/minion/p15-37B/association_tests/full_scatter.png', 'Genome Wide Chi Square Results')
association_scatter(df, 'X:/volume2/noam/minion/p15-37B/association_tests/scatter_300_600.png', 'Positions 300 to 600', start_pos=300, end_pos=600)
association_scatter(df, 'X:/volume2/noam/minion/p15-37B/association_tests/scatter_1400_1700.png', 'Positions 1400 to 1700', start_pos=1400, end_pos=1700)
association_scatter(df, 'X:/volume2/noam/minion/p15-37B/association_tests/scatter_1700_2000.png', 'Positions 1700 to 2000', start_pos=1700, end_pos=2000)
association_scatter(df, 'X:/volume2/noam/minion/p15-37B/association_tests/scatter_2900_3200.png', 'Positions 2900 to 3200', start_pos=2900, end_pos=3200)


df = pd.read_csv('X:/volume2/noam/minion/p15-37A/association_tests/association_results.csv')
association_scatter(df, 'X:/volume2/noam/minion/p15-37A/association_tests/full_scatter.png', 'Genome Wide Chi Square Results')
association_scatter(df, 'X:/volume2/noam/minion/p15-37A/association_tests/scatter_500_700.png', 'Positions 500 to 700', start_pos=500, end_pos=700)
association_scatter(df, 'X:/volume2/noam/minion/p15-37A/association_tests/scatter_900_1200.png', 'Positions 900 to 1200', start_pos=900, end_pos=1200)
association_scatter(df, 'X:/volume2/noam/minion/p15-37A/association_tests/scatter_1600_1800.png', 'Positions 1600 to 1800', start_pos=1600, end_pos=1800)

df = pd.read_csv('X:/volume2/noam/minion/p1-37A/association_tests/association_results.csv')
association_scatter(df, 'X:/volume2/noam/minion/p1-37A/association_tests/full_scatter.png', 'Genome Wide Chi Square Results')


def scaled_back_full_genome_heatmap(df,output_png, title, scale=10):
    '''
    Plot a heatmap for the whole genome, but scale it down. Use one value for
    every x by x square of positions, where x is defined by scale (default 10).
    '''
    df = df.drop_duplicates()
    df['pos1_scaled'] = df.pos1 // scale * scale
    df['pos2_scaled'] = df.pos2 // scale * scale
    df['scaled_chi2_max'] = df.groupby(['pos1_scaled', 'pos2_scaled']).chi2.transform('max')
    #df['scaled_chi2_max'] = df.groupby(['pos1_scaled', 'pos2_scaled']).zscore.transform('max')
    df = df[['pos1_scaled', 'pos2_scaled', 'scaled_chi2_max']].drop_duplicates()
    df['scaled_chi2_max_log'] = np.log(df.scaled_chi2_max)
    
    df = df[~((df.pos1_scaled == df.pos2_scaled) | ((df.pos2_scaled - df.pos1_scaled).abs() <= scale))]
    
    pivot = df.drop_duplicates().pivot(index='pos1_scaled', columns='pos2_scaled', values='scaled_chi2_max_log')
    fig, ax = plt.subplots(figsize=(20, 20))
    sns.heatmap(pivot, ax=ax, cbar_kws={'label': 'Ln(chi squared statistic)'}, xticklabels=20, yticklabels=20, square=True)
    ax.figure.axes[-1].yaxis.label.set_size(24)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.set_ylabel('Position in Genome', fontsize=24)
    ax.set_xlabel('Position in Genome', fontsize=24)
    ax.tick_params(labelsize=20)
    ax.set_title(title + '\n', fontsize=30)
    fig.savefig(output_png, dpi=800, bbox_inches='tight')
    return

df = pd.read_csv('X:/volume2/noam/minion/p15-37B/association_tests/association_results.csv')
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p15-37B/association_tests/full_heatmap.png', 'p15B', scale=20)
df = pd.read_csv('X:/volume2/noam/minion/p15-37A/association_tests/association_results.csv')
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p15-37A/association_tests/full_heatmap.png', 'p15A', scale=20)
df = pd.read_csv('X:/volume2/noam/minion/p1-37A/association_tests/association_results.csv')
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p1-37A/association_tests/full_heatmap.png', 'p1A', scale=20)

################ Ranked 1 errors
g = pd.read_csv('X:/volume2/noam/minion/p15-37B/pipeline_q0_b60/fastq.freqs', '\t')
g2 = pd.read_csv('X:/volume2/noam/passages/freqs/p15-37B.freqs', '\t')
pd.merge(g[(g.Rank == 0) & (g.Base != g.Ref) & (g.Ref != '-')], g2, how='left', on= ['Pos', 'Ref', 'Base']).sort_values('Freq_y')

#### variants that exceeded 10%
g = pd.read_csv('X:/volume2/noam/minion/p15-37B/pipeline_q0_b60/fastq.freqs', '\t')
g = pd.read_csv('X:/volume2/noam/passages/freqs/p15-37B.freqs', '\t')
g = pd.read_csv('X:/volume2/noam/passages/freqs/p15-37A.freqs', '\t')
g = pd.read_csv('X:/volume2/noam/minion/p15-37A/pipeline_q0_b60/fastq.freqs', sep='\t')
g = pd.read_csv('X:/volume2/noam/passages/freqs/p1-37A.freqs', '\t')
g = pd.read_csv('X:/volume2/noam/minion/p1-37A/pipeline_q0_b60/FAK41568.freqs', '\t')

len(g[(g.Base != g.Ref) & (g.Ref != '-') & (g.Freq >= 0.1) & ~(g.Pos.isin(n))])
len(g[(g.Base != g.Ref) & (g.Ref != '-') & (g.Freq >= 0.1) & ~(g.Pos.isin(n))].Pos.drop_duplicates())

len(g[(g.Base != g.Ref) & (g.Ref != '-') & (g.Freq >= 0.01) & ~(g.Pos.isin(n))])
len(g[(g.Base != g.Ref) & (g.Ref != '-') & (g.Freq >= 0.01) & ~(g.Pos.isin(n))].Pos.drop_duplicates())

########## version 3 ################

### associations ######

def scaled_back_full_genome_heatmap(df,output_png, title, scale=10):
    '''
    Plot a heatmap for the whole genome, but scale it down. Use one value for
    every x by x square of positions, where x is defined by scale (default 10).
    '''
    df = df.drop_duplicates()
    df = df[(df.pos1 - df.pos2).abs() > 15]
    df['pos1_scaled'] = df.pos1 // scale * scale
    df['pos2_scaled'] = df.pos2 // scale * scale
    df['scaled_chi2_max'] = df.groupby(['pos1_scaled', 'pos2_scaled']).chi2.transform('max')
    #df['scaled_chi2_max'] = df.groupby(['pos1_scaled', 'pos2_scaled']).zscore.transform('max')
    df = df[['pos1_scaled', 'pos2_scaled', 'scaled_chi2_max']].drop_duplicates()
    df['scaled_chi2_max_log'] = np.log(df.scaled_chi2_max)
    
    #df = df[~((df.pos1_scaled == df.pos2_scaled) | ((df.pos2_scaled - df.pos1_scaled).abs() <= scale))]
    #df = df[~(df.pos1_scaled == df.pos2_scaled)]
    pivot = df.drop_duplicates().pivot(index='pos1_scaled', columns='pos2_scaled', values='scaled_chi2_max_log')
    fig, ax = plt.subplots(figsize=(12, 12))
    #sns.heatmap(pivot, ax=ax, cbar_kws={'label': 'Ln(chi squared statistic)'}, xticklabels=40, yticklabels=40, square=True, vmin=1, vmax=7)
    sns.heatmap(pivot, ax=ax, cbar_kws={'label': 'Ln(chi squared statistic)'}, xticklabels=40, yticklabels=40, square=True, vmin=1, vmax=9.5)
    ax.figure.axes[-1].yaxis.label.set_size(24)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=24)
    ax.set_ylabel('Position in Genome', fontsize=30)
    ax.set_xlabel('Position in Genome', fontsize=30)
    ax.tick_params(labelsize=24)
    ax.set_title(title + '\n', fontsize=34)
    fig.savefig(output_png, dpi=800, bbox_inches='tight')
    return

df = pd.read_csv('X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.csv')
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.full_heatmap.png', 'p15B', scale=20)
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.full_heatmap_30.png', 'p15-37B', scale=30)
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.full_heatmap_40.png', 'p15-37B', scale=40)
df = pd.read_csv('X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.csv')
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.full_heatmap.png', 'p15A', scale=20)
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.full_heatmap_30.png', 'p15-37A', scale=30)
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.full_heatmap_40.png', 'p15-37A', scale=40)
df = pd.read_csv('X:/volume2/noam/minion/p1-37A/association_tests_only_long_reads/association_results_30_3539.csv')
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p1-37A/association_tests_only_long_reads/association_results_30_3539.full_heatmap.png', 'control', scale=20)
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p1-37A/association_tests_only_long_reads/association_results_30_3539.full_heatmap_30.png', 'p1-37A', scale=30)
scaled_back_full_genome_heatmap(df, 'X:/volume2/noam/minion/p1-37A/association_tests_only_long_reads/association_results_30_3539.full_heatmap_40.png', 'p1-37A', scale=40)



def scaled_back_full_genome_heatmap_modified_zscore(df,output_png, title, scale=10):
    '''
    Plot a heatmap for the whole genome, but scale it down. Use one value for
    every x by x square of positions, where x is defined by scale (default 10).
    '''
    df = df.drop_duplicates()
    df = df[(df.pos1 - df.pos2).abs() > 15]
    df['pos1_scaled'] = df.pos1 // scale * scale
    df['pos2_scaled'] = df.pos2 // scale * scale
    df['scaled_chi2_max'] = df.groupby(['pos1_scaled', 'pos2_scaled']).modified_zscore.transform('max')
    #df['scaled_chi2_max'] = df.groupby(['pos1_scaled', 'pos2_scaled']).zscore.transform('max')
    df = df[['pos1_scaled', 'pos2_scaled', 'scaled_chi2_max']].drop_duplicates()
    df['scaled_chi2_max_log'] = np.log(df.scaled_chi2_max)
    
    #df = df[~((df.pos1_scaled == df.pos2_scaled) | ((df.pos2_scaled - df.pos1_scaled).abs() <= scale))]
    #df = df[~(df.pos1_scaled == df.pos2_scaled)]
    pivot = df.drop_duplicates().pivot(index='pos1_scaled', columns='pos2_scaled', values='scaled_chi2_max_log')
    fig, ax = plt.subplots(figsize=(20, 20))
    #sns.heatmap(pivot, ax=ax, cbar_kws={'label': 'Ln(normalized chi square statistic)'}, xticklabels=20, yticklabels=20, square=True, vmin=1, vmax=7)
    sns.heatmap(pivot, ax=ax, cbar_kws={'label': 'Ln(normalized chi square statistic)'}, xticklabels=20, yticklabels=20, square=True, vmin=1, vmax=9)
    ax.figure.axes[-1].yaxis.label.set_size(24)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.set_ylabel('Position in Genome', fontsize=24)
    ax.set_xlabel('Position in Genome', fontsize=24)
    ax.tick_params(labelsize=20)
    ax.set_title(title + '\n', fontsize=30)
    fig.savefig(output_png, dpi=800, bbox_inches='tight')
    return

df = pd.read_csv('X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv')
scaled_back_full_genome_heatmap_modified_zscore(df, 'X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.modified_ztest.full_heatmap.png', 'p15-37B', scale=20)
df = pd.read_csv('X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv')
scaled_back_full_genome_heatmap_modified_zscore(df, 'X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.full_heatmap.png', 'p15-37A', scale=20)
df = pd.read_csv('X:/volume2/noam/minion/p1-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv')
scaled_back_full_genome_heatmap_modified_zscore(df, 'X:/volume2/noam/minion/p1-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.full_heatmap.png', 'p1-37A', scale=20)




def association_scatter(df, out_png, title, proximity_limit=15, start_pos=False, end_pos=False):
    '''
    For every position, plot all the chi square values, except for values of tests 
    between that position and positions [proximity_limit] bases away or closer (default 15).
    Can use start_pos and end_pos to zoom in on part of the genome.
    '''
    df = df[(df.pos1 - df.pos2).abs() > proximity_limit]
    if start_pos and end_pos:
        df = df[df.pos1.isin(range(start_pos, end_pos))]
    plot = df.plot(x='pos1', y='modified_zscore', kind='scatter', figsize=(10, 2))
    plot.set_ylabel('Normalized chi square statistic', fontsize=14)
    plot.set_xlabel('\nPosition in genome (bp)', fontsize=14)
    #plot.set_xlabel('', fontsize=14)
    plot.set_title(title, fontsize=16, x=0.5, y=1.05)
    fig = plot.get_figure()
    fig.savefig(out_png, dpi=800, bbox_inches='tight')
    return

df = pd.read_csv('X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv')
association_scatter(df, 'X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.modified_ztest.full_scatter.png', 'Genome Wide Chi Square Results')
association_scatter(df, 'X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.modified_ztest.scatter_1400_2000.png', 'Positions 1400 to 2000', start_pos=1400, end_pos=2000)


df = pd.read_csv('X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv')
association_scatter(df, 'X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.full_scatter.png', 'Genome Wide Chi Square Results')

df = pd.read_csv('X:/volume2/noam/minion/p1-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv')
association_scatter(df, 'X:/volume2/noam/minion/p1-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.full_scatter.png', 'Control')


df = pd.read_csv('X:/volume2/noam/minion/control/association_test_only_long_reads/association_results.modified_ztest.csv')
association_scatter(df, 'X:/volume2/noam/minion/control/association_test_only_long_reads/association_results.modified_ztest.full_scatter.png', 'Genome Wide Chi Square Results')



## association scatter with colors

def association_scatter(df, modified_zscore_peaks, out_png, title, gold_standard_mutation_list, cutoff, proximity_limit=15, start_pos=False, end_pos=False):
    '''
    For every position, plot all the chi square values, except for values of tests 
    between that position and positions [proximity_limit] bases away or closer (default 15).
    Can use start_pos and end_pos to zoom in on part of the genome.
    '''
    df = df[(df.pos1 - df.pos2).abs() > proximity_limit]
    df = pd.merge(df, modified_zscore_peaks[['pos1', 'pos2', 'is_peak']], how='left', on=['pos1', 'pos2'])
    if start_pos and end_pos:
        df = df[df.pos1.isin(range(start_pos, end_pos))]
    false_positives = df[(df.is_peak == True) & ~(df.pos1.isin(gold_standard_mutation_list))]
    true_positives = df[(df.is_peak == True) & (df.pos1.isin(gold_standard_mutation_list))]
    other_associations = df[df.is_peak != True]
    fig, axes = plt.subplots(nrows=1, ncols=1)
    for dots, color in zip([other_associations, true_positives, false_positives], ['#075E94', '#13CB12', 'red']):
        plot = dots.plot(x='pos1', y='modified_zscore', kind='scatter', figsize=(10, 2), marker='o', color=color, ax=axes)
    axes.axhline(y=cutoff, color='r', linestyle='--')
    plot.set_ylabel('Normalized chi square statistic', fontsize=14)
    plot.set_xlabel('\nPosition in genome (bp)', fontsize=14)
    #plot.set_xlabel('', fontsize=14)
    plot.set_title(title + '\n', fontsize=16, x=0.5, y=1.05)
    fig = plot.get_figure()
    fig.savefig(out_png, dpi=800, bbox_inches='tight')
    return


joined = pd.read_csv('X:/volume2/noam/passages/201902/all_freqs.csv')
joined['Full_mutation'] = joined.Ref + joined.Pos.astype(str) + joined.Base
joined = joined[joined.Time == 15]
joined = joined[joined.Degree == 37]
mutations_37A = joined[(joined.Replica == 'A') & (joined.Time == 15) & (joined.Base != joined.Ref) & (joined.Ref != '-') & ~(joined.Pos.isin(n)) & (joined.Freq >= 0.01)].sort_values('Freq', ascending=False).Pos.tolist()
mutations_37B = joined[(joined.Replica == 'B') & (joined.Time == 15) & (joined.Base != joined.Ref) & (joined.Ref != '-') & ~(joined.Pos.isin(n)) & (joined.Freq >= 0.01)].sort_values('Freq', ascending=False).Pos.tolist()

peaks_37B = pd.read_csv('X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.modified_ztest.peaks_z114.csv')
peaks_37A = pd.read_csv('X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.peaks_z114.csv')

df = pd.read_csv('X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv')
association_scatter(df, peaks_37B,'X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.modified_ztest.full_scatter_colors.png', 'p15B', mutations_37B, 114)
#association_scatter(df, peaks_37B,'X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.modified_ztest.scatter_1400_2000_colors.png', 'p15-37B - Positions 1400 to 2000', mutations_37B, 114, start_pos=1400, end_pos=2000,)


df = pd.read_csv('X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv')
association_scatter(df, peaks_37A, 'X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.full_scatter_colors.png', 'p15A', mutations_37A, 114)


## enolase

def association_scatter(df, out_png, title, proximity_limit=15, start_pos=False, end_pos=False):
    '''
    For every position, plot all the chi square values, except for values of tests 
    between that position and positions [proximity_limit] bases away or closer (default 15).
    Can use start_pos and end_pos to zoom in on part of the genome.
    '''
    df = df[(df.pos1 - df.pos2).abs() > proximity_limit]
    if start_pos and end_pos:
        df = df[df.pos1.isin(range(start_pos, end_pos))]
    fig, axes = plt.subplots(nrows=1, ncols=1)
    for dots, color in zip([df[(df.pos1.isin([363, 606])) &  (df.pos2.isin([363, 606]))], df[~((df.pos1.isin([363, 606])) &  (df.pos2.isin([363, 606])))]], ['#13CB12', '#940786']):
        plot = dots.plot(x='pos1', y='modified_zscore', kind='scatter', figsize=(10, 2), marker='o', color=color, ax=axes)
    plot.set_ylabel('Normalized chi square statistic', fontsize=14)
    plot.set_xlabel('\nPosition in genome (bp)', fontsize=14)
    #plot.set_xlabel('', fontsize=14)
    plot.set_title(title, fontsize=16, x=0.5, y=1.05)
    fig = plot.get_figure()
    fig.savefig(out_png, dpi=800, bbox_inches='tight')
    return

df = pd.read_csv('X:/volume2/noam/minion/control/association_test_only_long_reads/association_results.modified_ztest.csv')
association_scatter(df, 'X:/volume2/noam/minion/control/association_test_only_long_reads/association_results.modified_ztest.full_scatter_color.png', 'Yeast enolase mRNA')


###########3

def association_scatter_zika(df, out_png, title, proximity_limit=15, start_pos=False, end_pos=False):
    '''
    For every position, plot all the chi square values, except for values of tests 
    between that position and positions [proximity_limit] bases away or closer (default 15).
    Can use start_pos and end_pos to zoom in on part of the genome.
    '''
    df = df[(df.pos1 - df.pos2).abs() > proximity_limit]
    if start_pos and end_pos:
        df = df[df.pos1.isin(range(start_pos, end_pos))]
    fig, axes = plt.subplots(nrows=1, ncols=1)
    for dots, color in zip([df[(df.pos1.isin([1309, 1315, 1414, 1447, 1471, 1627])) &  (df.pos2.isin([1309, 1315, 1414, 1447, 1471, 1627]))], df[~((df.pos1.isin([1309, 1315, 1414, 1447, 1471, 1627])) &  (df.pos2.isin([1309, 1315, 1414, 1447, 1471, 1627])))]], ['#13CB12', '#940786']):
        plot = dots.plot(x='pos1', y='modified_zscore', kind='scatter', figsize=(10, 2), marker='o', color=color, ax=axes)
    plot.set_ylabel('Normalized chi square statistic', fontsize=14)
    plot.set_xlabel('\nPosition in genome (bp)', fontsize=14)
    #plot.set_xlabel('', fontsize=14)
    plot.set_title(title, fontsize=16, x=0.5, y=1.05)
    fig = plot.get_figure()
    fig.savefig(out_png, dpi=800, bbox_inches='tight')
    return

df = pd.read_csv('Z:/volume1/noam/zika_minion/zika_1_association_1229_1695/association_results.1259-1665.mztest_4peaks.csv')
association_scatter_zika(df, 'Z:/volume1/noam/zika_minion/zika_1_association_1229_1695/association_results.1259-1665.mztest_4peaks.full_scatter_color.png', 'Zika')

################################ mutation frequenices graph


COLORS = {'T1764.0-':'#F50202', 'A1664.0G':'#F49D09', 'A535.0G':'#5EC0D2', 'T1440.0C':'#F1F87A', 'T1440.0G':'#C4A0E4', 'A1443.0G':'#FE1D1D', 'A1611.0G':'#327CFD', 'C1724.0T':'#8FD95A', 'A1744.0G':'#FBB3DD', 'G1560.0A':'#A7F0A2', 'G1906.0A':'#A3A3A3', 'C3358.0T':'#26451C', 'G3114.0A':'#B37A42', 'A1770.0G':'#7603BA', 'G2310.0A':'#033E86', 'A2626.0G':'#8FD95A', 'C3299.0T':'#211785', 'C1718.0T':'#DFC236', 'T862.0C':'#880E05', 'A2790.0T':'#DF36C6', 'G1736.0A':'#CFFD2F', 'C1549.0T':'#2CA403', 'G531.0A':'#972FFE', 'C1050.0T':'#13B908'}
#n = [17, 18, 19, 20, 21, 22, 3462, 3463, 3524, 3542, 3543, 3544, 3545, 3546, 3547, 3548, 3564, 3566]
n = list(range(30)) + list(range(3540, 3570))
 

joined = pd.read_csv('X:/volume2/noam/passages/201902/all_freqs.csv')
joined['Full_mutation'] = joined.Ref + joined.Pos.astype(int).astype(str) + joined.Base
#
joined = joined[joined.Time == 15]
joined = joined[joined.Degree == 37]
joined['Sample'] = 'Illumina-p' + joined.Time.astype(str) + joined.Replica

mutations = joined[(joined.Time == 15) & (joined.Base != joined.Ref) & (joined.Ref != '-') & ~(joined.Pos.isin(n)) & (joined.Freq >= 0.1)].sort_values('Freq', ascending=False).Full_mutation.drop_duplicates().tolist()
joined = joined[joined.Full_mutation.isin(mutations)].sort_values('Freq', ascending=False)


fig, axes = plt.subplots(nrows=1, ncols=2)
for sample, a in zip(['Illumina-p15A', 'Illumina-p15B'], axes):
    sns.stripplot(x='Sample', y='Freq', hue='Full_mutation', data=joined[joined.Sample==sample].sort_values('Freq', ascending=False), jitter=True, size=10, palette={i.replace('.0', ''):COLORS[i] for i in COLORS}, ax=a)
    a.set_ylabel('Frequency')
    a.set_xlabel('')
    a.minorticks_on()
    a.grid(which='major', alpha=0.7, axis='y')
    a.get_legend().remove()
    a.set_ylim(-0.02,0.6)
a.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
#fig.suptitle('Mutations After Evolution of 15 Passages')
fig.set_size_inches(5, 3)
fig.subplots_adjust(wspace=0.45)
fig.savefig('X:/volume2/noam/passages/201902/only_p15_for_minion5.png', bbox_inches='tight', dpi=800)





############ roc curves
df1 = pd.read_csv('X:/volume2/noam/minion/p15-37B/roc_curve_naive_frequncy.png.csv')
df1['Sample'] = 'p15B'
df1['method'] = 'Frequencies'
df2 = pd.read_csv('X:/volume2/noam/minion/p15-37A/roc_curve_naive_frequncy.png.csv')
df2['Sample'] = 'p15A'
df2['method'] = 'Frequencies'
df3 = pd.read_csv('X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/roc_curve.png.csv')
df3['Sample'] = 'p15A'
df3['method'] = 'AssociVar'
df4 = pd.read_csv('X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/roc_curve.png.csv')
df4['Sample'] = 'p15B'
df4['method'] = 'AssociVar'

df = pd.concat([df1, df2, df3, df4])


fig, axes = plt.subplots(nrows=1, ncols=2)
colors={'10%':'#028ABC', '5%':'orange', '1%':'green'}
for sample, a in zip(['p15A', 'p15B'], axes):
    df_sample = df[df.Sample == sample]
    for c in [0.1, 0.05, 0.01]:
        df_sample[(df_sample.method == 'AssociVar') & (df_sample.gold_standard_cutoff == c)].plot(x='fpr', y='tpr', label=str(int(c*100)) + '% (AssociVar)' , ax=a, color=colors[str(int(c*100)) + '%']) 
    for c in [0.1, 0.05, 0.01]:
        df_sample[(df_sample.method == 'Frequencies') & (df_sample.gold_standard_cutoff == c)].plot(x='fpr', y='tpr', label=str(int(c*100)) + '% (naive)' , ax=a, linestyle='--', color=colors[str(int(c*100)) + '%']) 
        a.set_xlabel('False Positive Rate', fontsize = 14)
        a.set_ylabel('True Positive Rate', fontsize = 14)
        a.set_title(sample, fontsize=16, x=0.5, y=1.05)
        a.get_legend().remove()
a.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)        
fig.set_size_inches(11, 4)
fig.subplots_adjust(wspace=0.3)
fig.savefig('X:/volume2/noam/minion/associvar_naive_roc_curve.png', bbox_inches='tight', dpi=800)
    


############################### Tombo

import matplotlib.pyplot as plt

### enolase distributions

a = pd.read_csv('Z:/volume1/noam/minion/tombo/enolase-sample/enolase.alt_model.5mC.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
fig, ax = plt.subplots(nrows=1, ncols=1)
a.plot(x='position', y='dampened_fraction', kind='scatter', ax=ax, marker='.')
ax.set_ylabel('Read fraction modified')
ax.set_xlabel('Position')
ax.set_title('5mC Modified Base Detection')
fig.set_size_inches(4,3)
fig.savefig('Z:/volume1/noam/minion/tombo/plots/enolase_5mc_scatter.png', bbox_inches='tight', dpi=800)


fig, ax = plt.subplots(nrows=1, ncols=1)
sns.boxplot(a.dampened_fraction, ax=ax)
ax.set_xlabel('Read fraction modified')
ax.set_title('5mC Modified Base Detection')
fig.set_size_inches(4,3)
fig.savefig('Z:/volume1/noam/minion/tombo/plots/enolase_5mc_boxplot.png', bbox_inches='tight', dpi=800)



a = pd.read_csv('Z:/volume1/noam/minion/tombo/enolase-sample/enolase.de_novo.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
fig, ax = plt.subplots(nrows=1, ncols=1)
a.plot(x='position', y='dampened_fraction', kind='scatter', ax=ax, marker='.')
ax.set_ylabel('Read fraction modified')
ax.set_xlabel('Position')
ax.set_title('De Novo Modified Base Detection')
fig.set_size_inches(4,3)
fig.savefig('Z:/volume1/noam/minion/tombo/plots/enolase_de_novo_scatter.png', bbox_inches='tight', dpi=800)


fig, ax = plt.subplots(nrows=1, ncols=1)
sns.boxplot(a.dampened_fraction, ax=ax)
ax.set_xlabel('Read fraction modified')
ax.set_title('De Novo Modified Base Detection')
fig.set_size_inches(4,3)
fig.savefig('Z:/volume1/noam/minion/tombo/plots/enolase_de_novo_boxplot.png', bbox_inches='tight', dpi=800)



### distributions
p1 = pd.read_csv('Z:/volume1/noam/minion/tombo/p1-37A-sample/p1-37A.de_novo.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
p1['sample'] = 'p1A'
p15a = pd.read_csv('Z:/volume1/noam/minion/tombo/p15-37A-sample/p15-37A.de_novo.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
p15a['sample'] = 'p15A'
p15b = pd.read_csv('Z:/volume1/noam/minion/tombo/p15-37B-sample/p15-37B.de_novo.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
p15b['sample'] = 'p15B'
a = pd.read_csv('Z:/volume1/noam/minion/tombo/enolase-sample/enolase.de_novo.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
a['sample'] = 'enolase'
fig, ax = plt.subplots(nrows=1, ncols=1)
sns.boxplot(x='dampened_fraction', y='sample', data=pd.concat([p1, p15a, p15b, a]), ax=ax)
ax.set_xlabel('Read fraction modified')
ax.set_title('De Novo Modified Base Detection')
fig.set_size_inches(4,3)
fig.savefig('Z:/volume1/noam/minion/tombo/plots/de_novo_distributions.png', bbox_inches='tight', dpi=800)


p1 = pd.read_csv('Z:/volume1/noam/minion/tombo/p1-37A-sample/p1-37A.alt_model.5mC.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
p1['sample'] = 'p1A'
p15a = pd.read_csv('Z:/volume1/noam/minion/tombo/p15-37A-sample/p15-37A.alt_model.5mC.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
p15a['sample'] = 'p15A'
p15b = pd.read_csv('Z:/volume1/noam/minion/tombo/p15-37B-sample/p15-37B.alt_model.5mC.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
p15b['sample'] = 'p15B'
a = pd.read_csv('Z:/volume1/noam/minion/tombo/enolase-sample/enolase.alt_model.5mC.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
a['sample'] = 'enolase'
fig, ax = plt.subplots(nrows=1, ncols=1)
sns.boxplot(x='dampened_fraction', y='sample', data=pd.concat([p1, p15a, p15b, a]), ax=ax)
ax.set_xlabel('Read fraction modified')
ax.set_title('5mC Modified Base Detection')
fig.set_size_inches(4,3)
fig.savefig('Z:/volume1/noam/minion/tombo/plots/5mc_distributions.png', bbox_inches='tight', dpi=800)



### correlation associvar to tombo
def correlation_tombo_zscores(tombo_dampened_fractions_file, associvar_modified_zscores_file, miseq_freqs, out_path, title):
    f = pd.read_csv(miseq_freqs, '\t')
    remove_positions = f[(f.Base != f.Ref) & (f.Ref != '-') & ~(f.Pos.isin(n)) & (f.Freq >= 0.01)].Pos.tolist()
    tombo = pd.read_csv(tombo_dampened_fractions_file, ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
    z = pd.read_csv(associvar_modified_zscores_file)
    z_max = pd.concat([z[['pos1', 'modified_zscore']].rename(columns={'pos1':'position'}), z[['pos2', 'modified_zscore']].rename(columns={'pos2':'position'})]).groupby('position').modified_zscore.max().reset_index()
    df = pd.merge(z_max, tombo,  on='position')
    df = df[~(df.position.isin(remove_positions))]
    
    fig, ax = plt.subplots(nrows=1, ncols=1)
    sns.regplot(x='dampened_fraction', y='modified_zscore', data=df, line_kws={"color": "red"}, ax=ax)
    ax.set_xlabel('Read fraction modified')
    ax.set_ylabel('Normalized chi square statistic')
    ax.set_title(title)
    fig.set_size_inches(4,3)
    fig.savefig(out_path, bbox_inches='tight', dpi=800)
    return

correlation_tombo_zscores('Z:/volume1/noam/minion/tombo/p1-37A-sample/p1-37A.alt_model.5mC.browser.dampened_fraction_modified_reads.plus.wig', 'X:/volume2/noam/minion/p1-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv', 'X:/volume2/noam/passages/201907/freqs/p1-37A.freqs', 'Z:/volume1/noam/minion/tombo/plots/p1A_5mc_correlation_associvar_vs_tombo.png', 'p1A - 5mC')
correlation_tombo_zscores('Z:/volume1/noam/minion/tombo/p1-37A-sample/p1-37A.de_novo.browser.dampened_fraction_modified_reads.plus.wig', 'X:/volume2/noam/minion/p1-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv', 'X:/volume2/noam/passages/201907/freqs/p1-37A.freqs', 'Z:/volume1/noam/minion/tombo/plots/p1A_de_novo_correlation_associvar_vs_tombo.png', 'p1A - De Novo')
correlation_tombo_zscores('Z:/volume1/noam/minion/tombo/p15-37A-sample/p15-37A.alt_model.5mC.browser.dampened_fraction_modified_reads.plus.wig', 'X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv', 'X:/volume2/noam/passages/201907/freqs/p15-37A.freqs', 'Z:/volume1/noam/minion/tombo/plots/p15A_5mc_correlation_associvar_vs_tombo.png', 'p15A - 5mC')
correlation_tombo_zscores('Z:/volume1/noam/minion/tombo/p15-37A-sample/p15-37A.de_novo.browser.dampened_fraction_modified_reads.plus.wig', 'X:/volume2/noam/minion/p15-37A/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv', 'X:/volume2/noam/passages/201907/freqs/p15-37A.freqs', 'Z:/volume1/noam/minion/tombo/plots/p15A_de_novo_correlation_associvar_vs_tombo.png', 'p15A - De Novo')
correlation_tombo_zscores('Z:/volume1/noam/minion/tombo/p15-37B-sample/p15-37B.alt_model.5mC.browser.dampened_fraction_modified_reads.plus.wig', 'X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv', 'X:/volume2/noam/passages/201907/freqs/p15-37B.freqs', 'Z:/volume1/noam/minion/tombo/plots/p15B_5mc_correlation_associvar_vs_tombo.png', 'p15B - 5mC')
correlation_tombo_zscores('Z:/volume1/noam/minion/tombo/p15-37B-sample/p15-37B.de_novo.browser.dampened_fraction_modified_reads.plus.wig', 'X:/volume2/noam/minion/p15-37B/association_tests_only_long_reads/association_results_30_3539.modified_ztest.csv', 'X:/volume2/noam/passages/201907/freqs/p15-37B.freqs', 'Z:/volume1/noam/minion/tombo/plots/p15B_de_novo_correlation_associvar_vs_tombo.png', 'p15B - De Novo')


## correlation between three tombos

joined = pd.read_csv('X:/volume2/noam/passages/201902/all_freqs.csv')
joined['Full_mutation'] = joined.Ref + joined.Pos.astype(str) + joined.Base
joined = joined[joined.Degree == 37]
mutations_37A = joined[(joined.Replica == 'A') & (joined.Time == 15) & (joined.Base != joined.Ref) & (joined.Ref != '-') & ~(joined.Pos.isin(n)) & (joined.Freq >= 0.01)].Pos.tolist()
mutations_37B = joined[(joined.Replica == 'B') & (joined.Time == 15) & (joined.Base != joined.Ref) & (joined.Ref != '-') & ~(joined.Pos.isin(n)) & (joined.Freq >= 0.01)].Pos.tolist()
mutations_p1 = joined[(joined.Replica == 'A') & (joined.Time == 1) & (joined.Base != joined.Ref) & (joined.Ref != '-') & ~(joined.Pos.isin(n)) & (joined.Freq >= 0.01)].Pos.tolist()



a = pd.read_csv('Z:/volume1/noam/minion/tombo/p1-37A-sample/p1-37A.de_novo.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
b = pd.read_csv('Z:/volume1/noam/minion/tombo/p15-37A-sample/p15-37A.de_novo.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
c = pd.read_csv('Z:/volume1/noam/minion/tombo/p15-37B-sample/p15-37B.de_novo.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])

d = pd.merge(a, (pd.merge(b, c, on='position', suffixes=('_p15A', '_p15B'))), on='position')
d = d[~(d.position.isin(mutations_37A + mutations_37B + mutations_p1)) & ~(d.position.isin(n))]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(d['dampened_fraction'], d['dampened_fraction_p15A'], d['dampened_fraction_p15B'], marker='.', color='orange')
ax.azim = 190
ax.elev = 20
ax.set_xlabel('p1A')
ax.set_ylabel('p15A')
ax.set_zlabel('p15B')
ax.set_title('De Novo Modified Base Detection')
fig.savefig('Z:/volume1/noam/minion/tombo/plots/de_novo_compare_MS2.png', bbox_inches='tight', dpi=800)

a = pd.read_csv('Z:/volume1/noam/minion/tombo/p1-37A-sample/p1-37A.alt_model.5mC.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
b = pd.read_csv('Z:/volume1/noam/minion/tombo/p15-37A-sample/p15-37A.alt_model.5mC.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])
c = pd.read_csv('Z:/volume1/noam/minion/tombo/p15-37B-sample/p15-37B.alt_model.5mC.browser.dampened_fraction_modified_reads.plus.wig', ' ', skiprows=2, header=None, names=['position', 'dampened_fraction'])

d = pd.merge(a, (pd.merge(b, c, on='position', suffixes=('_p15A', '_p15B'))), on='position')
d = d[~(d.position.isin(mutations_37A + mutations_37B + mutations_p1)) & ~(d.position.isin(n))]
print(len(d))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(d['dampened_fraction'], d['dampened_fraction_p15A'], d['dampened_fraction_p15B'], marker='.', color='orange')
ax.azim = 190
ax.elev = 20
ax.set_xlabel('p1A')
ax.set_ylabel('p15A')
ax.set_zlabel('p15B')
ax.set_title('5mC Modified Base Detection')
fig.savefig('Z:/volume1/noam/minion/tombo/plots/5mC_compare_MS2.png', bbox_inches='tight', dpi=800)

