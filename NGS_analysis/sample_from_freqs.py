"""
This script simulates a bottleneck effect on a population represented by a 
freqs file. Before sampling, it normalizes the population using an important 
assumption: we only have two allelles per position! 
Does not deal with insertions.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def normalize_two_alleles(freqs_path):
    '''
    keeps only the two highest frequency alleles per position and recalculates
    their frequency accordingly
    '''
    df = pd.read_csv(freqs_path, '\t')
    df = df[df.Rank.isin([0,1])]
    df = df[df.Ref != '-']
    df['pos_sum'] = df.groupby('Pos').Freq.transform(sum)
    df['Freq'] = df.Freq / df.pos_sum
    df.Read_count = df.pos_sum * df.Read_count
    df.drop(columns=['pos_sum'], inplace=True)
    return df

def sample_freqs(freqs_path, sample_size):
    '''
    Gets a path to a freqs file, and return a df after sampling n allelles for 
    each position, where n equals sample size.
    '''
    df = normalize_two_alleles(freqs_path)
    major_allele = df[df.Rank == 0].copy()
    major_allele['Sampled_freq'] = major_allele['Freq'].apply(lambda x: np.random.binomial(sample_size, x) / sample_size)
    minor_allele = df[df.Rank == 1].copy()
    minor_allele = pd.merge(minor_allele, major_allele[['Pos', 'Sampled_freq']], on='Pos')
    minor_allele['Sampled_freq'] = 1 - minor_allele['Sampled_freq']
    return pd.concat([major_allele, minor_allele]).sort_values('Pos')



#
#dfs = []
#for i in range(1000):
#    dfs.append(sample_freqs('X:/volume2/noam/loop_genomics/sample_1_pipeline/contig.freqs', 50))
#df = pd.concat(dfs)
#
#df['sampling_diff'] = df.Freq - df.Sampled_freq
#
#fig, ax = plt.subplots(nrows=1, ncols=1)
#sns.regplot(data=df, x='Freq', y='Sampled_freq', ax=ax, color='purple')
#slope, intercept, r_value, p_value, std_err = stats.linregress(df['Freq'], df['Sampled_freq'])
#ax.set_title("y=%fx+%f, r2=%f" % (slope,intercept, r_value**2))
#
#fig, ax = plt.subplots(nrows=1, ncols=1)
#sns.regplot(data=df[(df.Freq < 0.5)], x='Freq', y='Sampled_freq', ax=ax, color='purple')
#slope, intercept, r_value, p_value, std_err = stats.linregress(df[(df.Freq < 0.5)]['Freq'], df[(df.Freq < 0.5)]['Sampled_freq'])
#ax.set_title("y=%fx+%f, r2=%f" % (slope,intercept, r_value**2))
#
#fig, ax = plt.subplots(nrows=1, ncols=1)
#sns.regplot(data=df[(df.Freq > 0.1) & (df.Freq < 0.5)], x='Freq', y='Sampled_freq', ax=ax, color='purple')
#slope, intercept, r_value, p_value, std_err = stats.linregress(df[(df.Freq > 0.1) & (df.Freq < 0.5)]['Freq'], df[(df.Freq > 0.1) & (df.Freq < 0.5)]['Sampled_freq'])
#ax.set_title("y=%fx+%f, r2=%f" % (slope,intercept, r_value**2))