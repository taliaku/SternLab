"""
This script simulates a bottleneck effect on a population represented by a 
freqs file. Before sampling, it normalized the population using an important 
assumption: we only have two allelles per position! also removes insertions.
"""
import pandas as pd
import numpy as np

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
    df = normalize_two_alleles(freqs_path)
    major_allele = df[df.Rank == 0].copy()
    major_allele['Sampled_freq'] = major_allele['Freq'].apply(lambda x: np.random.binomial(sample_size, x) / sample_size)
    minor_allele = df[df.Rank == 1].copy()
    minor_allele = pd.merge(minor_allele, major_allele[['Pos', 'Sampled_freq']], on='Pos')
    minor_allele['Sampled_freq'] = 1 - minor_allele['Sampled_freq']
    return pd.concat([major_allele, minor_allele]).sort_values('Pos')

def simulate_and_get_statistics():