
from scipy.stats import poisson
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import itertools
import argparse



CONSTANT_BACTERIA_COUNT = 10**10
CONSTANT_WT_COUNT = 10**10

PAYOFF_MATRIX = {'wt':{'wt':1, 'del':0, 'syn':1},
                 'del':{'wt':3, 'del':0, 'syn':0.1}, 
                 'syn':{'wt':3, 'del':0, 'syn':0.1}}
DEL_INITIAL_COUNT = 10**5
SYN_INITIAL_COUNT = 10**5

PASSAGES = 25
ITERATIONS_WITHIN_PASSAGES = 1

MAX_INFECTIONS_PER_CELL = 1000

def simulate(payoff_matrix=PAYOFF_MATRIX, del_initial_count=DEL_INITIAL_COUNT, syn_initial_count=SYN_INITIAL_COUNT, passages=PASSAGES, iterations_within_passages=ITERATIONS_WITHIN_PASSAGES):
    passages_list = []
    wt_frequencies = []
    del_frequencies = []
    syn_frequencies = []
    wt_counts = []
    del_counts = []
    syn_counts = []
    
    infection_frequencies_list = []
    viral_counts = {'wt':CONSTANT_WT_COUNT, 'del':del_initial_count, 'syn':syn_initial_count}
    for passage in range(1, (int(passages * iterations_within_passages) + 1)):
        # calculate frequencies for each infection types (poison)
        infection_frequencies = get_infection_type_frequencies(CONSTANT_BACTERIA_COUNT, viral_counts)
        #if sum(infection_frequencies.values()) != 1:
        #    print('problem with infection frequencies calculation')
        #    print(sum(infection_frequencies.values()))
        # calculate frequencies after passage
        virus_frequencies = get_viral_frequencies(infection_frequencies, payoff_matrix)
        # calculate virus numbers after passage
        if virus_frequencies['wt'] == 0:
            break
        for strain in virus_frequencies:    
            viral_counts[strain] = virus_frequencies[strain] * CONSTANT_WT_COUNT / virus_frequencies['wt']
        passages_list.append(passage)
        wt_frequencies.append(virus_frequencies['wt'])
        del_frequencies.append(virus_frequencies['del'])
        syn_frequencies.append(virus_frequencies['syn'])
        wt_counts.append(viral_counts['wt'])
        del_counts.append(viral_counts['del'])
        syn_counts.append(viral_counts['syn'])
        infection_frequencies_list.append(infection_frequencies)
    data = pd.DataFrame({'passage':passages_list, 'wt_frequency':wt_frequencies, 'del_frequency':del_frequencies, 'syn_frequency':syn_frequencies, 'wt_count':wt_counts, 'del_count':del_counts, 'syn_count':syn_counts, 'infection_frequencies':infection_frequencies_list})
    data.passage = data.passage.astype(float) / iterations_within_passages
    return data     
                

def get_infection_type_frequencies(bacteria_count, viral_counts_dict):
    infection_type_dict = {}
    # no infections
    f = 1
    for i in viral_counts_dict:
        f *= (poisson.cdf(0 , viral_counts_dict[i]/bacteria_count))
    infection_type_dict['none'] = f
    # alone
    for i in viral_counts_dict: # only i in cell
        f = (poisson.cdf(MAX_INFECTIONS_PER_CELL, viral_counts_dict[i]/bacteria_count) - poisson.cdf(0, viral_counts_dict[i]/bacteria_count))
        for j in viral_counts_dict:
            if j != i:
                f *= poisson.cdf(0 , viral_counts_dict[j]/bacteria_count)
        infection_type_dict[i] = f
    # two virus types
    for i in viral_counts_dict: # only i not in cell
        f = poisson.cdf(0 , viral_counts_dict[i]/bacteria_count)
        for j in viral_counts_dict:
            if j != i:
                f *= (poisson.cdf(MAX_INFECTIONS_PER_CELL , viral_counts_dict[j]/bacteria_count) - poisson.cdf(0 , viral_counts_dict[j]/bacteria_count))
        infection_type_dict['_'.join([k for k in viral_counts_dict if k != i])] = f
    # all virus types
    f = 1
    for i in viral_counts_dict:
        f *= (poisson.cdf(MAX_INFECTIONS_PER_CELL , viral_counts_dict[i]/bacteria_count) - poisson.cdf(0 , viral_counts_dict[i]/bacteria_count))
    infection_type_dict['all_virus_types'] = f
    return infection_type_dict

def get_viral_frequencies(infection_type_dict, payoff_matrix):
    infection_type_dict2 = infection_type_dict.copy()
    prenormalized_frequencies = {}
    pairs = [key for key in infection_type_dict2.keys() if key.count('_') == 1]
    for pair in pairs:
        infection_type_dict2['_'.join(pair.split('_')[::-1])] = infection_type_dict2[pair]
    for i in payoff_matrix:
        f = (payoff_matrix[i][i] * infection_type_dict2[i]) # alone
        for j in payoff_matrix: # two types per cell
            if j != i:
                f+= (payoff_matrix[i][j] * infection_type_dict2[i + '_' + j])
        # all types
        for dominant in payoff_matrix:
            for less_dom in payoff_matrix:
                if dominant != i and less_dom != i and dominant != less_dom:
                    if infection_type_dict2[i + '_' + dominant] + infection_type_dict2[i + '_' + less_dom] > 0:
                        f += ((infection_type_dict2[i + '_' + dominant] / (infection_type_dict2[i + '_' + dominant] + infection_type_dict2[i + '_' + less_dom])) * payoff_matrix[i][dominant] * infection_type_dict2['all_virus_types'])
        prenormalized_frequencies[i] = f
        # normalize
        normalized_frequencies = {}
        for i in prenormalized_frequencies:
            normalized_frequencies[i] = prenormalized_frequencies[i] / sum(prenormalized_frequencies.values())
    return normalized_frequencies

def plot_cheaters(data, out_path=False, params_text=False):
    fig, ax = plt.subplots(nrows=1, ncols=1)    
    data.plot(x='passage', y='wt_frequency', color='grey', ax=ax)
    data.plot(x='passage', y='del_frequency', color='red', ax=ax)
    data.plot(x='passage', y='syn_frequency', color='orange', ax=ax)
    if params_text:
        ax.text(0.8, 0.8, params_text)
    if out_path:
        fig.savefig(out_path)
    fig.show()        
    #plt.close('all')