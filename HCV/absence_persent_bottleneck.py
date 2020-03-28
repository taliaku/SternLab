# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 09:48:30 2020

@author: Noam
"""
import pandas as pd




def calculate_likelihood(pairs, nb, cutoff):
    # pairs: a list with tuples of the order (donor_freq, recepient_freq)
    present_pairs = [p for p in pairs if p[1] >= cutoff]
    absent_pairs = [p for p in pairs if p[1] < cutoff]
    likelihood = 1
    for pair in absent_pairs:
        likelihood *= ((1 - pair[0])**nb)
    for pair in present_pairs:
        likelihood *= (1 - ((1 - pair[0])**nb))
    return likelihood

def find_maximum_likelihood(pairs, cutoff=0.01):
    likelihoods = []
    for nb in range(200):
        likelihood = calculate_likelihood(pairs, nb, cutoff)
        likelihoods.append((nb, likelihood))
    likelihoods = pd.DataFrame(likelihoods, columns=['nb', 'likelihood'])
    #likelihoods.plot(x='nb', y='likelihood', kind='scatter')
    print('nb with max likelihood: ' + str(likelihoods.sort_values('likelihood', ascending=False).head(1).nb.max()))
    return likelihoods

infection_order = {'HCV-P8':1, 'HCV-P1':2, 'HCV-P11':3, 'HCV-P7':5, 'HCV-P4':6, 'HCV-P3':7, 'HCV-P5':8, 'HCV-P2':9, 'HCV-P6':10, 'HCV-P9':11, 'HCV-P10':12, 'HCV-PS2':0}
infection_order = pd.DataFrame.from_dict(infection_order, orient='index').reset_index().rename(columns={'index':'Sample', 0:'infection_order'})

all_results = []
for sample in infection_order.Sample.tolist():
    pairs = [tuple(x) for x in pd.read_csv('X:/volume2/noam/hcv/BB_bottleneck_output/7/' + sample + '.csv', '\t', header=None)[[0,1]].to_numpy()]
    a = find_maximum_likelihood(pairs, cutoff=0.04)
    nb_max = a.sort_values('likelihood', ascending=False).head(1).nb.max()
    all_results.append((sample, nb_max))
all_results = pd.DataFrame(all_results, columns=['Sample', 'bottleneck'])
all_results = pd.merge(all_results, infection_order, on=['Sample'], how='left')
all_results = all_results[~(all_results.infection_order.isna()) & (all_results.Sample != 'HCV-PS2')].sort_values('infection_order')
fig, ax = plt.subplots(nrows=1, ncols=1)
slope, intercept, r_value, p_value, std_err = stats.linregress(all_results.infection_order, all_results.bottleneck)
sns.regplot(all_results.infection_order, all_results.bottleneck, ax = ax)
ax.set_title("y=%fx+%f, r=%f, p=%f" % (slope,intercept, r_value, p_value))
