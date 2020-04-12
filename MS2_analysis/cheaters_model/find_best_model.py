
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from tqdm import tqdm
import numpy as np
import itertools


import sys, os
sys.path.append('/sternadi/home/volume2/noam/SternLab/MS2_analysis/cheaters_model/')
sys.path.append('X:/volume2/noam/SternLab/MS2_analysis/cheaters_model')
from model_triple_payoffs import simulate, plot_cheaters

#
########################################
######### find best parameters #########
########################################
#
## organize 37B real data for fit
#df = pd.read_csv('X:/volume2/noam/passages/201909/all_freqs.csv')
##df = pd.read_csv('/sternadi/home/volume2/noam/passages/201909/all_freqs.csv')
#df = df[~(df.Time.isin([11,12,14]))]
#real_data = df[(df.Full_mutation.isin(['T1764.0-', 'A1664.0G'])) & (df.Replica == 'B')].pivot_table(values='Freq', index=['Replica', 'Time'], columns='Full_mutation').reset_index()
#real_data.rename(columns={'Time':'passage', 'A1664.0G':'syn_frequency', 'T1764.0-':'del_frequency'}, inplace=True)
#real_data.columns.name = None
#real_data = real_data[['passage', 'syn_frequency', 'del_frequency']]
#real_data['wt_frequency'] = 1 - real_data.syn_frequency - real_data.del_frequency
#real_data['passage'] = real_data.passage.astype(float)
#real_data.to_csv('X:/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv', index=False)
#
## organize 37A real data for fit
#
#df = pd.read_csv('X:/volume2/noam/passages/201909/all_freqs.csv')
##df = pd.read_csv('/sternadi/home/volume2/noam/passages/201909/all_freqs.csv')
#df = df[~(df.Time.isin([11,12,14]))]
#real_data = df[(df.Full_mutation.isin(['T1764.0-', 'A1664.0G'])) & (df.Replica == 'A')].pivot_table(values='Freq', index=['Replica', 'Time'], columns='Full_mutation').reset_index()
#real_data.rename(columns={'Time':'passage', 'A1664.0G':'syn_frequency', 'T1764.0-':'del_frequency'}, inplace=True)
#real_data.columns.name = None
#real_data = real_data[['passage', 'syn_frequency', 'del_frequency']]
#real_data['wt_frequency'] = 1 - real_data.syn_frequency - real_data.del_frequency
#real_data['passage'] = real_data.passage.astype(float)
#real_data.to_csv('X:/volume2/noam/cheater_model/37A_realdata_until_p23_for_fit.csv', index=False)
##
#
#
#
#
#
## organize 37B real data for fit - no p18,20,22
#df = pd.read_csv('X:/volume2/noam/passages/201909/all_freqs.csv')
##df = pd.read_csv('/sternadi/home/volume2/noam/passages/201909/all_freqs.csv')
#df = df[~(df.Time.isin([11,12,14,18,20,22]))]
#real_data = df[(df.Full_mutation.isin(['T1764.0-', 'A1664.0G'])) & (df.Replica == 'B')].pivot_table(values='Freq', index=['Replica', 'Time'], columns='Full_mutation').reset_index()
#real_data.rename(columns={'Time':'passage', 'A1664.0G':'syn_frequency', 'T1764.0-':'del_frequency'}, inplace=True)
#real_data.columns.name = None
#real_data = real_data[['passage', 'syn_frequency', 'del_frequency']]
#real_data['wt_frequency'] = 1 - real_data.syn_frequency - real_data.del_frequency
#real_data['passage'] = real_data.passage.astype(float)
#real_data.to_csv('X:/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit_no_18_20_22.csv', index=False)
#
## organize 37A real data for fit - no p18,20,22
#df = pd.read_csv('X:/volume2/noam/passages/201909/all_freqs.csv')
##df = pd.read_csv('/sternadi/home/volume2/noam/passages/201909/all_freqs.csv')
#df = df[~(df.Time.isin([11,12,14,18,20,22]))]
#real_data = df[(df.Full_mutation.isin(['T1764.0-', 'A1664.0G'])) & (df.Replica == 'B')].pivot_table(values='Freq', index=['Replica', 'Time'], columns='Full_mutation').reset_index()
#real_data.rename(columns={'Time':'passage', 'A1664.0G':'syn_frequency', 'T1764.0-':'del_frequency'}, inplace=True)
#real_data.columns.name = None
#real_data = real_data[['passage', 'syn_frequency', 'del_frequency']]
#real_data['wt_frequency'] = 1 - real_data.syn_frequency - real_data.del_frequency
#real_data['passage'] = real_data.passage.astype(float)
#real_data.to_csv('X:/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit_no_18_20_22.csv', index=False)
#
#
##real_data = pd.read_excel('/sternadi/home/volume2/noam/passages/201908_w_2019_passages/new_2019/WT_1764_1664.xlsx')
###real_data = pd.read_excel('X:/volume2/noam/passages/201908_w_2019_passages/new_2019/WT_1764_1664.xlsx')
##real_data = real_data[real_data.Replica == 'B']
##real_data = real_data[~real_data.Time.isin([11, 12, 14])]
##real_data = real_data[['Time', 'A1664.0G frequency', 'T1764.0- frequency', 'WT frequency']]
##real_data.rename(columns={'Time':'passage', 'A1664.0G frequency':'syn_frequency', 'T1764.0- frequency':'del_frequency', 'WT frequency':'wt_frequency'}, inplace=True)
##real_data['passage'] = real_data.passage.astype(float)
##
##
##def find_best_parameters_euclidean(real_data):
##    payoff_options = [0, 0.1, 0.25, 0.5, 1, 2, 2.5, 3]
##    #payoff_options = [0, 0.1,0.5,1,2,]
##    payoff_matrix = itertools.product(payoff_options, repeat=7)
##    del_initial_count = [10**5]
##    syn_initial_count = [10**5]
##    passages = [25]
##    iterations_per_passage = [1,2]
##    results = []
##    for p in tqdm(payoff_matrix):
##        #payoff = {'wt':{'wt':1, 'del':p[0], 'syn':p[1]},
##        #                 'del':{'wt':p[2], 'del':0, 'syn':p[3]}, 
##        #                 'syn':{'wt':p[4], 'del':p[5], 'syn':p[6]}}
##        for d in del_initial_count:
##            for s in syn_initial_count:
##                for pa in passages:
##                    for iteration in iterations_per_passage:
##                        #results.append((payoff, p[0], p[1], p[2], p[3], p[4], p[5], p[6], d, s, pa, iteration))
##                        results.append((p[0], p[1], p[2], p[3], p[4], p[5], p[6], d, s, pa, iteration))
##    #results = pd.DataFrame(results, columns=['payoff_matrix', 'wt_with_del_p', 'wt_with_syn_p', 'del_with_wt_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'syn_with_syn_p', 'del_initial_count', 'syn_initial_count', 'passages', 'iterations_within_passage'])
##    results = pd.DataFrame(results, columns=['wt_with_del_p', 'wt_with_syn_p', 'del_with_wt_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'syn_with_syn_p', 'del_initial_count', 'syn_initial_count', 'passages', 'iterations_within_passage'])
##    results['iterator'] = results.index % 100
##    return results
##
##results = find_best_parameters_euclidean(real_data)
##results.to_csv('X:/volume2/noam/cheater_model/p23/iteration_job_ids_p23.csv', index=False)
##
##
##
##
##
##results = pd.read_csv('X:/volume2/noam/cheater_model/iteration_job_ids_2.results.csv')
###for row in results.sort_values('diff').head(20).itertuples():
###for i, row in enumerate(results[(results.wt_final_passage == 20) & (results.del_with_wt_p > 1) & (results.syn_with_wt_p > 1) & (results.wt_with_del_p < 1) & (results.syn_with_syn_p > 0)].sort_values('diff').head(20).itertuples()):
##for row in results[(results.wt_final_passage == 20) & (results.iterations_within_passage == 1) & (results.syn_with_syn_p < 1) & (results.syn_with_syn_p > 0)].sort_values('diff').head(10).itertuples():
###for row in results[(results.syn_with_syn_p < 1) & (results.syn_with_syn_p > 0) & (results.wt_final_passage == 20)& (results['diff'] > 0.3) & (results['diff'] < 0.35)].sort_values('diff').head(5).itertuples():
##    payoff_matrix={'wt':{'wt':1, 'del':row.wt_with_del_p, 'syn':row.wt_with_syn_p},
##                         'del':{'wt':row.del_with_wt_p, 'del':0, 'syn':row.del_with_syn_p}, 
##                         'syn':{'wt':row.syn_with_wt_p, 'del':row.syn_with_del_p, 'syn':row.syn_with_syn_p}}
##    simulation = simulate(payoff_matrix=payoff_matrix, del_initial_count=row.del_initial_count, syn_initial_count=row.syn_initial_count, passages=row.passages, iterations_within_passages=row.iterations_within_passage)
##    #if simulation[~simulation.isna().any(axis=1)].passage.max() >= 15:
##    #plot_cheaters(simulation, out_path='X:/volume2/noam/cheater_model/iteration_job_ids_2_figs/' + str(i), params_text='\n'.join([str(i) for i in payoff_matrix.items()]) + '\ninitial_del_count: ' + str(row.del_initial_count) + '\ninitial_syn_count: ' + str(row.syn_initial_count) + '\niterations_per_passage: ' + str(row.iterations_within_passage))
##    #plot_cheaters(simulation, out_path='X:/volume2/noam/cheater_model/iteration_job_ids_2_figs/' + str(i), params_text='\n'.join([str(i) for i in payoff_matrix.items()]) + '\niterations_per_passage: ' + str(row.iterations_within_passage) + '\ndiff: ' + str(row.diff))
##    plot_cheaters(simulation, out_path=False, params_text='\n'.join([str(i) for i in payoff_matrix.items()]) + '\niterations_per_passage: ' + str(row.iterations_within_passage) + '\ndiff: ' + str(row.diff))
##
##
##
#
#
#
############################## plot best models
#results = pd.read_csv('X:/volume2/noam/cheater_model/iteration_job_ids_2.results.csv')
#
#
#df = pd.read_csv('X:/volume2/noam/passages/201909/all_freqs.csv')
#df = df[~(df.Time.isin([11,12,14]))]
#real_data = df[df.Full_mutation.isin(['T1764.0-', 'A1664.0G'])].pivot_table(values='Freq', index=['Replica', 'Time'], columns='Full_mutation').reset_index()
#real_data.rename(columns={'Time':'passage', 'A1664.0G':'syn_frequency', 'T1764.0-':'del_frequency'}, inplace=True)
#real_data['passage'] = real_data.passage.astype(float)
#
#
#def plot_model_against_real2(model_data, real_data, out_path=False, replica=False):
#    fig, ax = plt.subplots(nrows=1, ncols=1)
#    model_data = model_data[model_data.passage <= 23]
#    model_data.plot(x='passage', y='del_frequency', color='red', ax=ax, label='$\Delta$1764 simulation' , linestyle='--')
#    model_data.plot(x='passage', y='syn_frequency', color='orange', ax=ax, label='A1664G simulation' , linestyle='--')
#    if replica:
#        real_data[real_data.Replica == replica].plot(x='passage', y='del_frequency', color='red', ax=ax, label='$\Delta$1764 line ' + replica)
#        real_data[real_data.Replica == replica].plot(x='passage', y='syn_frequency', color='orange', ax=ax, label='A1664G line ' + replica)    
#    ax.set_xticks([0,5,10,15,20,])
#    fig.set_size_inches(6.5,3)
#    ax.set_xlabel('Passage', fontsize=14)
#    ax.set_ylabel('Mutation Frequency', fontsize=14)
#    if out_path:
#        fig.savefig(out_path, dpi=800, bbox_inches='tight')
#    fig.show()   
#
#
## 2 cheater model
#for row in results[(results.wt_final_passage == 20) & (results.iterations_within_passage == 1) & (results.syn_with_syn_p < 1) & (results.syn_with_syn_p > 0)].sort_values('diff').head(1).itertuples():
##for row in results[(results.syn_with_syn_p < 1) & (results.syn_with_syn_p > 0) & (results.wt_final_passage == 20)& (results['diff'] > 0.3) & (results['diff'] < 0.35)].sort_values('diff').head(5).itertuples():
#    payoff_matrix={'wt':{'wt':1, 'del':row.wt_with_del_p, 'syn':row.wt_with_syn_p},
#                         'del':{'wt':row.del_with_wt_p, 'del':0, 'syn':row.del_with_syn_p}, 
#                         'syn':{'wt':row.syn_with_wt_p, 'del':row.syn_with_del_p, 'syn':row.syn_with_syn_p}}
#    #simulation = simulate(payoff_matrix=payoff_matrix, del_initial_count=row.del_initial_count, syn_initial_count=row.syn_initial_count, passages=row.passages, iterations_within_passages=row.iterations_within_passage)
#    simulation = simulate(payoff_matrix=payoff_matrix, del_initial_count=row.del_initial_count, syn_initial_count=row.syn_initial_count, passages=25, iterations_within_passages=row.iterations_within_passage)
#    #if simulation[~simulation.isna().any(axis=1)].passage.max() >= 15:
#    #plot_cheaters(simulation, out_path='X:/volume2/noam/cheater_model/iteration_job_ids_2_figs/' + str(i), params_text='\n'.join([str(i) for i in payoff_matrix.items()]) + '\ninitial_del_count: ' + str(row.del_initial_count) + '\ninitial_syn_count: ' + str(row.syn_initial_count) + '\niterations_per_passage: ' + str(row.iterations_within_passage))
#    #plot_cheaters(simulation, out_path='X:/volume2/noam/cheater_model/iteration_job_ids_2_figs/' + str(i), params_text='\n'.join([str(i) for i in payoff_matrix.items()]) + '\niterations_per_passage: ' + str(row.iterations_within_passage) + '\ndiff: ' + str(row.diff))
#    plot_cheaters(simulation, out_path=False, params_text='\n'.join([str(i) for i in payoff_matrix.items()]) + '\niterations_per_passage: ' + str(row.iterations_within_passage) + '\ndiff: ' + str(row.diff))
#    plot_model_against_real2(simulation, real_data, 'X:/volume2/noam/cheater_model/best_2_cheater_lineA.png', replica='A')
#    plot_model_against_real2(simulation, real_data, 'X:/volume2/noam/cheater_model/best_2_cheater_lineB.png', replica='B')
#
#
## 1 cheater model 
#
#def plot_model_against_real1(model_data, real_data, out_path=False):
#    fig, ax = plt.subplots(nrows=1, ncols=1)
#    model_data = model_data[model_data.passage <= 23]
#    model_data.plot(x='passage', y='del_frequency', color='red', ax=ax, label='$\Delta$1764 simulation', linestyle='--')
#    real_data[real_data.Replica =='B'].plot(x='passage', y='del_frequency', color='red', ax=ax , label='$\Delta$1764 line B')
#    real_data[real_data.Replica =='A'].plot(x='passage', y='del_frequency', color='#6D0C09', ax=ax, label='$\Delta$1764 line A')
#    fig.set_size_inches(6.5,3)
#    ax.set_xticks([0,5,10,15,20,])
#    ax.set_xlabel('Passage', fontsize=14)
#    ax.set_ylabel('Mutation Frequency', fontsize=14)
#    if out_path:
#        fig.savefig(out_path, dpi=800, bbox_inches='tight')
#    fig.show()  
#    
#p = {'wt':{'wt':1, 'del':1, 'syn':1},
#                 'del':{'wt':2, 'del':0, 'syn':0}, 
#                 'syn':{'wt':0, 'del':0, 'syn':0}}
#
#simulation = simulate(payoff_matrix=p, del_initial_count=10**5, syn_initial_count=0, passages=23, iterations_within_passages=2)
#
#plot_model_against_real1(simulation, real_data, 'X:/volume2/noam/cheater_model/best_1_cheater.png')
#
#
#######
## ineefective frequencies
#simulation.infection_frequencies.apply(lambda x: x['del'] + x['syn'] + x['del_syn'])
#
#
#
#
#
######################## triple infections #########
######################## hill climbing #############
#
#
##### initial filtering
#a = pd.read_csv('X:/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/all_results.csv')
## mean killing ratio = mean killing at del peak / mean killing at p1
#
#
#for row in a[(a.wt_final_passage == 24) & (a.mean_killing_ratio < 0.9) & (a.triple_wt <= 1)].sort_values('diff').head(20).itertuples():
#    payoff_matrix={'wt':{'wt':1, 'del':row.wt_with_del_p, 'syn':row.wt_with_syn_p},
#                         'del':{'wt':row.del_with_wt_p, 'del':0, 'syn':row.del_with_syn_p}, 
#                         'syn':{'wt':row.syn_with_wt_p, 'del':row.syn_with_del_p, 'syn':row.syn_with_syn_p}}
#    triple_payoff_matrix = {'wt':row.triple_wt, 'del':row.triple_del, 'syn':row.triple_syn}
#    simulation = simulate(payoff_matrix=payoff_matrix, triple_payoff_matrix=triple_payoff_matrix, iterations_within_passages=row.iterations_within_passage)
#    plot_cheaters(simulation, out_path=False, params_text='\n'.join([str(i) for i in payoff_matrix.items()]) + '\ntriple_payoff:\n' + str(triple_payoff_matrix) + '\niterations_per_passage: ' + str(row.iterations_within_passage) + '\ndiff: ' + str(row.diff))
#    
#    
#### more optimization
#a = pd.read_csv('X:/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/all_results_optimized_mean_killing_0.9.csv')
#a['iterations_within_passage'] = 1
#a = a.drop_duplicates()
#
#for row in a[(a.wt_final_passage == 24) & (a.mean_killing_ratio < 0.9)].sort_values('diff').head(10).itertuples():
##for row in a[(a.wt_final_passage == 24) & (a.mean_killing_ratio < 0.9) & (a.triple_wt <= 1)].sort_values('diff').head(20).itertuples():
#    print(row.syn_with_del_p)
#    payoff_matrix={'wt':{'wt':1, 'del':row.wt_with_del_p, 'syn':row.wt_with_syn_p},
#                         'del':{'wt':row.del_with_wt_p, 'del':0, 'syn':row.del_with_syn_p}, 
#                         'syn':{'wt':row.syn_with_wt_p, 'del':row.syn_with_del_p, 'syn':row.syn_with_syn_p}}
#    triple_payoff_matrix = {'wt':row.triple_wt, 'del':row.triple_del, 'syn':row.triple_syn}
#    simulation = simulate(payoff_matrix=payoff_matrix, triple_payoff_matrix=triple_payoff_matrix, iterations_within_passages=row.iterations_within_passage)
#    plot_cheaters(simulation, out_path=False, params_text='\n'.join([str(i) for i in payoff_matrix.items()]) + '\ntriple_payoff:\n' + str(triple_payoff_matrix) + '\niterations_per_passage: ' + str(row.iterations_within_passage) + '\ndiff: ' + str(row.diff) + '\nkilling ratio: ' + str(row.mean_killing_ratio))
#    
#
#
## take top 1% for further analysis
#a = pd.read_csv('X:/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/all_results_optimized_mean_killing_0.9.csv')
#a['iterations_within_passage'] = 1
#a = a.drop_duplicates()
#b = a[(a.wt_final_passage == 24) & (a.mean_killing_ratio < 0.9) & (a.syn_with_syn_p > 0) & (a.syn_with_syn_p < 1)]
#b = b[b['diff'] < b['diff'].quantile(0.01)]
##profile = b.profile_report(title='Pandas Profiling Report')
##profile.to_file(output_file='X:/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/all_results_optimized.html')
#
#c = b.copy()
#c['wt_with_wt_p'] = 1
#c['del_with_del_p'] = 0
#fig, axes = plt.subplots(nrows=4, ncols=3)
#fig.set_size_inches(8,8)
#fig.tight_layout()
#axes = axes.flatten()
#fig.subplots_adjust(hspace=0.5)
#fig.text(0.5, -0.02, 'Payoff', ha='center', fontsize = 14)
#fig.text(-0.02, 0.5, 'Frequency', va='center', rotation='vertical', fontsize = 14)
#for i, column in enumerate(['wt_with_wt_p', 'wt_with_del_p', 'wt_with_syn_p', 'del_with_wt_p', 'del_with_del_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'syn_with_syn_p', 'triple_wt', 'triple_del', 'triple_syn']):
#    c[column].plot.hist(ax=axes[i])
#    axes[i].set_xlabel(None)
#    axes[i].set_ylabel(None)
#    axes[i].set_xlim(0,3.8)
#    if 'with' in column:
#        axes[i].set_title(' when infecting with '.join(column.replace('_p', '').split('_with_')).replace('wt', 'WT').replace('del', '$\Delta$1764').replace('syn', 'A1664G'), fontsize=10)
#    if 'triple' in column:
#        axes[i].set_title(column.replace('triple_', '').replace('wt', 'WT').replace('del', '$\Delta$1764').replace('syn', 'A1664G') + ' during triple infection', fontsize=10)
#fig.savefig('X:/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/all_results_optimized_mean_killing_0.9_top_0.01/payoffs.png', dpi=800, bbox_inches='tight')
#
#
#fig, ax = plt.subplots(nrows=1, ncols=1)
#fig.set_size_inches(3,2)
#(b.del_with_wt_p / b.syn_with_wt_p).plot.hist(ax=ax)
#ax.set_xlabel('Ratio')
#ax.set_title('$\Delta$1764 when infecting with WT / \nA1664G when infecting with WT')
#fig.savefig('X:/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/all_results_optimized_mean_killing_0.9_top_0.01/ratio_del_wt_to_syn_wt.png', dpi=800, bbox_inches='tight')
#
#
#fig, ax = plt.subplots(nrows=1, ncols=1)
#fig.set_size_inches(3,2)
##((b.syn_with_del_p + b.triple_syn) - (b.del_with_syn_p + b.triple_del)).plot.hist()
#((b.del_with_syn_p + b.triple_del) / (b.syn_with_del_p + b.triple_syn)).plot.hist(ax=ax)
#ax.set_xlabel('Ratio')
#ax.set_title('$\Delta$1764 when infecting with A1664G + $\Delta$1764 during triple infection \n/ \nA1664G when infecting with $\Delta$1764 + A1664G during triple infection')
#fig.savefig('X:/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/all_results_optimized_mean_killing_0.9_top_0.01/ratio_del_syn.png', dpi=800, bbox_inches='tight')
#
#
#(b.del_with_wt_p / b.syn_with_wt_p).plot.hist()
#sns.regplot(x=b['triple_syn'], y=b['syn_with_del_p'])
#(b['triple_syn'] / b['syn_with_del_p']).plot.hist()
#sns.regplot(x=b['triple_wt'], y=b['wt_with_syn_p'])
#sns.regplot(x=b['del_with_wt_p'], y=b['syn_with_wt_p'])
#sns.regplot(x=b['triple_wt'], y=b['wt_with_syn_p'])
#sns.regplot(x=b['syn_with_syn_p'], y=b['wt_with_syn_p'])
#sns.regplot(x=b['triple_syn'], y=b['syn_with_del_p'])
#sns.regplot(x=b['triple_del'], y=b['del_with_syn_p'])
#(b.wt_with_del_p / b.wt_with_syn_p).plot.hist()
#
########################### cloud graph for best simulations
#simulations = []
#for row in tqdm(b.itertuples()):
#    payoff_matrix={'wt':{'wt':1, 'del':row.wt_with_del_p, 'syn':row.wt_with_syn_p},
#                         'del':{'wt':row.del_with_wt_p, 'del':0, 'syn':row.del_with_syn_p}, 
#                         'syn':{'wt':row.syn_with_wt_p, 'del':row.syn_with_del_p, 'syn':row.syn_with_syn_p}}
#    triple_payoff_matrix = {'wt':row.triple_wt, 'del':row.triple_del, 'syn':row.triple_syn}
#    simulation = simulate(payoff_matrix=payoff_matrix, triple_payoff_matrix=triple_payoff_matrix, iterations_within_passages=row.iterations_within_passage)
#    simulations.append(simulation)
#s = pd.concat(simulations)
#
#real_data = pd.read_csv('X:/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv')
#s = s[s.passage.isin(real_data.passage.tolist())]
#
#fig, ax = plt.subplots(nrows=1, ncols=1)
#ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().syn_frequency, s.groupby('passage').max().reset_index().syn_frequency, color='orange', alpha=0.5)   
#ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().del_frequency, s.groupby('passage').max().reset_index().del_frequency, color='red', alpha=0.5) 
##ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().wt_frequency, s.groupby('passage').max().reset_index().wt_frequency, color='grey', alpha=0.5)     
##real_data.plot(x='passage', y='wt_frequency', color='grey', ax=ax)
#real_data.plot(x='passage', y='del_frequency', color='red', ax=ax)
#real_data.plot(x='passage', y='syn_frequency', color='orange', ax=ax)
#
#
#simulations = []
#for row in tqdm(b[b['diff'] < b['diff'].quantile(0.01)].itertuples()):
#    payoff_matrix={'wt':{'wt':1, 'del':row.wt_with_del_p, 'syn':row.wt_with_syn_p},
#                         'del':{'wt':row.del_with_wt_p, 'del':0, 'syn':row.del_with_syn_p}, 
#                         'syn':{'wt':row.syn_with_wt_p, 'del':row.syn_with_del_p, 'syn':row.syn_with_syn_p}}
#    triple_payoff_matrix = {'wt':row.triple_wt, 'del':row.triple_del, 'syn':row.triple_syn}
#    simulation = simulate(payoff_matrix=payoff_matrix, triple_payoff_matrix=triple_payoff_matrix, iterations_within_passages=row.iterations_within_passage)
#    simulations.append(simulation)
#s = pd.concat(simulations)
#
#real_data = pd.read_csv('X:/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv')
#s = s[s.passage.isin(real_data.passage.tolist())]
#
#fig, ax = plt.subplots(nrows=1, ncols=1)
#ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().syn_frequency, s.groupby('passage').max().reset_index().syn_frequency, color='orange', alpha=0.5)   
#ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().del_frequency, s.groupby('passage').max().reset_index().del_frequency, color='red', alpha=0.5) 
##ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().wt_frequency, s.groupby('passage').max().reset_index().wt_frequency, color='grey', alpha=0.5)     
##real_data.plot(x='passage', y='wt_frequency', color='grey', ax=ax)
#real_data.plot(x='passage', y='del_frequency', color='red', ax=ax)
#real_data.plot(x='passage', y='syn_frequency', color='orange', ax=ax)
#
#
#simulations = []
#for row in tqdm(b[b['diff'] < b['diff'].quantile(0.1)].itertuples()):
#    payoff_matrix={'wt':{'wt':1, 'del':row.wt_with_del_p, 'syn':row.wt_with_syn_p},
#                         'del':{'wt':row.del_with_wt_p, 'del':0, 'syn':row.del_with_syn_p}, 
#                         'syn':{'wt':row.syn_with_wt_p, 'del':row.syn_with_del_p, 'syn':row.syn_with_syn_p}}
#    triple_payoff_matrix = {'wt':row.triple_wt, 'del':row.triple_del, 'syn':row.triple_syn}
#    simulation = simulate(payoff_matrix=payoff_matrix, triple_payoff_matrix=triple_payoff_matrix, iterations_within_passages=row.iterations_within_passage)
#    simulations.append(simulation)
#s = pd.concat(simulations)
#
#real_data = pd.read_csv('X:/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv')
#s = s[s.passage.isin(real_data.passage.tolist())]
#
#plt.style.use('default')
#fig, ax = plt.subplots(nrows=1, ncols=1)
#ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().syn_frequency, s.groupby('passage').max().reset_index().syn_frequency, color='orange', alpha=0.5)   
#ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().del_frequency, s.groupby('passage').max().reset_index().del_frequency, color='red', alpha=0.5) 
##ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().wt_frequency, s.groupby('passage').max().reset_index().wt_frequency, color='grey', alpha=0.5)     
##real_data.plot(x='passage', y='wt_frequency', color='grey', ax=ax)
#real_data.plot(x='passage', y='del_frequency', color='red', ax=ax,  label='$\Delta$1764 simulation')
#real_data.plot(x='passage', y='syn_frequency', color='orange', ax=ax, label='A1664G simulation')
#ax.set_xticks([0,5,10,15,20])
#fig.set_size_inches(6,4)
#


####################################################
################ ABC SMC ###########################
####################################################

### open best simulations
def get_best_parameters(astroabc_output, n):
    data = pd.read_csv(astroabc_output, '\t', index_col=False)
    data = data.tail(n)
    data['wt_with_wt_p'] = 1
    data['del_with_del_p'] = 0
    data = data.rename(columns={
     'param#0 ':'wt_with_del_p',
     ' param#1 ':'wt_with_syn_p',
     ' param#2 ':'del_with_wt_p',
     ' param#3 ':'del_with_syn_p',
     ' param#4 ':'syn_with_wt_p',
     ' param#5 ':'syn_with_del_p',
     ' param#6 ':'syn_with_syn_p',
     ' param#7 ':'triple_wt',
     ' param#8 ':'triple_del',
     ' param#9 ':'triple_syn',
     ' param#10 ':'starting_n_syn'})
    return data

# multi panel graph

def prior_posterior_graph(data, outpath):
    fig, axes = plt.subplots(nrows=4, ncols=3)
    fig.set_size_inches(8,8)
    fig.tight_layout()
    axes = axes.flatten()
    fig.subplots_adjust(hspace=0.5)
    fig.text(0.5, -0.02, 'Payoff', ha='center', fontsize = 14)
    fig.text(-0.02, 0.5, 'Frequency', va='center', rotation='vertical', fontsize = 14)
    for i, column in enumerate(['wt_with_wt_p', 'wt_with_del_p', 'wt_with_syn_p', 'del_with_wt_p', 'del_with_del_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'syn_with_syn_p', 'triple_wt', 'triple_del', 'triple_syn']):
        if column in ['wt_with_del_p', 'wt_with_syn_p', 'syn_with_syn_p', 'triple_wt']:
            axes[i].hist(np.linspace(0,1,len(data)), alpha=0.2, color='#267BB8')
        elif column in ['del_with_wt_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'triple_del', 'triple_syn']:
            axes[i].hist(np.linspace(0,4,len(data)), alpha=0.2, color='#267BB8')
        
        if column not in ['del_with_del_p', 'wt_with_wt_p']:
            data[column].plot.hist(ax=axes[i], color='#267BB8')
            axes[i].set_xlabel(None)
            axes[i].set_ylabel(None)
            axes[i].set_xlim(0,3.8)
        else:   
            axes[i].set_visible(False)
                        
        if 'with' in column and column != 'syn_with_syn_p':
            axes[i].set_title(' when infecting with '.join(column.replace('_p', '').split('_with_')).replace('wt', 'WT').replace('del', '$\Delta$1764').replace('syn', 'A1664G'), fontsize=10)
        elif 'triple' in column:
            axes[i].set_title(column.replace('triple_', '').replace('wt', 'WT').replace('del', '$\Delta$1764').replace('syn', 'A1664G') + ' during triple infection', fontsize=10)
        elif column == 'syn_with_syn_p':
            axes[i].set_title('A1664G alone', fontsize=10)
        axes[i].grid(which='major', alpha=0.2, linestyle='-')
    fig.text(0.18, 0.88, 'WT alone\nfixed at 1', ha='center', va='center', fontsize=10)
    fig.text(0.5, 0.65, '$\Delta$1764 alone\nfixed at 0', ha='center', va='center', fontsize=10)
    prior_patch = mpatches.Patch(color='#267BB8', label='Prior', alpha=0.2)
    posterior_patch = mpatches.Patch(color='#267BB8', label='Posterior')
    fig.legend(handles=[posterior_patch, prior_patch], bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    fig.savefig(outpath, dpi=800, bbox_inches='tight')
    return
    
# cloud graph for best simulations

def cloud_graph(data, real_data, outpath):
    real_data = pd.read_csv(real_data)
    simulations = []
    for row in tqdm(data.itertuples()):
        payoff_matrix={'wt':{'wt':1, 'del':row.wt_with_del_p, 'syn':row.wt_with_syn_p},
                             'del':{'wt':row.del_with_wt_p, 'del':0, 'syn':row.del_with_syn_p}, 
                             'syn':{'wt':row.syn_with_wt_p, 'del':row.syn_with_del_p, 'syn':row.syn_with_syn_p}}
        triple_payoff_matrix = {'wt':row.triple_wt, 'del':row.triple_del, 'syn':row.triple_syn}
        simulation = simulate(payoff_matrix=payoff_matrix, triple_payoff_matrix=triple_payoff_matrix)
        simulations.append(simulation)
    s = pd.concat(simulations)
    s = s[s.passage.isin(real_data.passage.tolist())]
    fig, ax = plt.subplots(nrows=1, ncols=1)
    fig.set_size_inches(5.5,3)
    #plt.style.use('classic')
    #fig.patch.set_facecolor('xkcd:white')
    ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().syn_frequency, s.groupby('passage').max().reset_index().syn_frequency, color='#F49D09',alpha=0.35)   
    ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().del_frequency, s.groupby('passage').max().reset_index().del_frequency, color='#F50202', alpha=0.35) 
    real_data.plot(x='passage', y='del_frequency', color='#F50202', ax=ax, label=r'$\Delta$1764    $\bigtriangleup$', linewidth=1.5)
    real_data.plot(x='passage', y='syn_frequency', color='#F49D09', ax=ax,  label=r'A1664G  $\bigcirc$', linewidth=1.5)
    ax.set_xlabel('Passage', fontsize = 14)
    ax.set_ylabel('Mutation Frequency', fontsize=14)
    ax.set_xticks([0,5,10,15,20,])
    ax.grid(which='major', alpha=0.3, linestyle='-')
    
    handles, labels = ax.get_legend_handles_labels()
    patch = mpatches.Patch(color='grey', label='Best fit\nsimulations', alpha=0.3)
    handles.append(patch) 
    ax.legend(handles=handles, loc='upper left')
    fig.savefig(outpath, dpi=800, bbox_inches='tight')
    return 
#
#prior_posterior_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.n10000.txt', 10000), 'X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.n10000.prior_posterior.png')
#prior_posterior_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.txt', 5000), 'X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.prior_posterior.png')
#prior_posterior_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.pertubation1.txt', 5000), 'X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.pertubation1.prior_posterior.png')
#
#prior_posterior_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.n10000.txt', 10000), 'X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.n10000.prior_posterior.png')
#prior_posterior_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.txt', 5000), 'X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.prior_posterior.png')
#prior_posterior_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.pertubation1.txt', 5000), 'X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.pertubation1.prior_posterior.png')


#cloud_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.n10000.txt', 10000), 'X:/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv', 'X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.n10000.prior_posterior.png')
#cloud_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.txt', 10000), 'X:/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv', 'X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.prior_posterior.png')
#cloud_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.pertubation1.txt', 10000), 'X:/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv', 'X:/volume2/noam/cheater_model/abc_smc/37B_all_passages.pertubation1.prior_posterior.png')
#
#cloud_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.n10000.txt', 10000), 'X:/volume2/noam/cheater_model/37A_realdata_until_p23_for_fit.csv', 'X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.n10000.prior_posterior.png')
#cloud_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.txt', 10000), 'X:/volume2/noam/cheater_model/37A_realdata_until_p23_for_fit.csv', 'X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.prior_posterior.png')
#cloud_graph(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.pertubation1.txt', 10000), 'X:/volume2/noam/cheater_model/37A_realdata_until_p23_for_fit.csv', 'X:/volume2/noam/cheater_model/abc_smc/37A_all_passages.pertubation1.prior_posterior.png')

#cloud_graph(get_best_parameters('/sternadi/home/volume2/noam/cheater_model/abc_smc/37B_all_passages.n10000.txt', 10000), '/sternadi/home/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/37B_all_passages.n10000.cloud.png')
#cloud_graph(get_best_parameters('/sternadi/home/volume2/noam/cheater_model/abc_smc/37B_all_passages.txt', 5000), '/sternadi/home/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/37B_all_passages.cloud.png')
#cloud_graph(get_best_parameters('/sternadi/home/volume2/noam/cheater_model/abc_smc/37B_all_passages.pertubation1.txt', 5000), '/sternadi/home/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/37B_all_passages.pertubation1.cloud.png')
#
#cloud_graph(get_best_parameters('/sternadi/home/volume2/noam/cheater_model/abc_smc/37A_all_passages.n10000.txt', 10000), '/sternadi/home/volume2/noam/cheater_model/37A_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/37A_all_passages.n10000.cloud.png')
#cloud_graph(get_best_parameters('/sternadi/home/volume2/noam/cheater_model/abc_smc/37A_all_passages.txt', 5000), '/sternadi/home/volume2/noam/cheater_model/37A_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/37A_all_passages.cloud.png')
#cloud_graph(get_best_parameters('/sternadi/home/volume2/noam/cheater_model/abc_smc/37A_all_passages.pertubation1.txt', 5000), '/sternadi/home/volume2/noam/cheater_model/37A_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/37A_all_passages.pertubation1.cloud.png')


#cloud_graph(get_best_parameters('/sternadi/home/volume2/noam/cheater_model/abc_smc/new_lim_only_wt_with_del/37B_abc.txt', 10000), '/sternadi/home/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/new_lim_only_wt_with_del/37B_abc.txt.cloud.png')
#cloud_graph(get_best_parameters('/sternadi/home/volume2/noam/cheater_model/abc_smc/new_lim_only_wt_with_del/37A_abc.txt', 10000), '/sternadi/home/volume2/noam/cheater_model/37A_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/new_lim_only_wt_with_del/37A_abc.txt.cloud.png')


#def cloud_graph_starting_n(data, real_data, outpath):
#    real_data = pd.read_csv(real_data)
#    simulations = []
#    for row in tqdm(data.itertuples()):
#        payoff_matrix={'wt':{'wt':1, 'del':row.wt_with_del_p, 'syn':row.wt_with_syn_p},
#                             'del':{'wt':row.del_with_wt_p, 'del':0, 'syn':row.del_with_syn_p}, 
#                             'syn':{'wt':row.syn_with_wt_p, 'del':row.syn_with_del_p, 'syn':row.syn_with_syn_p}}
#        triple_payoff_matrix = {'wt':row.triple_wt, 'del':row.triple_del, 'syn':row.triple_syn}
#        simulation = simulate(payoff_matrix=payoff_matrix, triple_payoff_matrix=triple_payoff_matrix, syn_initial_count=row.starting_n_syn)
#        simulations.append(simulation)
#    s = pd.concat(simulations)
#    s = s[s.passage.isin(real_data.passage.tolist())]
#    fig, ax = plt.subplots(nrows=1, ncols=1)
#    fig.set_size_inches(5,3)
#    #plt.style.use('classic')
#    #fig.patch.set_facecolor('xkcd:white')
#    ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().syn_frequency, s.groupby('passage').max().reset_index().syn_frequency, color='#F49D09',alpha=0.35)   
#    ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().del_frequency, s.groupby('passage').max().reset_index().del_frequency, color='#F50202', alpha=0.35) 
#    real_data.plot(x='passage', y='del_frequency', color='#F50202', ax=ax, label=r'$\Delta$1764    $\bigtriangleup$', linewidth=1.5)
#    real_data.plot(x='passage', y='syn_frequency', color='#F49D09', ax=ax,  label=r'A1664G  $\bigcirc$', linewidth=1.5)
#    ax.set_xlabel('Passage', fontsize = 14)
#    ax.set_ylabel('Mutation Frequency', fontsize=14)
#    ax.set_xticks([0,5,10,15,20,])
#    ax.grid(which='major', alpha=0.2, linestyle='-')
#    fig.savefig(outpath, dpi=800, bbox_inches='tight')
#    return 
#
#cloud_graph_starting_n(get_best_parameters('/sternadi/home/volume2/noam/cheater_model/abc_smc/new_lim_only_wt_with_del/37A_abc_starting_n.txt', 10000), '/sternadi/home/volume2/noam/cheater_model/37A_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/new_lim_only_wt_with_del/37A_abc_starting_n.txt.cloud.png')


### single cheater
def cloud_graph_single_cheater(data, real_data, outpath):
    real_data = pd.read_csv(real_data)
    real_data = real_data[real_data.passage <= 15]
    simulations = []
    data = pd.read_csv(data, '\t', index_col=False)
    data = data.tail(10000)
    #data = data.tail(100)
    data = data.rename(columns={'param#0 ':'wt_with_del_p', ' param#1 ':'del_with_wt_p'})
    for row in tqdm(data.itertuples()):
        payoff_matrix={'wt':{'wt':1, 'del':row.wt_with_del_p, 'syn':0},
                             'del':{'wt':row.del_with_wt_p, 'del':0, 'syn':0}, 
                             'syn':{'wt':0, 'del':0, 'syn':0}}
        triple_payoff_matrix = {'wt':0, 'del':0, 'syn':0}
        simulation = simulate(payoff_matrix=payoff_matrix, triple_payoff_matrix=triple_payoff_matrix, syn_initial_count=0, iterations_within_passages=2)
        simulations.append(simulation)
    s = pd.concat(simulations)
    s = s[s.passage.isin(real_data.passage.tolist())]
    fig, ax = plt.subplots(nrows=1, ncols=1)
    fig.set_size_inches(5.5,3)
    #plt.style.use('classic')
    #fig.patch.set_facecolor('xkcd:white')
    #ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().syn_frequency, s.groupby('passage').max().reset_index().syn_frequency, color='#F49D09',alpha=0.35)   
    ax.fill_between(s.groupby('passage').min().reset_index().passage, s.groupby('passage').min().reset_index().del_frequency, s.groupby('passage').max().reset_index().del_frequency, color='#F50202', alpha=0.35) 
    real_data.plot(x='passage', y='del_frequency', color='#F50202', ax=ax, label=r'$\Delta1764  \bigtriangleup$', linewidth=1.5)
    #real_data.plot(x='passage', y='syn_frequency', color='#F49D09', ax=ax,  label='A1664G', linewidth=1.5)
    ax.set_xlabel('Passage', fontsize = 14)
    ax.set_ylabel('Mutation Frequency', fontsize=14)
    ax.set_xticks([0,5,10,15])
    ax.grid(which='major', alpha=0.3, linestyle='-')
    
    handles, labels = ax.get_legend_handles_labels()
    patch = mpatches.Patch(color='grey', label='Best fit\nsimulations', alpha=0.3)
    handles.append(patch) 
    ax.legend(handles=handles, loc='upper left')
    fig.savefig(outpath, dpi=800, bbox_inches='tight')
    return

#cloud_graph_single_cheater('X:/volume2/noam/cheater_model/abc_smc/37B_single_cheater.txt', 'X:/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv', 'X:/volume2/noam/cheater_model/abc_smc/37B_single_cheater.txt.cloud.png')
#cloud_graph_single_cheater('/sternadi/home/volume2/noam/cheater_model/abc_smc/37B_single_cheater.txt', '/sternadi/home/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/37B_single_cheater.txt.cloud.png')
#cloud_graph_single_cheater('/sternadi/home/volume2/noam/cheater_model/abc_smc/37A_single_cheater.txt', '/sternadi/home/volume2/noam/cheater_model/37A_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/37A_single_cheater.txt.cloud.png')

cloud_graph_single_cheater('/sternadi/home/volume2/noam/cheater_model/abc_smc/new_single_cheater_2_iterations_per_passage/37B', '/sternadi/home/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/new_single_cheater_2_iterations_per_passage/37B_cloud.png')
cloud_graph_single_cheater('/sternadi/home/volume2/noam/cheater_model/abc_smc/new_single_cheater_2_iterations_per_passage/37A', '/sternadi/home/volume2/noam/cheater_model/37A_realdata_until_p23_for_fit.csv', '/sternadi/home/volume2/noam/cheater_model/abc_smc/new_single_cheater_2_iterations_per_passage/37A_cloud.png')



def fix_symbols(string):
    return string.replace('\Delta1764', r'\bigtriangleup').replace('A1664G', r'\bigcirc')


## single cheater prior posterior graph
def prior_posterior_graph_single_cheater(astroabc_output, n, outpath):
    data = pd.read_csv(astroabc_output, '\t', index_col=False)
    data = data.tail(n)
    data['wt_with_wt_p'] = 1
    data['del_with_del_p'] = 0
    data = data.rename(columns={
     'param#0 ':'wt_with_del_p',
     ' param#1 ':'del_with_wt_p',})
    fig, axes = plt.subplots(nrows=2, ncols=2)
    fig.set_size_inches(4,2.75)
    fig.tight_layout()
    axes = axes.flatten()
    
    temp_fig, temp_ax = plt.subplots(nrows=1, ncols=1)
    
    #fig.subplots_adjust(hspace=0.5)
    for i, column in enumerate(['wt_with_wt_p', 'wt_with_del_p','del_with_wt_p', 'del_with_del_p']):
        if column:
            if column in ['wt_with_del_p']:
                a,b,c = temp_ax.hist(data[column], bins=np.linspace(0, 4, 21))
                axes[i].hist(np.linspace(0,1,max(a)/2), alpha=0.4, color='#C5C6C6', bins=1)
            elif column in ['del_with_wt_p']:
                a,b,c = temp_ax.hist(data[column], bins=np.linspace(0, 4, 21))
                axes[i].hist(np.linspace(0,4,max(a)/2), alpha=0.4, color='#C5C6C6', bins=1)
            if column not in ['del_with_del_p', 'wt_with_wt_p']:
                axes[i].set_xticks([0,1,2,3,4])
                axes[i].set_xlim(-0.1,4.1)
                data[column].plot.hist(ax=axes[i], color='#267BB8', bins=np.linspace(0, 4, 21))
                axes[i].set_xlabel(None)
                axes[i].set_ylabel(None)
                axes[i].set_yticks([])
                #axes[i].axvline(data[column].median(), color='black', linestyle='dashed', linewidth=1)
            else:    
                axes[i].set_visible(False)
            if 'with' in column:
                title = 'W$_\mathrm{' + '|'.join(column.replace('_p', '').split('_with_')).replace('wt', 'WT').replace('del', '\Delta1764').replace('syn', 'A1664G') + '}$'
                axes[i].set_title(fix_symbols(title), fontsize=12)
            axes[i].grid(which='major', alpha=0.3, linestyle='-') 
    prior_patch = mpatches.Patch(color='#C5C6C6', label='Prior', alpha=0.4)
    posterior_patch = mpatches.Patch(color='#267BB8', label='Posterior')
    axes[2].legend(handles=[posterior_patch, prior_patch], bbox_to_anchor=(1.1,-0.7), loc='upper center', borderaxespad=0., ncol=2)
    fig.text(0.26, 0.8, fix_symbols('W$_\mathrm{WT|WT}$\nfixed at 1'), ha='center', va='center', fontsize=12)
    fig.text(0.76, 0.3, fix_symbols('W$_\mathrm{\Delta1764|\Delta1764}$\nfixed at 0'), ha='center', va='center', fontsize=12)
    text1 = fig.text(0.5, -0.12, 'Fitness\n', ha='center', fontsize = 14)
    text2 = fig.text(0, 0.5, 'Frequency', va='center', rotation='vertical', fontsize = 14)
    fig.savefig(outpath, dpi=800, bbox_inches='tight')
    return

#prior_posterior_graph_single_cheater('X:/volume2/noam/cheater_model/abc_smc/37B_single_cheater.txt', 10000, 'X:/volume2/noam/cheater_model/abc_smc/37B_single_cheater.txt.prior_posterior.png')
#prior_posterior_graph_single_cheater('X:/volume2/noam/cheater_model/abc_smc/37A_single_cheater.txt', 10000, 'X:/volume2/noam/cheater_model/abc_smc/37A_single_cheater.txt.prior_posterior.png')

prior_posterior_graph_single_cheater('X:/volume2/noam/cheater_model/abc_smc/new_single_cheater_2_iterations_per_passage/37B', 10000, 'X:/volume2/noam/cheater_model/abc_smc/new_single_cheater_2_iterations_per_passage/37B.prior_posterior.png')
prior_posterior_graph_single_cheater('X:/volume2/noam/cheater_model/abc_smc/new_single_cheater_2_iterations_per_passage/37A', 10000, 'X:/volume2/noam/cheater_model/abc_smc/new_single_cheater_2_iterations_per_passage/37A.prior_posterior.png')
