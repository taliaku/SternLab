import pandas as pd
from tqdm import tqdm
import itertools

PAYOFF_OPTIONS = (0, 0.25, 0.5, 1, 2, 3)
ITERATIONS_PER_PASSAGE = [1]

def find_best_parameters_euclidean():
    payoff_matrix = itertools.product(PAYOFF_OPTIONS, repeat=10)
    results = []
    for p in tqdm(payoff_matrix):
        for iteration in ITERATIONS_PER_PASSAGE:
            results.append(p + (iteration,))
    results = pd.DataFrame(results, columns=['wt_with_del_p', 'wt_with_syn_p', 'del_with_wt_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'syn_with_syn_p', 'triple_wt', 'triple_del', 'triple_syn', 'iterations_within_passage'])
    #results = results[(results.del_with_wt_p >=2) & (results.syn_with_wt_p >= 2)]
    results = results[(results.syn_with_syn_p < 1) & (results.syn_with_syn_p > 0)]
    results['iterator'] = results.index % 100
    return results

results = find_best_parameters_euclidean()
for i in tqdm(results.iterator.unique()):
    results[results.iterator == i].to_csv('/sternadi/home/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/iteration_job_ids/' + str(i) + '.csv', index=False)
#results.to_csv('/sternadi/home/volume2/noam/cheater_model/triple_infections/iteration_job_ids.csv', index=False)

 

##### for optimization, top 1% options
# with mean killing ratio under 0.9
results = pd.read_csv('X:/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/all_results.csv')
results = results[(results.mean_killing_ratio < 0.9) & (results.wt_final_passage == 24)]
diff_cutoff = results['diff'].quantile(0.01)
results_to_optimize = results[results['diff'] < diff_cutoff]
results_to_optimize = results_to_optimize.drop(columns=['Unnamed: 0'])
results_to_optimize = results_to_optimize.reset_index(drop=True)
results_to_optimize['iterator'] = results_to_optimize.index % 100
for i in tqdm(results_to_optimize.iterator.unique()):
    results_to_optimize[results_to_optimize.iterator == i].to_csv('X:/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/optimize_iteration_job_ids_mean_killing_0.9/' + str(i) + '.csv', index=False)
    
# no mean killing ratio filter
results = pd.read_csv('X:/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/all_results.csv')
results = results[(results.wt_final_passage == 24)]
diff_cutoff = results['diff'].quantile(0.01)
results_to_optimize = results[results['diff'] < diff_cutoff]
results_to_optimize = results_to_optimize.drop(columns=['Unnamed: 0'])
results_to_optimize = results_to_optimize.reset_index(drop=True)
results_to_optimize['iterator'] = results_to_optimize.index % 100
for i in tqdm(results_to_optimize.iterator.unique()):
    results_to_optimize[results_to_optimize.iterator == i].to_csv('X:/volume2/noam/cheater_model/triple_infections/no_passage_18_20_22/optimize_iteration_job_ids/' + str(i) + '.csv', index=False)
    
    
    
#PAYOFF_OPTIONS = (0, 0.1, 0.25, 0.5, 1, 1.5, 2, 2.5, 3)
#PARAM_NUM = 10
#PAYOFF_OPTIONS = (0, 0.1, 0.25, 0.5)
#PARAM_NUM = 3
#def index_to_params(index):
#    params = []
#    for i in range(PARAM_NUM):    
#        params.append(PAYOFF_OPTIONS[index % len(PAYOFF_OPTIONS)])
#        index = index // len(PAYOFF_OPTIONS)
#    return params
    