
import pandas as pd
from tqdm import tqdm
import argparse
import sys
import numpy as np
import os
sys.path.append('/sternadi/home/volume2/noam/scripts/cheaters_model/')
sys.path.append('X:/volume2/noam/scripts/cheaters_model')
from model_triple_payoffs import simulate

# columns for fixation order
FIXATION = ['wt_with_del_p', 'wt_with_syn_p', 'del_with_wt_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'syn_with_syn_p', 'triple_wt', 'triple_del', 'triple_syn']

    
def get_graph_similarity(real_data, model_data):
    # inner merge, only keeps values for passges in both datasets
    df = pd.merge(real_data, model_data[real_data.columns.tolist()], on='passage', suffixes=['_real','_model'])
    df['syn_diff'] = (df.syn_frequency_model - df.syn_frequency_real).abs()
    df['del_diff'] = (df.del_frequency_model - df.del_frequency_real).abs()
    df['wt_diff'] = (df.wt_frequency_model - df.wt_frequency_real).abs()
    df['total_diff'] = df.del_diff + df.wt_diff + df.syn_diff
    return df.total_diff.mean()

def calculate_row_diff(row, real_data):
    max_passage = real_data.passage.max()
    payoff_matrix={'wt':{'wt':1, 'del':row.wt_with_del_p, 'syn':row.wt_with_syn_p},
               'del':{'wt':row.del_with_wt_p, 'del':0, 'syn':row.del_with_syn_p}, 
               'syn':{'wt':row.syn_with_wt_p, 'del':row.syn_with_del_p, 'syn':row.syn_with_syn_p}}
    triple_payoff_matrix = {'wt':row.triple_wt, 'del':row.triple_del, 'syn':row.triple_syn}
    model_data = simulate(payoff_matrix=payoff_matrix, triple_payoff_matrix=triple_payoff_matrix, passages=max_passage + 1, iterations_within_passages=row.iterations_within_passage)
    diff = get_graph_similarity(real_data, model_data)
    wt_final_passage = model_data[~model_data.isna().any(axis=1)].passage.max()
    mean_killing_ratio = model_data[model_data.passage == model_data[model_data.del_frequency == model_data.del_frequency.max()].passage.min()].mean_killing.min() / model_data[model_data.passage == 1].mean_killing.min()
    return [diff, wt_final_passage, mean_killing_ratio]

def row_optimize(row, real_data, out_file):
    row = row.to_frame().T
    dfs = []
    for c in FIXATION:
        current_value = row[c].min()
        print(c)
        print(current_value)
        new = pd.DataFrame(np.linspace(current_value - 0.5, current_value + 0.5, num=11), columns=[c])
        new = new[new[c] >= 0]
        old = row.drop(columns=[c])
        new['a'] = 1
        old['a'] = 1
        df = pd.merge(old, new, on='a').drop(columns=['a', 'iterator'])
        df['diff'], df['wt_final_passage'], df['mean_killing_ratio'] = zip(*df.apply(calculate_row_diff, real_data=real_data, axis=1))
        dfs.append(df)
        # fixate best values for current column
        #row[c] = df[df['diff'] == df[(df.wt_final_passage == 24) & (df.mean_killing_ratio < 0.9)]['diff'].min()][c].min()
        row[c] = df[df['diff'] == df[(df.wt_final_passage == 24)]['diff'].min()][c].min()
    df = pd.concat(dfs, ignore_index=True, sort=False)
    df = df[FIXATION + ['diff', 'wt_final_passage','mean_killing_ratio']]
    if not os.path.isfile(out_file): 
        df.to_csv(out_file, header=True, index=False, mode='a')
    else:
        df.to_csv(out_file, header=False, index=False, mode='a')
    return 


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--out_directory", type=str, required=True)
    parser.add_argument('-i', '--input_iterations_df', type=str, required=True)
    parser.add_argument('-j', '--job_id', type=int, required=True)
    parser.add_argument('-r', '--real_data', type=str, required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    #real_data = pd.read_csv('/sternadi/home/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv')
    real_data = pd.read_csv(args.real_data)
    iterations_df = pd.read_csv(args.input_iterations_df + '/' + str(args.job_id) + '.csv')
    #iterations_df = iterations_df[iterations_df.iterator == args.job_id]
    tqdm.pandas()
    iterations_df.progress_apply(row_optimize, real_data=real_data, out_file=args.out_directory + '/' + str(args.job_id) + '.csv', axis=1)