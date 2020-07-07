import pandas as pd
import argparse
import sys
import numpy as np
sys.path.append('/sternadi/home/volume2/noam/scripts/cheaters_model/')
#sys.path.append('X:/volume2/noam/scripts/cheaters_model')
from model_triple_payoffs import simulate
#sys.path.insert(0, 'X:/volume2/noam/cheater_model/astroABC')
sys.path.insert(0, '/sternadi/home/volume2/noam/cheater_model/astroABC')
import astroabc


PARAMETERS = ['wt_with_del_p', 'del_with_wt_p',]

PRIORS =  [                    ('uniform', [0,1]), 
           ('uniform', [0,4]),                     ]


### calculate distance between real data and model data
def calculate_particle_diff(real_data, model_data):
    if type(model_data) != pd.core.frame.DataFrame: # impossibly high diff when parameters not in prior
        return np.inf
    # inner merge, only keeps values for passges in both datasets
    df = pd.merge(real_data, model_data[real_data.columns.tolist()], on='passage', suffixes=['_real','_model'])
    #df['syn_diff'] = (df.syn_frequency_model - df.syn_frequency_real).abs()
    df['del_diff'] = (df.del_frequency_model - df.del_frequency_real).abs()
    #df['total_diff'] = df.del_diff + df.syn_diff
    diff = df.del_diff.mean()
    wt_final_passage = model_data[~model_data.isna().any(axis=1)].passage.max()
    if wt_final_passage < MAX_PASSAGE: # impossibly high diff when WT crashes
        return np.inf
    return diff


# simulate data for a parameter set
def simulation(parameter_list):
    if all([i[0] >= i[1][0] and i[0] <= i[1][1] for i in  list(zip(parameter_list, [i[1] for i in PRIORS]))]):
        payoff_matrix={'wt':{'wt':1, 'del':parameter_list[0], 'syn':0},
                   'del':{'wt':parameter_list[1], 'del':0, 'syn':0}, 
                   'syn':{'wt':0, 'del':0, 'syn':0}}
        triple_payoff_matrix = {'wt':1, 'del':1, 'syn':1}
        model_data = simulate(payoff_matrix=payoff_matrix, triple_payoff_matrix=triple_payoff_matrix, passages=MAX_PASSAGE, syn_initial_count=0, iterations_within_passages=2)
        return model_data
    else:
        return np.inf


#### run abc smc        
############  
# problem - particles outside priors (for exmaple including negatives) are still included in the N count, although weighted to 0, so we end up using less than N parameters at each level t
# solution by giving infinity rho results when parameters in particles are outside priors
#############

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_path", type=str, help="a path to an output directory, optional", required=True)
    parser.add_argument("-i", "--real_data_file", type=str, help="input directory with .seq files",required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
        
    #REAL_DATA = pd.read_csv('/sternadi/home/volume2/noam/cheater_model/37B_realdata_until_p23_for_fit.csv')
    REAL_DATA = pd.read_csv(args.real_data_file)
    REAL_DATA = REAL_DATA[REAL_DATA.passage <=15]
    MAX_PASSAGE = REAL_DATA.passage.max()
    prop={'dfunc':calculate_particle_diff, 
          'outfile':args.output_path, 
          'tol_type':'exp', 
          'verbose':1 , 
          'pert_kernel':2, # multivariate perturbation based on local covariance
          'mp':True} 
    sampler = astroabc.ABC_class(2,10000,REAL_DATA,[0.3,0.03],15,PRIORS,**prop)
    sampler.sample(simulation)
    



