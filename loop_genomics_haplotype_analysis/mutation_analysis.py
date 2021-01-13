import glob

import pandas as pd
pd.set_option('display.max_rows', 500)
# pd.set_option('display.width', 3000)

import seaborn as sns
sns.set_context("poster")
import matplotlib.pyplot as plt

from RG_HIVC_analysis.diversity_calculation import pis_calc
from loop_genomics_haplotype_analysis.summarize_visualize_strain_results import env_samples, shafer_root_dir

def get_freqs(freq_filename = 'freqs_strain_data_all_samples_v6_with_new_pipelines.csv', filter = True, remove_strain_data = True):
    freqs_all = pd.read_csv('{}/results_merged/{}'.format(shafer_root_dir, freq_filename))
    print('data loaded')

    # select env samples only
    freqs_all = freqs_all[freqs_all['sample'].isin(env_samples)]

    # remove strain data
    if remove_strain_data:
        freqs_all = freqs_all.drop(columns=['strain'])
        freqs_all = freqs_all.drop_duplicates()

    # validation
    print('len(freqs_all): {}'.format(len(freqs_all)))
    print(freqs_all.groupby('sample').count()['Pos'])

    # visualize un-alignment between pipelines
    # print(freqs_all.groupby(['sample', 'source']).max()['Pos'])

    # filter
    # ** some duplication here but ok
    if filter:
        freqs_all['Base_count'] = freqs_all['Read_count'] * freqs_all['Freq']
        freqs_all = freqs_all[(freqs_all['Rank'] != 0)
                                   & (freqs_all['Prob'] > 0.8)
                                   & (freqs_all['Base_count'] > 20)
                                   ]
        print('len(freqs_all): {}'.format(len(freqs_all)))
        print(freqs_all.groupby('sample').count()['Pos'])
    return freqs_all

def cd4_diversity_regression(plot= False):
    # get freqs
    freqs = get_freqs(filter=False,
                      remove_strain_data = False)

    # side-analysis- coverage stats
    print(freqs.groupby('source').agg(['mean', 'median'])['Read_count'])

    # custom filtering
    freqs['Base_count'] = freqs['Read_count'] * freqs['Freq']
    freqs = freqs[(freqs['Prob'] > 0.8)
                  & (freqs['Base_count'] > 20)
                  ]
    # filtering validation
    print('len(freqs_all): {}'.format(len(freqs)))

    # filter interesting samples
    freqs['sample_source'] = freqs['sample'] + '_' + freqs['source']
    freqs = freqs[freqs['sample_source'].isin(
        ['env9_loop', 'env10_loop', 'env12_loop', 'env7_loop', 'env11_loop', 'env15_loop', 'env2_accungs',
         'env6_loop', 'env8_loop', 'env4_loop', 'env3_loop', 'env5_loop', 'env13_loop', 'env14_loop', 'env1_accungs'])]

    # get diversity per sample
    # TODO- verify pi calc
    pi_rates_by_sample = pis_calc(data=freqs, pivot_cols=['sample'])
    pi_rates_by_sample = pi_rates_by_sample.groupby('sample').mean().reset_index()
    # print(pi_rates_by_sample)

    # get CD4 levels
    cd4_counts = pd.DataFrame(data={
        'sample': ['env9', 'env11', 'env2', 'env4', 'env15', 'env10', 'env7', 'env12', 'env13', 'env1', 'env8', 'env5',
                   'env6', 'env14', 'env3'],
        'cd4_count': [842.0, 456.0, 401.0, 393.0, 377.0, 322.0, 308.0, 216.0, 129.0, 90.0, 33.0, 23.0, 23.0, 22.0,
                      3.0]})
    # print(cd4_counts)

    # merge by sample
    plot_data = pd.merge(pi_rates_by_sample, cd4_counts, on='sample', how='inner')
    plot_data = plot_data.rename(columns={'Pi': 'Pi diversity', 'cd4_count': 'CD4 count'})
    print(plot_data)
    plot_data.to_csv('{}/results_merged/CD4_diversity_rates.csv'.format(shafer_root_dir))

    # plot
    if plot:
        g = sns.lmplot(x='CD4 count',
                       y='Pi diversity',
                       height=10,
                       data=plot_data);
        plt.xlim(0, 900)
        plt.ylim(0, 0.14)
        plt.savefig(fname='{}/results_merged/CD4_diversity_regression.pdf'.format(shafer_root_dir))
        plt.show()


def get_pi_rates():
    run_name = 'new_pipeline_loop_v6_blast_2500'
    pipeline_dir = '/sternadi/nobackup/volume1/HIVB_shafer_accungs_seq/%s' % run_name
    freqs = glob.glob('%s/env*/freqs.tsv' % pipeline_dir)

    for freq in freqs:
        # TODO- imp
        # read
        # add run_name column
        # add sample column
        # transform to old format


    # filter
    # TODO- imp & move to general freqs
    freqs['Base_count'] = freqs['Read_count'] * freqs['Freq']
    freqs = freqs[(freqs['Prob'] > 0.8)
                  & (freqs['Base_count'] > 20)
                  ]
    filtered_freq_df = apply_pi_related_filters(freq_df, freq_threshold= 0.001, min_read_count= 100)

    # get pi rates
    # TODO- imp & move to general freqs
    pi_rates_by_sample = pis_calc(data=freqs, pivot_cols=['sample'])
    pi_rates_by_sample = pi_rates_by_sample.groupby('sample').mean().reset_index()

if __name__ == '__main__':
    cd4_diversity_regression()