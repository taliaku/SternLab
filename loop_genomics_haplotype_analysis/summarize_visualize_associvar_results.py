import glob
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 3000)

import seaborn as sns
# sns.set_context("poster")

import matplotlib.pyplot as plt

shafer_root_dir = '/Volumes/STERNADILABHOME$/volume1/shared/analysis/HIV_shafer'

env_samples = ['env10', 'env11', 'env12', 'env13', 'env14', 'env15', 'env3', 'env4', 'env5', 'env6','env7', 'env8', 'env9']
# env_samples = ['env10']

def summarize_associvar_strain_analysis():
    associvar_root_dir = '{}/associvar'.format(shafer_root_dir)

    text_summaries = []
    for sample in env_samples:
        print('Sample: {}'.format(sample))

        # get 3 contigs
        # TODO- maybe replace with less filtered reliable.csv
        contigs = glob.glob('{}/{}/associvar_strain_{}_contig_*_freq0.1.csv.reliable.csv'.format(associvar_root_dir, sample, sample))
        print(contigs)

        # filter top 3 haplotypes only (or above some read threshold?)
        contigs_df = [pd.read_csv(x).head() for x in contigs]

        # add strain identifier
        for i, contig in enumerate(contigs_df):
            # print(contig)
            contig.reset_index(inplace=True)
            contig.rename(columns={contig.columns[0]: "strain_id"}, inplace=True)
            contig['contig'] = chr(i+65)
            contig['strain_id'] = contig['strain_id'].astype(str)
            print(contig)

        # concat to one table per sample
        sample_strains_df = pd.concat(contigs_df)
        sample_strains_df['strain'] = sample_strains_df['strain'].apply(str.replace, args = (" ", ""))

        # text summary (row\column per strain)
        # freq + read_count - given
        # length + first\last positions
        # TODO- handle WT lines
        sample_strains_df['first_pos'] = sample_strains_df['strain'].apply(lambda x: x.split(',')[0][1:-1])
        print(sample_strains_df['first_pos'])
        sample_strains_df['last_pos'] = sample_strains_df['strain'].apply(lambda x: x.split(',')[-1][1:-1])
        print(sample_strains_df['last_pos'])
        # TODO- add length
        # sample_strains_df['length'] = sample_strains_df['last_pos'] - sample_strains_df['last_pos']

        # # mut count (founder)
        sample_strains_df['mutation_count'] = sample_strains_df['strain'].apply(lambda x: len(x.split(',')))

        # G->A count (APOBEC3)
        def get_ga_count(strain_string):
            filter_ga = [x for x in strain_string.split(',') if (x[0] == 'G' and x[-1] == 'A')]
            return len(filter_ga)

        sample_strains_df['GA_count'] = sample_strains_df['strain'].apply(get_ga_count)
        sample_strains_df['GA_freq'] = sample_strains_df['GA_count'] / sample_strains_df['mutation_count']

        print(sample_strains_df)

        sample_strains_df['sample'] = sample
        text_summaries.append(sample_strains_df)

    # combined summary (all env_samples)
    text_summaries_df = pd.concat(text_summaries)
    # cosmetics
    col_order = ['sample', 'contig', 'strain_id', 'strain', 'read_count', 'frequency', 'first_pos', 'last_pos',
             'mutation_count', 'GA_count', 'GA_freq']
    text_summaries_df = text_summaries_df[col_order]
    text_summaries_df['strain'] = text_summaries_df['strain'].apply(str.replace, args=(",", ", "))

    # export
    print(text_summaries_df)
    text_summaries_df.to_csv('{}/strain_summary_v6.csv'.format(associvar_root_dir), index=False)


def plot_associvar_accungs_strains():
    # load associvar strain summary
    associvar_rawdata = pd.read_csv('{}/associvar/strain_summary_v5.csv'.format(shafer_root_dir)).reset_index()
    # load accungs strain summary
    accungs_summarized_data = pd.read_csv('{}/accungs_haplotype_inference/mutation_type_summaries/global_summary_env_v3.csv'.format(shafer_root_dir))

    # plot 1 - stretches "strength" (according to mutation count)
    # data template: sample_id, stretch_id, first_pos, last_pos, freq, mutation_count
    associvar_plot1_data = associvar_rawdata[['sample', 'contig', 'strain_id', 'first_pos', 'last_pos', 'frequency', 'mutation_count']]
    associvar_plot1_data['strain_id'] = associvar_plot1_data['contig'] + associvar_plot1_data['strain_id'].astype(str)
    associvar_plot1_data['source'] = 'loop'

    accungs_plot1_data = accungs_summarized_data[['Stretch','total','length','first_pos','last_pos','meandist']].rename(columns={'total': 'mutation_count', 'meandist': 'frequency'})
    accungs_plot1_data[['sample','strain_id']] = accungs_plot1_data['Stretch'].str.split('_',expand=True)
    accungs_plot1_data = accungs_plot1_data[~accungs_plot1_data['sample'].isin(['env1','env2'])]
    accungs_plot1_data['source'] = 'AccuNGS'

    plot1_data = pd.concat([associvar_plot1_data, accungs_plot1_data], join="inner")
    print(plot1_data.head())

    # plot
    plot1_data = plot1_data.melt(id_vars=['source', 'sample', 'strain_id', 'frequency', 'mutation_count'], var_name='type', value_name='pos')
    plot1_data['sample_source'] = plot1_data['sample'] + '_' + plot1_data['source']
    plot1_data.sort_values(by= ['sample_source'], inplace= True)
    print(plot1_data.head())

    # emphasize mutation count (+additional amplification for loop mut. count)
    plot1_data['mutation_count'] = np.where(plot1_data['source'] == 'loop',
                                            plot1_data['mutation_count'] * 20,
                                            plot1_data['mutation_count'] * 5)

    g = sns.relplot(x='pos',
                    y='frequency',
                    col='sample_source',
                    hue='strain_id',
                    size= 'mutation_count',
                    col_wrap=2,
                    kind='line',
                    legend=False,
                    data=plot1_data)

    plt.yscale('log')
    # plt.show()
    plt.savefig(fname= '{}/results_merged/haplotypes_comparison_with_mutation_count_v1.pdf'.format(shafer_root_dir))


    # plot 2 - stretches APOBEC
    # data template: sample_id, stretch_id, pos, mutation, freq
    # TODO- transform both data


if __name__ == '__main__':
    # summarize_associvar_strain_analysis()

    plot_associvar_accungs_strains()