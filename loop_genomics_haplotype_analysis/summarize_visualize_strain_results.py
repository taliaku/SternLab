import glob
import pandas as pd
import numpy as np

from freqs_utilities import change_ref_to_consensus

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 3000)

import seaborn as sns
# sns.set_context("poster")

import matplotlib.pyplot as plt

shafer_root_dir = '/Volumes/STERNADILABHOME$/volume1/shared/analysis/HIV_shafer'
accungs_root_dir = '/Volumes/STERNADILABHOME$/volume2/ita/haplotypes'

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
        # length + first\last positions
        # TODO- handle WT lines - first\last pos + length
        sample_strains_df['first_pos'] = sample_strains_df['strain'].apply(lambda x: x.split(',')[0][1:-1])
        print(sample_strains_df['first_pos'])
        sample_strains_df['last_pos'] = sample_strains_df['strain'].apply(lambda x: x.split(',')[-1][1:-1])
        print(sample_strains_df['last_pos'])
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


def plot_compare_associvar_accungs_strains():
    # load associvar strain summary
    associvar_rawdata = pd.read_csv('{}/associvar/strain_summary_v5.csv'.format(shafer_root_dir)).reset_index()
    # load accungs strain summary
    accungs_summarized_data = pd.read_csv('{}/accungs_haplotype_inference/mutation_type_analysis/global_summary_env_v3.csv'.format(shafer_root_dir))

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


def plot_compare_loop_accungs_mutations():
    #  0.fetch freq files
    loop_freq_files = glob.glob('{}/loop_genomics_pipeline/envs_output/*/joined.freqs'.format(shafer_root_dir))
    accungs_freq_files_env = glob.glob('{}/envs_output/*/joined.freqs'.format(accungs_root_dir))
    accungs_freq_files_gag = glob.glob('{}/gags_output/*/joined.freqs'.format(accungs_root_dir))

    files_all = loop_freq_files + accungs_freq_files_env + accungs_freq_files_gag
    print(len(files_all))

    # 1. freq-> dataset
    dfs_all = []
    for file in files_all:
        print(file)
        df = pd.read_csv(file, sep='\t')

        # remove insertions + fix ref
        df = df[df['Pos'] == np.round(df['Pos'])]  # remove insertion
        df = change_ref_to_consensus(df)

        df['sample'] = file.split('/')[-2]
        if file in loop_freq_files:
            df['source'] = 'loop'
        else:
            df['source'] = 'accungs'

        dfs_all.append(df)

    freqs_all = pd.concat(dfs_all)
    print(freqs_all.head())

    # 2. add mutation column according to reference
    freqs_all['mutation'] = freqs_all['Ref'] + '->' + freqs_all['Base']

    # 3. add strain data
    print('add strain data')
    # get data
    strain_df = get_strain_data()
    strain_df.to_csv('{}/results_merged/strain_data.csv'.format(shafer_root_dir), index=False)
    # merge by sample + pos
    tmp = pd.merge(freqs_all, strain_df,
                   how='left',
                   on=['source', 'sample', 'Pos'])
    tmp['strain'] = tmp['strain'].fillna('WT_or_other')

    freqs_all = tmp

    # 4. add syn\non-syn data
    # TODO-imp
    # get ORFs
    # produce consensus codons \ mutated codons
    # TODO- mark adjescent mutations- additional possible codon inference
    #  (mutated-codon-per-position, mutated-codon-with-neighbour-BPs)
    # -> translate to consensus \ mutated AA
    # -> infer syn\non-syn

    # 5. plot
    freqs_all['sample_serial_id'] = freqs_all['sample'].apply(lambda x: x[3:]).astype(int)
    freqs_all['sample_source'] = freqs_all['sample'] + '_' + freqs_all['source']
    freqs_all.sort_values(by=['sample_serial_id','sample_source'], inplace=True)
    print(freqs_all.head())

    # filters for plot
    # base count- 100
    # transitions only
    freqs_all['Base_count'] = freqs_all['Read_count'] * freqs_all['Freq']
    freqs_all_filtered = freqs_all[(freqs_all['Rank'] != 0)
                                    & (freqs_all['Prob'] > 0.8)
                                    & (freqs_all['mutation'].isin(['G->A', 'A->G', 'G->G', 'A->A', 'C->T', 'T->C', 'C->C', 'T->T']))
                                    & (freqs_all['Base_count'] > 100)]


    print(len(freqs_all_filtered))

    strain_names = freqs_all_filtered['strain'].unique().tolist()
    palette = dict(zip(strain_names, sns.color_palette(n_colors=len(strain_names))))
    palette.update({"WT_or_other": "gray"})

    # leave selected samples
    freqs_all_filtered = freqs_all_filtered[~freqs_all_filtered['sample_serial_id'].isin([1,2,4,7])]

    if True:
        g = sns.relplot(x='Pos',
                        y='Freq',
                        col='sample_source',
                        hue='strain',
                        # hue_order=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                        #            'WT_or_other'],
                        palette= palette,
                        style='mutation',
                        col_wrap=3,
                        data=freqs_all_filtered)

        plt.yscale('log')
        plt.ylim(1e-04, 1)

        # extract
        # plt.show()
        plt.savefig(fname='{}/results_merged/base_mutations_comparison_v7.pdf'.format(shafer_root_dir))


def get_strain_data():
    strain_df = []

    # get accungs strain data
    accungs_strain_files = glob.glob(
        '{}/accungs_haplotype_inference/mutation_type_analysis/stretch_data_with_mutation_type/env*_mutations_by_stretch_env_v1.csv'.format(
            shafer_root_dir))

    for file in accungs_strain_files:
        sample = (file.split('/')[-1]).split('_')[0]
        print(sample)
        df = pd.read_csv(file)

        df = df[['Stretch', 'Pos']]
        df['sample'] = sample
        df['source'] = 'accungs'

        # convert to minimal values
        sample_strain_reindex = df['Stretch'].unique().tolist()
        df['Stretch'] = df['Stretch'].apply(lambda x: sample_strain_reindex.index(x))

        strain_df.append(df)

    # get loop strain data
    loop_strain_data = pd.read_csv(
        '/Volumes/STERNADILABHOME$/volume1/shared/analysis/HIV_shafer/associvar/strain_summary_v5.csv')

    loop_strain_data['strain_id'] = loop_strain_data['contig'] + loop_strain_data['strain_id'].astype(str)
    # convert to minimal values
    strain_id_reindex_dict = {"A0": 0,"A1": 1,"A2": 2,"A3": 3,"A4": 4,
                              "B0": 5,"B1": 6,"B2": 7,"B3": 8,"B4": 9,
                              "C0": 10,"C1": 11,"C2": 12,"C3": 13,"C4": 14
                              }
    loop_strain_data['strain_id'] = loop_strain_data['strain_id'].apply(lambda x: strain_id_reindex_dict.get(x))

    loop_strain_data = loop_strain_data[['sample', 'strain_id', 'strain']]

    # Removing WT strains- irrelevant in this plot #TODO- consult & verify
    loop_strain_data = loop_strain_data[loop_strain_data['strain'] != 'WT']

    loop_strain_data['mutation'] = loop_strain_data['strain'].apply(lambda x: x.replace(" ", "").split(','))
    loop_strain_data = loop_strain_data.explode('mutation')
    loop_strain_data['Pos'] = loop_strain_data['mutation'].apply(lambda x: int(float(x[1:-1])))
    loop_strain_data = loop_strain_data.rename(columns={"strain_id": "Stretch"})
    loop_strain_data['source'] = 'loop'
    loop_strain_data = loop_strain_data[['Stretch', 'Pos', 'sample', 'source']]

    # combine data
    strain_df.append(loop_strain_data)
    strain_df = pd.concat(strain_df)

    # cosmetics
    strain_df = strain_df[['source', 'sample', 'Stretch', 'Pos']]
    strain_df = strain_df.rename(columns={"Stretch": "strain"})

    return strain_df


if __name__ == '__main__':
    # summarize_associvar_strain_analysis()
    # plot_compare_associvar_accungs_strains()

    plot_compare_loop_accungs_mutations()