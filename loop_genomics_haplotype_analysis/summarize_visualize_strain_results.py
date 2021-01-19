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
shafer_dir_nobackup = '/Volumes/STERNADILABTEMP$/volume1/HIVB_shafer_accungs_seq'

shafer_root_dir_power9 = '/sternadi/home/volume1/shared/analysis/HIV_shafer'
accungs_root_dir_power9 = '/sternadi/home/volume2/ita/haplotypes'

env_samples = ['env1', 'env2', 'env3', 'env4', 'env5', 'env6','env7', 'env8', 'env9', 'env10', 'env11', 'env12', 'env13', 'env14', 'env15']
# env_samples = ['env10']

def summarize_associvar_strain_analysis():
    associvar_root_dir = '{}/associvar'.format(shafer_root_dir)

    text_summaries = []
    for sample in env_samples:
        print('Sample: {}'.format(sample))

        # get 2 contigs
        contigA = glob.glob('{}/{}/associvar_strain_{}_contigA_freq0.1.csv.reliable.csv'.format(associvar_root_dir, sample, sample))[0]
        contigB = glob.glob('{}/{}/associvar_strain_{}_contigB_freq0.1.csv.reliable.csv'.format(associvar_root_dir, sample, sample))[0]
        contigs = contigA + contigB
        print(contigs)

        # filter top 5 haplotypes only (or above some read threshold?)
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
    text_summaries_df.to_csv('{}/strain_summary_v7.csv'.format(associvar_root_dir), index=False)


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


def plot_compare_loop_accungs_mutations(run_name,
                                        plot_gag = True,
                                        loop_only = False,
                                        new_pipe_runs = False,
                                        generate_freq_data = True,
                                        generate_strain_data = True,
                                        add_strain_data = True,
                                        plot = True):

    freqs_all = get_freq_data(run_name, new_pipe_runs, generate_freq_data, generate_strain_data, add_strain_data)

    # 4. add syn\non-syn data
    # TODO-imp
    # get ORFs
    # produce consensus codons \ mutated codons
    # TODO- mark adjescent mutations- additional possible codon inference
    #  (mutated-codon-per-position, mutated-codon-with-neighbour-BPs)
    # -> translate to consensus \ mutated AA
    # -> infer syn\non-syn

    # 5. plot
    if plot:
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
                                        & ( (freqs_all['mutation'].isin(['G->A', 'A->G', 'G->G', 'A->A', 'C->T', 'T->C', 'C->C', 'T->T']))
                                          | (freqs_all['indel'].isin(['insertion','deletion'])) )
                                        & (freqs_all['Base_count'] > 100)]
        print(len(freqs_all_filtered))

        # leave selected samples
        # if not plot_gag:
        #     print('filter gag')
        #     print(len(freqs_all_filtered))
        #     freqs_all_filtered = freqs_all_filtered[freqs_all_filtered['sample'].apply(lambda x: x[:3]) == 'gag']
        #     print(len(freqs_all_filtered))
        #

        # freqs_all_filtered = freqs_all_filtered[~freqs_all_filtered['sample_serial_id'].isin([1,2,4,7])]
        freqs_all_filtered = freqs_all_filtered[(freqs_all_filtered['source'] == 'loop') | (freqs_all_filtered['sample_serial_id'].isin([1,2]))]
        # freqs_all_filtered = freqs_all_filtered[freqs_all_filtered['sample_source'] == 'env15_loop']
        print(freqs_all_filtered['sample_source'].unique().tolist())

        # strain coloring
        strain_names = freqs_all_filtered['strain'].unique().tolist()
        print(strain_names)
        palette = dict(zip(strain_names, sns.color_palette(n_colors=len(strain_names))))
        palette.update({"WT_or_other": "gray", 5: "black"})
        # palette.update({"WT_or_other": "gray"})
        # palette.update({"WT_or_other": "gray",
        #                 5: "black",
        #                 12: "black"})

        g = sns.relplot(x='Pos',
                        y='Freq',
                        col='sample_source',
                        # col_order= ['env15_loop'],
                        # col_order= ['env9_loop', 'env10_loop', 'env12_loop', 'env7_loop'],# low diversity
                        # col_order= ['env11_loop', 'env15_loop', 'env2_accungs'], # multiple foudners
                        # col_order= ['env6_loop', 'env8_loop', 'env4_loop',
                        #             'env3_loop', 'env5_loop', 'env13_loop', 'env14_loop', 'env1_accungs'], # high diversity
                        col_order= ['env9_loop', 'env10_loop', 'env12_loop', 'env7_loop',# low diversity
                                    'env11_loop', 'env15_loop', 'env2_accungs', '', # multiple foudners
                                    'env6_loop', 'env8_loop', 'env4_loop', # high diversity, few strains
                                    'env3_loop', 'env5_loop', 'env13_loop', 'env14_loop', 'env1_accungs'], # high diversity, many strains
                        hue='strain',
                        # hue_order=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,'WT_or_other',
                        #            14],
                        palette= palette,
                        # style='mutation',
                        # style_order= ['A->G', 'G->A', 'T->C', 'C->T'],
                        style='indel',
                        style_order= ['insertion', 'deletion', 'mutation'],
                        col_wrap=4,
                        data=freqs_all_filtered)

        plt.yscale('log')
        plt.ylim(1e-04, 1)

        # extract
        # plt.show()
        plt.savefig(fname='{}/results_merged/mutations_strains_{}.pdf'.format(shafer_root_dir, run_name))


def get_freq_data(run_name,
                    new_pipe_runs = False,
                    generate_freq_data = True,
                    generate_strain_data = True,
                    add_strain_data = True
                  ):

    if generate_freq_data:
        #  0.fetch freq files
        files_all = []
        loop_freq_files = glob.glob('{}/loop_genomics_pipeline/envs_output/*/joined.freqs'.format(shafer_root_dir))
        files_all.extend(loop_freq_files)
        accungs_freq_files_env = glob.glob('{}/envs_output/*/joined.freqs'.format(accungs_root_dir))
        files_all.extend(accungs_freq_files_env)
        accungs_freq_files_gag = glob.glob('{}/gags_output/*/joined.freqs'.format(accungs_root_dir))
        files_all.extend(accungs_freq_files_gag)

        if new_pipe_runs:
            # TODO- can run new accungs for env1/2 & gag samples too
            accungs_new = glob.glob('{}/new_pipeline_accungs_v2_15_iter/env*/freqs.tsv'.format(shafer_dir_nobackup))
            files_all.extend(accungs_new)
            loop_new_q10 = glob.glob('{}/new_pipeline_loop_v2_q10/env*/freqs.tsv'.format(shafer_root_dir))
            files_all.extend(loop_new_q10)

            loop_new_long_reads = glob.glob(
                '{}/new_pipeline_loop_v3_long_reads/env*/freqs.tsv'.format(shafer_dir_nobackup))
            # files_all.extend(loop_new_long_reads)

        print(len(files_all))

        # 1. freq-> dataset (old format)
        dfs_all = []
        for file in files_all:
            print(file)
            df = pd.read_csv(file, sep='\t')

            if file in loop_freq_files:
                df['source'] = 'loop'
                df = fix_ref_temp(df)
            elif file in accungs_freq_files_env + accungs_freq_files_gag:
                df['source'] = 'accungs'
                df = fix_ref_temp(df)
            elif (new_pipe_runs and file in loop_new_q10 + loop_new_long_reads + accungs_new):
                df.rename(columns={'ref_pos': 'Pos', 'read_base': 'Base', 'ref_base': 'Ref',
                                   'frequency': 'Freq', 'base_rank': 'Rank', 'coverage': 'Read_count',
                                   'probability': 'Prob'},
                          inplace=True)
                if file in loop_new_q10:
                    df['source'] = 'loop_new_q10'
                elif file in loop_new_long_reads:
                    df['source'] = 'loop_new_long_reads'
                elif file in accungs_new:
                    df['source'] = 'accungs_new'
                else:
                    raise Exception("file belongs nowhere: {}".format(file))

            df['sample'] = file.split('/')[-2]

            dfs_all.append(df)

        freqs_all = pd.concat(dfs_all)

        # validations
        print(freqs_all.head())
        print(freqs_all.groupby('source').count()['Pos'])
        print(freqs_all.groupby('sample').count()['Pos'])

        # 2. add mutation column according to reference
        freqs_all['mutation'] = freqs_all['Ref'] + '->' + freqs_all['Base']
        freqs_all['indel'] = freqs_all['Ref'] + '->' + freqs_all['Base']
        freqs_all['indel'] = np.select(
            [freqs_all['Ref'] == '-', freqs_all['Base'] == '-'],
            ['insertion', 'deletion'],
            default='mutation'
        )

        # 3. add strain data
        if add_strain_data:
            print('add strain data')
            # get data
            strain_df = get_strain_data(generate_data=generate_strain_data)
            strain_df = strain_df.drop(columns=['mutation'])
            # merge by sample + pos
            # TODO- merge on Rank1 only?
            tmp = pd.merge(freqs_all, strain_df,
                           how='left',
                           on=['source', 'sample', 'Pos'])
            tmp['strain'] = tmp['strain'].fillna('WT_or_other')
            freqs_all = tmp

        # export
        freqs_all.to_csv('{}/results_merged/freqs_strain_data_all_samples_{}.csv'.format(shafer_root_dir, run_name),
                         index=False)
    else:
        freqs_all = pd.read_csv(
            '{}/results_merged/freqs_strain_data_all_samples_{}.csv'.format(shafer_root_dir, run_name))
    return freqs_all


def fix_ref_temp(df):
    # remove insertions + fix ref
    # TODO- keep insertions (fix change_ref_to_consensus to handle data with insertions)
    #  or keep dont run this on new pipeline? just check status?
    df = df[df['Pos'] == np.round(df['Pos'])]  # remove insertion
    df = change_ref_to_consensus(df)
    return df


def get_strain_data(generate_data = True):
    if generate_data:
        strain_df = []

        # get accungs strain data
        accungs_strain_files = glob.glob(
            '{}/accungs_haplotype_inference/mutation_type_analysis/stretch_data_with_mutation_type/env*_mutations_by_stretch_env_v5.csv'.format(
                shafer_root_dir))

        for file in accungs_strain_files:
            sample = (file.split('/')[-1]).split('_')[0]
            print(sample)
            df = pd.read_csv(file)

            df = df[['Stretch', 'Pos', 'mutation']]
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
        loop_strain_data['mutation'] = loop_strain_data['mutation'].apply(lambda x: x[0]+x[-1])
        loop_strain_data = loop_strain_data.rename(columns={"strain_id": "Stretch"})
        loop_strain_data['source'] = 'loop'
        loop_strain_data = loop_strain_data[['Stretch', 'Pos', 'mutation', 'sample', 'source']]

        # combine data
        strain_df.append(loop_strain_data)
        strain_df = pd.concat(strain_df)

        # cosmetics
        strain_df = strain_df[['source', 'sample', 'Stretch', 'Pos', 'mutation']]
        strain_df = strain_df.rename(columns={"Stretch": "strain"})
        strain_df['strain'] = strain_df['strain'].astype(int)

        # export
        strain_df.to_csv('{}/results_merged/strain_data_v4_param2.csv'.format(shafer_root_dir), index=False)
    else:
        strain_df = pd.read_csv('{}/results_merged/strain_data_v4_param2.csv'.format(shafer_root_dir))

    return strain_df

def blending_strains_analysis():
    strain_df = pd.read_csv('{}/results_merged/strain_data_v3_param2.csv'.format(shafer_root_dir))
    strain_df['source_strain'] = strain_df['source'] + '_' + strain_df['strain'].astype(str)
    # print(strain_df.head())

    # for sample in strain_df.sample.unique():
    for sample in env_samples:
        print(sample)

        df = strain_df[strain_df['sample'] == sample]
        all_strains = df['source_strain'].unique().tolist()
        loop_strains = list(filter(lambda k: 'loop' in k, all_strains))
        accungs_strains = list(filter(lambda k: 'accungs' in k, all_strains))

        df = df.merge(df, on='Pos')
        df = pd.crosstab(df['source_strain_x'], df['source_strain_y'])

        df = df.drop(loop_strains, axis=1) # drop columns
        df = df.drop(accungs_strains, axis=0) # drop rows
        # print(df)
        df.to_csv('/Users/omer/Documents/Lab/HIV_shafer/strain_intersect/{}_strain_intersect.csv'.format(sample))


if __name__ == '__main__':
    # summarize_associvar_strain_analysis()
    # plot_compare_associvar_accungs_strains()

    # plot_compare_loop_accungs_mutations()

    # blending_strains_analysis()

    # plot_compare_loop_accungs_mutations(generate_strain_data= False, generate_freq_data= False)

    # plot_compare_loop_accungs_mutations(run_name='v7_no_indels_for_seq_comparison', new_pipe_runs= True, add_strain_data= False, plot=False)

    plot_compare_loop_accungs_mutations(run_name='v8_new_for_indels', new_pipe_runs=True,
                                        generate_strain_data= False, plot=False)
