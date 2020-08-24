import glob
from pathlib import Path

import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

from FITS.create_fits_input import split_by_mutation
from freqs_utilities import change_ref_to_consensus
from pbs_runners import fits_runner


def generate_fits_input():
    run_folder = 'orig_high'

    # orig patients
    # patients = ['12796', '13003', '15664', '16207', '17339', '19937', '22097', '22763', '22828', '23271', '26892','28545', '28841', '29219', '29447', '31254', '34253', '47939']
    # patients = ['22097', '22763', '22828', '23271', '26892','28545', '28841', '29219', '29447', '31254', '34253', '47939']
    # patients = ['13003', '15664', '16207', '22097', '22763', '22828', '26892', '29447', '31254', '47939']  # high_q
    # patients = ['15664', '16207', '22097', '22763', '22828', '29447', '31254', '47939']  # high_q

    # control patients
    # patients = ['24277', '6773', '26755', '4956', '7965', '8992', '22992', '1689', '13694', '4845', '15687', '8670', '14201', '324', '7878', '4179']

    # zn patients
    patients = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']
    patients = ['p7', 'p6', 'p5', 'p4', 'p3', 'p2', 'p1']

    # get freq data from unified freq
    # unified_freq_path = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/{}/unified_freqs_filtered_verbose.csv'.format(run_folder)
    unified_freq_path = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_rawdata/freq_files/unified_freqs_filtered_verbose.csv'
    df = pd.read_csv(unified_freq_path)

    a = df[(df['Ref'] == 'M')]
    # filter indels
    df = df[(df['Base'] != '-') & (df['Ref'] != '-')]
    df = df.astype({"Pos": int})

    # filtering to selected patients
    df = df[df['ind_id'].isin(patients)]

    # run per patient
    # for patient in ['15664']:
    for patient in df.ind_id.unique():
        print('Generating fits input files for patient: {}'.format(patient))
        df_patient = df[(df['ind_id'] == patient)]

        # filtering to syn positions
        # syn_pos_files_dir = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/MR/syn_pos_by_ZN_no_entropy_filter/'
        # syn_pos_file_path_per_patient = syn_pos_files_dir + ('mutation_rate_positions_orig_high_v5_no_entropy_filter_%s_0.txt') % patient
        syn_pos_file_path_per_patient = '/Users/omer/PycharmProjects/HIV_fitness_landscape/data/mutation_rates/mutation_rate_positions_0.01.txt'
        synonymous_positions = get_syn_pos_from_file(syn_pos_file_path_per_patient)
        df_patient = df_patient[df_patient['Pos'].isin(synonymous_positions)]

        # convert years_since_infection to Gen{0,1,2..}
        df_patient['Gen'] = df_patient['years_since_infection'].apply(lambda ysi: int((ysi*365)/2))

        # force ref=founder consensus in all timepoints
        df_patient_t = force_founder_con_as_ref_across_timepoints(df_patient)
        df_patient = df_patient_t

        # create FITS input files
        # output_path = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/fits_inputs/{}/'.format(run_folder)
        output_path = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_rawdata/fits_inputs/'
        patient_out_path = output_path + '{}/'.format(patient)
        Path(patient_out_path).mkdir(exist_ok=True, parents=False)

        # truncate Gen>350
        # gen_threshold = 350
        # if patient in [13003, 15664, 22097, 22763, 22828, 26892, 29447]:
        #     df_patient = df_patient[df_patient['Gen'] < gen_threshold]

        # generate params file (executed here to ensure same filtering as input file)
        generate_params_file = False
        if generate_params_file:
            generate_fits_params_file(df_patient, patient)

        ### per patient
        per_patient = False
        if per_patient:
            split_by_mutation(df_patient, output_path, ID= patient)

        ### per position
        # for pos in [1685]:
        for pos in tqdm(df_patient.Pos.unique()):
            print('pos: {}'.format(pos))
            df_pos = df_patient[(df_patient['Pos'] == pos)]
            split_by_mutation(df_pos, output_path + str(patient) + '/', ID='{}_{}'.format(patient, pos))

        ## alternative to for loop? TODO- examine this. or simply give it 2 hours to run.
        # incompletly_implemented_alternative(df, df_patient)


def incompletly_implemented_alternative(df, df_patient):
    mutations = ['AA', 'AG', 'GG', 'GA', 'CC', 'CT', 'TT', 'TC']
    for i in range(0, len(mutations), 2):
        print("Starting splitting mutation type {}....\n".format(mutations[i + 1]))
        mut_df = df[(df['Mutation'] == mutations[i]) | (df['Mutation'] == mutations[i + 1])]

        # change bases to numeric value 0 for wt and 1 for mutant
        mut_df['Base'] = 1 - (mut_df['Base'] == mut_df['Ref']).astype(int)
        mut_df['Ref'] = 0
        mut_df.reset_index(drop=True, inplace=True)

        # DEBUG
        # for p in range(0, len(mut_df['Base'].values), 2):
        #     if mut_df['Base'].values[p] != 0 or mut_df['Base'].values[p+1] != 1:
        #         print(mut_df.ix[p])

        count = Counter(mut_df['Base'].values)
        if count[0] != count[1]:
            raise Exception('Error: incompatible fits format. Base to index failed\n')

        # change the frequencies so the sum of each mutant wt will be 1
        # TODO- normalize per (pos+gen) using merge
        # mut_df = sum_2_one(mut_df)

        # skip if no data
        if len(mut_df) == 0:
            # print('**no data**')
            continue

        # remove unnecessary columns and change the order of the columns
        mut_df = mut_df.drop(['Prob', 'Rank', 'ind_id', 'sample_id', 'years_since_infection', 'Mutation'], axis=1)
        mut_df = mut_df[['Gen', 'Base', 'Freq', 'Pos']]
        mut_df = mut_df.rename(columns={'Gen': 'gen', 'Base': 'allele', 'Freq': 'freq', 'Pos': 'pos'})

        # save the result into a file
        # filename = 'FITS_input_file_{}_{}.txt'.format('_{}_{}'.format(patient, pos), mutations[i + 1])  # mutation type is at i + 1
        #
        # # print("Done parsing mutation type {}. Saving..\n".format(mutations[i+1]))
        # for pos in df_patient.Pos.unique():
        #     print('pos: {}'.format(pos))
        #     df_pos = df_patient[(df_patient['Pos'] == pos)]
        #     if len(mut_df) != 0:
        #         mut_df.to_csv(os.path.join(output_path + str(patient) + '/', filename), sep='\t', encoding='utf-8', index=False)


def force_founder_con_as_ref_across_timepoints(df_patient):
    # get first available reference for every position in patient
    founder_consensus = df_patient.sort_values(by=['Pos', 'Gen', 'Rank'])
    founder_consensus = founder_consensus.drop_duplicates("Pos")
    founder_consensus = founder_consensus[['Pos', 'Base']]

    df_patient_t = pd.merge(df_patient, founder_consensus, on='Pos', how='left', suffixes=('', '_r'))
    df_patient_t.Ref = df_patient_t.Base_r
    df_patient_t = df_patient_t.drop(columns=['Base_r'])

    # verifications-
    # 1. sanity- same length
    if len(df_patient) != len(df_patient_t):
        raise BaseException('row count should not change')
    # 2. verify identical ref in first timepoint in both df's
    df_founder_ref = df_patient[df_patient['Gen'] == 0]
    df_t_founder_ref = df_patient_t[df_patient_t['Gen'] == 0]
    merge_not_precise = df_founder_ref.set_index(df_founder_ref['Pos']).join(df_t_founder_ref.set_index(df_t_founder_ref['Pos']),
                                                                             rsuffix='_r')
    founder_ref_corruption = merge_not_precise[merge_not_precise['Ref'] != merge_not_precise['Ref_r']]
    if len(founder_ref_corruption) != 0:
        raise BaseException('original founder ref changed')

    # 3. verify identical refs in all timepoints in new df
    tmp = (df_patient_t[['Pos', 'Ref']]).drop_duplicates()
    tmp2 = tmp.groupby(['Pos']).count()['Ref']
    positions_with_discrapencies = tmp2[tmp2 != 1]
    if len(positions_with_discrapencies) != 0:
        raise BaseException('ref not aligned across all timepoints')

    return df_patient_t


def force_con_as_ref_per_timepoint(df_patient):
    dfs = []
    for timepoint in df_patient.Gen.unique():
        df_timepoint = df_patient[df_patient['Gen'] == timepoint]
        # TODO- remove indels before change_ref_to_consensus()
        df_timepoint = change_ref_to_consensus(df_timepoint)
        dfs.append(df_timepoint)
    df_patient = pd.concat(dfs)
    return df_patient


def generate_fits_params_file(df_patient, p):
    highest_gen = df_patient['Gen'].max()
    example_param = "/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_rawdata/fits_inputs/mr_params_p1_example.txt"
    patient_param = "/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_rawdata/fits_inputs/mr_params_{}.txt".format(p)
    with open(example_param, "rt") as fin:
        with open(patient_param, "wt") as fout:
            for line in fin:
                fout.write(line.replace('$$$$$', str(highest_gen)))


def get_syn_pos_from_file(file):
    with open(file) as f:
        syn_pos = f.readlines()
    syn_pos = [int(x.strip()) for x in syn_pos]
    return syn_pos


###############################


def run_fits():
    input_files_orig_high=  '/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/fits/input_files/orig_high/'

    ## per patient
    # patients_file = glob.glob(input_files_orig_high+'FITS_input_file*')
    # for file in patients_file:
    #     patient_id= int(file.split('_')[-3])
    #     params_file = input_files_orig_high + 'mr_params_{}.txt'.format(patient_id)
    #
    #     print(params_file)
    #     print(file)
    #     fits_runner(1, file, params_file, job_name='FITS_{}'.format(patient_id),
    #                 posterior_file=file+'.posterior', summary_file=file+'.summary')

    ## per pos
    # patients = ['12796', '13003', '15664', '16207', '17339', '19937', '22097', '22763', '22828', '23271', '26892','28545', '28841', '29219', '29447', '31254', '34253', '47939']
    patients = ['47939', '34253', '31254', '29447', '29219', '28841', '28545', '26892', '23271', '22828', '22763'] # opposite order
    # patients = ['13003', '15664', '22097', '22763', '22828', '26892', '29447'] # relevant after truncation
    # patients = ['12796']

    for p in patients:
        # params_file = input_files_orig_high + 'mr_params_{}.txt'.format(p)
        params_file = input_files_orig_high + 'mr_params_bottleneck/mr_params_{}.txt'.format(p)

        # using pbs array per mut
        for mut in ['GA', 'AG', 'CT', 'TC']:
            input_file = input_files_orig_high + '%s/FITS_input_file_%s_$PBS_ARRAY_INDEX\\_%s.txt' % (p,p,mut)

            # use truncated files
            # input_file = input_file + '.trunc'

            print(params_file)
            print(input_file)
            fits_runner(1, input_file, params_file,
                        alias='FITS_{}_{}'.format(p,mut),
                        posterior_file=input_file+'.posterior.bottleneck',
                        summary_file=input_file+'.summary.bottleneck',
                        queue= 'short',
                        batch= 9062)


###############################


def custom_summary_2_csv_biallelic(summary_dir, filenames_pattern = 'FITS*.summary', out=None, patient='', inverse_direction = False):
    """
    this method creates a summary csv file for all results in summary dir
    :param summary_dir: a directory containing all summary files
    :param out: output file path
    :return: a data frame of all results
    """

    files = glob.glob(summary_dir + '/' + filenames_pattern)

    dfs = []
    for f in tqdm(files):
        # print(str(f))
        pos = int(f.split('_')[-2])
        mt = f.split('_')[-1].split('.')[0]

        with open(f, 'r') as o:
            lines = o.readlines()
        if inverse_direction:
            line = [l for l in lines if '1     0     ' in l][0]
        else:
            line = [l for l in lines if '0     1     ' in l][0]

        rate = line.split()[2]
        if '*' in rate:
            rate = float(rate.split('*')[-1])
            significance = 'non-significant'
        else:
            rate = float(rate)
            significance = 'significant'

        df = pd.DataFrame({'Patient':patient, 'Pos':pos, 'Mutation':mt, 'MR':rate, 'significance':significance}, index=[0])

        dfs.append(df)

    final = pd.concat(dfs)
    final = final.sort_values(by=['Patient', 'Pos', 'Mutation'])

    if out!= None:
        final.to_csv(out, index=False)

    return final

def fits_outputs_to_csv():

    # run_name = 'orig_high'
    run_name = 'zn'

    fits_dir = '/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/fits/input_files/%s/' % run_name
    outputs_dir = fits_dir + 'results_analysis/'
    Path(outputs_dir).mkdir(exist_ok=True, parents=False)

    # patients = ['12796', '13003', '15664', '16207', '17339', '19937', '22097', '22763', '22828', '23271', '26892', '28545', '28841', '29219', '29447', '31254', '34253', '47939']
    # patients = ['13003', '15664', '22097', '22763', '22828', '26892', '29447'] # relevant after truncation
    # patients = ['26892']
    patients = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']
    # patients = ['12796', '13003', '15664', '16207', '17339', '19937', '22097', '22763', '22828', '47939'] # relevant for bottleneck

    filenames_pattern = 'FITS*.summary'

    ### fits outputs to CSV
    dfs = []
    for patient in patients:
        patient_summary_files_dir = fits_dir + '{}/'.format(patient)
        patient_mr_summary = custom_summary_2_csv_biallelic(patient_summary_files_dir,
                                                            out=outputs_dir + 'fits_mr_summary_{}.csv'.format(patient),
                                                            filenames_pattern=filenames_pattern,
                                                            patient=patient)
        dfs.append(patient_mr_summary)
    final = pd.concat(dfs)
    final.to_csv(outputs_dir + 'fits_mr_summary_' + run_name + '.csv', index=False)


###############################


def analysis():
    run_name = 'orig_high'
    # run_name = 'orig_high-bottleneck'
    # run_name = 'zn'
    fits_summaries_dir = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/MR/fits_runs_summaries/%s/' % run_name
    fits_summary_unified = pd.read_csv(fits_summaries_dir + 'fits_mr_summary_{}.csv'.format(run_name))

    ## stats by patient & mut
    # significant only
    fits_summary_unified = fits_summary_unified[fits_summary_unified['significance'] == 'significant']

    # optional filtering low entropy positions
    filter_low_entropy = True

    if filter_low_entropy:
        orig_high = True

        if orig_high:
            orig_high_pos_after_entropy = []
            for patient in fits_summary_unified.Patient.unique():
                syn_entropy_pos_per_patient = get_syn_pos_from_file('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/MR/syn_pos_by_ZN_with_entropy_filter/mutation_rate_positions_orig_high_v4_%s_0.3.txt' % patient)
                orig_high_pos_after_entropy = orig_high_pos_after_entropy + syn_entropy_pos_per_patient

            orig_high_pos_after_entropy = list(set(orig_high_pos_after_entropy))
            orig_high_pos_after_entropy.sort()
            print(len(orig_high_pos_after_entropy))
            pos_after_entropy = orig_high_pos_after_entropy
        else:
            zn_pos_file_with_entropy = '/Users/omer/PycharmProjects/HIV_fitness_landscape/data/mutation_rates/mutation_rate_positions_0.3.txt'
            pos_after_entropy = get_syn_pos_from_file(zn_pos_file_with_entropy)

        print(len(fits_summary_unified))
        fits_summary_unified = fits_summary_unified[fits_summary_unified['Pos'].isin(pos_after_entropy)]
        print(len(fits_summary_unified))

    stats = fits_summary_unified.groupby(['Patient', 'Mutation'])['MR'].agg(['mean'])

    stats_table_path = fits_summaries_dir + 'fits_mr_summary_{}_stats_entropy.csv'.format(run_name)
    stats.to_csv(stats_table_path)
    stats = pd.read_csv(stats_table_path)  # TODO- replace with unstack or whatever needed
    stats = stats.pivot(index='Mutation', columns='Patient', values=['mean'])
    stats = stats.T
    stats['Avg'] = stats.mean(numeric_only=True, axis=1)
    stats = stats.sort_values(by=['Avg'], ascending=False)
    stats.to_csv(stats_table_path)  # TODO- replace with unstack or whatever needed
    print(stats)

    # plot distribution by patient & mut
    dist_plots = True
    if dist_plots:
        # v1
        # for p in fits_summary_unified.Patient.unique():
        #     print(p)
        #     for mut in fits_summary_unified.Mutation.unique():
        #         print(mut)
        #         ax = sns.distplot(fits_summary_unified[(fits_summary_unified.Mutation == mut) & (fits_summary_unified.Patient == p)]['MR'])
        #         plot_header = mut
        #         ax.set_title(plot_header)
        #         ax.set_xscale('log')
        #         # plt.show()
        #         plt.savefig(fits_summaries_dir + str(p) + '_' + str(plot_header) + '.png')
        #         plt.cla()
        # v2
        g = sns.FacetGrid(fits_summary_unified, row="Patient", col="Mutation", margin_titles=True)
        g.set(xscale="log")
        bins = np.linspace(0, 5e-3, 20)
        g.map(plt.hist, 'MR', bins=bins)

        # plt.show()
        mr_dist_plot_path = fits_summaries_dir + 'mr_dist_plot_{}_entropy.png'.format(run_name)
        plt.savefig(mr_dist_plot_path)



def get_gens():
    run_folder = 'orig_high'
    df = pd.read_csv(
        '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/{}/unified_freqs_filtered_verbose.csv'.format(
            run_folder))

    df['Gen'] = df['years_since_infection'].apply(lambda ysi: int((ysi * 365) / 2))
    df = df[df['Gen']> 350]
    df2 = df.groupby(['ind_id', 'sample_id'])['Gen'].agg(['max', 'min'])
    print(df2)


if __name__ == "__main__":
    # generate_fits_input()
    # run_fits()
    # fits_outputs_to_csv()
    analysis()


    ## Hacks etc.
    # (trunc quick analysis)
    # post_analysis_reduced()
    # get_gens()

    # outputs_dir = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/MR/SUMMARY/fits_runs_summaries/trunc/'
    # stats = pd.read_csv('%s/fits_mr_summary_orig_high_stats.csv' % outputs_dir)  # TODO- replace with unstack or whatever needed
    # stats = stats.pivot(index='Mutation', columns='Patient', values=['median', 'mean'])
    # stats.to_csv('%s/fits_mr_summary_orig_high_stats.csv' % outputs_dir)  # TODO- replace with unstack or whatever needed


