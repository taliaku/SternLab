import glob
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

from FITS.create_fits_input import split_by_mutation
from freqs_utilities import change_ref_to_consensus
from pbs_runners import fits_runner


def generate_fits_input():
    run_folder = 'orig_high'

    # orig patients
    patients = ['12796', '13003', '15664', '16207', '17339', '19937', '22097', '22763', '22828', '23271', '26892','28545', '28841', '29219', '29447', '31254', '34253', '47939']
    # patients = ['13003', '15664', '16207', '22097', '22763', '22828', '26892', '29447', '31254', '47939']  # high_q
    # patients = ['15664', '16207', '22097', '22763', '22828', '29447', '31254', '47939']  # high_q
    # patients = ['22097', '22763', '22828', '23271', '26892','28545', '28841', '29219', '29447', '31254', '34253', '47939']

    # control patients
    # patients = ['24277', '6773', '26755', '4956', '7965', '8992', '22992', '1689', '13694', '4845', '15687', '8670', '14201', '324', '7878', '4179']

    # get freq data
    df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/{}/unified_freqs_filtered_verbose.csv'.format(run_folder))

    # filter indels
    df = df[(df['Base'] != '-') & (df['Ref'] != '-')]
    df = df.astype({"Pos": int})


    # filtering to selected patients
    df = df[df['ind_id'].isin(patients)]

    # run per patient
    for patient in df.ind_id.unique():
        print('Generating fits input files for patient: {}'.format(patient))
        df_patient = df[(df['ind_id'] == patient)]

        # filtering to syn positions
        pos_file_with_entropy= '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/mr_summary_files/syn_pos_by_ZN_with_entropy_filter/mutation_rate_positions_orig_high_v4_%s_0.3.txt' % patient
        pos_file_no_entropy= '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/mr_summary_files/syn_pos_by_ZN_no_entropy_filter/mutation_rate_positions_orig_high_v4_no_entropy_filter_%s_0.txt' % patient
        synonymous_selected_positions = get_syn_pos_from_file(pos_file_no_entropy)
        df_patient = df_patient[df_patient['Pos'].isin(synonymous_selected_positions)]

        # convert years_since_infection to Gen{0,1,2..}
        df_patient['Gen'] = df_patient['years_since_infection'].apply(lambda ysi: int((ysi*365)/2))

        # TODO- force founder consensus in all timepoints? needs to fit split_by_mutation() code
        # force ref=consensus in each timepoint (rank 0 == ref)
        dfs = []
        for timepoint in df_patient.Gen.unique():
            df_timepoint = df_patient[df_patient['Gen'] == timepoint]
            df_timepoint = change_ref_to_consensus(df_timepoint)
            dfs.append(df_timepoint)
        df_patient = pd.concat(dfs)

        # create FITS input files
        output_path = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/{}/fits/'.format(run_folder)

        ### per patient
        # split_by_mutation(df_patient, output_path, ID= patient)
        ### per position
        # for pos in df_patient.Pos.unique():
        #     print('pos: {}'.format(pos))
        #     df_pos = df_patient[(df_patient['Pos'] == pos)]
        #     split_by_mutation(df_pos, output_path + str(patient) + '/', ID='no_entropy_{}_{}'.format(patient, pos))

        ### Alternative
        # fits_input = freq_2_fits(df, filter='transition', out=r'/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/{}/fits_input_unified.csv')

        ### generate params file - HACK- is based on filtering above
        generate_fits_params_file(df_patient, patient)


def generate_fits_params_file(df_patient, patient):
    highest_gen = df_patient['Gen'].max()
    example_param = "/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/orig_high/fits/mr_params_26892_example.txt"
    patient_param = "/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/orig_high/fits/mr_params_{}.txt".format(
        patient)
    with open(example_param, "rt") as fin:
        with open(patient_param, "wt") as fout:
            for line in fin:
                fout.write(line.replace('$$$$$', str(highest_gen)))


def get_syn_pos_from_file(file):
    with open(file) as f:
        syn_pos = f.readlines()
    syn_pos = [int(x.strip()) for x in syn_pos]
    return syn_pos

def run_fits():
    input_files_orig_high=  '/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/fits/input_files/orig_high_qual/'

    ## per patient
    # patients_file = glob.glob(input_files_orig_high+'FITS_input_file*')
    # for file in patients_file:
    #     patient_id= int(file.split('_')[-3])
    #     params_file = input_files_orig_high + 'mr_params_{}.txt'.format(patient_id)
    #
    #     print(params_file)
    #     print(file)
    #     fits_runner(1, file, params_file, alias='FITS_{}'.format(patient_id),
    #                 posterior_file=file+'.posterior', summary_file=file+'.summary')

    ## per pos
    patients = ['12796', '13003', '15664', '16207', '17339', '19937', '22097', '22763', '22828', '23271', '26892','28545', '28841', '29219', '29447', '31254', '34253', '47939']
    # patients = ['12796']

    for p in patients:
        params_file = input_files_orig_high + 'mr_params_{}.txt'.format(p)

        # using pbs array per mut
        for mut in ['GA', 'AG', 'CT', 'TC']:
            input_file = input_files_orig_high + '%s/FITS_input_file_no_entropy_%s_$PBS_ARRAY_INDEX\\_%s.txt' % (p,p,mut)

            #filter Gen > 400
            df = pd.read_csv(input_file, sep='\t')
            df = df[df['Gen'] < 400]
            input_file = input_file + '.trunc'
            df.to_csv(input_file, sep='\t', encoding='utf-8', index=False)

            print(params_file)
            print(input_file)
            fits_runner(1, input_file, params_file,
                        alias='FITS_{}_{}'.format(p,mut),
                        posterior_file=input_file+'.posterior',
                        summary_file=input_file+'.summary',
                        batch= 9000)


def post_analysis():
    # entropy_position_files_dir = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/mr_summary_files/syn_pos_by_ZN_with_entropy_filter/'
    entropy_position_files_dir = '/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/fits/syn_pos_by_ZN_with_entropy_filter/'

    # fits_dir = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/orig_high/fits/try/'
    fits_dir = '/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/fits/input_files/orig_high_qual/'
    outputs_dir = fits_dir + 'results_analysis/'

    patients = ['12796', '13003', '15664', '16207', '17339', '19937', '22097', '22763', '22828', '23271', '26892','28545', '28841', '29219', '29447', '31254', '34253', '47939']
    # patients = ['26892']

    dfs=[]
    for patient in patients:
        patient_summary_files_dir = '%s/%s/' % (fits_dir, patient)
        patient_mr_summary = custom_summary_2_csv_biallelic(patient_summary_files_dir,
                                                     out= outputs_dir + 'fits_mr_summary_%s.csv' % patient,
                                                     patient=patient,
                                                     filter_entropy=True,
                                                     syn_position_files_dir= entropy_position_files_dir,
                                                     inverse_direction= False)
        dfs.append(patient_mr_summary)

    final = pd.concat(dfs)
    final.to_csv('%s/fits_mr_summary_orig_high.csv' % outputs_dir, index=False)

    # stats by patient & mut
    final = final[final['significance'] == 'significant']
    stats = final.groupby(['Patient', 'Mutation'])['MR'].agg(['median', 'mean'])
    stats = stats.pivot(index='Mutation', columns='Patient', values=['median', 'mean'])
    stats.to_csv('%s/fits_mr_summary_orig_high_stats.csv' % outputs_dir)
    print(stats)


    # plot distribution by patient & mut
    for p in final.Patient.unique():
        for mut in final.Mutation.unique():
            print(mut)
            ax = sns.distplot(final[(final.Mutation == mut) & (final.Patient == p)]['MR'])
            plot_header = mut
            ax.set_title(plot_header)
            # plt.show()
            plt.savefig(outputs_dir + str(p) + '_' + str(plot_header) + '.png')
            plt.cla()

def custom_summary_2_csv_biallelic(summary_dir, out=None, patient='', filter_entropy= True, syn_position_files_dir = None, inverse_direction = False):
    """
    this method creates a summary csv file for all results in summary dir
    :param summary_dir: a directory containing all summary files
    :param out: output file path
    :return: a data frame of all results
    """

    files = glob.glob(summary_dir + '/FITS*summary')

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

    # filtering to entropy positions
    if filter_entropy:
        pos_file_with_entropy = '%smutation_rate_positions_orig_high_v4_%s_0.3.txt' % (syn_position_files_dir, patient)
        synonymous_selected_positions = get_syn_pos_from_file(pos_file_with_entropy)
        final = final[final['Pos'].isin(synonymous_selected_positions)]

    if out!= None:
        final.to_csv(out, index=False)

    return final

if __name__ == "__main__":
    # generate_fits_input()
    run_fits()
    # post_analysis()

