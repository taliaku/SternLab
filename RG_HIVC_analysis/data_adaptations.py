import glob
import os

import numpy as np
import pandas as pd

from RG_HIVC_analysis.constants import excluded_samples


def create_unified_samples_to_patient_and_dsi():
    extension = 'tsv'
    all_filenames = [i for i in glob.glob('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_input/tables_control/samples_*.{}'.format(extension))]

    # combine all files in the list
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])
    # export to csv
    combined_csv.to_csv("samples_to_patient_and_dsi_control.csv", index=False, encoding='utf-8-sig')


def generate_unified_filtered_verbose_freqs_df(min_read_count = 100, freq_threshold=1e-03):
    freqs_folder_path = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/control_low/'
    samples_data_file = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/output_tables/samples_to_patient_and_dsi_control.csv'
    freq_files_with_muts = glob.glob(f'{freqs_folder_path}/*.freqs')

    freq_dfs_with_id = []
    samples_to_patient_and_dates = pd.read_csv(samples_data_file, sep='\t')
    samples_to_patient_and_dates.set_index('id', inplace= True)

    for file in freq_files_with_muts:
        sample_id = os.path.splitext(os.path.basename(file))[0]

        # filtering samples
        if sample_id in excluded_samples:
            print('Excluded sample: {} - Skipping'.format(sample_id))
            continue
        print('Handling sample: ' + sample_id)
        freq_df = pd.read_csv(file, sep='\t')

        # filtering freq file
        # TODO- additional filters?
        freq_df = freq_df[freq_df["Read_count"] > min_read_count]
        freq_df['Freq'] = np.where(freq_df['Freq'] >= freq_threshold, freq_df['Freq'], 0)
        # freq_df = freq_df[freq_df['Mutation_type'] != 'consensus']

        # adding id & dsi columns
        patient_id = samples_to_patient_and_dates.loc[f'{sample_id}', 'patient']
        days_since_infection = samples_to_patient_and_dates.loc[f'{sample_id}', 'days since infection']
        freq_df['ind_id'] = patient_id
        freq_df['sample_id'] = sample_id
        freq_df['years_since_infection'] = str(np.round((days_since_infection) / float(365),2))

        # if patient_id == 29447:
        # print(freq_df)
        freq_dfs_with_id.append(freq_df)


    unified_freq_df_with_ids_ysi = pd.concat(freq_dfs_with_id)
    print(unified_freq_df_with_ids_ysi.shape)
    unified_freq_df_with_ids_ysi.to_csv( f'{freqs_folder_path}/unified_freqs_filtered_verbose.csv', index=False)

    # coverage checks
    # cov_filtered = unified_freq_df_with_ids_ysi[unified_freq_df_with_ids_ysi["Read_count"] > min_read_count]
    # pos_filtered = unified_freq_df_with_ids_ysi[unified_freq_df_with_ids_ysi["Pos"] > 3000]
    # pos_cov_filtered = unified_freq_df_with_ids_ysi[(unified_freq_df_with_ids_ysi["Pos"] > 3000) & (unified_freq_df_with_ids_ysi["Read_count"] > min_read_count)]
    # print(unified_freq_df_with_ids_ysi.shape)
    # print(cov_filtered.shape)
    # print(pos_filtered.shape)
    # print(pos_cov_filtered.shape)
    return unified_freq_df_with_ids_ysi


if __name__ == "__main__":
    # create_unified_samples_to_patient_and_dsi()
    generate_unified_filtered_verbose_freqs_df()
