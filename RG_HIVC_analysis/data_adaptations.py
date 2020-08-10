import glob
import os

import numpy as np
import pandas as pd
from Bio.Seq import Seq

from RG_HIVC_analysis import constants
from RG_HIVC_analysis.constants import orig_excluded_samples, get_ET86_region
from freqs_utilities import change_ref_to_consensus


def create_unified_samples_to_patient_and_dsi():
    extension = 'tsv'
    all_filenames = [i for i in glob.glob('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_input/tables_control/samples_*.{}'.format(extension))]

    # combine all files in the list
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])
    # combined_csv['timepoint'] = combined_csv.groupby('patient').cumcount() # TODO- handle column seperation to make this work
    # export to csv
    combined_csv.to_csv("samples_to_patient_and_dsi_control_timepoint.csv", index=False, encoding='utf-8-sig')

def add_codon_aa_info(freq_df):
    df_next_base = freq_df[freq_df['Rank'] == 0][['Pos', 'Base']]
    df_next_base['Pos'] = df_next_base['Pos'] - 1
    new_freq_df = freq_df.merge(df_next_base, on='Pos', how='left', sort=False,suffixes=('', '_next'))
    df_next2 = freq_df[freq_df['Rank'] == 0][['Pos', 'Base']]
    df_next2['Pos'] = df_next2['Pos'] - 2
    new_freq_df = new_freq_df.merge(df_next2, on='Pos', how='left', sort=False,suffixes=('', '_next2'))

    new_freq_df["codon_next"] = new_freq_df["Ref"] + new_freq_df["Base_next"] + new_freq_df["Base_next2"] # TODO- Ref or Base?
    # new_freq_df["codon_mid"] = new_freq_df["context"] # TODO
    # new_freq_df["codon_prev"] = new_freq_df["prev2bases"] + new_freq_df["consensus"] # TODO

    # Set reading frame and translation for HIV-pol: Locate TGG as "start codon" after position 3000
    data_hivpol = new_freq_df[new_freq_df["virus"] == "HIV-pol"]
    data_nonhivpol = new_freq_df[~(new_freq_df["virus"] == "HIV-pol")]

    start_codons = data_hivpol[(data_hivpol["codon_next"] == "TGG") & (data_hivpol["Pos"] >= 2990)].groupby("control")[
        "Pos"].aggregate(min)
    data_hivpol = data_hivpol.join(start_codons, on=["control"], how="left", rsuffix='_start_codon')
    data_hivpol = data_hivpol[(data_hivpol["Pos"] >= data_hivpol["Pos_start_codon"]) & (data_hivpol["Pos"] <= 8000)]
    data_hivpol["frame"] = (data_hivpol["Pos"] - data_hivpol["Pos_start_codon"]) % 3

    data_hivpol["consensus_codon"] = np.where(data_hivpol["frame"] == 0, data_hivpol["codon_next"],
                                              np.where(data_hivpol["frame"] == 1, data_hivpol["codon_mid"],
                                                       data_hivpol["codon_prev"]))
    data_hivpol["mutated_codon"] = np.where(data_hivpol["frame"] == 0, data_hivpol["Base"] + data_hivpol["next2bases"],
                                            np.where(data_hivpol["frame"] == 1,
                                                     data_hivpol["prevBase"] + data_hivpol["Base"] + data_hivpol[
                                                         "nextBase"], data_hivpol["prev2bases"] + data_hivpol["Base"]))

    data_hivpol["consensus_aa"] = data_hivpol["consensus_codon"].apply(
        lambda x: Seq(x).translate()[0] if "-" not in x else "-")
    data_hivpol["mutated_aa"] = data_hivpol["mutated_codon"].apply(
        lambda x: Seq(x).translate()[0] if "-" not in x else "-")
    data_hivpol["mutated_type"] = np.where(data_hivpol["mutated_aa"] == "*", "STOP",
                                           np.where(data_hivpol["mutated_aa"] == data_hivpol["consensus_aa"], "SYN",
                                                    np.where(data_hivpol["mutated_aa"] == "-", "INDEL", "NONSYN")))
    data_hivpol = data_hivpol[data_hivpol["Pos"] < 5020]  # ignore everything from the overlap until the end
    data_hivpol["orf"] = "pol"

    new_freq_df = pd.concat([data_hivpol, data_nonhivpol])
    return new_freq_df


def generate_unified_filtered_verbose_freqs_df(min_read_count = constants.coverage_threshold, freq_threshold=constants.freq_threshold, prob_threshold=constants.prob_threshold):
    freqs_folder_path = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/orig_high/'
    # freqs_folder_path = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_rawdata/freq_files/'

    samples_info_file = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/info_tables/samples_to_patient_and_dsi_timepoint.tsv'
    samples_to_patient_and_dates = pd.read_csv(samples_info_file, sep='\t')
    samples_to_patient_and_dates.set_index('id', inplace= True)

    freq_files_with_muts = glob.glob(f'{freqs_folder_path}/*.freqs')
    freq_dfs_with_id = []
    for file in freq_files_with_muts:
        sample_id = os.path.splitext(os.path.basename(file))[0]

        # skip excluded samples
        # TODO- is this needed if I filter by position?
        if sample_id in orig_excluded_samples:
            print('Excluded sample: {} - Skipping'.format(sample_id))
            continue

        # get df
        print('Handling sample: ' + sample_id)
        freq_df = pd.read_csv(file, sep='\t')

        # filter freq file
        tmp_df = freq_df[freq_df["Read_count"] > min_read_count]
        print('Coverage filter: {}% of current freqs file'.format((1-len(tmp_df)/len(freq_df))*100))
        freq_df = tmp_df
        freq_df['Freq'] = np.where(freq_df['Freq'] >= freq_threshold, freq_df['Freq'], 0)

        # TODO- verify prob threshold
        # tmp_df = freq_df[freq_df["Prob"] >= prob_threshold]
        # print('Filtered Prob: {}% of current freqs file'.format(len(tmp_df) * 100.0 / len(freq_df)))
        # freq_df = tmp_df
        # Optional- filter Read_count*Freq = 1 (only ~13 rows per freq file)

        # filter indels
        before_indel_filter_len = len(freq_df)
        print('insertion ratio: {}%'.format((len(freq_df[(freq_df['Ref'] == '-')])/len(freq_df)) * 100))
        # print('deletions ratio: {}%'.format((len(freq_df[(freq_df['Base'] == '-')])/len(freq_df)) * 100))

        # TODO- notice- might make freqs per positions not sum to 1
        freq_df = freq_df[(freq_df['Base'] != '-') & (freq_df['Ref'] != '-')]
        freq_df = freq_df.astype({"Pos": int})

        after_indel_filter_len = len(freq_df)
        removed_indels_ratio = (1- (after_indel_filter_len/before_indel_filter_len))
        print('removed_indels_ratio: {}%\n'.format(removed_indels_ratio * 100))
        if removed_indels_ratio > 0.5:
            raise BaseException('too many changes')


        # force ref=consensus in each sample (rank 0 == ref)
        freq_df = change_ref_to_consensus(freq_df)


        # add patient & sample data
        # Could be done much simpler by merge. But easier for ZN data this way
        patient_id = samples_to_patient_and_dates.loc[f'{sample_id}', 'patient']
        days_since_infection = samples_to_patient_and_dates.loc[f'{sample_id}', 'days_since_infection']
        timepoint = samples_to_patient_and_dates.loc[f'{sample_id}', 'timepoint']

        # for zn data: (should also be done for mine)
        # patient_id = sample_id.split("_")[0]
        # days_since_infection = int(sample_id.split("_")[1])

        freq_df['ind_id'] = patient_id
        freq_df['sample_id'] = sample_id
        freq_df['years_since_infection'] = str(np.round((days_since_infection) / float(365),2))
        freq_df['timepoint'] = timepoint

        # add annotation info
        freq_df['region'] = freq_df.apply(lambda row: get_ET86_region(row['Pos']), axis=1)

        # adding AA & codon info - TODO
        # freq_df = add_codon_aa_info(freq_df)

        freq_dfs_with_id.append(freq_df)


    unified_freq_with_additional_fields = pd.concat(freq_dfs_with_id)
    print(unified_freq_with_additional_fields.shape)
    unified_freq_with_additional_fields = unified_freq_with_additional_fields.sort_values(by=['ind_id', 'years_since_infection', 'Pos', 'Rank']) # sorting by Pos, Rank before YSI can be relevant sometimes
    unified_freq_with_additional_fields.to_csv(f'{freqs_folder_path}/unified_freqs_filtered_verbose.csv', index=False)



if __name__ == "__main__":
    # create_unified_samples_to_patient_and_dsi()
    generate_unified_filtered_verbose_freqs_df()
