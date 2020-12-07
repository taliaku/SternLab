import pandas as pd

def create_chosen_mutation_from_freqs(freq_file, output_dir, run_name, freq_threshold = 0.001, prob_threshold = 0.85,
                                      coverage_threshold = 1000, range = (0, 100000)):

    freqs = pd.read_csv(freq_file, sep='\t')

    invalid_ref_positions = set(freqs[(freqs['Rank'] == 0) & (freqs['Ref'] != freqs['Base']) & (freqs['Ref'] != '-')]['Pos'].tolist())
    print('len(invalid_ref_positions): {}'.format(len(invalid_ref_positions)))

    chosen_mutations = freqs[(freqs['Rank'] != 0) # all minor variants
                           & (freqs['Freq'] > freq_threshold)
                           & (freqs['Prob'] > prob_threshold)
                           & (freqs['Ref'] != '-') # remove insertions
                           & ( ~(freqs['Pos'].isin(invalid_ref_positions)) )
                           & (freqs['Pos'] >= range[0]) & (freqs['Pos'] < range[1])
                           & (freqs['Read_count'] > coverage_threshold)]

    print(len(chosen_mutations))

    chosen_mutations['variant'] = chosen_mutations['Ref'] + chosen_mutations['Pos'].astype(str) + chosen_mutations['Base']

    # extract
    mutations_filename = 'chosen_mutations_{}'.format(run_name)
    chosen_mutations['variant'].to_csv('{}/{}.csv'.format(output_dir, mutations_filename), header= True, index=False)
    # chosen_mutations.to_csv('{}/{}_with_details.csv'.format(output_dir, mutations_filename), index=False)


def main_hiv_shafer_loop():
    loop_freq_files_root = '/Volumes/STERNADILABHOME$/volume1/shared/analysis/HIV_shafer/loop_genomics_pipeline/envs_output'
    output_root = '/Volumes/STERNADILABHOME$/volume1/shared/analysis/HIV_shafer/associvar'
    samples = ['env10', 'env11', 'env12', 'env13', 'env14', 'env15', 'env3', 'env4', 'env5', 'env6',
               'env7', 'env8', 'env9']
    for sample in samples:
        print('Sample: {}'.format(sample))
        freq_file = '{}/{}/pipeline_3/s{}.freqs'.format(loop_freq_files_root, sample, sample[3:])
        output_dir = '{}/{}'.format(output_root, sample)

        create_chosen_mutation_from_freqs(freq_file,
                                          output_dir,
                                          run_name = 'contigB_freq0.1',
                                          freq_threshold = 0.001,
                                          range = (1000, 3051))


if __name__ == '__main__':

    main_hiv_shafer_loop()