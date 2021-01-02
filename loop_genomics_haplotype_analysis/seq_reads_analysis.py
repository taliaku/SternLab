import glob

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def visualize_read_mapping():
    dual_reads_all = []
    for i in range(3,16):
        try:
            # get all reads
            sample = 'env{}'.format(i)
            print('sample: {}'.format(sample))
            blast_file = '/Volumes/STERNADILABHOME$/volume2/ita/haplotypes/envs_output/{}/joined.blast'.format(sample)
            blast_df = pd.read_csv(blast_file,
                                   sep='\t',
                                   names=['read_id', 'start_pos_ref', 'end_pos_ref', 'start_pos_read', 'end_pos_read',
                                          'direction', 'match_length', 'mutation_string'])
            print('plotting..')
            # length distribution
            # p = blast_df.hist(bins=100, column=['length'])
            # plt.show()

            # get DUAL reads
            a = blast_df.groupby('read_id').count()
            b = a[a['direction'] > 2]
            c = blast_df[blast_df['read_id'].isin(b.index.to_list())]
            c.to_csv('/Volumes/STERNADILABHOME$/volume1/shared/analysis/HIV_shafer/accungs_haplotype_inference/dual_reads_analysis/{}.blast.filtered.dual.csv'.format(sample), index=False)
            print('DUAL ratio: {}'.format(len(b)/len(a)*100))

            # distributions
            p = c.hist(bins=100, column=['start_pos_ref'])
            plt.show()

            # visualize mapping?
            # how?

            # get SHORT reads
            # start_pos distribution
            # visualize mapping
            c['sample'] = sample
            dual_reads_all.append(c)
        except:
            print('{} not working'.format(i))

    pd.concat(dual_reads_all).to_csv('/Volumes/STERNADILABHOME$/volume1/shared/analysis/HIV_shafer/accungs_haplotype_inference/dual_reads_analysis/all.blast.filtered.dual.csv', index=False)


def read_distribution():
    # get tables
    loop_metadata_files = glob.glob('/Users/omer/Documents/Lab/HIV_shafer/loop_genomics_seq_results/262*_Agentek_Stern_4_amplicon_*/262*_output/262*_sample_s*/262*_sample_s*_stats_trimmed.csv')
    print(len(loop_metadata_files))

    # to df
    reads_metadata_dfs = [pd.read_csv(s) for s in loop_metadata_files]
    reads_metadata = pd.concat(reads_metadata_dfs)
    reads_metadata['sample'] = reads_metadata['molecule_id'].apply(lambda x: x.split('_')[2])
    reads_metadata['sample_idx'] = reads_metadata['sample'].apply(lambda x: int(x[1:]))
    reads_metadata.sort_values(by='sample_idx', inplace=True)

    validation = reads_metadata.groupby('sample').count()
    print(validation)

    # hist
    # reads_metadata.groupby('sample')['contig1_len'].hist(bins=50, alpha=0.8, legend=True)
    reads_metadata.hist(column='contig1_len',
                        by= 'sample_idx',
                        figsize=(20,10),
                        sharex= True,
                        bins=50)
    # plt.savefig(fname='{}/results_merged/accungs_loop_freq_dist_agg.pdf'.format(shafer_root_dir))
    plt.show()

if __name__ == '__main__':
    # visualize_read_mapping()
    read_distribution()