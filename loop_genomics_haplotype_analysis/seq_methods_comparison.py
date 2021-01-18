import glob
from pathlib import Path

import pandas as pd

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 3000)
import seaborn as sns
# sns.set_context("poster")
import matplotlib.pyplot as plt

def run_new_pipeline_on_server():
    import sys
    ITAMAR_REPO_PATH = '/sternadi/home/volume2/ita/sternlab-public'
    sys.path.insert(1, ITAMAR_REPO_PATH)
    from pbs_multi_runner import multi_runner

    ref_path = '/sternadi/home/volume2/ita/haplotypes/references/HXB2_sequence_env_5853-8904.fasta'

    input_path_accungs = '/sternadi/datasets/volume2/MY060520202/FASTQ_Generation_2020-05-08_19_30_42Z-246623382/env{}_*'
    # output_path_accungs = '/sternadi/home/volume1/shared/analysis/HIV_shafer/new_pipeline_accungs_v2/env{}'
    output_path_accungs = '/sternadi/nobackup/volume1/HIVB_shafer_accungs_seq/new_pipeline_accungs_v3_no_indel/env{}'

    # q10
    # input_path_loop = '/sternadi/datasets/volume2/HIV_shafer_loop_genomics/fastq_only/s{}'
    # output_path_loop = '/sternadi/home/volume1/shared/analysis/HIV_shafer/new_pipeline_loop_v4_q10_no_indel/env{}'
    # long reads
    # input_path_loop = '/sternadi/datasets/volume2/HIV_shafer_loop_genomics/fastq_long_only/s{}'
    # output_path_loop = '/sternadi/nobackup/volume1/HIVB_shafer_accungs_seq/new_pipeline_loop_v5_long_reads_no_indel/env{}'
    # blast_mapped_2500 + q10
    input_path_loop = '/sternadi/datasets/volume2/HIV_shafer_loop_genomics/fastq_mapped_2500_only/s{}'
    output_path_loop = '/sternadi/nobackup/volume1/HIVB_shafer_accungs_seq/new_pipeline_loop_v6_blast_2500/env{}'

    # create param list
    accungs_run_param_list = []
    loop_run_param_list = []

    for sample_id in range(3,16):
        print('sample: env{}'.format(sample_id))
        # accungs
        input_dir_glob = input_path_accungs.format(sample_id)
        input_dir = glob.glob(input_dir_glob)[0]
        output_dir = output_path_accungs.format(sample_id)
        Path(input_dir).mkdir(parents=True, exist_ok=True)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        print(input_dir)
        print(output_dir)

        sample_dict_accungs = {'input_dir': input_dir,
                                'output_dir': output_dir,
                                'reference_file': ref_path,
                                # "with_indels": "N", #TODO- notice
                                "overlapping_reads": "Y",
                                "calculate_haplotypes": "N", # TODO- run Y on nobackup, when there's space
                                "max_basecall_iterations": 10}
        accungs_run_param_list.append(sample_dict_accungs)

        # loop
        input_dir = input_path_loop.format(sample_id)
        output_dir = output_path_loop.format(sample_id)
        Path(input_dir).mkdir(parents=True, exist_ok=True)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        print(input_dir)
        print(output_dir)

        sample_dict_loop = {'input_dir': input_dir,
                            'output_dir': output_dir,
                            'reference_file': ref_path,
                            "with_indels": "Y",
                            "overlapping_reads": "N",
                            "quality_threshold": "10",
                            "max_basecall_iterations": 10}
        loop_run_param_list.append(sample_dict_loop)

    print(len(accungs_run_param_list))
    print(len(loop_run_param_list))

    # run
    # multi_runner(accungs_run_param_list)
    multi_runner(loop_run_param_list)


# Note from sheri- dont take these results too seriously
def accungs_loop_muts_correlation(run_name,
                                  accungs_source_name = 'accungs',
                                  loop_source_name = 'loop',
                                  plot_correlation = True,
                                  freq_boxplot = False,
                                  freq_dist_all_samples = False):

    from loop_genomics_haplotype_analysis.summarize_visualize_strain_results import env_samples, shafer_root_dir, shafer_dir_nobackup
    from loop_genomics_haplotype_analysis.mutation_analysis import get_freqs

    print('* accungs_loop_muts_correlation *')
    # get filtered freqs
    freqs = get_freqs(freq_filename= 'freqs_strain_data_all_samples_v7_no_indels_for_seq_comparison.csv', remove_strain_data= False) # original dataset doesnt include strains

    # 1. plot per-mut correlation
    if plot_correlation:
        # print(freqs['source'].unique().tolist())
        freqs_accungs = freqs[freqs['source'] == accungs_source_name]
        freqs_loop = freqs[freqs['source'] == loop_source_name]

        # align positions from 2 sources
        align_freqs_positions_according_to_consensus_alignment()

        print(freqs_accungs.groupby('sample').count()['Pos'])
        print(freqs_loop.groupby('sample').count()['Pos'])
        cross_freqs = pd.merge(freqs_accungs,
                               freqs_loop,
                               how='inner',
                               suffixes= ('_accungs','_loop'),
                               on=['sample', 'Pos', 'Rank'],
                               sort= True
                            )

        print(cross_freqs.groupby('sample').count()['Pos'])

        # run on high & low freqs separatly
        freq_threshold = 0.05
        cross_freqs_high = cross_freqs[(cross_freqs['Freq_accungs'] >= freq_threshold) & (cross_freqs['Freq_loop'] >= freq_threshold)]
        cross_freqs_low = cross_freqs[(cross_freqs['Freq_accungs'] < freq_threshold) & (cross_freqs['Freq_loop'] < freq_threshold)]
        print('len(cross_freqs_high): {}'.format(len(cross_freqs_high)))
        print('len(cross_freqs_low): {}'.format(len(cross_freqs_low)))

        for name, freq_set, axis in [('high',cross_freqs_high, (1e-02,0.7)), ('low',cross_freqs_low, (1e-04,0.07))]:
            print(name)
            g = sns.lmplot(x="Freq_accungs",
                       y="Freq_loop",
                       col="sample",
                       col_order=env_samples,
                       col_wrap=5,
                       height=9,
                       data=freq_set);

            plt.xlim(axis[0], axis[1])
            plt.ylim(axis[0], axis[1])
            # plt.xscale('log')
            # plt.yscale('log')

            # g.set_title(name)
            # add diagonal
            for ax in g.axes.flat:
                ax.plot((0, 1), (0, 1), ls='--', color='0.7')

            plt.savefig(fname='{}/results_merged/accungs_loop_mutation_correlation_linear_{}_{}.pdf'.format(shafer_root_dir, name, run_name))
            plt.show()

    # 2. boxplot freq distribution
    if freq_boxplot:
        print('boxplot')
        #  1. regular
        # v1
        boxplot(freqs, '{}/results_merged/accungs_loop_freq_medians.pdf')
        # v2
        # mutations_box_plots_itamar(freqs, '', True)

        #  2. non-ox too (CA\GT)
        print(len(freqs))
        freqs_nonox = freqs[~freqs['mutation'].isin(['C->A', 'G->T'])]
        print(len(freqs_nonox))
        boxplot(freqs, '{}/results_merged/accungs_loop_freq_medians_nonox.pdf')
        #  3. high\low freqs?
        #  4. syn only - requires HXB2 coords

    # 3. mut-type count avg. (per sample, per seq-type)
    # I did count per strand, can also do per seq type (accungs is more but we'll see how much)
    # TODO- find correct agg method\column to agg on
    freqs_transitions = freqs[freqs['mutation'].isin(['G->A', 'C->T', 'A->G', 'T->C'])]
    print(freqs_transitions.groupby(['source', 'mutation']).agg(['count', 'mean'])['Freq'])
    print(freqs_transitions.groupby(['sample', 'mutation']).agg(['count', 'mean'])['Freq'])

    # 4. mut-freq distribution, per seq-type
    if freq_dist_all_samples:
        freqs.groupby('source')['Freq'].hist(bins=50, log=True, alpha=0.8, legend = True)
        plt.savefig(fname='{}/results_merged/accungs_loop_freq_dist_agg.pdf'.format(shafer_root_dir))
        plt.show()


def align_freqs_positions_according_to_consensus_alignment():
    pass
    # TODO- imp (start with one sample)
    # for sample in ['env12']:
    #     sample_all_sources = freqs[freqs['sample'] == sample]
    #
    #     # produce cons alignment?
    #     # TODO- alignment headache
    #     #  option?-
    #     pairwise2.align.localms(seq, target,10,-5,-10,-3)
    #
    #     # get cons alignment
    #     # TODO- fasta_to_df small-headache
    #     two_cons_alignment = fasta_to_df(
    #         '{}/env_consensuses_alignments/{}_accungs_loop_long.fasta'.format(shafer_dir_nobackup, sample))
    #     two_cons_alignment = two_cons_alignment.rename(columns={"Pos": "Aligned_pos", "Ref": "Rhino_ref"})
    #
    #     # add synched pos to df
    #     sample_accungs = two_cons_alignment[two_cons_alignment.Rhino_ref != "-"]
    #     max_pos = len(sample_accungs.Aligned_pos)
    #     pos_range = list(range(1, max_pos + 1))
    #     sample_accungs["Rhino_pos"] = pos_range


def boxplot(freqs_filtered, title):

    sns.boxplot(x='sample',
                y='Freq',
                hue='source',
                showmeans=True,
                meanprops={"markersize": 10, "marker": "o", "markerfacecolor": "red",
                           "markeredgecolor": "red"},
                data=freqs_filtered,
                )
    plt.ylim(1e-04, 1)
    plt.yscale('log')
    plt.savefig(fname=title.format(shafer_root_dir))
    plt.show()


# def mutations_box_plots_itamar(data, save_fig_path=None):
#     fig, axes = plt.subplots(figsize=(60, 15), ncols=5, nrows=1)
#     # plt.suptitle(title, fontsize=32)
#     plt.subplots_adjust(wspace=0.4, hspace=0.3)
#     daplot = sns.boxplot(data=data, x='sample', y='Freq', hue='source', ax=axes[i], showmeans=True,
#                          meanprops={"markersize":10, "marker":"o" ,"markerfacecolor":"red", "markeredgecolor":"red"})
#     daplot.set_yscale('log')
#     if save_fig_path:
#         fig.savefig(save_fig_path)
#     return daplot


# TODO- not good for me- concats seqs in fasta file
def fasta_to_df(fasta_file):
    with open(fasta_file, "r") as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                continue
            line = line.rstrip("\n")
            seq += line

    pos_dict = {}
    for i in range(1,len(seq)):
        # TODO- can be done without loop
        pos_dict[i] = seq[i-1]

    output_df = pd.DataFrame(pos_dict.items(), columns=["Pos", "Ref"])
    return output_df

def fasta_to_df2(tmp_cirseq_dir):
    fasta_files = glob.glob(tmp_cirseq_dir + "*.fasta")
    records = {}
    for file in fasta_files:
        records = parse_fasta(file)
        df = pd.DataFrame.from_dict(records, orient='index')
        df.index.name = 'id'
        df.columns = ['seq']
        return df


if __name__ == '__main__':
    run_new_pipeline_on_server()

    # accungs_loop_muts_correlation(run_name='v2_rerun', accungs_source_name= 'accungs', loop_source_name='loop')
    # accungs_loop_muts_correlation(run_name='v3_q10_no_indel', accungs_source_name= 'accungs_new_no_indel', loop_source_name='loop_new_q10_no_indel')
    # accungs_loop_muts_correlation(run_name='v4_long_reads_no_indel', accungs_source_name= 'accungs_new_no_indel', loop_source_name='loop_new_long_reads_no_indel')