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
    output_path_accungs = '/sternadi/nobackup/volume1/HIVB_shafer_accungs_seq/new_pipeline_accungs_v2_15_iter/env{}'

    # q10
    # input_path_loop = '/sternadi/datasets/volume2/HIV_shafer_loop_genomics/fastq_only/s{}'
    # output_path_loop = '/sternadi/home/volume1/shared/analysis/HIV_shafer/new_pipeline_loop_v2_q10/env{}'
    # long reads
    input_path_loop = '/sternadi/datasets/volume2/HIV_shafer_loop_genomics/fastq_long_only/s{}'
    output_path_loop = '/sternadi/nobackup/volume1/HIVB_shafer_accungs_seq/new_pipeline_loop_v3_long_reads/env{}'

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
                                "with_indels": "Y",
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
                            # "quality_threshold": "10",  # TODO- trying other Q scores
                            "max_basecall_iterations": 10}
        loop_run_param_list.append(sample_dict_loop)

    print(len(accungs_run_param_list))
    print(len(loop_run_param_list))

    # run
    # multi_runner(accungs_run_param_list)
    multi_runner(loop_run_param_list)


# TODO- from sheri- dont take these results too seriously
def accungs_loop_muts_correlation():
    from loop_genomics_haplotype_analysis.summarize_visualize_strain_results import env_samples, shafer_root_dir
    from loop_genomics_haplotype_analysis.mutation_analysis import get_freqs

    print('* accungs_loop_muts_correlation *')
    # get filtered freqs
    freqs_filtered = get_freqs()

    # 1. plot per-mut correlation
    plot_correlation = True
    if plot_correlation:
        freqs_accungs = freqs_filtered[freqs_filtered['source'] == 'accungs']
        freqs_loop = freqs_filtered[freqs_filtered['source'] == 'loop']
        print('len(freqs_accungs): {}'.format(len(freqs_accungs)))
        print('len(freqs_loop): {}'.format(len(freqs_loop)))
        print(freqs_accungs.groupby('sample').count()['Pos'])
        print(freqs_loop.groupby('sample').count()['Pos'])

        # align 2 cons coordinates
        # no need- all files are without insertions.

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
                       # height=3,
                       data=freq_set);

            plt.xlim(axis[0], axis[1])
            plt.ylim(axis[0], axis[1])
            # plt.xscale('log')
            # plt.yscale('log')
            # g.set_title(name)
            plt.savefig(fname='{}/results_merged/accungs_loop_mutation_correlation_linear_{}_v2.pdf'.format(shafer_root_dir, name))
            plt.show()

    # 2. boxplot freq distribution
    freq_boxplot = False
    if freq_boxplot:
        print('boxplot')
        #  1. regular
        # v1
        boxplot(freqs_filtered, '{}/results_merged/accungs_loop_freq_medians.pdf')
        # v2
        # mutations_box_plots_itamar(freqs_filtered, '', True)

        #  2. non-ox too (CA\GT)
        print(len(freqs_filtered))
        freqs_nonox = freqs_filtered[~freqs_filtered['mutation'].isin(['C->A', 'G->T'])]
        print(len(freqs_nonox))
        boxplot(freqs_filtered, '{}/results_merged/accungs_loop_freq_medians_nonox.pdf')
        #  3. high\low freqs?
        #  4. syn only - requires HXB2 coords

    # 3. mut-type count avg. (per sample, per seq-type)
    # I did count per strand, can also do per seq type (accungs is more but we'll see how much)
    # TODO- imp
    # freqs_filtered = freqs_filtered.pivot(index='mutation', columns='source', values='mean')
    freqs_filtered = freqs_filtered.groupby(index='mutation', columns='source', values='mean')

    # 4. mut-freq distribution, per seq-type
    # TODO- imp
    freq_dist_all_samples = True
    if freq_dist_all_samples:
        freqs_filtered.groupby('source')['Freq'].hist(bins=50, log=True, alpha=0.8, legend = True)
        plt.savefig(fname='{}/results_merged/accungs_loop_freq_dist_agg.pdf'.format(shafer_root_dir))
        plt.show()


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


if __name__ == '__main__':
    run_new_pipeline_on_server()

    # accungs_loop_muts_correlation()